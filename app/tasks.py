print(">>> loading app.tasks …")

import os
import time
import subprocess
import uuid
import math
import numpy as np
import boto3
from redis import Redis
from rq import Queue, get_current_job

# ─── ECS INTEGRATION ───────────────────────────────────────────────────────────

ecs = boto3.client(
    "ecs",
    region_name=os.getenv("AWS_REGION", "us-east-2")
)

# These must be set as env vars in your corigami-web (or worker) task definition:
CLUSTER         = os.getenv("ECS_CLUSTER")             # e.g. "corigami-gpu-cluster"
CPU_PROVIDER    = os.getenv("CPU_CAP_PROVIDER")        # e.g. "corigami-cpu-capacity-provider"
TASK_DEF        = os.getenv("CPU_TASK_DEF")            # e.g. "corigami-cpu-rq-worker:1"
SUBNETS         = os.getenv("ECS_SUBNETS", "").split(",")
SECURITY_GROUPS = os.getenv("ECS_SGS", "").split(",")

def run_ecs_prediction_task(payload_bytes):
    """
    Launch an ECS EC2 task on your CPU capacity provider,
    passing the raw payload bytes in an environment var.
    Returns the launched taskArn.
    """
    resp = ecs.run_task(
        cluster=CLUSTER,
        launchType="EC2",
        capacityProviderStrategy=[{
            "capacityProvider": CPU_PROVIDER,
            "weight": 1
        }],
        taskDefinition=TASK_DEF,
        networkConfiguration={
            "awsvpcConfiguration": {
                "subnets": SUBNETS,
                "securityGroups": SECURITY_GROUPS,
                "assignPublicIp": "DISABLED"
            }
        },
        overrides={
            "containerOverrides": [{
                "name": "worker",
                "environment": [
                    {"name": "REDIS_URL", "value": os.getenv("REDIS_URL")},
                    {"name": "PAYLOAD_BYTES", "value": payload_bytes}
                ]
            }]
        }
    )
    return resp["tasks"][0]["taskArn"]

def wait_for_ecs_task(task_arn, timeout=1800):
    """
    Poll until the ECS task stops. Raises on non-zero exit code or timeout.
    """
    deadline = time.time() + timeout
    while time.time() < deadline:
        desc = ecs.describe_tasks(cluster=CLUSTER, tasks=[task_arn])
        task = desc["tasks"][0]
        status = task["lastStatus"]
        if status == "STOPPED":
            code = task["containers"][0].get("exitCode")
            if code == 0:
                return True
            else:
                raise RuntimeError(f"ECS task failed (exitCode={code})")
        time.sleep(2)
    raise TimeoutError("Timed out waiting for ECS task to finish")

# ─── RQ WORKER SETUP ────────────────────────────────────────────────────────────
from redis import Redis
from rq    import Queue

# Prefer the Valkey endpoint from the task-definition env-vars.
# Falls back to localhost only when you’re running it on a laptop.
redis_url = (
    os.getenv("RQ_REDIS_URL")      # primary
    or os.getenv("REDIS_URL")      # secondary
    or "redis://localhost:6379/0"  # dev fallback
)

redis_conn = Redis.from_url(redis_url)
q = Queue("default", connection=redis_conn)

def run_cmd(cmd, env=None, cwd=None, timeout=None):
    """
    Execute a subprocess, streaming its output to the worker logs,
    while respecting RQ job cancellation and an optional timeout.
    """
    print(f"[RQ-Worker] Executing: {' '.join(cmd)}")
    process = subprocess.Popen(cmd, env=env, cwd=cwd)
    start = time.time()

    while True:
        ret = process.poll()
        if ret is not None:
            break

        job = get_current_job()
        if job and job.get_status() == "canceled":
            process.kill()
            raise RuntimeError("Job was canceled by the user.")

        if timeout and (time.time() - start) > timeout:
            process.kill()
            raise RuntimeError("Timed out: " + " ".join(cmd))

        time.sleep(0.2)

    if process.wait() != 0:
        raise RuntimeError(f"Command failed (exit code {process.returncode}): {' '.join(cmd)}")

# ─── TASK 1: run_editing ────────────────────────────────────────────────────────

def run_editing_task(
    editing_script,
    region_chr,
    region_start,
    model_path,
    seq_dir_for_prediction,
    atac_bw_for_model,
    del_start,
    del_width,
    output_folder,
    ctcf_bw_for_model=None,
    env=None
):
    cmd = [
        "python", editing_script,
        "--chr", region_chr,
        "--start", str(region_start),
        "--model", model_path,
        "--seq", seq_dir_for_prediction,
        "--atac", atac_bw_for_model,
        "--del-start", str(del_start),
        "--del-width", str(del_width),
        "--out", output_folder,
    ]
    if ctcf_bw_for_model:
        cmd += ["--ctcf", ctcf_bw_for_model]

    run_cmd(cmd, env=env)
    return True

# ─── TASK 2: run_prediction ───────────────────────────────────────────────────

def run_prediction_task(
    prediction_script,
    region_chr,
    region_start,
    region_end,
    model_path,
    seq_dir_for_prediction,
    atac_bw_for_model,
    output_folder,
    ctcf_bw_for_model=None,
    env=None
):
    """
    Picks between local RQ invocation or ECS based on
    CORIGAMI_USE_ECS_env var set in your web container.
    """
    use_ecs = os.getenv("CORIGAMI_USE_ECS", "false").lower() == "true"
    if use_ecs:
        # Prepare the raw bytes payload
        from corigami.inference.utils import inference_utils as infer
        seq_arr, ctcf_arr, atac_arr = infer.load_region(
            region_chr,
            region_start,
            seq_dir_for_prediction,
            ctcf_bw_for_model,
            atac_bw_for_model
        )
        payload = np.stack([seq_arr, ctcf_arr, atac_arr], axis=0).astype("float32").tobytes()

        # Run on ECS and wait
        task_arn = run_ecs_prediction_task(payload)
        wait_for_ecs_task(task_arn)

        # Once done, result.npy will be written to your EFS mount at user_data
        return True

    # Fallback to local subprocess + RQ
    print(f"[RQ-Worker] Sending start={region_start}, end={region_end} to prediction.py")
    cmd = [
        "python", prediction_script,
        "--chr", region_chr,
        "--start", str(region_start),
        "--end",   str(region_end),
        "--model", model_path,
        "--seq",   seq_dir_for_prediction,
        "--atac",  atac_bw_for_model,
        "--out",   output_folder,
    ]
    if ctcf_bw_for_model:
        cmd += ["--ctcf", ctcf_bw_for_model]

    run_cmd(cmd, env=env)
    return True

# ─── TASK 4: run_screening ────────────────────────────────────────────────────
def _ensure_interior_peaks(peaks_file, chrom, screen_start, screen_end, window_bp, outdir):
    """
    Keep only peaks whose midpoint is >= window_bp/2 away from both edges so a full window fits.
    If nothing survives, seed a few interior peaks so screening can run.
    Returns the path to the filtered peaks file.
    """
    try:
        if not (peaks_file and os.path.exists(peaks_file)):
            return peaks_file

        os.makedirs(outdir, exist_ok=True)
        out_path = os.path.join(outdir, "auto_peaks_interior.narrowPeak")

        min_mid = int(screen_start + window_bp // 2)
        max_mid = int(screen_end   - window_bp // 2)

        kept = 0
        with open(peaks_file) as fin, open(out_path, "w") as fout:
            for line in fin:
                if not line.strip():
                    continue
                parts = line.split("\t")
                if parts[0] != chrom:
                    continue
                try:
                    start = int(parts[1]); end = int(parts[2])
                except ValueError:
                    continue
                mid = (start + end) // 2
                if min_mid <= mid <= max_mid:
                    fout.write(line)
                    kept += 1

        if kept == 0:
            # seed a few interior peaks to force at least some windows
            step = 400_000
            seed_half = 150
            with open(out_path, "w") as fout:
                pos = min_mid + 200_000
                while pos <= max_mid - 200_000:
                    s = max(screen_start, pos - seed_half)
                    e = min(screen_end,   pos + seed_half)
                    fout.write(f"{chrom}\t{s}\t{e}\tseed\t1000\t.\n")
                    pos += step
            print(f"[worker] no interior peaks; seeded synthetic peaks in {out_path}")
        else:
            print(f"[worker] interior peaks kept: {kept} -> {out_path}")

        return out_path
    except Exception as e:
        print(f"[worker] _ensure_interior_peaks error: {e}; using original {peaks_file}")
        return peaks_file 

def run_screening_task(
    screening_script,
    region_chr,
    screen_start,
    screen_end,
    model_path,
    seq_dir,
    atac_bw_path,
    ctcf_bw_path,
    peaks_file,
    output_dir,
    perturb_width,
    env=None
):
    if not ctcf_bw_path or not os.path.isfile(ctcf_bw_path):
        raise ValueError(f"CTCF file is required: {ctcf_bw_path}")
    if not peaks_file or not os.path.isfile(peaks_file):
        raise ValueError(f"Peaks file is required: {peaks_file}")

    cmd = [
        "python", screening_script,
        "--no-server",
        "--chr", region_chr,
        "--screen-start", str(screen_start),
        "--screen-end",   str(screen_end),
        "--model", model_path,
        "--seq",   seq_dir,
        "--atac",  atac_bw_path,
        "--out",   output_dir,
        "--perturb-width", str(perturb_width),
        "--plot-impact-score",
        "--save-pred",
        "--save-perturbation",
        "--save-diff",
        "--save-bedgraph",
        "--ctcf",       ctcf_bw_path,
        "--peaks-file", peaks_file
    ]

    run_cmd(cmd, env=env)
    return True

def preprocessing_task(form_dict):
    """
    Enqueued by the API.  Runs all heavy utils code inside the worker.
    """
    from app.preprocessing import preprocessing
    result = preprocessing(form_dict, use_session=False)
    return result   # RQ serialises the dict; the API can fetch job.result

# ─── TASK X:  worker_predict_and_norm  ────────────────────────────────────────
"""
   One-shot helper that performs every heavy step inside the worker:

      • prepare_inputs()  – normalisation, CTCF prediction, MACS2 peaks …
      • run_prediction_and_render()  – model inference + chart configs

   All resulting artefacts (plot configs, axis config, etc.) are stashed
   in job.meta so the /api/job/<id>/html route can render them later.
"""
import uuid, os, json, numpy as np
from app.api import prepare_inputs    

def worker_predict_and_norm(form_json: dict, outdir: str) -> bool:
    """
    Phase‑1 job.
      • heavy preprocessing + fast 2 Mb prediction
      • stores plot configs in parent‑job meta
      • keeps the screening_job_id (written by /api/predict) intact
    """
    import shutil
    from rq import get_current_job
    from app.api   import prepare_inputs
    from app.routes import run_prediction_and_render

    # ── fresh /output folder ─────────────────────────────────────────────
    if os.path.exists(outdir):
        for p in os.listdir(outdir):
            fp = os.path.join(outdir, p)
            shutil.rmtree(fp) if os.path.isdir(fp) else os.unlink(fp)
    else:
        os.makedirs(outdir, exist_ok=True)

    # ── heavy preprocessing + prediction ────────────────────────────────
    run_args, axis_cfg, gene_cfg = prepare_inputs(form_json, outdir)
    result   = run_prediction_and_render(**run_args)

    # ── update job.meta without losing the child job id ─────────────────
    job            = get_current_job()
    existing_child = job.meta.get("screening_job_id")        # preserve if present

    job.meta.update(
        **result,                         # hi_c_config, atac_config, …
        custom_axis_config = axis_cfg,
        screening_mode     = (run_args["ds_option"] == "screening"),
        run_args           = run_args,
        form_json          = form_json,
        outdir             = outdir,
    )

    if existing_child:                    # re‑insert after bulk update
        job.meta["screening_job_id"] = existing_child

    job.save_meta()
    return True


def worker_run_screening(parent_job_id: str, **kwargs) -> bool:
    """
    Phase‑2 job (runs after worker_predict_and_norm).

    • reads paths & args from the parent job’s meta
    • executes corigami screening
    • writes the finished Highcharts config into **both**
      – parent.meta  (so /api/job/<id>/html can embed it immediately)
      – child.meta   (so /api/run_screening can fetch it directly)
    """
    import json, os
    from rq import get_current_job
    from rq.job import Job
    from app.tasks import run_screening_task
    from app.utils import generate_peaks_from_bigwig_macs2

    # ── look up parent job & its stored data ────────────────────────────
    parent = Job.fetch(parent_job_id, connection=get_current_job().connection)
    meta   = parent.meta

    run_args  = meta["run_args"]
    form_json = meta["form_json"]
    outdir    = meta["outdir"]

    # ── ensure a peaks file exists ──────────────────────────────────────
    peaks_file = form_json.get("peaks_file_path")
    if not peaks_file or peaks_file.lower() == "none":
        peaks_file, _ = generate_peaks_from_bigwig_macs2(
            bw_path = run_args["atac_bw_for_model"],
            chrom   = run_args["region_chr"],
            start   = run_args["region_start"],
            end     = run_args["region_end"],
            outdir  = outdir,
        )
        # NEW: for chrCHIM, keep only interior peaks so windows pass the edge filter
    if run_args["region_chr"] == "chrCHIM" and peaks_file:
        # use the model window size (2,097,152 bp), not the screen span
        peaks_file = _ensure_interior_peaks(
            peaks_file = peaks_file,
            chrom = run_args["region_chr"],
            screen_start = run_args["region_start"],
            screen_end   = run_args["region_end"],
            window_bp = 2_097_152,
            outdir = os.path.join(outdir, "screening"),
        )
        print(f"[worker] screening will use interior peaks file: {peaks_file}")

    # ── run the screening script ────────────────────────────────────────
    screening_script = os.path.join(
        "C.Origami", "src", "corigami", "inference", "screening.py"
    )
    run_screening_task(
        screening_script,
        run_args["region_chr"],
        run_args["region_start"],
        run_args["region_end"],
        run_args["model_path"],
        run_args["seq_dir"],
        run_args["atac_bw_for_model"],
        run_args["ctcf_bw_for_model"],
        peaks_file,
        outdir,
        perturb_width = int(form_json.get("perturb_width", 2000)),
        env = os.environ.copy() | {"PYTHONPATH": "C.Origami/src"},
    )

    # ── build the Highcharts config ─────────────────────────────────────
    scr_json = os.path.join(outdir, "screening", "screening_results.json")
    with open(scr_json) as fh:
        scr = json.load(fh)

    if "screening_config" not in scr:
        x = scr["window_midpoints_mb"]
        y = scr["impact_scores"]
        scr["screening_config"] = {
            "chart": {"type": "line", "height": 100},
            "xAxis": {
                # use EXACTLY the same domain the other charts got
                "min": run_args["region_start"] / 1e6,
                "max": run_args["region_end"]   / 1e6
            },
            "yAxis": {},
            "series": [{
                "name": "Screening score",
                "data": [[xi, yi] for xi, yi in zip(x, y)],
                "color": "#9FC0DE"
            }]
        }

    # ── write into BOTH metas so nothing can clobber it ─────────────────
    parent.meta["screening_config"] = scr["screening_config"]
    parent.save_meta()

    job = get_current_job()                # <- this child job
    job.meta["screening_config"] = scr["screening_config"]
    job.save_meta()

    return True
