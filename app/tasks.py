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
