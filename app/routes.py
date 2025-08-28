# routes.py

import os
import subprocess
import json
import time
import uuid
import numpy as np
from PIL import Image
from scipy.ndimage import rotate
from flask import Blueprint, render_template, request, url_for, jsonify, session, current_app
import re
import shutil


from app.utils import (
    get_upload_folder,
    get_user_upload_folder,
    get_user_output_folder,
    cleanup_session_folders,
    generate_peaks_from_bigwig_macs2,
    get_bigwig_signal,
    save_uploaded_file,
    normalize_bigwig,
    predict_ctcf,
    predict_peaks,
    prepare_chimeric_gene_track_config,
    prepare_gene_track_config,
    WINDOW_WIDTH,
    clear_folder,
    assemble_chimeric_arrays,
    write_chimeric_bigwig,
    extract_bigwig_region_array,
    assemble_chimeric_fasta
)
from app.tasks import q
###############################################################################
# Blueprint & Globals
###############################################################################
main = Blueprint('main', __name__)
test_mode = False

@main.before_request
def before_request():
    # print("BEFORE REQUEST: session_id =", session.get('session_id'))
    # print("BEFORE REQUEST: current_job_id =", session.get('current_job_id'))
    cleanup_session_folders()


###############################################################################
# Standard Plot Config Preparation
###############################################################################
def prepare_plot_configs(
    hi_c_matrix,
    region_chr,
    region_start,
    region_end,
    ds_option,
    ctcf_bw_for_model,
    atac_bw_for_model,
    del_start=None,
    del_width=None
):
    """
    Renders the final charts for standard (non-chimeric) or single-chrom
    predictions, or simplified for 'chrCHIM'. `raw_atac_path` is used
    for the ATAC track, `ctcf_bw_for_model` for the CTCF track.
    """
    import numpy as np
    from scipy.ndimage import rotate
    print(f"prepare_plot_configs received: atac_bw_for_model={atac_bw_for_model}")


    # 1) Remove rows and columns that are entirely zero BEFORE rotation
    # nonzero_row_mask = ~np.all(hi_c_matrix == 0, axis=1)
    # hi_c_matrix = hi_c_matrix[nonzero_row_mask, :]
    # nonzero_col_mask = ~np.all(hi_c_matrix == 0, axis=0)
    # hi_c_matrix = hi_c_matrix[:, nonzero_col_mask]

    # 2) Rotate the Hi-C matrix by 45 degrees, filling with zeros
    hi_c_matrix = rotate(
        hi_c_matrix,
        angle=45,
        reshape=True,
        mode='constant',  # constant fill
        cval=0
    )

    # 3) Remove rows/columns that are "near-zero" after rotation
    threshold = 1e-6
    row_mask_after = ~np.all(hi_c_matrix < threshold, axis=1)
    hi_c_matrix = hi_c_matrix[row_mask_after, :]
    col_mask_after = ~np.all(hi_c_matrix < threshold, axis=0)
    hi_c_matrix = hi_c_matrix[:, col_mask_after]

    # 4) Keep only the "top" (array-index) half of the rotated matrix
    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[: num_rows // 2, :]

    # 5) Remove the first 22 rows (the "top" 22 pixels in array terms)
    if (region_end - region_start) > WINDOW_WIDTH:
        print(f'ACCORDING to prepare_plot_configs, region_end - region_start = {region_end - region_start} and WINDOW_WIDTH = {WINDOW_WIDTH}')
        if hi_c_matrix.shape[0] > 28:
            hi_c_matrix = hi_c_matrix[28:, :]

    # Build the final data array for the Hi-C plot
    n_rows, n_cols = hi_c_matrix.shape
    hi_c_data = []
    for i in range(n_rows):
        for j in range(n_cols):
            hi_c_data.append([j, i, float(hi_c_matrix[i, j])])

    # Convert region to Mb
    x_start_mb = region_start / 1e6
    x_end_mb = region_end / 1e6

    hi_c_height = 170 if (region_end - region_start) > WINDOW_WIDTH else 200
    # Hi-C chart config
    hi_c_chart_config = {
        "chart": {"height": hi_c_height},
        "xAxis": {"min": x_start_mb, "max": x_end_mb, "title": ""},
        "yAxis": {"min": 0, "max": 100},
        "colorAxis": {
            "min": 0,
            "max": float(np.nanmax(hi_c_matrix)) if hi_c_matrix.size else 0.0,
        },
        "series": [{"name": "Hi-C Signal", "data": hi_c_data}],
        "n_cols": n_cols,
        "hideXAxis": True
    }

    # CTCF track
    ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_for_model, region_chr, region_start, region_end)
    ctcf_positions_mb = [p / 1e6 for p in ctcf_positions]
    ctcf_chart_config = {
        "chart": {"height": 100},
        "xAxis": {"min": x_start_mb, "max": x_end_mb, "title": ""},
        "yAxis": {},
        "series": [{
            "name": "CTCF Signal",
            "data": [[x, y] for x, y in zip(ctcf_positions_mb, ctcf_values)],
            "color": "#779ECC"
        }]
    }

    # ATAC track
    atac_positions, atac_values = get_bigwig_signal(atac_bw_for_model, region_chr, region_start, region_end)
    atac_positions_mb = [p / 1e6 for p in atac_positions]
    atac_chart_config = {
        "chart": {"height": 100},
        "xAxis": {"min": x_start_mb, "max": x_end_mb, "title": ""},
        "yAxis": {},
        "series": [{
            "name": "ATAC Signal",
            "data": [[x, y] for x, y in zip(atac_positions_mb, atac_values)],
            "color": "#FF985A"
        }]
    }

    # If 'deletion' => hide x-axis on Hi-C and add axisBreak to CTCF/ATAC
    if ds_option == "deletion" and del_start is not None and del_width is not None:
        hi_c_chart_config["hideXAxis"] = True
        del_start_mb = del_start / 1e6
        del_end_mb = (del_start + del_width) / 1e6

        for chart in (ctcf_chart_config, atac_chart_config):
            chart["xAxis"] = {
                "axisBreak": True,
                "min": x_start_mb,
                "max": x_end_mb,
                "leftDomain": [x_start_mb, del_start_mb],
                "rightDomain": [del_end_mb, x_end_mb],
                "title": ""
            }
    print(f"prepare_plot_configs: atac_bw={atac_bw_for_model}, ctcf_bw={ctcf_bw_for_model}")
    return hi_c_chart_config, ctcf_chart_config, atac_chart_config, None, {}


###############################################################################
# NEW HELPER: Unified "Prediction + Render" function
###############################################################################

def run_prediction_and_render(
    region_chr,
    region_start,
    region_end,
    ds_option,
    model_path,
    genome,
    seq_dir,
    atac_bw_for_model,
    ctcf_bw_for_model,
    output_folder,
    del_start=None,
    del_width=None,
    screening_requested=False # for chimeric mode to pass chrCHIM.fa.gz folder
):
    """
    1. Runs editing.py (if ds_option='deletion') or prediction.py otherwise.
    2. Waits for result.npy, loads => hi_c_matrix.
    3. Builds chart configs.
    4. Optionally sets screening config.
    5. Returns dict with chart configs in JSON form.
    """
    print(f"run_prediction_and_render received: atac_bw_for_model={atac_bw_for_model}")

    BASE_DIR        = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "..", "C.Origami", "src")
    env             = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH

    script_path       = os.path.join(PYTHON_SRC_PATH, "corigami", "inference")
    hi_c_matrix_path  = os.path.join(output_folder, "result.npy")

    # Decide which seq_dir to us

    # 1) Launch editing.py or prediction.py (as an RQ job)
    from app.tasks import run_editing_task, run_prediction_task

    if ds_option == "deletion" and del_start is not None and del_width is not None:
        editing_script = os.path.join(script_path, "editing.py")
        run_editing_task(
            editing_script,
            region_chr,
            region_start,
            model_path,
            seq_dir,
            atac_bw_for_model,
            del_start,
            del_width,
            output_folder,
            ctcf_bw_for_model=ctcf_bw_for_model,
            env=env,
        )
    else:
        prediction_script = os.path.join(script_path, "prediction.py")
        run_prediction_task(
            prediction_script,
            region_chr,
            region_start,
            region_end,
            model_path,
            seq_dir,
            atac_bw_for_model,
            output_folder,
            ctcf_bw_for_model=ctcf_bw_for_model,
            env=env,
        )

    # By here, the job is done and 'result.npy' should be fully written.
    if not os.path.exists(hi_c_matrix_path):
        raise RuntimeError(f"Hi-C matrix not found at {hi_c_matrix_path}")

    hi_c_matrix = np.load(hi_c_matrix_path)

    if region_chr == "chrCHIM":
        atac_bw_for_model = os.path.join(output_folder, "chrCHIM_atac.bw")
        ctcf_bw_for_model = os.path.join(output_folder, "chrCHIM_ctcf.bw")
        print(f"[FORCE] chrCHIM plotting inputs → ATAC={atac_bw_for_model}  CTCF={ctcf_bw_for_model}")
    # 3) Build chart configs
    hi_c_config, ctcf_config, atac_config, _, _ = prepare_plot_configs(
        hi_c_matrix,
        region_chr,
        0 if region_chr=="chrCHIM" else region_start,
        region_end,
        ds_option,
        ctcf_bw_for_model,
        atac_bw_for_model,
        del_start,
        del_width
    )

    # 4) Build gene track config
    gene_track_config = prepare_gene_track_config(
        genome,
        region_chr,
        0 if region_chr=="chrCHIM" else region_start,
        region_end,
        ds_option,
        del_start,
        del_width
    )

    # 5) Screening config (placeholder)
    screening_config_json = None

        # ------------------------------------------------------------------
    # >>> NEW: stash everything in the current RQ job’s meta so that
    #          /api/job/<id>/html can render the template later.
    # ------------------------------------------------------------------
    from rq import get_current_job
    _job = get_current_job()
    if _job:
        _job.meta.update(
            hi_c_config      = hi_c_config,
            ctcf_config      = ctcf_config,
            atac_config      = atac_config,
            gene_track_config= gene_track_config,
            screening_config = screening_config_json,
        )
        _job.save_meta()

    return {
        "hi_c_config": hi_c_config,  # (a Python dict, not a string!)
        "ctcf_config": ctcf_config,
        "atac_config": atac_config,
        "gene_track_config": gene_track_config,
        "screening_config":  screening_config_json
    }

###############################################################################
# Main Page: GET => form, POST => run standard or chimeric logic
###############################################################################
# routes.py (excerpt) -- remove references to build_chimeric.py or build_chimeric_task

@main.route('/', methods=['GET', 'POST'])
def index():
    """
    Landing page (GET) or legacy synchronous workflow (POST).
    For POST we still run everything in‑process – useful for quick
    single‑machine setups without Redis – by re‑using exactly the same
    `prepare_inputs()` helper the JSON/RQ path uses.
    """
    output_folder = get_user_output_folder()

    # ------------------------------  GET  ------------------------------
    if request.method == 'GET':
        return render_template(
            "index.html",
            screening_mode=False,
            output_folder=output_folder
        )

    # ------------------------------  POST  -----------------------------
    # Make sure user‑specific output dir is fresh
    clear_folder(output_folder)
    os.makedirs(output_folder, exist_ok=True)

    # Consolidate all preprocessing / normalisation / chimeric logic
    run_args = prepare_inputs(request.form.to_dict(flat=True), output_folder)

    # Perform prediction + plot‑config building **in the web process**
    # (If you prefer the RQ path even here, simply enqueue instead.)
    configs = run_prediction_and_render(**run_args)

    # Extra configs used only for HTML rendering
    custom_axis_cfg, gene_track_cfg = build_axis_and_gene_configs(
        run_args, request.form   # tiny helper, omitted for brevity
    )

    # Decide whether the browser asked for just the partial fragment
    tmpl = (
        "plots_partial.html"
        if request.headers.get("X-Requested-With") == "XMLHttpRequest"
        else "index.html"
    )

    return render_template(
        tmpl,
        hi_c_config        = configs["hi_c_config"],
        ctcf_config        = configs["ctcf_config"],
        atac_config        = configs["atac_config"],
        screening_config   = configs["screening_config"],
        screening_mode     = run_args["screening_requested"],
        screening_params   = "{}",
        gene_track_config  = gene_track_cfg,
        user_output_folder = output_folder,
        custom_axis_config = custom_axis_cfg,
    )
###############################################################################
# Screening Route
###############################################################################
@main.route('/run_screening', methods=['GET'])
def run_screening_endpoint():
    from app.tasks import q, run_screening_task
    from app.utils import (
        generate_peaks_from_bigwig_macs2,
        get_user_output_folder
    )

    # 1) Retrieve final BigWigs and region info from session
    final_atac_bw  = session.get('final_atac_bw')
    final_ctcf_bw  = session.get('final_ctcf_bw')
    region_chr     = session.get('final_region_chr')
    region_start   = session.get('final_region_start')
    region_end     = session.get('final_region_end')
    seq_dir        = session.get('seq_dir')
    perturb_width  = session.get('perturb_width')
    model_path     = session.get('model_path')
    peaks_file     = session.get('peaks_file')
    print(f'screening starting with seq_dir: {seq_dir}')

    # If any are missing => user did not run the main route first
    if not final_atac_bw or not final_ctcf_bw or region_chr is None:
        return jsonify({"error": "No existing data from main run. Please run prediction first."}), 400

    # 2) Extract additional parameters for how to run the screening
    # e.g. perturb_width, model selection
    output_dir    = get_user_output_folder()


    print(f'run screening endpoint: peaks file: {peaks_file}')
    # 4) Now do the screening with final_atac_bw, final_ctcf_bw, peaks_file
    # (No re-assembly, no normalization, no chimeric logic needed)
    screening_script = os.path.join("C.Origami", "src", "corigami", "inference", "screening.py")
    env = os.environ.copy()
    env["PYTHONPATH"] = "C.Origami/src"

    job_screen = q.enqueue(
        run_screening_task,
        screening_script,
        region_chr,
        region_start,
        region_end,
        model_path,
        seq_dir,
        final_atac_bw,
        final_ctcf_bw,
        peaks_file,    # 9th
        output_dir,    # 10th
        perturb_width,
        env=env,
        job_timeout=2000
    )
    # 5) Wait for completion (or do a non-blocking approach)
    while not job_screen.is_finished and not job_screen.is_failed:
        time.sleep(1)
        job_screen.refresh()

    if job_screen.is_failed:
        return jsonify({"error": "Screening script failed"}), 500

    # 6) Load results, build chart config, respond
    results_file = os.path.join(output_dir, "screening", "screening_results.json")
    if not os.path.exists(results_file):
        return jsonify({"error": "No screening_results.json found."}), 500

    with open(results_file, "r") as f:
        results = json.load(f)

    # e.g. build chart config
    window_mb = results["window_midpoints_mb"]
    impact   = results["impact_scores"]
    series_data = [[x,y] for x,y in zip(window_mb, impact)]
    screening_chart = {
      "screening": True,
      "chart": {"type": "line", "height": 100},
      "xAxis": {"min": region_start/1e6, "max": region_end/1e6},
      "yAxis": {},
      "series": [{
          "name": "Screening Score",
          "data": series_data,
          "color": "#9FC0DE"
      }]
    }
    if session.get('chimeric_active'):
        print("screening route thinks chimeric is active")
        screening_chart["isChimeric"] = True
    else:
        print("screening route thinks chimeric is not active")
        screening_chart["isChimeric"] = False

    results["screening_config"] = screening_chart
    
    return jsonify(results)

###############################################################################
# File Upload
###############################################################################
@main.route('/upload_file', methods=['POST'])
def upload_file():
    file_type = request.args.get('file_type', '')
    if not file_type:
        return jsonify({"error": "Missing file_type parameter"}), 400

    upload_folder = get_upload_folder(file_type)

    if file_type == "atac":
        file_key = "atac_bw_file"
    elif file_type == "ctcf":
        file_key = "ctcf_bw_file"
    elif file_type == "peaks":
        file_key = "peaks_file"
    else:
        return jsonify({"error": "Unknown file_type"}), 400

    file_storage = request.files.get(file_key)
    if not file_storage:
        return jsonify({"error": "No file provided"}), 400

    saved_path = save_uploaded_file(file_storage, file_storage.filename, upload_folder)
    print(f"Uploaded {file_type} file saved to: {saved_path}")
    raw_name = os.path.basename(saved_path)
    display_name = raw_name.split("_", 1)[-1] if "_" in raw_name else raw_name
    return jsonify({"saved_path": saved_path, "display_name": display_name})

###############################################################################
# List Uploads
###############################################################################
@main.route('/list_uploads', methods=['GET'])
def list_uploads():
    file_type = request.args.get('file_type', '')
    folder = get_upload_folder(file_type) if file_type else get_user_upload_folder()
    files = []
    if os.path.exists(folder):
        for f in os.listdir(folder):
            full_path = os.path.join(folder, f)
            if os.path.isfile(full_path):
                display_name = f.split("_", 1)[-1] if "_" in f else f
                files.append({"value": full_path, "name": display_name})
    return jsonify(files)


###


@main.route("/health", methods=["GET"])
def health():
    return "OK", 200

@main.before_request
def before_request():
    if request.path == "/health":
        return           # skip heavy work for the ALB check
    cleanup_session_folders()


# ── AWS worker‑control endpoints (only if USE_AWS=1) ───────────
import os
from flask import jsonify, current_app
from botocore.exceptions import ClientError

USE_AWS = os.getenv("USE_AWS") == "1"

if USE_AWS:
    import boto3

    region = os.getenv("AWS_REGION", "us-east-2").replace("\u2011", "-")
    ecs = boto3.client("ecs", region_name=region)
    asg = boto3.client("autoscaling", region_name=region)

    CLUSTER = "corigami"
    SERVICE = "corigami-cpu-rq-worker-service"
    ASG     = "corigami-cpu-asg"

    @main.route("/api/worker-status")
    def worker_status():
        svc = ecs.describe_services(
            cluster=CLUSTER,
            services=[SERVICE]
        )["services"][0]

        # True if at least one task is already RUNNING
        running  = svc.get("runningCount", 0) > 0
        # True if at least one task is PENDING
        starting = svc.get("pendingCount", 0) > 0

        current_app.logger.debug(
            f"ECS svc desired={svc['desiredCount']} "
            f"pending={svc['pendingCount']} running={svc['runningCount']}"
        )
        return jsonify({"running": running, "starting": starting})

    @main.route("/api/start-worker", methods=["POST"])
    def start_worker():
        # Step 1: bump ECS desiredCount → this immediately sets pendingCount=1
        try:
            ecs.update_service(cluster=CLUSTER, service=SERVICE, desiredCount=1)
        except ClientError as e:
            current_app.logger.error(f"ECS update_service failed: {e}")
            return jsonify({"error": str(e)}), 500

        # Step 2: also ensure ASG can add capacity downstream
        try:
            asg.set_desired_capacity(
                AutoScalingGroupName=ASG,
                DesiredCapacity=1,
                HonorCooldown=True
            )
        except ClientError as e:
            err = e.response.get("Error", {})
            code = err.get("Code", "")
            if code != "ScalingActivityInProgress":
                current_app.logger.error(f"ASG set_desired_capacity failed: {e}")
                return jsonify({"error": str(e)}), 500

        return jsonify({"message": "starting"}), 200

else:
    # Dummy endpoints so the front‑end never 404s when AWS is off
    @main.route("/api/worker-status")
    def worker_status():
        return jsonify({"running": False, "starting": False})

    @main.route("/api/start-worker", methods=["POST"])
    def start_worker():
        return jsonify({"message": "AWS disabled"}), 501
 
