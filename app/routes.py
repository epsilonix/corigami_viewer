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

from app.utils import (
    get_upload_folder,
    get_user_upload_folder,
    get_user_output_folder,
    cleanup_session_folders,
    generate_peaks_from_bigwig_macs2,
    get_bigwig_signal,
    save_uploaded_file,
    normalize_atac,
    normalize_ctcf,
    predict_ctcf,
    predict_peaks,
    prepare_chimeric_gene_track_config,
    prepare_gene_track_config,
    WINDOW_WIDTH
)
from app.tasks import q
###############################################################################
# Blueprint & Globals
###############################################################################
main = Blueprint('main', __name__)
test_mode = False

@main.before_request
def before_request():
    """
    Cleanup session folders before each request.
    """
    cleanup_session_folders()

###############################################################################
# Standard Plot Config Preparation
###############################################################################
def prepare_plot_configs(hi_c_matrix, region_chr, region_start, region_end, ds_option,
                         ctcf_bw_for_model, raw_atac_path, del_start=None, del_width=None,
                         norm_atac=None, norm_ctcf=None):
    """
    Renders the final charts for standard (non-chimeric) or single-chrom
    predictions, or simplified for 'chrCHIM'. `raw_atac_path` is used
    for the ATAC track, `ctcf_bw_for_model` for the CTCF track.
    """
    from scipy.ndimage import rotate

    x_start_mb = region_start / 1e6
    x_end_mb   = region_end / 1e6

    # Rotate the Hi-C matrix
    hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)
    mask = np.any(hi_c_matrix > 0, axis=1)
    hi_c_matrix = hi_c_matrix[mask, :]
    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[: num_rows // 2, :]
    n_rows, n_cols = hi_c_matrix.shape

    hi_c_data = []
    for i in range(n_rows):
        for j in range(n_cols):
            hi_c_data.append([j, i, float(hi_c_matrix[i, j])])

    hi_c_chart_config = {
        "chart": {"height": 250},
        "xAxis": {"min": x_start_mb, "max": x_end_mb, "title": ""},
        "yAxis": {"min": 0, "max": 100},
        "colorAxis": {
            "min": 0,
            "max": float(np.nanmax(hi_c_matrix)),
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
        "yAxis": {"title": "CTCF Signal"},
        "series": [{
            "name": "CTCF Signal",
            "data": [[x, y] for x, y in zip(ctcf_positions_mb, ctcf_values)],
            "color": "blue"
        }]
    }

    # ATAC track
    atac_positions, atac_values = get_bigwig_signal(raw_atac_path, region_chr, region_start, region_end)
    atac_positions_mb = [p / 1e6 for p in atac_positions]
    atac_chart_config = {
        "chart": {"height": 100},
        "xAxis": {"min": x_start_mb, "max": x_end_mb, "title": ""},
        "yAxis": {"title": "ATAC Signal"},
        "series": [{
            "name": "ATAC Signal",
            "data": [[x, y] for x, y in zip(atac_positions_mb, atac_values)],
            "color": "green"
        }]
    }

    # If 'deletion' => hide x-axis on Hi-C and add axisBreak to CTCF/ATAC
    if ds_option == "deletion" and del_start is not None and del_width is not None:
        hi_c_chart_config["hideXAxis"] = True
        del_start_mb = del_start / 1e6
        del_end_mb   = (del_start + del_width) / 1e6

        for chart in (ctcf_chart_config, atac_chart_config):
            chart["xAxis"] = {
                "axisBreak": True,
                "leftDomain": [x_start_mb, del_start_mb],
                "rightDomain": [del_end_mb, x_end_mb],
                "title": ""
            }

    return hi_c_chart_config, ctcf_chart_config, atac_chart_config, None, {}

###############################################################################
# NEW HELPER: Unified "Prediction + Render" function
###############################################################################
##### RQ CHANGES #####
# We'll import RQ tasks & queue here
from app.tasks import q
from app.tasks import (
    run_editing_task,
    run_prediction_task,
    build_chimeric_task,
    run_screening_task
)

def run_prediction_and_render(
    region_chr,
    region_start,
    region_end,
    ds_option,
    model_path,
    genome,
    atac_bw_for_model,
    ctcf_bw_for_model,
    raw_atac_path,
    output_folder,
    del_start=None,
    del_width=None,
    norm_atac=None,
    norm_ctcf=None,
    screening_requested=False,
    chimeric_fa_dir=None  # for chimeric mode to pass chrCHIM.fa.gz folder
):
    """
    1. Runs editing.py (if ds_option='deletion') or prediction.py otherwise.
    2. Waits for result.npy, loads => hi_c_matrix.
    3. Builds chart configs.
    4. Optionally sets screening config.
    5. Returns dict with chart configs in JSON form.
    """
    BASE_DIR        = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "..", "C.Origami", "src")
    env             = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH

    script_path       = os.path.join(PYTHON_SRC_PATH, "corigami", "inference")
    hi_c_matrix_path  = os.path.join(output_folder, "result.npy")

    # Decide which seq_dir to use
    if region_chr == "chrCHIM" and chimeric_fa_dir:
        seq_dir_for_prediction = chimeric_fa_dir
    else:
        if genome == 'hg38':
            seq_dir_for_prediction = "./corigami_data/data/hg38/dna_sequence"
        elif genome == 'mm10':
            seq_dir_for_prediction = "./corigami_data/data/mm10/dna_sequence"
        else:
            raise ValueError(f"Unsupported genome: {genome}")

    # 1) Launch editing.py or prediction.py (as an RQ job)
    from app.tasks import q, run_editing_task, run_prediction_task

    if ds_option == "deletion" and del_start is not None and del_width is not None:
        editing_script = os.path.join(script_path, "editing.py")
        job = q.enqueue(
            run_editing_task,
            editing_script,
            region_chr,
            region_start,
            model_path,
            seq_dir_for_prediction,
            atac_bw_for_model,
            del_start,
            del_width,
            output_folder,
            ctcf_bw_for_model=ctcf_bw_for_model,
            env=env
        )
    else:
        prediction_script = os.path.join(script_path, "prediction.py")
        job = q.enqueue(
            run_prediction_task,
            prediction_script,
            region_chr,
            region_start,
            region_end,
            model_path,
            seq_dir_for_prediction,
            atac_bw_for_model,
            output_folder,
            ctcf_bw_for_model=ctcf_bw_for_model,
            env=env
        )

    # **Store the job ID in session so we can cancel if needed**
    session['current_job_id'] = job.id

    # 2) Wait for result (blocking) but also check for "canceled"
    start_wait = time.time()
    while not job.is_finished and not job.is_failed:
        job.refresh()
        # If user canceled => job.get_status() == 'canceled'
        if job.get_status() == 'canceled':
            raise RuntimeError("Job was canceled by the user.")

        if (time.time() - start_wait) >= 60:
            # If we haven't gotten result.npy after 60s, we can still keep waiting
            # or handle a timeout. For now, do nothing except break if you want a hard limit
            pass
        if os.path.exists(hi_c_matrix_path):
            break
        time.sleep(0.5)

    if job.is_failed:
        raise RuntimeError("Error in RQ job. Check worker logs.")

    if not os.path.exists(hi_c_matrix_path):
        raise RuntimeError(f"Hi-C matrix not found at {hi_c_matrix_path}")

    hi_c_matrix = np.load(hi_c_matrix_path)

    # 3) Build chart configs
    hi_c_config, ctcf_config, atac_config, _, _ = prepare_plot_configs(
        hi_c_matrix,
        region_chr,
        0 if region_chr=="chrCHIM" else region_start,
        region_end,
        ds_option,
        ctcf_bw_for_model,
        raw_atac_path,
        del_start,
        del_width,
        norm_atac,
        norm_ctcf
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
    if screening_requested:
        screening_config_json = json.dumps({"placeholder": "screening would happen"})

    return {
        "hi_c_config":       json.dumps(hi_c_config),
        "ctcf_config":       json.dumps(ctcf_config),
        "atac_config":       json.dumps(atac_config),
        "gene_track_config": json.dumps(gene_track_config),
        "screening_config":  screening_config_json
    }

###############################################################################
# Main Page: GET => form, POST => run standard or chimeric logic
###############################################################################
@main.route('/', methods=['GET', 'POST'])
def index():
    """
    Main route:
     - GET => show form
     - POST => decide standard vs chimeric
       * If ds_option="screening" and the user didn't provide a peaks file => call predict_peaks
       * Build gene track config in both standard + chimeric modes
       * Return partial or full page
    """
    if request.method == 'GET':
        user_output_folder = get_user_output_folder()
        return render_template("index.html", screening_mode=False, user_output_folder=user_output_folder)

    # 1) Check if second region is present => chimeric mode
    chimeric_active = (request.form.get('region_chr2', '').strip() != '')
    ds_option = request.form.get('ds_option', 'none')
    screening_requested = (ds_option == 'screening')

    # Common environment, etc.
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "..", "C.Origami", "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH

    output_folder = get_user_output_folder()
    os.makedirs(output_folder, exist_ok=True)

    # 2) Model selection
    region_model = request.form.get('model_select', 'V4')
    if region_model == 'V1':
        model_path = "corigami_data/model_weights/v1_jimin.ckpt"
    elif region_model == 'V2':
        model_path = "corigami_data/model_weights/v2_javier.ckpt"
    elif region_model == 'V3':
        model_path = "corigami_data/model_weights/v3_romane.ckpt"
    elif region_model == 'V4':
        model_path = "corigami_data/model_weights/v4_javier.ckpt"
    else:
        return "Invalid model selection."

    # 3) Genome
    genome = request.form.get('genome_select', 'hg38')
    if genome not in ['hg38', 'mm10']:
        return "Invalid genome selection."

    # 4) Distinguish standard vs chimeric
    if chimeric_active:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # CHIMERIC MODE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print("Running in CHIMERIC mode...")

        from app.utils import (
            normalize_atac, normalize_ctcf, predict_ctcf,
            prepare_chimeric_gene_track_config
        )

        # Region 1
        region_chr1 = request.form.get('region_chr1', '') or request.form.get('region_chr', '')
        region_start1 = int(request.form.get('region_start1', '') or request.form.get('region_start', ''))
        region_end1_str = request.form.get('region_end1', '') or request.form.get('region_end', '')
        if region_end1_str.strip():
            region_end1 = int(region_end1_str)
        else:
            region_end1 = region_start1 + WINDOW_WIDTH

        # Region 2
        region_chr2 = request.form.get('region_chr2', '').strip()
        region_start2 = int(request.form.get('region_start2', '').strip())
        region_end2_str = request.form.get('region_end2', '').strip()
        if region_end2_str:
            region_end2 = int(region_end2_str)
        else:
            region_end2 = region_start2 + WINDOW_WIDTH

        # ATAC / CTCF (raw or user-provided)
        raw_atac_path = request.form.get('atac_bw_path', '').strip()
        raw_ctcf_path = request.form.get('ctcf_bw_path', '').strip()
        if not raw_atac_path:
            return "No ATAC bigwig for chimeric mode."
        if not raw_ctcf_path or raw_ctcf_path == "none":
            return "No CTCF bigwig for chimeric mode."

        # Norm
        norm_atac = request.form.get('norm_atac')
        norm_ctcf = request.form.get('norm_ctcf')
        training_norm_selection = request.form.get('training_norm')

        # Normalize (like standard) => final_atac_bw, final_ctcf_bw
        if raw_ctcf_path and raw_ctcf_path != "none":
            if norm_atac == "none" and norm_ctcf == "none":
                if not training_norm_selection:
                    return "Please select training norm for chimeric."
                final_atac_bw = normalize_atac(raw_atac_path, region_chr1, region_start1, region_end1,
                                               training_norm_selection, output_folder)
                final_ctcf_bw = normalize_ctcf(raw_ctcf_path, region_chr1, region_start1, region_end1,
                                               training_norm_selection, output_folder)
            elif norm_atac != "none" and norm_ctcf == "none":
                final_atac_bw = raw_atac_path
                final_ctcf_bw = normalize_ctcf(raw_ctcf_path, region_chr1, region_start1, region_end1,
                                               norm_atac, output_folder)
            elif norm_atac == "none" and norm_ctcf != "none":
                final_ctcf_bw = raw_ctcf_path
                final_atac_bw = normalize_atac(raw_atac_path, region_chr1, region_start1, region_end1,
                                               norm_ctcf, output_folder)
            else:
                if norm_atac != norm_ctcf:
                    return "ATAC/CTCF must share the same normalization method (chimeric)."
                final_atac_bw = raw_atac_path
                final_ctcf_bw = raw_ctcf_path
        else:
            return "No CTCF bigwig? Not allowed for chimeric."

        # Optional flips
        first_reverse  = (request.form.get('first_reverse') == 'true')
        first_flip     = (request.form.get('first_flip') == 'true')
        second_reverse = (request.form.get('second_reverse') == 'true')
        second_flip    = (request.form.get('second_flip') == 'true')

        # Build chimeric
        build_chimeric_script = os.path.join(PYTHON_SRC_PATH, "corigami", "inference", "build_chimeric.py")

        ##### RQ CHANGES #####
        job = q.enqueue(
            build_chimeric_task,
            build_chimeric_script,
            region_chr1, region_start1, region_end1,
            region_chr2, region_start2, region_end2,
            model_path,
            f"./corigami_data/data/{genome}/dna_sequence",
            final_ctcf_bw,
            final_atac_bw,
            output_folder,
            first_reverse,
            first_flip,
            second_reverse,
            second_flip,
            env=env
        )
        # Block while job finishes
        while not job.is_finished and not job.is_failed:
            time.sleep(0.5)
            job.refresh()
        if job.is_failed:
            return jsonify({"error": "Chimeric build failed. Check RQ worker logs."}), 500

        try:
            chimera_dir, chim_info = job.result
            chim_name   = chim_info.get("CHIM_NAME", "chrCHIM")
            chim_len    = int(chim_info.get("CHIM_LENGTH","0"))
            chim_ctcf_bw= chim_info.get("CHIM_CTCF","")
            chim_atac_bw= chim_info.get("CHIM_ATAC","")
        except Exception as e:
            return jsonify({"error": f"Error building chimeric: {e}"}), 500

        # Next, do standard run_prediction_and_render
        from app.routes import run_prediction_and_render
        configs = run_prediction_and_render(
            region_chr=chim_name,
            region_start=0,
            region_end=chim_len,
            ds_option=ds_option,
            model_path=model_path,
            genome=genome,
            atac_bw_for_model=chim_atac_bw,
            ctcf_bw_for_model=chim_ctcf_bw,
            raw_atac_path=chim_atac_bw,
            output_folder=output_folder,
            del_start=None,
            del_width=None,
            norm_atac=norm_atac,
            norm_ctcf=norm_ctcf,
            screening_requested=screening_requested,
            chimeric_fa_dir=os.path.dirname(chim_info["CHIM_FASTA"]) if "CHIM_FASTA" in chim_info else None
        )

        # Also build a gene track for chimeric => combine gene coords from region1+2
        from app.utils import prepare_chimeric_gene_track_config
        if genome == 'hg38':
            annotation_file = "static/genes.gencode.v38.txt"
        else:
            annotation_file = "static/genes.gencode.M21.mm10.txt"

        gene_track_cfg = prepare_chimeric_gene_track_config(
            annotation_file,
            region_chr1, region_start1, region_end1,
            region_chr2, region_start2, region_end2
        )

        # Return partial or full
        if request.headers.get("X-Requested-With") == "XMLHttpRequest":
            return render_template(
                "plots_partial.html",
                hi_c_config=configs["hi_c_config"],
                ctcf_config=configs["ctcf_config"],
                atac_config=configs["atac_config"],
                screening_config=configs["screening_config"],
                screening_mode=screening_requested,
                screening_params="{}",
                norm_atac=norm_atac,
                norm_ctcf=norm_ctcf,
                gene_track_config=json.dumps(gene_track_cfg)
            )
        else:
            return render_template(
                "index.html",
                hi_c_config=configs["hi_c_config"],
                ctcf_config=configs["ctcf_config"],
                atac_config=configs["atac_config"],
                screening_config=configs["screening_config"],
                screening_mode=screening_requested,
                screening_params="{}",
                gene_track_config=json.dumps(gene_track_cfg),
                user_output_folder=output_folder
            )

    else:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # STANDARD MODE
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print("Running in STANDARD mode...")
        from app.utils import (
            normalize_atac, normalize_ctcf, predict_ctcf,
            predict_peaks, generate_peaks_from_bigwig_macs2,
            prepare_gene_track_config
        )

        # Region
        region_chr = request.form.get('region_chr')
        try:
            region_start = int(request.form.get('region_start'))
        except:
            return "Invalid start position."

        region_end_input = request.form.get('region_end', '').strip()
        if region_end_input:
            try:
                region_end = int(region_end_input)
            except:
                return "Invalid end position."
        else:
            region_end = region_start + WINDOW_WIDTH

        atac_bw_path = request.form.get('atac_bw_path','')
        ctcf_bw_path = (request.form.get('ctcf_bw_path','') or "").strip()
        peaks_file   = request.form.get('peaks_file_path', "").strip()
        norm_atac    = request.form.get('norm_atac')
        norm_ctcf    = request.form.get('norm_ctcf')
        training_norm= request.form.get('training_norm')

        raw_atac_path = atac_bw_path
        raw_ctcf_path = ctcf_bw_path if ctcf_bw_path != "none" else None

        # Normalization logic
        if raw_ctcf_path:
            # both files
            if norm_atac == "none" and norm_ctcf == "none":
                if not training_norm:
                    return "Please select training norm (standard)."
                final_atac_bw = normalize_atac(raw_atac_path, region_chr, region_start, region_end,
                                               training_norm, output_folder)
                final_ctcf_bw = normalize_ctcf(raw_ctcf_path, region_chr, region_start, region_end,
                                               training_norm, output_folder)
            elif norm_atac != "none" and norm_ctcf == "none":
                final_atac_bw = raw_atac_path
                final_ctcf_bw = normalize_ctcf(raw_ctcf_path, region_chr, region_start, region_end,
                                               norm_atac, output_folder)
            elif norm_atac == "none" and norm_ctcf != "none":
                final_ctcf_bw = raw_ctcf_path
                final_atac_bw = normalize_atac(raw_atac_path, region_chr, region_start, region_end,
                                               norm_ctcf, output_folder)
            else:
                if norm_atac != norm_ctcf:
                    return "ATAC/CTCF must share the same norm method."
                final_atac_bw = raw_atac_path
                final_ctcf_bw = raw_ctcf_path
        else:
            # No CTCF => predict
            if norm_atac == "log":
                return "If no CTCF, ATAC must be raw or minmax (log not allowed)."
            elif norm_atac == "minmax":
                final_atac_bw = raw_atac_path
            elif norm_atac == "raw":
                final_atac_bw = normalize_atac(raw_atac_path, region_chr, region_start, region_end,
                                               "minmax", output_folder)
            else:
                return "Unexpected ATAC norm if no CTCF."
            final_ctcf_bw = predict_ctcf(final_atac_bw, region_chr, region_start, region_end, output_folder)
            norm_ctcf = "minmax predicted"

        # If ds_option=deletion => read del start/width
        del_start, del_width = None, None
        if ds_option == "deletion":
            try:
                del_start = int(request.form.get('del_start','1500000'))
                del_width = int(request.form.get('del_width','500000'))
            except:
                return "Invalid deletion parameters."

        # ~~~ If in screening mode + user gave no peaks => auto-generate peaks ~~~
        if ds_option == "screening":
            if not peaks_file or peaks_file == "none":
                # Use the newly prepped final_atac_bw
                # The user wants MACS2-based peaks => so we do
                # predict_peaks( ) or generate_peaks_from_bigwig_macs2
                peaks_file = predict_peaks(final_atac_bw, region_chr, region_start, region_end, output_folder)
                print(f"[Standard] Auto-generated peaks => {peaks_file}")

        from app.routes import run_prediction_and_render
        configs = run_prediction_and_render(
            region_chr=region_chr,
            region_start=region_start,
            region_end=region_end,
            ds_option=ds_option,
            model_path=model_path,
            genome=genome,
            atac_bw_for_model=final_atac_bw,
            ctcf_bw_for_model=final_ctcf_bw,
            raw_atac_path=raw_atac_path,
            output_folder=output_folder,
            del_start=del_start,
            del_width=del_width,
            norm_atac=norm_atac,
            norm_ctcf=norm_ctcf,
            screening_requested=screening_requested,
            chimeric_fa_dir=None
        )

        # Build gene track for standard mode
        gene_track_cfg = prepare_gene_track_config(
            genome,
            region_chr,
            region_start,
            region_end,
            ds_option,
            del_start,
            del_width
        )

        if request.headers.get("X-Requested-With") == "XMLHttpRequest":
            return render_template(
                "plots_partial.html",
                hi_c_config=configs["hi_c_config"],
                ctcf_config=configs["ctcf_config"],
                atac_config=configs["atac_config"],
                screening_config=configs["screening_config"],
                screening_mode=screening_requested,
                screening_params="{}",
                norm_atac=norm_atac,
                norm_ctcf=norm_ctcf,
                gene_track_config=json.dumps(gene_track_cfg)
            )
        else:
            return render_template(
                "index.html",
                hi_c_config=configs["hi_c_config"],
                ctcf_config=configs["ctcf_config"],
                atac_config=configs["atac_config"],
                screening_config=configs["screening_config"],
                screening_mode=screening_requested,
                screening_params="{}",
                gene_track_config=json.dumps(gene_track_cfg),
                user_output_folder=output_folder
            )

###############################################################################
# Screening Route
###############################################################################
@main.route('/run_screening', methods=['GET'])
def run_screening_endpoint():
    """
    AJAX endpoint to run corigami's screening routine on a region,
    optionally auto-generating peaks from BigWig if none provided.
    If chimeric_active is true, we build a synthetic chromosome
    then override region parameters (chrCHIM, etc.) and proceed.
    Returns JSON including "screening_config" for the front-end plot.
    """

    from app.tasks import q, run_screening_task
    from app.utils import generate_peaks_from_bigwig_macs2, get_user_output_folder
    
    print("RUN_SCREENING_ENDPOINT CALLED", flush=True)
    current_app.logger.info("Entered run_screening_endpoint")
    params = request.args.to_dict()
    current_app.logger.info("Received screening parameters: %s", params)

    try:
        # ------------------------------------------------------
        # 1) Parse main screening parameters
        # ------------------------------------------------------
        from app.utils import generate_peaks_from_bigwig_macs2, get_user_output_folder

        chimeric_active = request.args.get('chimeric_active', 'false').lower() == 'true'

        region_chr = request.args.get('region_chr', 'chr2')
        screen_start = int(request.args.get('region_start', '500000'))
        region_end_str = request.args.get('region_end', '')
        if not region_end_str.strip():
            screen_end = screen_start + WINDOW_WIDTH
        else:
            screen_end = int(region_end_str)

        perturb_width = int(request.args.get('perturb_width', '1000'))
        step_size     = int(request.args.get('step_size', '1000'))

        BASE_DIR = os.path.dirname(current_app.root_path)
        output_dir = request.args.get('output_dir', get_user_output_folder())
        if not output_dir or not os.path.exists(output_dir):
            current_app.logger.error("Invalid or missing output directory: %s", output_dir)
            return jsonify({"error": "Invalid or missing output directory"}), 400

        # 2) Model
        region_model = request.args.get('model_select', 'V1')
        if region_model == 'V1':
            model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "v1_jimin.ckpt")
        elif region_model == 'V2':
            model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "v2_javier.ckpt")
        elif region_model == 'V3':
            model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "v3_romane.ckpt")
        elif region_model == 'V4':
            model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "v4_javier.ckpt")
        else:
            return jsonify({"error": "Invalid model selection"}), 400

        # 3) ATAC/CTCF
        default_atac_path = os.path.join(BASE_DIR, "corigami_data", "data", "hg38", "imr90", "genomic_features", "atac.bw")
        default_ctcf_path= os.path.join(BASE_DIR, "corigami_data", "data", "hg38", "imr90", "genomic_features", "ctcf_log2fc.bw")
        atac_bw_path = request.args.get('atac_bw_path', default_atac_path)
        ctcf_bw_path = request.args.get('ctcf_bw_path', default_ctcf_path)

        # 4) peaks_file
        peaks_file = request.args.get('peaks_file', "").strip()

        # ------------------------------------------------------
        # 5) If chimeric => Build Synthetic Chromosome
        # ------------------------------------------------------
        if chimeric_active:
            current_app.logger.info("SCREENING in chimeric mode => building chrCHIM first...")

            region_chr2     = request.args.get('region_chr2', 'chr9')
            region_start2   = int(request.args.get('region_start2', '1500000'))
            region_end2_str = request.args.get('region_end2', '')
            if region_end2_str.strip():
                region_end2 = int(region_end2_str)
            else:
                region_end2 = region_start2 + WINDOW_WIDTH

            genome_select = request.args.get('genome_select', 'hg38')
            seq_dir = os.path.join(BASE_DIR, "corigami_data", "data", genome_select, "dna_sequence")

            first_reverse  = (request.args.get('first_reverse') == 'true')
            first_flip     = (request.args.get('first_flip') == 'true')
            second_reverse = (request.args.get('second_reverse') == 'true')
            second_flip    = (request.args.get('second_flip') == 'true')

            PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
            env             = os.environ.copy()
            env["PYTHONPATH"] = PYTHON_SRC_PATH
            build_chimeric_script = os.path.join(PYTHON_SRC_PATH, "corigami", "inference", "build_chimeric.py")

            ##### RQ CHANGES #####
            from app.tasks import q, build_chimeric_task
            job_chim = q.enqueue(
                build_chimeric_task,
                build_chimeric_script,
                region_chr, screen_start, screen_end,
                region_chr2, region_start2, region_end2,
                model_path,
                seq_dir,
                ctcf_bw_path,
                atac_bw_path,
                output_dir,
                first_reverse,
                first_flip,
                second_reverse,
                second_flip,
                env=env
            )
            while not job_chim.is_finished and not job_chim.is_failed:
                time.sleep(0.5)
                job_chim.refresh()
            if job_chim.is_failed:
                raise ValueError("Could not build chimeric in screening. RQ job failed.")

            chimera_dir, chim_info = job_chim.result
            region_chr   = chim_info["CHIM_NAME"]  # "chrCHIM"
            screen_start = 0
            screen_end   = int(chim_info["CHIM_LENGTH"])
            atac_bw_path = chim_info["CHIM_ATAC"]
            ctcf_bw_path = chim_info["CHIM_CTCF"]

            current_app.logger.info("Now screening on %s:%d-%d", region_chr, screen_start, screen_end)

    except Exception as e:
        current_app.logger.exception("Invalid screening parameters or chimeric build error.")
        return jsonify({"error": "Invalid screening parameters or chimeric build error", "details": str(e)}), 400

    # ------------------------------------------------------
    # 6) Actual screening logic
    # ------------------------------------------------------
    try:
        # Possibly remove old screening_results
        results_file = os.path.join(output_dir, "screening", "screening_results.json")
        if os.path.exists(results_file):
            os.remove(results_file)

        # Auto-peaks if user did not specify peaks_file
        auto_peaks_file = None
        temp_bedgraph   = None
        if not peaks_file or peaks_file == "none":
            try:
                auto_peaks_file, temp_bedgraph = generate_peaks_from_bigwig_macs2(
                    bw_path=atac_bw_path,
                    chrom=region_chr,
                    start=screen_start,
                    end=screen_end,
                    outdir=output_dir
                )
                peaks_file = auto_peaks_file
            except subprocess.CalledProcessError as e:
                current_app.logger.exception("Failed to run bigWigToBedGraph or MACS2.")
                return jsonify({"error": "Failed to run bigWigToBedGraph or MACS2", "details": str(e)}), 500
            except Exception as e:
                current_app.logger.exception("Failed to generate peaks from BigWig")
                return jsonify({"error": "Failed to generate peaks from BigWig", "details": str(e)}), 500
        else:
            auto_peaks_file, temp_bedgraph = None, None

        # call screening.py
        PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
        env = os.environ.copy()
        env["PYTHONPATH"] = PYTHON_SRC_PATH
        screening_script = os.path.join(PYTHON_SRC_PATH, "corigami", "inference", "screening.py")
        seq_dir = os.path.join(BASE_DIR, "corigami_data", "data", "hg38", "dna_sequence")

        ##### RQ CHANGES #####
        job_screen = q.enqueue(
            run_screening_task,
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
            step_size,
            env=env
        )
        # block until done
        while not job_screen.is_finished and not job_screen.is_failed:
            time.sleep(1)
            job_screen.refresh()
        if job_screen.is_failed:
            current_app.logger.error("Screening command failed with RQ job id: %s", job_screen.id)
            return jsonify({"error": "Screening script failed with RQ worker"}), 500

        if not os.path.exists(results_file):
            current_app.logger.error("screening_results.json not found at %s", results_file)
            return jsonify({"error": "screening_results.json was not generated."}), 500

        # parse the result
        with open(results_file, "r") as f:
            results = json.load(f)

        window_midpoints_mb = results.get("window_midpoints_mb", [])
        impact_scores       = results.get("impact_scores", [])
        screening_series_data = sorted(
            [[x, y] for x, y in zip(window_midpoints_mb, impact_scores)],
            key=lambda d: d[0]
        )
        region_start_mb = screen_start / 1e6
        region_end_mb   = screen_end   / 1e6

        screening_chart_config = {
            "chart": {"type": "line", "height": 100},
            "xAxis": {
                "min": region_start_mb,
                "max": region_end_mb,
                "tickColor": "#000",
                "tickFont":  "10px sans-serif"
            },
            "yAxis": {
                "title": {"text": ""},
                "ticks": 3,
                "tickColor": "#000",
                "tickFont": "10px sans-serif"
            },
            "series": [{
                "name": "Screening Score",
                "data": screening_series_data,
                "color": "dodgerblue"
            }]
        }
        screening_config_json = json.dumps(screening_chart_config)
        results["screening_config"] = screening_config_json

        if auto_peaks_file and os.path.exists(auto_peaks_file):
            os.remove(auto_peaks_file)
        if temp_bedgraph and os.path.exists(temp_bedgraph):
            os.remove(temp_bedgraph)

        current_app.logger.info("Returning screening results with screening_config")
        return jsonify(results)

    except Exception as e:
        current_app.logger.exception("Exception in screening logic")
        return jsonify({"error": "Exception in screening logic", "details": str(e)}), 500

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

@main.route('/cancel_run', methods=['POST'])
def cancel_run():
    """
    Cancels the current RQ job if one is stored in session['current_job_id'].
    Returns JSON indicating success or error.
    """
    from app.tasks import q
    job_id = session.get('current_job_id')
    if not job_id:
        return jsonify({"error": "No job in progress"}), 400

    job = q.fetch_job(job_id)
    if not job:
        return jsonify({"error": f"No job found for ID: {job_id}"}), 404

    job.cancel()  # Mark the job as canceled
    return jsonify({"message": "Current run has been canceled."})
