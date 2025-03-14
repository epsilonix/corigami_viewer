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
    WINDOW_WIDTH,
    clear_folder
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
    run_screening_task
)

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
    raw_atac_path,
    output_folder,
    del_start=None,
    del_width=None,
    norm_atac=None,
    norm_ctcf=None,
    screening_requested=False # for chimeric mode to pass chrCHIM.fa.gz folder
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

    # Decide which seq_dir to us

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
            seq_dir,
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
            seq_dir,
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
# routes.py (excerpt) -- remove references to build_chimeric.py or build_chimeric_task

@main.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'GET':
        user_output_folder = get_user_output_folder()
        return render_template("index.html", screening_mode=False, user_output_folder=user_output_folder)

    chimeric_active = (request.form.get('region_chr2', '').strip() != '')
    ds_option = request.form.get('ds_option', 'none')
    screening_requested = (ds_option == 'screening')

    output_folder = get_user_output_folder()
    clear_folder(output_folder)
    os.makedirs(output_folder, exist_ok=True)

    # Model selection, genome, etc. (unchanged)
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

    # ... pick model_path ...
    genome = request.form.get('genome_select')
    session['final_genome'] = genome
    if chimeric_active:
        print("Running in CHIMERIC mode (inline) ...")
        from app.utils import (
            assemble_chimeric_arrays,
            write_chimeric_bigwig,
            normalize_atac,
            normalize_ctcf,
            prepare_chimeric_gene_track_config
        )

        # 1) Gather region1 + region2
        region_chr1 = request.form.get('region_chr1')
        region_start1 = int(request.form.get('region_start1'))
        region_end1 = region_start1 + WINDOW_WIDTH
        maybe_end1 = request.form.get('region_end1','').strip()
        if maybe_end1:
            region_end1 = int(maybe_end1)

        region_chr2 = request.form.get('region_chr2')
        region_start2 = int(request.form.get('region_start2'))
        region_end2 = region_start2 + WINDOW_WIDTH
        maybe_end2 = request.form.get('region_end2','').strip()
        if maybe_end2:
            region_end2 = int(maybe_end2)

        # 2) Raw ATAC / CTCF (no normalization yet)
        raw_atac_path = request.form.get('atac_bw_path','').strip()
        raw_ctcf_path = request.form.get('ctcf_bw_path','').strip()
        if not raw_atac_path:
            return "Chimeric mode requires an ATAC bigWig"
        if not raw_ctcf_path or raw_ctcf_path == "none":
            return "Chimeric mode requires a CTCF bigWig"

        # 3) Flips
        first_reverse  = (request.form.get('first_reverse') == 'true')
        first_flip     = (request.form.get('first_flip') == 'true')
        second_reverse = (request.form.get('second_reverse') == 'true')
        second_flip    = (request.form.get('second_flip') == 'true')

        # 4) Assemble arrays => "chrCHIM"
        chim_atac, chim_ctcf = assemble_chimeric_arrays(
            raw_atac_path,
            raw_ctcf_path,
            region_chr1, region_start1, region_end1,
            region_chr2, region_start2, region_end2,
            first_reverse, first_flip,
            second_reverse, second_flip
        )
        chim_len = len(chim_atac)  # same as len(chim_ctcf)

        # 5) Write “chrCHIM_atac.bw” + “chrCHIM_ctcf.bw”
        chim_atac_bw, chim_ctcf_bw = write_chimeric_bigwig(
            chim_atac, chim_ctcf, output_folder, chim_len, chim_name="chrCHIM"
        )
        print(f"Chimeric chromosome length => {chim_len} bp")

        # 5a) create synthetic FASTA
        # We'll assume you keep your reference FASTAs in e.g. "corigami_data/data/hg38/dna_sequence/chr1.fa.gz", etc.
        # We'll attempt to locate them based on region_chr1, region_chr2:
        # This is just an example. You can adapt how you find them:
        fa1_path = os.path.join(f"./corigami_data/data/{genome}/dna_sequence", f"{region_chr1}.fa.gz")
        fa2_path = os.path.join(f"./corigami_data/data/{genome}/dna_sequence", f"{region_chr2}.fa.gz")

        from app.utils import assemble_chimeric_fasta
        chimera_files_dir = os.path.join(output_folder, "chimera_files")
        os.makedirs(chimera_files_dir, exist_ok=True)

        chim_fa_gz = assemble_chimeric_fasta(
            fa1_path, region_start1, region_end1, 
            do_reverse1=first_reverse, do_flip1=first_flip,
            fa2_path=fa2_path, chr2_start=region_start2, chr2_end=region_end2,
            do_reverse2=second_reverse, do_flip2=second_flip,
            output_folder=chimera_files_dir,
            chim_name="chrCHIM"
        )

        seq_dir = chimera_files_dir
        session['seq_dir'] = seq_dir
        
        peaks_file = request.form.get('peaks_file_path', '').strip()
        if peaks_file:
            session['final_peaks_file'] = peaks_file


        # 6) Normalization if requested (chimeric mode) -- updated to mimic standard mode normalization logic
        norm_atac_method = request.form.get('norm_atac')   # e.g. 'none', 'minmax', 'log'
        norm_ctcf_method = request.form.get('norm_ctcf')     # e.g. 'none', 'minmax', etc.
        training_norm = request.form.get('training_norm')    # e.g. 'minmax' or another training norm method

        # Both chim_atac_bw and chim_ctcf_bw should be available in chimeric mode.
        if chim_ctcf_bw:
            if norm_atac_method == "none" and norm_ctcf_method == "none":
                # Require a training norm method when neither normalization is explicitly chosen.
                if not training_norm:
                    return "Please select training norm (standard)."
                final_atac_bw = normalize_atac(chim_atac_bw, "chrCHIM", 0, chim_len, training_norm, output_folder)
                final_ctcf_bw = normalize_ctcf(chim_ctcf_bw, "chrCHIM", 0, chim_len, training_norm, output_folder)
            elif norm_atac_method != "none" and norm_ctcf_method == "none":
                # Use the ATAC norm method for CTCF if only ATAC's norm method is provided.
                final_atac_bw = chim_atac_bw
                final_ctcf_bw = normalize_ctcf(chim_ctcf_bw, "chrCHIM", 0, chim_len, norm_atac_method, output_folder)
            elif norm_atac_method == "none" and norm_ctcf_method != "none":
                # Use the CTCF norm method for ATAC if only CTCF's norm method is provided.
                final_ctcf_bw = chim_ctcf_bw
                final_atac_bw = normalize_atac(chim_atac_bw, "chrCHIM", 0, chim_len, norm_ctcf_method, output_folder)
            else:
                # If both are provided, they must match.
                if norm_atac_method != norm_ctcf_method:
                    return "ATAC/CTCF must share the same norm method."
                final_atac_bw = chim_atac_bw
                final_ctcf_bw = chim_ctcf_bw
        else:
            # Fallback in the unlikely event that chim_ctcf_bw is missing.
            final_atac_bw = chim_atac_bw
            final_ctcf_bw = None


        # 7) run_prediction_and_render
        from app.routes import run_prediction_and_render
        configs = run_prediction_and_render(
            region_chr="chrCHIM",
            region_start=0,
            region_end=chim_len,
            ds_option=ds_option,
            model_path=model_path,
            genome=genome,
            atac_bw_for_model=final_atac_bw,
            ctcf_bw_for_model=final_ctcf_bw,
            raw_atac_path=final_atac_bw,
            output_folder=output_folder,
            norm_atac=norm_atac_method,
            norm_ctcf=norm_ctcf_method,
            screening_requested=screening_requested,
            seq_dir=seq_dir  # <-- points to where chrCHIM.fa.gz resides
        )

        session['final_atac_bw'] = final_atac_bw
        session['final_ctcf_bw'] = final_ctcf_bw
        session['final_region_chr'] = "chrCHIM"
        session['final_region_start'] = 0
        session['final_region_end'] = chim_len

        # 8) Build a gene track for “chrCHIM”
        if genome == 'hg38':
            annotation_file = "static/genes.gencode.v38.txt"
        else:
            annotation_file = "static/genes.gencode.M21.mm10.txt"

        offset = (region_end1 - region_start1)
        gene_track_cfg = prepare_chimeric_gene_track_config(
            annotation_file,
            region_chr1, region_start1, region_end1,
            region_chr2, region_start2, region_end2
        )

        # 9) Return whichever template
        if request.headers.get("X-Requested-With") == "XMLHttpRequest":
            return render_template(
                "plots_partial.html",
                hi_c_config=configs["hi_c_config"],
                ctcf_config=configs["ctcf_config"],
                atac_config=configs["atac_config"],
                screening_config=configs["screening_config"],
                screening_mode=screening_requested,
                screening_params="{}",
                norm_atac=norm_atac_method,
                norm_ctcf=norm_ctcf_method,
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
        region_chr1 = request.form.get('region_chr1')
        try:
            region_start1 = int(request.form.get('region_start1'))
        except:
            return "Invalid start position."

        region_end1 = request.form.get('region_end1', '').strip()
        if region_end1:
            try:
                region_end1 = int(region_end1)
            except:
                return "Invalid end position."
        else:
            region_end1 = region_start1 + WINDOW_WIDTH

        if genome == 'hg38':
            seq_dir = "./corigami_data/data/hg38/dna_sequence"
        elif genome == 'mm10':
            seq_dir = "./corigami_data/data/mm10/dna_sequence"
        session['seq_dir'] = seq_dir


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
                final_atac_bw = normalize_atac(raw_atac_path, region_chr1, region_start1, region_end1,
                                               training_norm, output_folder)
                final_ctcf_bw = normalize_ctcf(raw_ctcf_path, region_chr1, region_start1, region_end1,
                                               training_norm, output_folder)
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
                final_atac_bw = normalize_atac(raw_atac_path, region_chr1, region_start1, region_end1,
                                               "minmax", output_folder)
            else:
                return "Unexpected ATAC norm if no CTCF."
            final_ctcf_bw = predict_ctcf(final_atac_bw, region_chr1, region_start1, region_end1, output_folder)
            norm_ctcf = "minmax predicted"

        # If ds_option=deletion => read del start/width
        del_start, del_width = None, None
        if ds_option == "deletion":
            try:
                del_start = int(request.form.get('del_start','1500000'))
                del_width = int(request.form.get('del_width','500000'))
            except:
                return "Invalid deletion parameters."

        from app.routes import run_prediction_and_render
        configs = run_prediction_and_render(
            region_chr=region_chr1,
            region_start=region_start1,
            region_end=region_end1,
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
            seq_dir=seq_dir
        )

        session['final_atac_bw'] = final_atac_bw
        session['final_ctcf_bw'] = final_ctcf_bw
        session['final_region_chr'] = region_chr1
        session['final_region_start'] = region_start1
        session['final_region_end'] = region_end1

        # Build gene track for standard mode
        gene_track_cfg = prepare_gene_track_config(
            genome,
            region_chr1,
            region_start1,
            region_end1,
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
    

    print(f'screening starting with seq_dir: {seq_dir}')


    # If any are missing => user did not run the main route first
    if not final_atac_bw or not final_ctcf_bw or region_chr is None:
        return jsonify({"error": "No existing data from main run. Please run prediction first."}), 400

    # 2) Extract additional parameters for how to run the screening
    # e.g. perturb_width, step_size, model selection
    perturb_width = int(request.args.get('perturb_width', '1000'))
    step_size     = int(request.args.get('step_size', '1000'))
    model_select  = request.args.get('model_select')
    output_dir    = get_user_output_folder()
    # ... map model_select to model_path ...

    if model_select == 'V1':
        model_path = "corigami_data/model_weights/v1_jimin.ckpt"
    elif model_select == 'V2':
        model_path = "corigami_data/model_weights/v2_javier.ckpt"
    elif model_select == 'V3':
        model_path = "corigami_data/model_weights/v3_romane.ckpt"
    elif model_select == 'V4':
        model_path = "corigami_data/model_weights/v4_javier.ckpt"

    peaks_file = session.get('final_peaks_file', '').strip()
    print(f'run screening endpoint: peaks file: {peaks_file}')
    if not peaks_file:
        # If user didn’t select one in the form, or you want to handle the fallback
        # If you want to force the user to pick a peaks file, you can raise an error instead.
        peaks_file, temp_bedgraph = generate_peaks_from_bigwig_macs2(
            bw_path=final_atac_bw,
            chrom=region_chr,
            start=region_start,
            end=region_end,
            outdir=output_dir
        )
        print(f'auto-generated peaks file: {peaks_file}')
    print(f'run screening endpoint: peaks file: {peaks_file}')
    # 4) Now do the screening with final_atac_bw, final_ctcf_bw, peaks_file
    # (No re-assembly, no normalization, no chimeric logic needed)
    screening_script = os.path.join("C.Origami", "src", "corigami", "inference", "screening.py")
    env = os.environ.copy()
    env["PYTHONPATH"] = "C.Origami/src"

    # Determine the DNA sequence directory for screening:
    seq_dir = session.get('seq_dir')

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
        step_size,
        env=env
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
      "chart": {"type": "line", "height": 100},
      "xAxis": {"min": region_start/1e6, "max": region_end/1e6},
      "yAxis": {"title": {"text": ""}},
      "series": [{
          "name": "Screening Score",
          "data": series_data,
          "color": "dodgerblue"
      }]
    }
    results["screening_config"] = json.dumps(screening_chart)
    
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
