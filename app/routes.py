test_mode = False

import os
import subprocess
import json
import time
import uuid
import numpy as np
from PIL import Image
from scipy.ndimage import rotate
from flask import Blueprint, render_template, request, url_for, jsonify, session, current_app


from app.utils import (
    get_upload_folder,
    get_user_upload_folder,
    get_user_output_folder,
    cleanup_session_folders,
    generate_peaks_from_bigwig_macs2,
    get_bigwig_signal,
    save_uploaded_file,
    normalize_file,
    WINDOW_WIDTH
)

main = Blueprint('main', __name__)

@main.before_request
def before_request():
    cleanup_session_folders()

###############################################################################
# AJAX File Upload Endpoint
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

    if file_type == "atac":
        saved_path = save_uploaded_file(file_storage, "atac.bw", upload_folder)
    elif file_type == "ctcf":
        saved_path = save_uploaded_file(file_storage, "ctcf.bw", upload_folder)
    elif file_type == "peaks":
        saved_path = save_uploaded_file(file_storage, "peaks.narrowPeak", upload_folder)

    print(f"Uploaded {file_type} file saved to: {saved_path}")
    display_name = os.path.basename(saved_path).split("_", 1)[-1]
    return jsonify({"saved_path": saved_path, "display_name": display_name})

###############################################################################
# Main Page: Handle Form Submission or Render
###############################################################################

def prepare_plot_configs(hi_c_matrix, region_chr, region_start, region_end, ds_option,
                         ctcf_bw_for_model, raw_atac_path, del_start=None, del_width=None, norm_atac=None, norm_ctcf=None):
    import numpy as np
    from scipy.ndimage import rotate

    x_start_mb = region_start / 1e6
    x_end_mb   = region_end / 1e6

    hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)

    mask = np.any(hi_c_matrix > 0, axis=1)
    hi_c_matrix = hi_c_matrix[mask, :]

    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[:num_rows // 2, :]
    n_rows, n_cols = hi_c_matrix.shape

    hi_c_data = []
    for i in range(n_rows):
        for j in range(n_cols):
            hi_c_data.append([j, i, float(hi_c_matrix[i, j])])
    
    hi_c_chart_config = {
        "chart": {"height": 250},
        "xAxis": {
            "min": x_start_mb,
            "max": x_end_mb,
        },
        "yAxis": {"min": 0, "max": 100},
        "colorAxis": {
            "min": 0,
            "max": float(np.nanmax(hi_c_matrix)),
            "minColor": "#ffffff",
            "maxColor": "#7f0000"
        },
        "series": [{
            "name": "Hi-C Signal",
            "data": hi_c_data
        }],
        "n_cols": n_cols
    }

    ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_for_model, region_chr, region_start, region_end)
    atac_positions, atac_values = get_bigwig_signal(raw_atac_path, region_chr, region_start, region_end)

    ctcf_positions_mb = [p / 1e6 for p in ctcf_positions]
    atac_positions_mb = [p / 1e6 for p in atac_positions]

    ctcf_data = [[x, y] for x, y in zip(ctcf_positions_mb, ctcf_values)]
    atac_data = [[x, y] for x, y in zip(atac_positions_mb, atac_values)]

    if norm_ctcf is None or norm_ctcf == "none":
        ctcf_title = "raw"
    elif "predicted with maxATAC" in norm_ctcf:
        ctcf_title = norm_ctcf
    else:
        ctcf_title = norm_ctcf + " normalized"

    ctcf_chart_config = {
        "chart": {"height": 100},
        "title": {
            "text": f"CTCF signal ({ctcf_title})",
            "style": {"fontSize": "9px"}  # Adjust the font size here
        },
        "xAxis": {"min": x_start_mb, "max": x_end_mb},
        "yAxis": {"title": "CTCF Signal"},
        "series": [{"name": "CTCF Signal", "data": ctcf_data, "color": "blue"}]
    }
    atac_chart_config = {
        "chart": {"height": 100},
        "title": {
            "text": f"ATAC signal ({'raw' if (norm_atac is None or norm_atac == 'none') else norm_atac + ' normalized'})",
            "style": {"fontSize": "9px"}  # Adjust the font size here
        },
        "xAxis": {"min": x_start_mb, "max": x_end_mb},
        "yAxis": {"title": "ATAC Signal"},
        "series": [{"name": "ATAC Signal", "data": atac_data, "color": "green"}]
    }

    screening_chart_config = None
    screening_params = {}
    if ds_option == "screening":
        try:
            perturb_width = int(request.form.get('perturb_width', '1000'))
            step_size = int(request.form.get('step_size', '1000'))
        except ValueError:
            perturb_width, step_size = 1000, 1000
        screening_params = {
            "region_chr": region_chr,
            "region_start": region_start,
            "perturb_width": perturb_width,
            "step_size": step_size,
            "atac_bw_path": ctcf_bw_for_model,
            "ctcf_bw_path": ctcf_bw_for_model,
            "peaks_file": request.form.get('peaks_file_path', "").strip(),
            "output_dir": get_user_output_folder()
        }
        screening_chart_config = {
            "chart": {"height": 250},
            "yAxis": {"title": "Impact score"},
            "series": [{"name": "Impact score", "data": [], "color": "dodgerblue"}]
        }

    return hi_c_chart_config, ctcf_chart_config, atac_chart_config, screening_chart_config, screening_params

@main.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'GET':
        user_output_folder = get_user_output_folder()
        return render_template("index.html", screening_mode=False, user_output_folder=user_output_folder)

    print("=== Starting new submission ===")
    atac_upload_folder = get_upload_folder("atac")
    ctcf_upload_folder = get_upload_folder("ctcf")
    peaks_upload_folder = get_upload_folder("peaks")

    output_folder = get_user_output_folder()
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "..", "C.Origami", "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH

    region_model = request.form.get('model_select')
    region_chr = request.form.get('region_chr')
    try:
        region_start = int(request.form.get('region_start'))
    except ValueError:
        return "Invalid start position. Please enter an integer value."
    DEFAULT_WINDOW = 2097152
    MAX_WINDOW = 20971520
    if request.form.get("allow_wider_windows"):
        try:
            region_end = int(request.form.get("region_end"))
        except ValueError:
            return "Invalid end position. Please enter an integer value."
        if region_end > region_start + MAX_WINDOW:
            return f"End position must be no greater than {region_start + MAX_WINDOW}."
    else:
        region_end = region_start + DEFAULT_WINDOW

    print(f"Window: {region_chr}:{region_start}-{region_end}")


    atac_bw_path = request.form.get('atac_bw_path')
    print("ATAC file path (from dropdown):", os.path.abspath(atac_bw_path))
    ctcf_bw_path = (request.form.get('ctcf_bw_path') or "").strip()
    peaks_file = request.form.get('peaks_file_path', "").strip()

    norm_atac = request.form.get('norm_atac')
    norm_ctcf = request.form.get('norm_ctcf')
    training_norm_selection = request.form.get('training_norm')
    print(f"ATAC normalization: {norm_atac}")
    print(f"CTCF normalization: {norm_ctcf}")
    print(f"Training normalization selection: {training_norm_selection}")

    raw_atac_path = atac_bw_path
    raw_ctcf_path = ctcf_bw_path if ctcf_bw_path and ctcf_bw_path != "none" else None

    if raw_ctcf_path:
        if norm_atac != "none" or norm_ctcf != "none":
            if norm_atac != "none" and norm_ctcf != "none" and norm_atac != norm_ctcf:
                return "Uploaded files must use the same normalization method."
            training_norm = norm_atac if norm_atac != "none" else norm_ctcf
            if norm_atac == "none":
                normalized_atac = normalize_file(raw_atac_path, region_chr, region_start, region_end, training_norm, output_folder, "normalized_atac")
            else:
                normalized_atac = raw_atac_path
            if norm_ctcf == "none":
                normalized_ctcf = normalize_file(raw_ctcf_path, region_chr, region_start, region_end, training_norm, output_folder, "normalized_ctcf")
            else:
                normalized_ctcf = raw_ctcf_path
        else:
            if not training_norm_selection:
                return "Please select the normalization method used during training."
            training_norm = training_norm_selection
            normalized_atac = normalize_file(raw_atac_path, region_chr, region_start, region_end, training_norm, output_folder, "normalized_atac")
            normalized_ctcf = normalize_file(raw_ctcf_path, region_chr, region_start, region_end, training_norm, output_folder, "normalized_ctcf")
    else:
        training_norm = "minmax"
        if norm_atac != "none":
            return "Non-normalized ATAC signal is required to predict CTCF when no CTCF file is provided."
        normalized_atac = raw_atac_path
        normalized_ctcf = None

    atac_bw_for_model = normalized_atac

    if not raw_ctcf_path:
        print("No CTCF file provided; generating CTCF file using raw ATAC signal...")
        roi_file = os.path.join(output_folder, "temp_roi.bed")
        with open(roi_file, "w") as f:
            f.write(f"{region_chr}\t{region_start}\t{region_end}\n")
        ctcf_generated = os.path.join(output_folder, "predicted_ctcf.bw")
        generate_cmd = [
            "maxatac", "predict", "--tf", "CTCF",
            "--signal", atac_bw_for_model,
            "--bed", roi_file,
            "--out", output_folder,
            "--name", "predicted_ctcf"
        ]
        print("Running maxATAC for CTCF generation:", " ".join(generate_cmd))
        # Create a temporary environment for this subprocess call.
        temp_env = env.copy()
        temp_env["PATH"] = "/Users/everett/anaconda3/bin:" + temp_env["PATH"]
        try:
            subprocess.run(generate_cmd, check=True, env=temp_env, capture_output=True, text=True)
            normalized_predicted_ctcf = normalize_file(ctcf_generated, region_chr, region_start, region_end, "minmax", output_folder, "normalized_predicted_ctcf")
            ctcf_bw_for_model = normalized_predicted_ctcf
            # Set norm_ctcf to a special string indicating it was predicted.
            norm_ctcf = "minmax normalized, predicted with maxATAC"
            print("Predicted CTCF file generated and normalized:", ctcf_bw_for_model)
        except subprocess.CalledProcessError as e:
            return f"Error generating CTCF file: {e.stderr}"
        finally:
            if os.path.exists(roi_file):
                os.remove(roi_file)
    else:
        ctcf_bw_for_model = normalized_ctcf



    model_path = "corigami_data/model_weights/atac_ctcf_model.ckpt"
    ds_option = request.form.get('ds_option', 'none')
    print("Selected ds_option:", ds_option)

    if ds_option == "deletion":
        print("Running deletion process...")
        try:
            del_start = int(request.form.get('del_start', '1500000'))
            del_width = int(request.form.get('del_width', '500000'))
        except ValueError:
            return "Invalid deletion parameters."
        script_path = os.path.join(PYTHON_SRC_PATH, "corigami/inference/editing.py")
        cmd = [
            "python", script_path,
            "--chr", region_chr,
            "--start", str(region_start),
            "--model", model_path,
            "--seq", "./corigami_data/data/hg38/dna_sequence",
            "--atac", atac_bw_for_model,
            "--del-start", str(del_start),
            "--del-width", str(del_width),
            "--out", output_folder,
            "--ctcf", ctcf_bw_for_model
        ]
        print("Command:", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
            print("Editing script completed.")
        except subprocess.CalledProcessError as e:
            return f"Error running editing script: {e.stderr}"
        hi_c_matrix_path = os.path.join(output_folder, "deletion", "npy", f"{region_chr}_{region_start}_del_{del_start}_{del_width}_padding_zero.npy")
    else:
        print("Running standard prediction...")
        script_path = os.path.join(BASE_DIR, "..", "C.Origami", "src", "corigami/inference/prediction.py")
        cmd = [
            "python", script_path,
            "--chr", region_chr,
            "--start", str(region_start),
            "--model", model_path,
            "--seq", "./corigami_data/data/hg38/dna_sequence",
            "--atac", atac_bw_for_model,
            "--out", output_folder,
            "--ctcf", ctcf_bw_for_model
        ]
        if region_end > region_start + DEFAULT_WINDOW:
            cmd.extend(["--full_end", str(region_end)])
        print("Command:", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, env=env, text=True)
            print("Prediction script completed.")
        except subprocess.CalledProcessError as e:
            return f"Error running prediction script: {e.stderr}"
        if region_end > region_start + 2097152:
            # Consensus prediction was run, adjust the expected file path accordingly.
            hi_c_matrix_path = os.path.join(output_folder, f"{region_chr}_{region_start}_{region_end}_consensus.npy")
        else:
            hi_c_matrix_path = os.path.join(output_folder, "prediction", "npy", f"{region_chr}_{region_start}.npy")


    print("Waiting for Hi-C matrix file...")
    timeout = 60
    start_time_wait = time.time()
    while not os.path.exists(hi_c_matrix_path) and (time.time() - start_time_wait) < timeout:
        time.sleep(0.5)
    if not os.path.exists(hi_c_matrix_path):
        return f"Error: Hi-C matrix not found at {hi_c_matrix_path}"
    try:
        hi_c_matrix = np.load(hi_c_matrix_path)
        print("Hi-C matrix loaded.")
    except Exception as e:
        return f"Error loading Hi-C matrix: {e}"
    try:
        os.remove(hi_c_matrix_path)
        print("Temporary Hi-C matrix file removed.")
    except Exception as e:
        print(f"Warning: could not remove {hi_c_matrix_path} - {e}")

    if ds_option == "deletion":
        pass
    else:
        del_start, del_width = None, None

    hi_c_chart_config, ctcf_chart_config, atac_chart_config, screening_chart_config, screening_params = \
        prepare_plot_configs(
            hi_c_matrix, region_chr, region_start, region_end, ds_option,
            ctcf_bw_for_model, raw_atac_path, del_start, del_width, norm_atac, norm_ctcf
        )

    hi_c_config_json = json.dumps(hi_c_chart_config)
    ctcf_config_json = json.dumps(ctcf_chart_config)
    atac_config_json = json.dumps(atac_chart_config)
    screening_config_json = json.dumps(screening_chart_config) if screening_chart_config else None
    screening_mode_flag = (ds_option == "screening")
    screening_params_json = json.dumps(screening_params) if screening_mode_flag else "{}"

    if test_mode:
        test_mode_data = {
            "hi_c_chart_config": hi_c_chart_config,
            "ctcf_chart_config": ctcf_chart_config,
            "atac_chart_config": atac_chart_config,
            "screening_chart_config": screening_chart_config,
            "screening_params": screening_params
        }
        try:
            with open("test_mode_output", "w") as f:
                json.dump(test_mode_data, f, indent=2)
            print("Test mode output saved to 'test_mode_output'")
        except Exception as e:
            print("Error saving test mode output:", e)

    if request.headers.get("X-Requested-With") == "XMLHttpRequest":
        return render_template(
            "plots_partial.html",
            hi_c_config=hi_c_config_json,
            ctcf_config=ctcf_config_json,
            atac_config=atac_config_json,
            screening_config=screening_config_json,
            screening_mode=screening_mode_flag,
            screening_params=screening_params_json
        )

    return render_template(
        "index.html",
        hi_c_config=hi_c_config_json,
        ctcf_config=ctcf_config_json,
        atac_config=atac_config_json,
        screening_config=screening_config_json,
        screening_mode=screening_mode_flag,
        screening_params=screening_params_json,
        user_output_folder=get_user_output_folder()
    )

###############################################################################
# Screening Route
###############################################################################
@main.route('/run_screening', methods=['GET'])
def run_screening_endpoint():
    """
    Endpoint that runs the corigami 'screening' workflow on a given region and
    returns the screening data in JSON form (including a D3 line chart config).
    Now prints output to terminal (no 'capture_output').
    Also creates a line-based config with a blank y-axis label.
    """
    print("RUN_SCREENING_ENDPOINT CALLED", flush=True)
    current_app.logger.info("Entered run_screening_endpoint")
    params = request.args.to_dict()
    current_app.logger.info("Received screening parameters: %s", params)
    
    try:
        region_chr = request.args.get('region_chr', 'chr2')
        screen_start = int(request.args.get('region_start', '500000'))
        screen_end = screen_start + WINDOW_WIDTH
        perturb_width = int(request.args.get('perturb_width', '1000'))
        step_size = int(request.args.get('step_size', '1000'))
        
        # Paths
        BASE_DIR = os.path.dirname(current_app.root_path)
        default_atac_path = os.path.join(
            BASE_DIR, "corigami_data", "data", "hg38", "imr90", "genomic_features", "atac.bw"
        )
        atac_bw_path = request.args.get('atac_bw_path', default_atac_path)

        # Output folder
        output_dir = request.args.get('output_dir', get_user_output_folder())
        peaks_file = request.args.get('peaks_file', "").strip()
        if not output_dir or not os.path.exists(output_dir):
            current_app.logger.error("Invalid or missing output directory for screening: %s", output_dir)
            return jsonify({"error": "Invalid or missing output directory for screening."}), 400

    except Exception as e:
        current_app.logger.exception("Invalid screening parameters")
        return jsonify({"error": "Invalid screening parameters", "details": str(e)}), 400

    # Model & CTCF paths
    model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "atac_ctcf_model.ckpt")
    default_ctcf_path = os.path.join(
        BASE_DIR, "corigami_data", "data", "hg38", "imr90", "genomic_features", "ctcf_log2fc.bw"
    )
    ctcf_bw_path = request.args.get('ctcf_bw_path', default_ctcf_path)

    # Check if a predicted CTCF file was generated
    if ctcf_bw_path == "none":
        predicted_ctcf_path = os.path.join(output_dir, "normalized_predicted_ctcf_minmax.bw")
        if os.path.exists(predicted_ctcf_path):
            ctcf_bw_path = predicted_ctcf_path
            current_app.logger.info("Using auto-generated predicted CTCF file for screening: %s", ctcf_bw_path)
        else:
            current_app.logger.info("No auto-generated predicted CTCF file found. Using 'none'.")

    # Final results
    results_file = os.path.join(output_dir, "screening", "screening_results.json")

    # Python environment & script
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH
    screening_script = os.path.join(PYTHON_SRC_PATH, "corigami", "inference", "screening.py")

    # Auto-generate peaks if none provided
    auto_peaks_file = None
    temp_bedgraph = None
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
            current_app.logger.exception("Failed to run bigWigToBedGraph or MACS2")
            return jsonify({"error": "Failed to run bigWigToBedGraph or MACS2.", "details": str(e)}), 500
        except Exception as e:
            current_app.logger.exception("Failed to generate peaks from BigWig")
            return jsonify({"error": "Failed to generate peaks from BigWig.", "details": str(e)}), 500

    # DNA sequence directory (chr*.fa.gz)
    seq_dir = os.path.join(BASE_DIR, "corigami_data", "data", "hg38", "dna_sequence")

    # Build command
    cmd = [
        "python", screening_script,
        "--no-server",
        "--chr", region_chr,
        "--screen-start", str(screen_start),
        "--screen-end", str(screen_end),
        "--model", model_path,
        "--seq", seq_dir,
        "--atac", atac_bw_path,
        "--out", output_dir,
        "--perturb-width", str(perturb_width),
        "--step-size", str(step_size),
        "--plot-impact-score",
        "--save-pred", "--save-perturbation", "--save-diff", "--save-bedgraph",
        "--ctcf", ctcf_bw_path,
        "--peaks-file", peaks_file
    ]

    current_app.logger.info("Running screening command: %s", " ".join(cmd))

    # --- PRINT OUTPUT TO TERMINAL ---
    # Instead of capture_output, we let stdout/stderr go directly to terminal.
    import sys
    try:
        result = subprocess.run(cmd, env=env, check=True, stdout=sys.stdout, stderr=sys.stderr)
    except subprocess.CalledProcessError as e:
        current_app.logger.exception("Screening script failed")
        return jsonify({"error": "Screening script failed.", "details": str(e)}), 500
    except subprocess.TimeoutExpired as e:
        current_app.logger.error("Screening script timed out")
        return jsonify({"error": "Screening script timed out."}), 500

    # Confirm the JSON results file exists
    if not os.path.exists(results_file):
        current_app.logger.error("Screening results JSON file was not generated at: %s", results_file)
        return jsonify({"error": "Screening results JSON file was not generated."}), 500

    # Parse results
    try:
        with open(results_file, "r") as f:
            results = json.load(f)
    except Exception as e:
        current_app.logger.exception("Failed to read screening results JSON")
        return jsonify({"error": "Failed to read screening results JSON.", "details": str(e)}), 500

    # Build a *line* chart config (like your other signals)
    try:
        window_midpoints_mb = results.get("window_midpoints_mb", [])
        impact_scores = results.get("impact_scores", [])

        # Pair up [x, y], then sort by x to ensure left-to-right
        screening_series_data = [[x, y] for x, y in zip(window_midpoints_mb, impact_scores)]
        screening_series_data.sort(key=lambda d: d[0])

        # Compute min & max from the midpoints
        min_val = min(window_midpoints_mb) if window_midpoints_mb else 0
        max_val = max(window_midpoints_mb) if window_midpoints_mb else 0

        # Build line chart config
        screening_chart_config = {
           "chart": {"type": "line", "height": 100},
           "title": {
                "text": "Impact score",
                "style": {"fontSize": "12px"}  # Adjust the font size here
            },
           "xAxis": {
               "min": min_val,
               "max": max_val
           },
           "yAxis": {"title": {"text": ""}},
           "series": [{
             "name": "Screening Score",
             "data": screening_series_data,
             "color": "dodgerblue"
           }]
        }
        screening_config_json = json.dumps(screening_chart_config)

    except Exception as e:
        current_app.logger.exception("Failed to generate screening plot from JSON")
        return jsonify({"error": "Failed to generate screening plot from JSON.", "details": str(e)}), 500

    # Clean up temp files if needed
    if auto_peaks_file and os.path.exists(auto_peaks_file):
        os.remove(auto_peaks_file)
    if temp_bedgraph and os.path.exists(temp_bedgraph):
        os.remove(temp_bedgraph)
    
    # Add chart config to results, return to client
    results["screening_config"] = screening_config_json
    current_app.logger.info("Returning screening results")
    return jsonify(results)
###############################################################################
# Additional Endpoint: List Files in Upload Folder
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
