test_mode = False

import os
import subprocess
import json
import time
import uuid
import numpy as np
from PIL import Image
from scipy.ndimage import rotate
# (Removed Bokeh imports)
from flask import Blueprint, render_template, request, url_for, jsonify, session, current_app

from app.utils import (
    get_user_folder,
    get_upload_folder,
    get_user_upload_folder,
    get_user_output_folder,
    cleanup_session_folders,
    generate_peaks_from_bigwig_macs2,
    create_aligned_figure,
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
                         ctcf_bw_for_model, raw_atac_path, del_start=None, del_width=None):
    """
    Prepares the Highcharts configurations for Hi-C, CTCF, and ATAC plots.
    
    Parameters:
        hi_c_matrix (np.array): The Hi-C matrix.
        region_chr (str): Chromosome.
        region_start (int): Start position.
        region_end (int): End position.
        ds_option (str): Option indicating deletion, screening, or standard prediction.
        ctcf_bw_for_model (str): Path to the (normalized or predicted) CTCF bigwig file.
        raw_atac_path (str): Path to the raw ATAC bigwig file.
        del_start (int, optional): Deletion start position (required if ds_option=='deletion').
        del_width (int, optional): Deletion width (required if ds_option=='deletion').
    
    Returns:
        A tuple with:
            - hi_c_chart_config (dict)
            - ctcf_chart_config (dict)
            - atac_chart_config (dict)
            - screening_chart_config (dict or None)
            - screening_params (dict)
    """
    # Convert genomic coordinates to Mb for plotting
    x_start_mb = region_start / 1e6
    x_end_mb = region_end / 1e6

    # Define transformation function for deletion mode.
    def transform_positions(positions, values, del_start, del_width):
        new_positions = []
        new_values = []
        for p, v in zip(positions, values):
            if p < del_start:
                new_positions.append(p)
                new_values.append(v)
            elif p >= del_start + del_width:
                new_positions.append(p - del_width)
                new_values.append(v)
        return new_positions, new_values

    if ds_option == "deletion" and (del_start is not None and del_width is not None):
        effective_region_end = region_end - del_width
    else:
        effective_region_end = region_end

    # Process Hi-C matrix: rotate and flip as before.
    hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)
    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[:num_rows // 2, :]
    hi_c_matrix = np.flipud(hi_c_matrix)

    # Prepare Hi-C heatmap data.
    n_rows, n_cols = hi_c_matrix.shape
    x_range_mb = x_end_mb - x_start_mb
    x_step = x_range_mb / n_cols
    y_range = 100  # Arbitrary y-range as before
    y_step = y_range / n_rows

    heatmap_data = []
    for i in range(n_rows):
        for j in range(n_cols):
            heatmap_data.append([j, i, float(hi_c_matrix[i, j])])
    x_categories = [round(x_start_mb + (j + 0.5) * x_step, 2) for j in range(n_cols)]
    y_categories = [round((i + 0.5) * y_step, 2) for i in range(n_rows)]

    hi_c_chart_config = {
        "chart": {"type": "heatmap", "height": 500},
        "xAxis": {
            "categories": x_categories,
            "title": {"text": "Genomic position (Mb)"}
        },
        "yAxis": {
            "categories": y_categories,
            "title": {"text": None},
            "reversed": False
        },
        "colorAxis": {
            "min": float(np.nanmin(hi_c_matrix)),
            "max": float(np.nanmax(hi_c_matrix)),
            "minColor": "#ffffff",
            "maxColor": "#7f0000"
        },
        "series": [{
            "name": "Hi-C Signal",
            "data": heatmap_data,
            "borderWidth": 0
        }],
        "tooltip": {
            "pointFormat": "<b>Value:</b> {point.value}"
        }
    }

    # Prepare CTCF Signal chart configuration.
    ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_for_model, region_chr, region_start, region_end)
    if ds_option == "deletion" and (del_start is not None and del_width is not None):
        ctcf_positions, ctcf_values = transform_positions(ctcf_positions, ctcf_values, del_start, del_width)
    ctcf_positions_mb = [p / 1e6 for p in ctcf_positions]
    ctcf_series_data = [[pos, val] for pos, val in zip(ctcf_positions_mb, ctcf_values)]
    ctcf_chart_config = {
        "chart": {"type": "line", "height": 250},
        "title": {"text": "CTCF Signal"},
        "xAxis": {"title": {"text": "Genomic position (Mb)"}, "min": x_start_mb, "max": x_end_mb},
        "yAxis": {"title": {"text": "CTCF Signal"}},
        "series": [{
            "name": "CTCF Signal",
            "data": ctcf_series_data,
            "color": "blue"
        }]
    }

    # Prepare ATAC Signal chart configuration.
    atac_positions, atac_values = get_bigwig_signal(raw_atac_path, region_chr, region_start, region_end)
    if ds_option == "deletion" and (del_start is not None and del_width is not None):
        atac_positions, atac_values = transform_positions(atac_positions, atac_values, del_start, del_width)
    atac_positions_mb = [p / 1e6 for p in atac_positions]
    atac_series_data = [[pos, val] for pos, val in zip(atac_positions_mb, atac_values)]
    atac_chart_config = {
        "chart": {"type": "line", "height": 250},
        "title": {"text": "ATAC Signal"},
        "xAxis": {"title": {"text": "Genomic position (Mb)"}, "min": x_start_mb, "max": x_end_mb},
        "yAxis": {"title": {"text": "ATAC Signal"}},
        "series": [{
            "name": "ATAC Signal",
            "data": atac_series_data,
            "color": "green"
        }]
    }

    # Prepare screening parameters and configuration if in screening mode.
    screening_chart_config = None
    screening_params = {}
    if ds_option == "screening":
        # These parameters could be adjusted as needed; here we read them from the request form.
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
            "atac_bw_path": ctcf_bw_for_model,  # Adjust as necessary
            "ctcf_bw_path": ctcf_bw_for_model,
            "peaks_file": request.form.get('peaks_file_path', "").strip(),
            "output_dir": get_user_output_folder()
        }
        screening_chart_config = {
            "chart": {"type": "column", "height": 250},
            "title": {"text": "Screening Impact Score"},
            "xAxis": {"title": {"text": "Genomic position (Mb)"}},
            "yAxis": {"title": {"text": "Impact Score"}},
            "series": [{
                "name": "Impact Score",
                "data": [],  # Placeholder for screening data
                "color": "dodgerblue"
            }]
        }

    return hi_c_chart_config, ctcf_chart_config, atac_chart_config, screening_chart_config, screening_params

###############################################################################
# Bokeh Plot Configurations
###############################################################################
import numpy as np
from scipy.ndimage import rotate

def prepare_plot_configs(hi_c_matrix, region_chr, region_start, region_end, ds_option,
                         ctcf_bw_for_model, raw_atac_path, del_start=None, del_width=None):
    # Convert genomic coordinates to Mb for plotting
    x_start_mb = region_start / 1e6
    x_end_mb = region_end / 1e6

    # Define transformation function for deletion mode.
    def transform_positions(positions, values, del_start, del_width):
        new_positions = []
        new_values = []
        for p, v in zip(positions, values):
            if p < del_start:
                new_positions.append(p)
                new_values.append(v)
            elif p >= del_start + del_width:
                new_positions.append(p - del_width)
                new_values.append(v)
        return new_positions, new_values

    if ds_option == "deletion" and (del_start is not None and del_width is not None):
        effective_region_end = region_end - del_width
    else:
        effective_region_end = region_end

    # Process Hi-C matrix: rotate and flip as before.
    hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)
    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[:num_rows // 2, :]
    hi_c_matrix = np.flipud(hi_c_matrix)

    # Prepare Hi-C heatmap data.
    n_rows, n_cols = hi_c_matrix.shape
    x_range_mb = x_end_mb - x_start_mb
    x_step = x_range_mb / n_cols
    y_range = 100  # Arbitrary y-range as before
    y_step = y_range / n_rows

    heatmap_data = []
    for i in range(n_rows):
        for j in range(n_cols):
            heatmap_data.append([j, i, float(hi_c_matrix[i, j])])

    x_categories = [round(x_start_mb + (j + 0.5) * x_step, 2) for j in range(n_cols)]
    y_categories = [round((i + 0.5) * y_step, 2) for i in range(n_rows)]

    # Hi-C chart configuration with requested modifications:
    hi_c_chart_config = {
        "chart": {
            "type": "heatmap",
            "height": 350
        },
        # Remove the chart title
        "title": {
            "text": None
        },
        "credits": {
            "enabled": False
        },
        "legend": {
            "enabled": False
        },
        "xAxis": {
            "categories": x_categories,
            "title": {"text": "Genomic position (Mb)"},
            # Force horizontal labels
            "labels": {
                "rotation": 0
            }
        },
        "yAxis": {
            "categories": y_categories,
            "title": {"text": None},
            "reversed": False,
            # Make the axis line/ticks/labels white
            "lineColor": "#ffffff",
            "tickColor": "#ffffff",
            "labels": {
                "style": {
                    "color": "#ffffff"
                }
            }
        },
        # Ensures 0 or below is colored white by flooring at 0
        "colorAxis": {
            "min": 0,
            "max": float(np.nanmax(hi_c_matrix)),
            "minColor": "#ffffff",
            "maxColor": "#7f0000"
        },
        "series": [{
            "name": "Hi-C Signal",
            "data": heatmap_data,
            "borderWidth": 0,
            # Disable hover by disabling mouse tracking on the series
            "tooltip": {
                "enabled": False
            }
        }],
        "tooltip": {
            # Globally disable tooltips
            "enabled": False
        }
    }

    # Prepare CTCF Signal chart configuration.
    ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_for_model, region_chr, region_start, region_end)
    if ds_option == "deletion" and (del_start is not None and del_width is not None):
        ctcf_positions, ctcf_values = transform_positions(ctcf_positions, ctcf_values, del_start, del_width)
    ctcf_positions_mb = [p / 1e6 for p in ctcf_positions]
    ctcf_series_data = [[pos, val] for pos, val in zip(ctcf_positions_mb, ctcf_values)]

    ctcf_chart_config = {
        "chart": {
            "type": "line",
            "height": 150
        },
        # Remove the chart title
        "title": {
            "text": None
        },
        "credits": {
            "enabled": False
        },
        "legend": {
            "enabled": False
        },
        "xAxis": {
            "title": {"text": "Genomic position (Mb)"},
            "min": x_start_mb,
            "max": x_end_mb
        },
        "yAxis": {
            "title": {"text": "CTCF Signal"}
        },
        "series": [{
            "name": "CTCF Signal",
            "data": ctcf_series_data,
            "color": "blue"
        }]
    }

    # Prepare ATAC Signal chart configuration.
    atac_positions, atac_values = get_bigwig_signal(raw_atac_path, region_chr, region_start, region_end)
    if ds_option == "deletion" and (del_start is not None and del_width is not None):
        atac_positions, atac_values = transform_positions(atac_positions, atac_values, del_start, del_width)
    atac_positions_mb = [p / 1e6 for p in atac_positions]
    atac_series_data = [[pos, val] for pos, val in zip(atac_positions_mb, atac_values)]

    atac_chart_config = {
        "chart": {
            "type": "line",
            "height": 150
        },
        # Remove the chart title
        "title": {
            "text": None
        },
        "credits": {
            "enabled": False
        },
        "legend": {
            "enabled": False
        },
        "xAxis": {
            "title": {"text": "Genomic position (Mb)"},
            "min": x_start_mb,
            "max": x_end_mb
        },
        "yAxis": {
            "title": {"text": "ATAC Signal"}
        },
        "series": [{
            "name": "ATAC Signal",
            "data": atac_series_data,
            "color": "green"
        }]
    }

    # Prepare screening parameters and configuration if in screening mode.
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
            "atac_bw_path": ctcf_bw_for_model,  # Adjust as necessary
            "ctcf_bw_path": ctcf_bw_for_model,
            "peaks_file": request.form.get('peaks_file_path', "").strip(),
            "output_dir": get_user_output_folder()
        }
        screening_chart_config = {
            "chart": {
                "type": "column",
                "height": 250
            },
            # Remove the chart title
            "title": {
                "text": None
            },
            "credits": {
                "enabled": False
            },
            "legend": {
                "enabled": False
            },
            "xAxis": {
                "title": {"text": "Genomic position (Mb)"}
            },
            "yAxis": {
                "title": {"text": "Impact Score"}
            },
            "series": [{
                "name": "Impact Score",
                "data": [],  # Placeholder for screening data
                "color": "dodgerblue"
            }]
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
    WINDOW_WIDTH = 500000  # Define your window width if not defined elsewhere.
    region_end = region_start + WINDOW_WIDTH
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
        try:
            subprocess.run(generate_cmd, check=True, env=env, capture_output=True, text=True)
            normalized_predicted_ctcf = normalize_file(ctcf_generated, region_chr, region_start, region_end, "minmax", output_folder, "normalized_predicted_ctcf")
            ctcf_bw_for_model = normalized_predicted_ctcf
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
        print("Command:", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, env=env, text=True)
            print("Prediction script completed.")
        except subprocess.CalledProcessError as e:
            return f"Error running prediction script: {e.stderr}"
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

    # For deletion mode, ensure that del_start and del_width are set.
    if ds_option == "deletion":
        pass
    else:
        del_start, del_width = None, None

    hi_c_chart_config, ctcf_chart_config, atac_chart_config, screening_chart_config, screening_params = \
        prepare_plot_configs(
            hi_c_matrix, region_chr, region_start, region_end, ds_option,
            ctcf_bw_for_model, raw_atac_path, del_start, del_width
        )

    hi_c_config_json = json.dumps(hi_c_chart_config)
    ctcf_config_json = json.dumps(ctcf_chart_config)
    atac_config_json = json.dumps(atac_chart_config)
    screening_config_json = json.dumps(screening_chart_config) if screening_chart_config else None
    screening_mode_flag = (ds_option == "screening")
    screening_params_json = json.dumps(screening_params) if screening_mode_flag else "{}"

    # Save the output configurations to a file if test_mode is enabled.
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
    params = request.args.to_dict()
    print("Received screening parameters:", params)

    try:
        region_chr = request.args.get('region_chr', 'chr2')
        screen_start = int(request.args.get('region_start', '500000'))
        screen_end = screen_start + WINDOW_WIDTH
        perturb_width = int(request.args.get('perturb_width', '1000'))
        step_size = int(request.args.get('step_size', '1000'))
        atac_bw_path = request.args.get('atac_bw_path', "./corigami_data/data/hg38/imr90/genomic_features/atac.bw")
        output_dir = request.args.get('output_dir', "")
        peaks_file = request.args.get('peaks_file', "").strip()
        if not output_dir or not os.path.exists(output_dir):
            return jsonify({"error": "Invalid or missing output directory for screening."}), 400
    except Exception as e:
        return jsonify({"error": "Invalid screening parameters", "details": str(e)}), 400

    model_path = "corigami_data/model_weights/atac_ctcf_model.ckpt"
    ctcf_bw_path = request.args.get('ctcf_bw_path', "./corigami_data/data/hg38/imr90/genomic_features/atac.bw")
    if ctcf_bw_path == "none":
        predicted_ctcf_path = os.path.join(output_dir, "normalized_predicted_ctcf_minmax.bw")
        if os.path.exists(predicted_ctcf_path):
            ctcf_bw_path = predicted_ctcf_path
            print("Using auto-generated predicted CTCF file for screening:", ctcf_bw_path)
        else:
            print("No auto-generated predicted CTCF file found. Using 'none'.")
    results_file = os.path.join(output_dir, "screening", "screening_results.json")
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "..", "C.Origami", "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH

    screening_script = os.path.join(PYTHON_SRC_PATH, "corigami/inference/screening.py")
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
            return jsonify({"error": "Failed to run bigWigToBedGraph or MACS2.", "details": str(e)}), 500
        except Exception as e:
            return jsonify({"error": "Failed to generate peaks from BigWig.", "details": str(e)}), 500

    cmd = [
        "python", screening_script, "--no-server",
        "--chr", region_chr,
        "--screen-start", str(screen_start),
        "--screen-end", str(screen_end),
        "--model", model_path,
        "--seq", "./corigami_data/data/hg38/dna_sequence",
        "--atac", atac_bw_path,
        "--out", output_dir,
        "--perturb-width", str(perturb_width),
        "--step-size", str(step_size),
        "--plot-impact-score",
        "--save-pred", "--save-perturbation", "--save-diff", "--save-bedgraph",
        "--ctcf", ctcf_bw_path,
        "--peaks-file", peaks_file
    ]
    try:
        subprocess.run(cmd, env=env, check=True)
    except subprocess.CalledProcessError as e:
        return jsonify({"error": "Screening script failed.", "details": str(e)}), 500
    except subprocess.TimeoutExpired as e:
        return jsonify({"error": "Screening script timed out."}), 500

    if not os.path.exists(results_file):
        return jsonify({"error": "Screening results JSON file was not generated."}), 500

    try:
        with open(results_file, "r") as f:
            results = json.load(f)
    except Exception as e:
        return jsonify({"error": "Failed to read screening results JSON.", "details": str(e)}), 500

    try:
        window_midpoints_mb = results.get("window_midpoints_mb", [])
        impact_scores = results.get("impact_scores", [])
        screening_series_data = [[x, y] for x, y in zip(window_midpoints_mb, impact_scores)]
        screening_chart_config = {
           "chart": {"type": "column", "height": 250},
           "title": {"text": "Screening Impact Score"},
           "xAxis": {"title": {"text": "Genomic position (Mb)"}},
           "yAxis": {"title": {"text": "Impact Score"}},
           "series": [{
             "name": "Impact Score",
             "data": screening_series_data,
             "color": "dodgerblue"
           }]
        }
        screening_config_json = json.dumps(screening_chart_config)
    except Exception as e:
        return jsonify({"error": "Failed to generate screening plot from JSON.", "details": str(e)}), 500

    if auto_peaks_file and os.path.exists(auto_peaks_file):
        os.remove(auto_peaks_file)
    if temp_bedgraph and os.path.exists(temp_bedgraph):
        os.remove(temp_bedgraph)
    
    results["screening_config"] = screening_config_json
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