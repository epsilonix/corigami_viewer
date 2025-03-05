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
test_mode = False

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

    saved_path = save_uploaded_file(file_storage, file_storage.filename, upload_folder)
    print(f"Uploaded {file_type} file saved to: {saved_path}")
    raw_name = os.path.basename(saved_path)
    display_name = raw_name.split("_", 1)[-1] if "_" in raw_name else raw_name
    return jsonify({"saved_path": saved_path, "display_name": display_name})

###############################################################################
# Main Page: Handle Form Submission or Render
###############################################################################
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
    if region_model == 'V1':
        model_path = "corigami_data/model_weights/v1_jimin.ckpt"
    elif region_model == 'V2':
        model_path = "corigami_data/model_weights/v2_javier.ckpt"
    elif region_model == 'V3':
        model_path = "corigami_data/model_weights/v3_romane.ckpt"
    else:
        return "Invalid model selection."
    
    genome = request.form.get('genome_select', 'hg38')
    if genome == 'hg38':
        seq_dir = "./corigami_data/data/hg38/dna_sequence"
    elif genome == 'mm10':
        seq_dir = "./corigami_data/data/mm10/dna_sequence"
    else:
        return "Invalid genome selection."

    region_chr = request.form.get('region_chr')
    try:
        region_start = int(request.form.get('region_start'))
    except ValueError:
        return "Invalid start position. Please enter an integer value."

    DEFAULT_WINDOW = WINDOW_WIDTH
    region_end_input = request.form.get("region_end")
    if region_end_input:
        try:
            region_end = int(region_end_input)
        except ValueError:
            return "Invalid end position. Please enter an integer value."
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
            normalized_atac = raw_atac_path if norm_atac != "none" else normalize_file(raw_atac_path, region_chr, region_start, region_end,
                                                                                          training_norm, output_folder, "normalized_atac")
            normalized_ctcf = raw_ctcf_path if norm_ctcf != "none" else normalize_file(raw_ctcf_path, region_chr, region_start, region_end,
                                                                                        training_norm, output_folder, "normalized_ctcf")
        else:
            if not training_norm_selection:
                return "Please select the normalization method used during training."
            training_norm = training_norm_selection
            normalized_atac = normalize_file(raw_atac_path, region_chr, region_start, region_end,
                                             training_norm, output_folder, "normalized_atac")
            normalized_ctcf = normalize_file(raw_ctcf_path, region_chr, region_start, region_end,
                                             training_norm, output_folder, "normalized_ctcf")
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
            normalized_predicted_ctcf = normalize_file(ctcf_generated, region_chr, region_start, region_end,
                                                       "minmax", output_folder, "normalized_predicted_ctcf")
            ctcf_bw_for_model = normalized_predicted_ctcf
            norm_ctcf = "minmax normalized, predicted with maxATAC"
            print("Predicted CTCF file generated and normalized:", ctcf_bw_for_model)
        except subprocess.CalledProcessError as e:
            return f"Error generating CTCF file: {e.stderr}"
        finally:
            if os.path.exists(roi_file):
                os.remove(roi_file)
    else:
        ctcf_bw_for_model = normalized_ctcf

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
            "--seq", seq_dir,
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
        hi_c_matrix_path = os.path.join(
            output_folder, "deletion", "npy",
            f"{region_chr}_{region_start}_del_{del_start}_{del_width}_padding_zero.npy"
        )
    else:
        del_start, del_width = None, None
        print("Running standard prediction...")
        script_path = os.path.join(BASE_DIR, "..", "C.Origami", "src", "corigami/inference/prediction.py")
        cmd = [
            "python", script_path,
            "--chr", region_chr,
            "--start", str(region_start),
            "--model", model_path,
            "--seq", seq_dir,
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
        if region_end > region_start + WINDOW_WIDTH:
            hi_c_matrix_path = os.path.join(
                output_folder, f"{region_chr}_{region_start}_{region_end}_consensus.npy"
            )
        else:
            hi_c_matrix_path = os.path.join(
                output_folder, "prediction", "npy", f"{region_chr}_{region_start}.npy"
            )

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

    from app.routes import prepare_plot_configs, prepare_gene_track_config
    hi_c_chart_config, ctcf_chart_config, atac_chart_config, screening_chart_config, screening_params = \
        prepare_plot_configs(
            hi_c_matrix,
            region_chr,
            region_start,
            region_end,
            ds_option,
            ctcf_bw_for_model,
            raw_atac_path,
            del_start,
            del_width,
            norm_atac,
            norm_ctcf
        )

    gene_track_config = prepare_gene_track_config(
        genome,
        region_chr,
        region_start,
        region_end,
        ds_option=ds_option,
        del_start=del_start,
        del_width=del_width
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
            screening_params=screening_params_json,
            norm_atac=norm_atac,
            norm_ctcf=norm_ctcf,
            gene_track_config=json.dumps(gene_track_config)
        )

    return render_template(
        "index.html",
        hi_c_config=hi_c_config_json,
        ctcf_config=ctcf_config_json,
        atac_config=atac_config_json,
        screening_config=screening_config_json,
        screening_mode=screening_mode_flag,
        screening_params=screening_params_json,
        gene_track_config=json.dumps(gene_track_config),
        user_output_folder=get_user_output_folder()
    )

###############################################################################
# Screening Route
###############################################################################
@main.route('/run_screening', methods=['GET'])
def run_screening_endpoint():
    print("RUN_SCREENING_ENDPOINT CALLED", flush=True)
    current_app.logger.info("Entered run_screening_endpoint")
    params = request.args.to_dict()
    current_app.logger.info("Received screening parameters: %s", params)
    
    try:
        region_chr = request.args.get('region_chr', 'chr2')
        screen_start = int(request.args.get('region_start', '500000'))
        screen_end = int(request.args.get('region_end', screen_start + WINDOW_WIDTH))
        perturb_width = int(request.args.get('perturb_width', '1000'))
        step_size = int(request.args.get('step_size', '1000'))
        
        BASE_DIR = os.path.dirname(current_app.root_path)
        default_atac_path = os.path.join(
            BASE_DIR, "corigami_data", "data", "hg38", "imr90", "genomic_features", "atac.bw"
        )
        atac_bw_path = request.args.get('atac_bw_path', default_atac_path)
        output_dir = request.args.get('output_dir', get_user_output_folder())
        peaks_file = request.args.get('peaks_file', "").strip()
        if not output_dir or not os.path.exists(output_dir):
            current_app.logger.error("Invalid or missing output directory for screening: %s", output_dir)
            return jsonify({"error": "Invalid or missing output directory for screening."}), 400

    except Exception as e:
        current_app.logger.exception("Invalid screening parameters")
        return jsonify({"error": "Invalid screening parameters", "details": str(e)}), 400

    region_model = request.args.get('model_select', 'V1')
    if region_model == 'V1':
        model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "v1_jimin.ckpt")
    elif region_model == 'V2':
        model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "v2_javier.ckpt")
    elif region_model == 'V3':
        model_path = os.path.join(BASE_DIR, "corigami_data", "model_weights", "v3_romane.ckpt")
    else:
        return jsonify({"error": "Invalid model selection"}), 400

    default_ctcf_path = os.path.join(
        BASE_DIR, "corigami_data", "data", "hg38", "imr90", "genomic_features", "ctcf_log2fc.bw"
    )
    ctcf_bw_path = request.args.get('ctcf_bw_path', default_ctcf_path)

    if ctcf_bw_path == "none":
        predicted_ctcf_path = os.path.join(output_dir, "normalized_predicted_ctcf_minmax.bw")
        if os.path.exists(predicted_ctcf_path):
            ctcf_bw_path = predicted_ctcf_path
            current_app.logger.info("Using auto-generated predicted CTCF file for screening: %s", ctcf_bw_path)
        else:
            current_app.logger.info("No auto-generated predicted CTCF file found. Using 'none'.")

    results_file = os.path.join(output_dir, "screening", "screening_results.json")

    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH
    screening_script = os.path.join(PYTHON_SRC_PATH, "corigami", "inference", "screening.py")

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

    seq_dir = os.path.join(BASE_DIR, "corigami_data", "data", "hg38", "dna_sequence")

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
    try:
        process = subprocess.Popen(
            cmd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1
        )
        for line in iter(process.stdout.readline, ''):
            if line: print(line, end='')
        for line in iter(process.stderr.readline, ''):
            if line: print(line, end='')
        process.stdout.close()
        process.stderr.close()
        retcode = process.wait()
        if retcode != 0:
            current_app.logger.error("Screening command failed with return code %s", retcode)
            return jsonify({"error": "Screening script failed with return code {}".format(retcode)}), 500
    except Exception as e:
        current_app.logger.exception("Screening script encountered an error")
        return jsonify({"error": "Screening script encountered an error", "details": str(e)}), 500

    if not os.path.exists(results_file):
        current_app.logger.error("Screening results JSON file was not generated at: %s", results_file)
        return jsonify({"error": "Screening results JSON file was not generated."}), 500

    try:
        with open(results_file, "r") as f:
            results = json.load(f)
    except Exception as e:
        current_app.logger.exception("Failed to read screening results JSON")
        return jsonify({"error": "Failed to read screening results JSON.", "details": str(e)}), 500

    try:
        window_midpoints_mb = results.get("window_midpoints_mb", [])
        impact_scores = results.get("impact_scores", [])
        screening_series_data = sorted([[x, y] for x, y in zip(window_midpoints_mb, impact_scores)], key=lambda d: d[0])
        region_start_mb = screen_start / 1e6
        region_end_mb = screen_end / 1e6
        screening_chart_config = {
            "chart": {"type": "line", "height": 100},
            "xAxis": {
                "min": region_start_mb,
                "max": region_end_mb,
                "tickColor": "#000",
                "tickFont": "10px sans-serif"
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
    except Exception as e:
        current_app.logger.exception("Failed to generate screening plot from JSON")
        return jsonify({"error": "Failed to generate screening plot from JSON.", "details": str(e)}), 500

    if auto_peaks_file and os.path.exists(auto_peaks_file):
        os.remove(auto_peaks_file)
    if temp_bedgraph and os.path.exists(temp_bedgraph):
        os.remove(temp_bedgraph)
    
    results["screening_config"] = screening_config_json
    current_app.logger.info("Returning screening results")
    return jsonify(results)

###############################################################################
# List Files in Upload Folder Endpoint
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

def prepare_gene_track_config(genome, region_chr, region_start, region_end, ds_option=None, del_start=None, del_width=None):
    if genome == "hg38":
        annotation_file = "static/genes.gencode.v38.txt"
    elif genome == "mm10":
        annotation_file = "static/genes.gencode.M21.mm10.txt"
    else:
        annotation_file = ""
    
    gene_track_config = {
        "annotationFile": annotation_file,
        "region": {"chr": region_chr, "start": region_start, "end": region_end},
        "chart": {"height": 50}
    }
    
    if ds_option == "deletion" and del_start is not None and del_width is not None:
        gene_track_config["deletionStart"] = del_start
        gene_track_config["deletionEnd"] = del_start + del_width
        gene_track_config["xAxis"] = {
            "axisBreak": True,
            "min": region_start / 1e6,
            "max": region_end / 1e6,
            "leftDomain": [region_start / 1e6, del_start / 1e6],
            "rightDomain": [(del_start + del_width) / 1e6, region_end / 1e6],
            "title": ""
        }
    else:
        gene_track_config["xAxis"] = {
            "min": region_start / 1e6,
            "max": region_end / 1e6,
            "title": ""
        }
    
    return gene_track_config

def prepare_plot_configs(hi_c_matrix, region_chr, region_start, region_end, ds_option,
                         ctcf_bw_for_model, raw_atac_path, del_start=None, del_width=None,
                         norm_atac=None, norm_ctcf=None):
    from scipy.ndimage import rotate
    x_start_mb = region_start / 1e6
    x_end_mb   = region_end / 1e6

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
            "minColor": "#ffffff",
            "maxColor": "#7f0000"
        },
        "series": [{"name": "Hi-C Signal", "data": hi_c_data}],
        "n_cols": n_cols,
        "hideXAxis": False
    }

    ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_for_model, region_chr, region_start, region_end)
    atac_positions, atac_values = get_bigwig_signal(raw_atac_path, region_chr, region_start, region_end)
    ctcf_positions_mb = [p / 1e6 for p in ctcf_positions]
    atac_positions_mb = [p / 1e6 for p in atac_positions]

    const_xaxis = {"min": x_start_mb, "max": x_end_mb, "title": ""}
    ctcf_chart_config = {
        "chart": {"height": 100},
        "xAxis": const_xaxis.copy(),
        "yAxis": {"title": "CTCF Signal"},
        "series": [{"name": "CTCF Signal", "data": [[x, y] for x, y in zip(ctcf_positions_mb, ctcf_values)], "color": "blue"}]
    }
    atac_chart_config = {
        "chart": {"height": 100},
        "xAxis": const_xaxis.copy(),
        "yAxis": {"title": "ATAC Signal"},
        "series": [{"name": "ATAC Signal", "data": [[x, y] for x, y in zip(atac_positions_mb, atac_values)], "color": "green"}]
    }

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
