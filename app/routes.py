import os
import subprocess
import json
import time
import uuid
import numpy as np
from PIL import Image
from scipy.ndimage import rotate
from bokeh.palettes import Reds256

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

# Bokeh imports
from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.embed import components
from bokeh.models import LinearColorMapper, ColorBar, Span, Label

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
                normalized_atac = normalize_file(raw_atac_path, region_chr, region_start, region_end,
                                                 training_norm, output_folder, "normalized_atac")
            else:
                normalized_atac = raw_atac_path
            if norm_ctcf == "none":
                normalized_ctcf = normalize_file(raw_ctcf_path, region_chr, region_start, region_end,
                                                 training_norm, output_folder, "normalized_ctcf")
            else:
                normalized_ctcf = raw_ctcf_path
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
        hi_c_matrix_path = os.path.join(
            output_folder, "deletion", "npy",
            f"{region_chr}_{region_start}_del_{del_start}_{del_width}_padding_zero.npy"
        )
    else:
        print("Running standard prediction...")
        script_path = os.path.join(PYTHON_SRC_PATH, "corigami/inference/prediction.py")
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
        hi_c_matrix_path = os.path.join(
            output_folder, "prediction", "npy",
            f"{region_chr}_{region_start}.npy"
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

    # Define a helper function to transform positions for signal plots in deletion mode.
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
            # Skip positions within the deletion region.
        return new_positions, new_values

    # For signal plots, if deletion mode is active, transform the positions.
    if ds_option == "deletion":
        del_start = int(request.form.get('del_start', '1500000'))
        del_width = int(request.form.get('del_width', '500000'))
        # For the signal plots, we remove the deleted region.
        effective_region_end = region_end - del_width  # for signal plots only
        # We'll use transformed positions for CTCF and ATAC signals.
    else:
        effective_region_end = region_end

    # Process Hi-C matrix without modifying its content:
    hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)
    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[:num_rows // 2, :]
    hi_c_matrix = np.flipud(hi_c_matrix)

    # ----- Generate Bokeh Plots -----
    # For Hi-C plot, use the full region.
    x_start_mb = region_start / 1e6
    x_end_mb = region_end / 1e6

    p1 = figure(title="Hi-C Map", width=800, height=300,
                x_range=(x_start_mb, x_end_mb),
                tools="xpan,xwheel_zoom,reset", toolbar_location="above")
    mapper = LinearColorMapper(palette=list(reversed(Reds256)),
                               low=np.nanmin(hi_c_matrix),
                               high=np.nanmax(hi_c_matrix))
    p1.image(image=[hi_c_matrix], x=x_start_mb, y=0, dw=(x_end_mb - x_start_mb), dh=100, color_mapper=mapper)
    color_bar = ColorBar(color_mapper=mapper, location=(0,0))
    p1.add_layout(color_bar, 'right')
    p1.xaxis.axis_label = "Genomic position (Mb)"

    # Add deletion breakpoint annotation on Hi-C plot (do not modify hi-c content)
    if ds_option == "deletion":
        del_start_mb = del_start / 1e6
        del_end_mb = (del_start + del_width) / 1e6
        # Place a red vertical line at the deletion start (breakpoint)
        deletion_span = Span(location=del_start_mb, dimension='height', line_color='red', line_dash='dashed', line_width=2)
        p1.add_layout(deletion_span)
        deletion_label = Label(x=del_start_mb, y=-10, x_units='data', y_units='screen',
                               text=f"{del_start_mb:.2f}-{del_end_mb:.2f}", text_color='red',
                               text_font_size='10pt')
        p1.add_layout(deletion_label)

    # Plot 2: CTCF Signal
    ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_for_model, region_chr, region_start, region_end)
    if ds_option == "deletion":
        ctcf_positions, ctcf_values = transform_positions(ctcf_positions, ctcf_values, del_start, del_width)
    ctcf_positions_mb = [p/1e6 for p in ctcf_positions]
    p2 = figure(title="CTCF Signal", width=800, height=150,
                x_range=p1.x_range,
                tools="xpan,xwheel_zoom,reset", toolbar_location="above")
    p2.line(ctcf_positions_mb, ctcf_values, line_color="blue", line_width=2)
    p2.xaxis.axis_label = "Genomic position (Mb)"
    p2.yaxis.axis_label = "CTCF Signal"

    # Plot 3: ATAC Signal
    atac_positions, atac_values = get_bigwig_signal(raw_atac_path, region_chr, region_start, region_end)
    if ds_option == "deletion":
        atac_positions, atac_values = transform_positions(atac_positions, atac_values, del_start, del_width)
    atac_positions_mb = [p/1e6 for p in atac_positions]
    p3 = figure(title="ATAC Signal", width=800, height=150,
                x_range=p1.x_range,
                tools="xpan,xwheel_zoom,reset", toolbar_location="above")
    p3.line(atac_positions_mb, atac_values, line_color="green", line_width=2)
    p3.xaxis.axis_label = "Genomic position (Mb)"
    p3.yaxis.axis_label = "ATAC Signal"

    if ds_option == "screening":
        window_midpoints_mb = results.get("window_midpoints_mb", [])
        impact_scores = results.get("impact_scores", [])
        bar_width_mb = float(request.form.get('perturb_width', 1000)) / 1e6
        p4 = figure(title="Screening Impact Score", width=800, height=150,
                    x_range=p1.x_range,
                    tools="xpan,xwheel_zoom,reset", toolbar_location="above")
        p4.vbar(x=window_midpoints_mb, top=impact_scores, width=bar_width_mb, color="dodgerblue")
        p4.xaxis.axis_label = "Genomic position (Mb)"
        p4.yaxis.axis_label = "Impact Score"
        layout = column(p1, p2, p3, p4)
    else:
        layout = column(p1, p2, p3)

    script, div = components(layout)

    screening_mode = (ds_option == "screening")
    screening_params = {}
    if ds_option == "screening":
        try:
            perturb_width = int(request.form.get('perturb_width', '1000'))
            step_size = int(request.form.get('step_size', '1000'))
            screening_params = {
                "region_chr": region_chr,
                "region_start": region_start,
                "perturb_width": perturb_width,
                "step_size": step_size,
                "atac_bw_path": atac_bw_for_model,
                "ctcf_bw_path": ctcf_bw_for_model,
                "peaks_file": peaks_file,
                "output_dir": output_folder
            }
        except Exception as e:
            screening_params = {}
    screening_params_json = json.dumps(screening_params)

    if request.headers.get("X-Requested-With") == "XMLHttpRequest":
        return render_template(
            "plots_partial.html",
            bokeh_script=script,
            bokeh_div=div,
            screening_mode=screening_mode,
            screening_params=screening_params_json
        )

    return render_template(
        "index.html",
        bokeh_script=script,
        bokeh_div=div,
        screening_mode=screening_mode,
        screening_params=screening_params_json
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
        print("Output directory received:", output_dir)
        if not output_dir or not os.path.exists(output_dir):
            print("Output directory invalid or missing. os.path.exists(output_dir):", os.path.exists(output_dir))
            return jsonify({"error": "Invalid or missing output directory for screening."}), 400
    except Exception as e:
        print("Error parsing screening parameters:", e)
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

    print("Listing contents of output directory:", output_dir)
    try:
        print(os.listdir(output_dir))
    except Exception as e:
        print("Error listing output directory:", e)

    screening_script = os.path.join(PYTHON_SRC_PATH, "corigami/inference/screening.py")
    auto_peaks_file = None
    temp_bedgraph = None
    if not peaks_file or peaks_file == "none":
        print("No peaks file provided; generating via MACS2 from BigWig...")
        try:
            auto_peaks_file, temp_bedgraph = generate_peaks_from_bigwig_macs2(
                bw_path=atac_bw_path,
                chrom=region_chr,
                start=screen_start,
                end=screen_end,
                outdir=output_dir
            )
            peaks_file = auto_peaks_file
            print("Generated peaks file:", peaks_file)
        except subprocess.CalledProcessError as e:
            print("Error running bigWigToBedGraph/MACS2:", e)
            return jsonify({"error": "Failed to run bigWigToBedGraph or MACS2.", "details": str(e)}), 500
        except Exception as e:
            print("Error generating peaks from BigWig:", e)
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
    print("Running screening script with command:")
    print(" ".join(cmd))
    try:
        subprocess.run(cmd, env=env, check=True)
        print("Screening script completed.")
    except subprocess.CalledProcessError as e:
        print("Screening script failed with CalledProcessError:", e)
        return jsonify({"error": "Screening script failed.", "details": str(e)}), 500
    except subprocess.TimeoutExpired as e:
        print("Screening script timed out:", e)
        return jsonify({"error": "Screening script timed out."}), 500

    if not os.path.exists(results_file):
        print("Screening results JSON file not found at:", results_file)
        return jsonify({"error": "Screening results JSON file was not generated."}), 500

    try:
        with open(results_file, "r") as f:
            results = json.load(f)
        print("Screening results JSON loaded successfully.")
    except Exception as e:
        print("Error loading screening results JSON:", e)
        return jsonify({"error": "Failed to read screening results JSON.", "details": str(e)}), 500

    try:
        window_midpoints_mb = results.get("window_midpoints_mb", [])
        impact_scores = results.get("impact_scores", [])
        screen_start_mb = results.get("screen_start_mb", screen_start / 1e6)
        screen_end_mb = results.get("screen_end_mb", screen_end / 1e6)
        bar_width_mb = perturb_width / 1e6
        p_screen = figure(title="Screening Impact Score", width=800, height=150,
                          x_range=(screen_start_mb, screen_end_mb),
                          tools="xpan,xwheel_zoom,reset", toolbar_location="above")
        p_screen.vbar(x=window_midpoints_mb, top=impact_scores, width=bar_width_mb, color="dodgerblue")
        p_screen.xaxis.axis_label = "Genomic position (Mb)"
        p_screen.yaxis.axis_label = "Impact Score"
        screening_layout = column(p_screen)
        screening_script, screening_div = components(screening_layout)
    except Exception as e:
        print("Error generating screening plot from JSON:", e)
        return jsonify({"error": "Failed to generate screening plot from JSON.", "details": str(e)}), 500

    if auto_peaks_file and os.path.exists(auto_peaks_file):
        os.remove(auto_peaks_file)
        print("Removed auto-generated peaks file:", auto_peaks_file)
    if temp_bedgraph and os.path.exists(temp_bedgraph):
        os.remove(temp_bedgraph)
        print("Removed temporary bedGraph file:", temp_bedgraph)
    
    results["screening_image"] = screening_div
    results["screening_script"] = screening_script
    print("Returning results:", results)
    return jsonify(results)

###############################################################################
# Additional Endpoint: List Files in Upload Folder
###############################################################################
@main.route('/list_uploads', methods=['GET'])
def list_uploads():
    file_type = request.args.get('file_type', '')
    folder = get_upload_folder(file_type) if file_type else get_user_upload_folder()
    print(f"Listing files in folder: {folder}")
    files = []
    if os.path.exists(folder):
        for f in os.listdir(folder):
            full_path = os.path.join(folder, f)
            if os.path.isfile(full_path):
                display_name = f.split("_", 1)[-1] if "_" in f else f
                files.append({"value": full_path, "name": display_name})
    print(f"Files found: {files}")
    return jsonify(files)
