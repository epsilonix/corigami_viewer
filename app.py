import os
import subprocess
import json
import time
import uuid
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import pyBigWig
from flask import Flask, render_template, request, url_for, jsonify, session
from PIL import Image
from scipy.ndimage import rotate
import numpy.ma as ma
from matplotlib.ticker import FuncFormatter

###############################################################################
# Flask App and Config
###############################################################################
app = Flask(__name__)
app.secret_key = "YOUR_SECRET_KEY_CHANGE_THIS"

# Clean up session folders before each request
@app.before_request
def before_request():
    cleanup_session_folders()

# Base directory for session data (we use "user_data")
SESSION_BASE_FOLDER = os.path.join(os.getcwd(), 'user_data')
os.makedirs(SESSION_BASE_FOLDER, exist_ok=True)

# Time in seconds after which an entire session folder will be deleted (e.g., 24 hours)
MAX_SESSION_AGE = 86400

WINDOW_WIDTH = 2097152  # 2 Mb window

###############################################################################
# Folder Management and Cleanup
###############################################################################
def get_user_folder():
    """
    Returns a unique folder path for the current session, e.g.,
    user_data/{session_id}
    """
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
    session_id = session['session_id']
    user_folder = os.path.join(SESSION_BASE_FOLDER, session_id)
    os.makedirs(user_folder, exist_ok=True)
    return user_folder

def get_upload_folder(file_type):
    """
    Returns the path for uploads of a specific file type:
    user_data/{session_id}/uploads/{file_type}
    """
    user_folder = get_user_folder()
    folder = os.path.join(user_folder, "uploads", file_type)
    os.makedirs(folder, exist_ok=True)
    return folder

def get_user_upload_folder():
    """
    Returns the main uploads folder for the session:
    user_data/{session_id}/uploads
    """
    user_folder = get_user_folder()
    uploads_folder = os.path.join(user_folder, "uploads")
    os.makedirs(uploads_folder, exist_ok=True)
    return uploads_folder

def get_user_output_folder():
    """
    Returns the output folder for the session:
    user_data/{session_id}/output
    """
    user_folder = get_user_folder()
    output_folder = os.path.join(user_folder, "output")
    os.makedirs(output_folder, exist_ok=True)
    return output_folder

def cleanup_session_folders():
    now = time.time()
    for session_dir in os.listdir(SESSION_BASE_FOLDER):
        session_path = os.path.join(SESSION_BASE_FOLDER, session_dir)
        if os.path.isdir(session_path):
            age = now - os.path.getmtime(session_path)
            if age > MAX_SESSION_AGE:
                try:
                    for root, dirs, files in os.walk(session_path, topdown=False):
                        for f in files:
                            os.remove(os.path.join(root, f))
                        for d in dirs:
                            os.rmdir(os.path.join(root, d))
                    os.rmdir(session_path)
                    print(f"Removed old session folder: {session_path}")
                except Exception as e:
                    print(f"Error removing session folder {session_path}: {e}")

###############################################################################
# Peak Generation from BigWig using MACS2
###############################################################################
def generate_peaks_from_bigwig_macs2(bw_path, chrom, start, end, outdir):
    temp_bedgraph = os.path.join(outdir, "temp_region.bedGraph")
    auto_peaks = os.path.join(outdir, "auto_peaks.narrowPeak")
    convert_cmd = [
        "/Users/everett/anaconda3/bin/bigWigToBedGraph",
        bw_path,
        temp_bedgraph,
        f"-chrom={chrom}",
        f"-start={start}",
        f"-end={end}"
    ]
    print("Running:", " ".join(convert_cmd))
    result_convert = subprocess.run(
        convert_cmd,
        check=True,
        capture_output=True,
        text=True
    )
    print("bigWigToBedGraph stdout:", result_convert.stdout)
    print("bigWigToBedGraph stderr:", result_convert.stderr)
    macs2_cmd = [
        "macs2", "bdgpeakcall",
        "-i", temp_bedgraph,
        "-o", auto_peaks,
        "--cutoff", "2.0"
    ]
    print("Running:", " ".join(macs2_cmd))
    result_macs2 = subprocess.run(
        macs2_cmd,
        check=True,
        capture_output=True,
        text=True
    )
    print("macs2 bdgpeakcall stdout:", result_macs2.stdout)
    print("macs2 bdgpeakcall stderr:", result_macs2.stderr)
    return auto_peaks, temp_bedgraph

###############################################################################
# Plotting and BigWig Helper Functions
###############################################################################
def create_aligned_figure(figsize, left_margin, bottom_margin):
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=left_margin, right=0.95, top=0.95, bottom=bottom_margin)
    return fig, ax

def get_bigwig_signal(bw_path, chrom, start, end, bins=256):
    try:
        bw = pyBigWig.open(bw_path)
    except Exception as e:
        print(f"Error opening {bw_path}: {e}")
        return [], []
    bin_width = max(1, (end - start) // bins)
    values = []
    positions = []
    for i in range(bins):
        bin_start = start + i * bin_width
        bin_end = bin_start + bin_width
        if bin_end > end:
            bin_end = end
        avg = bw.stats(chrom, bin_start, bin_end, type="mean")[0]
        if avg is None:
            avg = 0
        values.append(avg)
        positions.append((bin_start + bin_end) / 2)
    bw.close()
    return positions, values

def save_uploaded_file(file_storage, default_filename, folder):
    """
    Save an uploaded file to the specified folder and return its full path.
    """
    if file_storage:
        unique_filename = f"{uuid.uuid4()}_{default_filename}"
        filepath = os.path.join(folder, unique_filename)
        file_storage.save(filepath)
        return filepath
    return None

def normalize_file(input_bw, chrom, start, end, method, output_folder, prefix):
    try:
        bw_in = pyBigWig.open(input_bw, "r")
        values = bw_in.values(chrom, start, end)
        chrom_lengths = bw_in.chroms()
        bw_in.close()
    except Exception as e:
        raise Exception(f"Error opening file {input_bw}: {e}")
    arr = np.array(values, dtype=np.float64)
    arr = np.nan_to_num(arr, nan=0)
    if method == "log":
        norm_arr = np.log(arr + 1)
    elif method == "minmax":
        arr_min = np.min(arr)
        arr_max = np.max(arr)
        range_val = arr_max - arr_min if arr_max != arr_min else 1
        norm_arr = (arr - arr_min) / range_val
    else:
        norm_arr = arr
    chr_length = chrom_lengths.get(chrom, end)
    output_path = os.path.join(output_folder, f"{prefix}_{method}.bw")
    try:
        bw_out = pyBigWig.open(output_path, "w")
        bw_out.addHeader([(chrom, chr_length)])
        bw_out.addEntries(chrom, start, values=norm_arr.tolist(), span=1, step=1)
        bw_out.close()
    except Exception as e:
        raise Exception(f"Error writing normalized file: {e}")
    print(f"Normalization ({method}) completed for {prefix}. Normalized file: {output_path}")
    return output_path

###############################################################################
# AJAX File Upload Endpoint
###############################################################################
@app.route('/upload_file', methods=['POST'])
def upload_file():
    # Expect a query parameter "file_type": "atac", "ctcf", or "peaks"
    file_type = request.args.get('file_type', '')
    if not file_type:
        return jsonify({"error": "Missing file_type parameter"}), 400

    # Get the proper upload folder: user_data/{session_id}/uploads/{file_type}
    upload_folder = get_upload_folder(file_type)
    
    # Determine the expected file key based on file_type
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

    # Save file based on type
    if file_type == "atac":
        saved_path = save_uploaded_file(file_storage, "atac.bw", upload_folder)
    elif file_type == "ctcf":
        saved_path = save_uploaded_file(file_storage, "ctcf.bw", upload_folder)
    elif file_type == "peaks":
        saved_path = save_uploaded_file(file_storage, "peaks.narrowPeak", upload_folder)

    print(f"Uploaded {file_type} file saved to: {saved_path}")
    # Return the saved path and a display-friendly name (remove UUID prefix)
    display_name = os.path.basename(saved_path).split("_", 1)[-1]
    return jsonify({"saved_path": saved_path, "display_name": display_name})

###############################################################################
# Main Page: Handle Form Submission or Render
###############################################################################
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'GET':
        user_output_folder = get_user_output_folder()
        return render_template("index.html", screening_mode=False, user_output_folder=user_output_folder)
    
    print("=== Starting new submission ===")
    # Get upload folders for each type (files should have been uploaded already via AJAX)
    atac_upload_folder = get_upload_folder("atac")
    ctcf_upload_folder = get_upload_folder("ctcf")
    peaks_upload_folder = get_upload_folder("peaks")

    output_folder = get_user_output_folder()  # user_data/{session_id}/output
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
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

    # Get file paths from the dropdowns (updated via AJAX upload)
    atac_bw_path = request.form.get('atac_bw_path')
    print("ATAC file path (from dropdown):", os.path.abspath(atac_bw_path))
    ctcf_bw_path = (request.form.get('ctcf_bw_path') or "").strip()
    peaks_file = request.form.get('peaks_file_path', "").strip()

    # Read normalization selections
    norm_atac = request.form.get('norm_atac')    # "none", "log", or "minmax"
    norm_ctcf = request.form.get('norm_ctcf')      # "none", "log", or "minmax"
    training_norm_selection = request.form.get('training_norm')  # "log" or "minmax"
    print(f"ATAC normalization: {norm_atac}")
    print(f"CTCF normalization: {norm_ctcf}")
    print(f"Training normalization selection: {training_norm_selection}")

    # Save raw file paths for plotting
    raw_atac_path = atac_bw_path
    raw_ctcf_path = ctcf_bw_path if ctcf_bw_path and ctcf_bw_path != "none" else None

    # Determine processing route based on whether a CTCF file was provided
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
            subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
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

    if ds_option == "deletion":
        col_avgs = np.mean(hi_c_matrix, axis=0)
        i = hi_c_matrix.shape[1] - 1
        while i >= 0 and col_avgs[i] < 0.1:
            i -= 1
        hi_c_matrix = hi_c_matrix[:, :i+1]
        row_avgs = np.mean(hi_c_matrix, axis=1)
        j = hi_c_matrix.shape[0] - 1
        while j >= 0 and row_avgs[j] < 0.1:
            j -= 1
        hi_c_matrix = hi_c_matrix[:j+1, :]
    hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)
    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[:num_rows // 2, :]

    session_id = session['session_id']
    hic_img_filename = f"hic_square_{session_id}.png"
    ctcf_img_filename = f"ctcf_signal_{session_id}.png"
    atac_img_filename = f"atac_signal_{session_id}.png"

    # Plot Hi-C map
    fig_square, ax_square = create_aligned_figure(figsize=(10, 10), left_margin=0.20, bottom_margin=0.2)
    ax_square.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{int(y)}"))
    for spine in ['top', 'right']:
        ax_square.spines[spine].set_visible(False)
    masked_matrix = ma.masked_where(hi_c_matrix == 0, hi_c_matrix)
    ax_square.imshow(masked_matrix, cmap='Reds', interpolation='nearest')
    ncols = hi_c_matrix.shape[1]
    if ds_option == "deletion":
        del_start_mb = del_start / 1e6
        del_end_mb = (del_start + del_width) / 1e6
        x_tick_mid = (((del_start_mb + del_end_mb) / 2 * 1e6 - region_start) / WINDOW_WIDTH) * (ncols - 1)
        xticks = [0, x_tick_mid, ncols - 1]
        tick_labels = [
            f"{region_start/1e6:.2f}",
            f"{del_start_mb:.2f}-{del_end_mb:.2f}",
            f"{(region_start+WINDOW_WIDTH)/1e6:.2f}"
        ]
        ax_square.set_xticks(xticks)
        ax_square.set_xticklabels(tick_labels)
        for tick in ax_square.xaxis.get_major_ticks():
            if tick.label1.get_text() == f"{del_start_mb:.2f}-{del_end_mb:.2f}":
                tick.label1.set_color("red")
    else:
        ax_square.set_xticks([0, ncols - 1])
        ax_square.xaxis.set_major_formatter(FuncFormatter(
            lambda x, pos: f"{region_start/1e6:.2f}" if x < 1 else f"{(region_start+WINDOW_WIDTH)/1e6:.2f}"
        ))
    ax_square.set_xlabel("Genomic position (Mb)", fontsize=8)
    for label in ax_square.get_yticklabels():
        label.set_color("white")
    hic_img_path = os.path.join("static", hic_img_filename)
    plt.savefig(hic_img_path, dpi=300, bbox_inches='tight')
    plt.close(fig_square)
    print("Hi-C map image saved to", hic_img_path)

    # Plot CTCF signal
    fig_ctcf, ax_ctcf = create_aligned_figure(figsize=(10, 2), left_margin=0.15, bottom_margin=0.15)
    ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_for_model, region_chr, region_start, region_end) if (raw_ctcf_path or ctcf_bw_for_model) else ([], [])
    ctcf_positions_mb = [p/1e6 for p in ctcf_positions]
    ax_ctcf.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.1f}"))
    ax_ctcf.plot(ctcf_positions_mb, ctcf_values, color='blue')
    ax_ctcf.set_xlim(region_start/1e6, region_end/1e6)
    ax_ctcf.set_ylabel("CTCF Signal")
    ax_ctcf.set_xlabel("")
    ax_ctcf.set_xticks([region_start/1e6, region_end/1e6])
    ctcf_img_path = os.path.join("static", ctcf_img_filename)
    plt.savefig(ctcf_img_path, dpi=300, bbox_inches='tight')
    plt.close(fig_ctcf)
    print("CTCF plot saved to", ctcf_img_path)

    # Plot ATAC signal
    fig_atac, ax_atac = create_aligned_figure(figsize=(10, 2), left_margin=0.15, bottom_margin=0.15)
    atac_positions, atac_values = get_bigwig_signal(raw_atac_path, region_chr, region_start, region_end)
    atac_positions_mb = [p/1e6 for p in atac_positions]
    ax_atac.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.1f}"))
    ax_atac.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{y:.1f}"))
    ax_atac.plot(atac_positions_mb, atac_values, color='green')
    ax_atac.set_xlim(region_start/1e6, region_end/1e6)
    ax_atac.set_ylabel("ATAC Signal")
    ax_atac.set_xlabel("")
    ax_atac.set_xticks([region_start/1e6, region_end/1e6])
    atac_img_path = os.path.join("static", atac_img_filename)
    plt.savefig(atac_img_path, dpi=300, bbox_inches='tight')
    plt.close(fig_atac)
    print("ATAC plot saved to", atac_img_path)

    hi_img = Image.open(hic_img_path)
    hi_width, _ = hi_img.size
    for img_path in [ctcf_img_path, atac_img_path]:
        sig_img = Image.open(img_path)
        sig_width, sig_height = sig_img.size
        if sig_width != hi_width:
            new_height = int(sig_height * hi_width / sig_width)
            sig_img = sig_img.resize((hi_width, new_height), resample=Image.Resampling.LANCZOS)
            sig_img.save(img_path)

    # Set up screening mode and parameters if selected
    screening_mode = False
    screening_params = {}
    if ds_option == "screening":
        try:
            perturb_width = int(request.form.get('perturb_width', '1000'))
            step_size = int(request.form.get('step_size', '1000'))
        except ValueError:
            return "Invalid perturb width or step size."
        screening_mode = True
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

    # If this is an AJAX request, return a partial template.
    if request.headers.get("X-Requested-With") == "XMLHttpRequest":
        if ds_option == "screening":
            return render_template(
                "plots_partial.html",
                hic_image=url_for('static', filename=hic_img_filename),
                ctcf_image=url_for('static', filename=ctcf_img_filename),
                atac_image=url_for('static', filename=atac_img_filename),
                screening_mode=True,
                screening_params=json.dumps(screening_params)
            )
        else:
            return render_template(
                "plots_partial.html",
                hic_image=url_for('static', filename=hic_img_filename),
                ctcf_image=url_for('static', filename=ctcf_img_filename),
                atac_image=url_for('static', filename=atac_img_filename)
            )

    # Otherwise, render the full page.
    return render_template(
        "index.html",
        hic_image=url_for('static', filename=hic_img_filename),
        ctcf_image=url_for('static', filename=ctcf_img_filename),
        atac_image=url_for('static', filename=atac_img_filename),
        screening_image="",
        screening_mode=screening_mode,
        screening_params=json.dumps(screening_params)
    )

###############################################################################
# Screening Route
###############################################################################

@app.route('/run_screening', methods=['GET'])
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
    results_file = os.path.join(output_dir, "screening", "screening_results.json")
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
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
        fig, ax = create_aligned_figure(figsize=(10, 2), left_margin=0.15, bottom_margin=0.15)
        ax.bar(window_midpoints_mb, impact_scores, width=bar_width_mb, align='center', color='dodgerblue')
        ax.set_xlim(screen_start_mb, screen_end_mb)
        ax.set_ylabel("Impact Score")
        ax.set_xticks([screen_start_mb, screen_end_mb])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.2f}"))
        session_id = session.get('session_id', 'unknown')
        screening_plot_filename = f"screening_plot_{session_id}.png"
        plot_path = os.path.join("static", screening_plot_filename)
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("Screening bar plot saved to", plot_path)
    except Exception as e:
        print("Error generating screening plot from JSON:", e)
        return jsonify({"error": "Failed to generate screening plot from JSON.", "details": str(e)}), 500

    if auto_peaks_file and os.path.exists(auto_peaks_file):
        os.remove(auto_peaks_file)
        print("Removed auto-generated peaks file:", auto_peaks_file)
    if temp_bedgraph and os.path.exists(temp_bedgraph):
        os.remove(temp_bedgraph)
        print("Removed temporary bedGraph file:", temp_bedgraph)
    
    results["screening_image"] = url_for('static', filename=screening_plot_filename)
    print("Returning results:", results)
    return jsonify(results)

###############################################################################
# Additional Endpoint: List Files in Upload Folder
###############################################################################
@app.route('/list_uploads', methods=['GET'])
def list_uploads():
    # This endpoint accepts a query parameter "file_type" (e.g., "atac", "ctcf", "peaks")
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

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    app.run(debug=True)
