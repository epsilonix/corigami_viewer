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
from scipy.signal import find_peaks
import sys

###############################################################################
# Flask App and Config
###############################################################################
app = Flask(__name__)
app.secret_key = "YOUR_SECRET_KEY_CHANGE_THIS"

# Clean up session folders before each request
@app.before_request
def before_request():
    cleanup_session_folders()

# Base directory for all session folders
SESSION_BASE_FOLDER = os.path.join(os.getcwd(), 'user_data')
os.makedirs(SESSION_BASE_FOLDER, exist_ok=True)

# Time in seconds after which an entire session folder will be deleted
MAX_SESSION_AGE = 86400  # e.g. 24 hours

WINDOW_WIDTH = 2097152  # 2 Mb window

###############################################################################
# Folder management and cleanup
###############################################################################
def get_user_folder():
    """
    Returns a unique folder path for the current session. Creates one if needed.
    """
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
    session_id = session['session_id']

    user_folder = os.path.join(SESSION_BASE_FOLDER, session_id)
    os.makedirs(user_folder, exist_ok=True)
    return user_folder

def get_user_upload_folder():
    """
    Returns the path to the 'uploads' subfolder for this session, creating it if necessary.
    """
    user_folder = get_user_folder()
    uploads_folder = os.path.join(user_folder, 'uploads')
    os.makedirs(uploads_folder, exist_ok=True)
    return uploads_folder

def get_user_output_folder():
    """
    Returns the path to the 'output' subfolder for this session, creating it if necessary.
    """
    user_folder = get_user_folder()
    output_folder = os.path.join(user_folder, 'output')
    os.makedirs(output_folder, exist_ok=True)
    return output_folder

def cleanup_session_folders():
    """
    Periodically remove session folders that haven't been modified in over MAX_SESSION_AGE seconds.
    Called before each request to keep the server from accumulating stale data.
    """
    now = time.time()
    for session_dir in os.listdir(SESSION_BASE_FOLDER):
        session_path = os.path.join(SESSION_BASE_FOLDER, session_dir)
        if os.path.isdir(session_path):
            age = now - os.path.getmtime(session_path)
            if age > MAX_SESSION_AGE:
                try:
                    # Recursively remove the old session directory
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
import os
import subprocess

def generate_peaks_from_bigwig_macs2(bw_path, chrom, start, end, outdir):
    """
    Naive fallback for generating peaks from a BigWig:
      1) Use bigWigToBedGraph to extract coverage over [chrom:start-end].
      2) Run 'macs2 bdgpeakcall' on that BedGraph to get a .narrowPeak.

    Returns (auto_peaks_path, temp_bedgraph_path).
    """
    temp_bedgraph = os.path.join(outdir, "temp_region.bedGraph")
    auto_peaks = os.path.join(outdir, "auto_peaks.narrowPeak")

    # Step 1: bigWigToBedGraph
    convert_cmd = [
        "bigWigToBedGraph",
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

    # Step 2: macs2 bdgpeakcall --cutoff 2.0 is arbitrary, adjust as needed
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
# Plotting and BigWig helper functions
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
    Save an uploaded file to the specified folder, returning the path to the saved file.
    """
    if file_storage:
        unique_filename = f"{uuid.uuid4()}_{default_filename}"
        filepath = os.path.join(folder, unique_filename)
        file_storage.save(filepath)
        return filepath
    return None

###############################################################################
# Main Page: Handle Form Submission or Render
###############################################################################
@app.route('/', methods=['GET', 'POST'])
def index():
    # If GET: Just show the form
    if request.method == 'GET':
        return render_template("index.html", screening_mode=False)

    # If POST: handle the form submission
    print("=== Starting new submission ===")

    upload_folder = get_user_upload_folder()
    output_folder = get_user_output_folder()

    # Prepare environment for C.Origami calls
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

    # Get file paths or presets
    atac_bw_path = request.form.get('atac_bw_path')
    print("Absolute ATAC file path:", os.path.abspath(atac_bw_path))
    ctcf_bw_path = (request.form.get('ctcf_bw_path') or "").strip()
    peaks_file = request.form.get('peaks_file_path', "").strip()

    # Handle file uploads
    uploaded_atac = request.files.get('atac_bw_file')
    uploaded_ctcf = request.files.get('ctcf_bw_file')
    uploaded_peaks = request.files.get('peaks_file')
    if uploaded_atac:
        print("Using uploaded ATAC file.")
        atac_bw_path = save_uploaded_file(uploaded_atac, "atac.bw", upload_folder)
    if uploaded_ctcf:
        print("Using uploaded CTCF file.")
        ctcf_bw_path = save_uploaded_file(uploaded_ctcf, "ctcf.bw", upload_folder)
    if uploaded_peaks:
        print("Using uploaded Peaks file.")
        peaks_file = save_uploaded_file(uploaded_peaks, "peaks.narrowPeak", upload_folder)

    norm_atac = request.form.get('norm_atac')
    norm_ctcf = request.form.get('norm_ctcf')
    print(f"ATAC normalization method: {norm_atac}")
    print(f"CTCF normalization method: {norm_ctcf}")

    # Normalize ATAC if required
    if norm_atac != "none":
        print(f"Normalizing ATAC file using {norm_atac} normalization.")
        try:
            bw_in = pyBigWig.open(atac_bw_path, "r")
            atac_values = bw_in.values(region_chr, region_start, region_end)
            chrom_lengths = bw_in.chroms()
            bw_in.close()
        except Exception as e:
            return f"Error opening ATAC file: {e}"

        atac_arr = np.array(atac_values, dtype=np.float64)
        atac_arr = np.nan_to_num(atac_arr, nan=0)
        if norm_atac == "log":
            norm_arr = np.log(atac_arr + 1)
        elif norm_atac == "minmax":
            arr_min = np.min(atac_arr)
            arr_max = np.max(atac_arr)
            range_val = arr_max - arr_min if arr_max != arr_min else 1
            norm_arr = (atac_arr - arr_min) / range_val
        else:
            norm_arr = atac_arr

        chr_length = chrom_lengths.get(region_chr, region_end)

        normalized_atac = os.path.join(output_folder, f"normalized_atac_{norm_atac}.bw")
        try:
            bw_out = pyBigWig.open(normalized_atac, "w")
            bw_out.addHeader([(region_chr, chr_length)])
            bw_out.addEntries(region_chr, region_start, values=norm_arr.tolist(), span=1, step=1)
            bw_out.close()
        except Exception as e:
            return f"Error writing normalized ATAC file: {e}"

        print("ATAC normalization completed. Normalized file:", normalized_atac)
        atac_bw_path = normalized_atac
    else:
        print("Skipping ATAC normalization; using raw ATAC file.")

    # Generate CTCF file if none provided
    if ctcf_bw_path == "none" or ctcf_bw_path == "":
        print(f"No CTCF file provided; generating CTCF file using {norm_ctcf} normalization...")
        roi_file = os.path.join(output_folder, "temp_roi.bed")
        with open(roi_file, "w") as f:
            f.write(f"{region_chr}\t{region_start}\t{region_end}\n")

        ctcf_generated = os.path.join(output_folder, "predicted_ctcf.bw")
        generate_cmd = [
            "maxatac", "predict", "--tf", "CTCF",
            "--signal", atac_bw_path,
            "--bed", roi_file,
            "--out", output_folder,
            "--name", "predicted_ctcf"
        ]
        print("Running maxATAC for CTCF generation:", " ".join(generate_cmd))
        try:
            subprocess.run(generate_cmd, check=True, env=env, capture_output=True, text=True)
            ctcf_bw_path = ctcf_generated
            print("CTCF file generated:", ctcf_bw_path)
        except subprocess.CalledProcessError as e:
            return f"Error generating CTCF file: {e.stderr}"
        finally:
            if os.path.exists(roi_file):
                os.remove(roi_file)

    model_path = "corigami_data/model_weights/atac_ctcf_model.ckpt"
    ds_option = request.form.get('ds_option', 'none')
    print("Selected ds_option:", ds_option)

    # Run prediction (or deletion) to generate the Hi-C matrix
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
            "--atac", atac_bw_path,
            "--del-start", str(del_start),
            "--del-width", str(del_width),
            "--out", output_folder,
            "--ctcf", ctcf_bw_path
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
            "--atac", atac_bw_path,
            "--out", output_folder,
            "--ctcf", ctcf_bw_path
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

    # Wait for the Hi-C matrix to appear
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

    # Remove the raw matrix file to save space
    try:
        os.remove(hi_c_matrix_path)
        print("Temporary Hi-C matrix file removed.")
    except Exception as e:
        print(f"Warning: could not remove {hi_c_matrix_path} - {e}")

    # If deletion, trim zero columns/rows
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

    # Rotate and truncate matrix for triangular map
    hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)
    num_rows = hi_c_matrix.shape[0]
    hi_c_matrix = hi_c_matrix[:num_rows // 2, :]

    # Generate unique image filenames for this session
    session_id = session['session_id']
    hic_img_filename = f"hic_square_{session_id}.png"
    ctcf_img_filename = f"ctcf_signal_{session_id}.png"
    atac_img_filename = f"atac_signal_{session_id}.png"

    # Generate Hi-C map image
    print("Generating Hi-C map image...")
    fig_square, ax_square = create_aligned_figure(figsize=(10,10), left_margin=0.20, bottom_margin=0.2)
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

    # Generate CTCF plot
    print("Generating CTCF plot...")
    fig_ctcf, ax_ctcf = create_aligned_figure(figsize=(10,2), left_margin=0.15, bottom_margin=0.15)
    ctcf_positions, ctcf_values = (get_bigwig_signal(ctcf_bw_path, region_chr, region_start, region_end)
                                   if ctcf_bw_path else ([], []))
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

    # Generate ATAC plot
    print("Generating ATAC plot...")
    fig_atac, ax_atac = create_aligned_figure(figsize=(10,2), left_margin=0.15, bottom_margin=0.15)
    atac_positions, atac_values = get_bigwig_signal(atac_bw_path, region_chr, region_start, region_end)
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

    # Resize signal plots to match Hi-C plot width
    hi_img = Image.open(hic_img_path)
    hi_width, _ = hi_img.size
    for img_path in [ctcf_img_path, atac_img_path]:
        sig_img = Image.open(img_path)
        sig_width, sig_height = sig_img.size
        if sig_width != hi_width:
            new_height = int(sig_height * hi_width / sig_width)
            sig_img = sig_img.resize((hi_width, new_height), resample=Image.Resampling.LANCZOS)
            sig_img.save(img_path)

    # If screening is selected, set a flag and pass parameters for AJAX call
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
            "atac_bw_path": atac_bw_path,
            "ctcf_bw_path": ctcf_bw_path,
            "peaks_file": peaks_file,  # might be "none" or an actual path
            "output_dir": output_folder
        }

    return render_template(
        "index.html",
        hic_image=url_for('static', filename=hic_img_filename),
        ctcf_image=url_for('static', filename=ctcf_img_filename),
        atac_image=url_for('static', filename=atac_img_filename),
        screening_image="",  # Will be updated via AJAX if screening_mode is True
        screening_mode=screening_mode,
        screening_params=json.dumps(screening_params)
    )

###############################################################################
# Screening Route
###############################################################################
@app.route('/run_screening', methods=['GET'])
def run_screening_endpoint():
    """
    Run the screening script asynchronously. If no peaks file is provided,
    generate one using MACS2 on the ATAC BigWig for the requested region.
    After screening, remove the auto peaks file so it won't affect future runs.
    """
    try:
        region_chr = request.args.get('region_chr', 'chr2')
        screen_start = int(request.args.get('region_start', '500000'))
        screen_end = screen_start + WINDOW_WIDTH
        perturb_width = int(request.args.get('perturb_width', '1000'))
        step_size = int(request.args.get('step_size', '1000'))
        atac_bw_path = request.args.get('atac_bw_path', "./corigami_data/data/hg38/imr90/genomic_features/atac.bw")
        output_dir = request.args.get('output_dir', "")  # the per-session output folder
        peaks_file = request.args.get('peaks_file', "").strip()

        if not output_dir or not os.path.exists(output_dir):
            return jsonify({"error": "Invalid or missing output directory for screening."}), 400
    except Exception as e:
        return jsonify({"error": "Invalid screening parameters", "details": str(e)}), 400

    model_path = "corigami_data/model_weights/atac_ctcf_model.ckpt"
    ctcf_bw_path = request.args.get('ctcf_bw_path', "./corigami_data/data/hg38/imr90/genomic_features/atac.bw")

    results_file = os.path.join(output_dir, "screening", "screening_results.json")

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH

    screening_script = os.path.join(PYTHON_SRC_PATH, "corigami/inference/screening.py")

    # If user didn't supply peaks file, generate from BigWig
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
        except subprocess.CalledProcessError as e:
            return jsonify({"error": "Failed to run bigWigToBedGraph or MACS2.", "details": str(e)}), 500
        except Exception as e:
            return jsonify({"error": "Failed to generate peaks from BigWig.", "details": str(e)}), 500

    # Build the screening command
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
        subprocess.run(cmd, env=env, check=True)  # This will block until completion
        print("Screening script completed.")
    except subprocess.CalledProcessError as e:
        return jsonify({"error": "Screening script failed.", "details": str(e)}), 500
    except subprocess.TimeoutExpired:
        return jsonify({"error": "Screening script timed out."}), 500

    # Check if screening results exist
    if not os.path.exists(results_file):
        return jsonify({"error": "Screening results JSON file was not generated."}), 500

    # Load screening results
    try:
        with open(results_file, "r") as f:
            results = json.load(f)
    except Exception as e:
        return jsonify({"error": "Failed to read screening results JSON.", "details": str(e)}), 500

    # Generate the screening bar plot
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
        return jsonify({"error": "Failed to generate screening plot from JSON.", "details": str(e)}), 500

    # Clean up auto-generated peaks so it won't affect future runs
    if auto_peaks_file and os.path.exists(auto_peaks_file):
        os.remove(auto_peaks_file)
        print(f"Removed auto-generated peaks file: {auto_peaks_file}")
    if temp_bedgraph and os.path.exists(temp_bedgraph):
        os.remove(temp_bedgraph)
        print(f"Removed temporary bedGraph: {temp_bedgraph}")

    # Return the path to the newly generated screening plot
    results["screening_image"] = url_for('static', filename=screening_plot_filename)
    return jsonify(results)

###############################################################################
# Additional Endpoint: List Files in Upload Folder
###############################################################################
@app.route('/list_uploads', methods=['GET'])
def list_uploads():
    """
    Return a JSON list of files in the current user's upload folder with the session ID removed.
    This lets the client populate dropdowns with the actual file names as they were uploaded.
    """
    user_folder = get_user_upload_folder()
    # List only files (not subfolders)
    files = []
    for f in os.listdir(user_folder):
        full_path = os.path.join(user_folder, f)
        if os.path.isfile(full_path):
            # Remove the session ID prefix (everything up to the first underscore)
            parts = f.split("_", 1)
            display_name = parts[1] if len(parts) == 2 else f
            files.append(display_name)
    return jsonify(files)

###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    app.run(debug=True)
