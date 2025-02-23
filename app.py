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
from scipy.ndimage import rotate  # For rotating the matrix
import numpy.ma as ma  # For masking array values
from matplotlib.ticker import FuncFormatter
from scipy.signal import find_peaks  # For peak calling in screening

app = Flask(__name__)
# Set your secret key for session security:
app.secret_key = "YOUR_SECRET_KEY_CHANGE_THIS"

# Base upload folder for all users
BASE_UPLOAD_FOLDER = os.path.join(os.getcwd(), 'uploads')
# Make sure the base folder exists
if not os.path.exists(BASE_UPLOAD_FOLDER):
    os.makedirs(BASE_UPLOAD_FOLDER)

DEFAULT_OUTPUT_DIR = "./output_new"
WINDOW_WIDTH = 2097152  # 2 Mb window

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
    bin_width = (end - start) // bins
    values = []
    positions = []
    for i in range(bins):
        bin_start = start + i * bin_width
        bin_end = bin_start + bin_width
        avg = bw.stats(chrom, bin_start, bin_end, type="mean")[0]
        if avg is None:
            avg = 0
        values.append(avg)
        positions.append((bin_start + bin_end) / 2)
    bw.close()
    return positions, values

def get_user_upload_folder():
    # Generate or retrieve the unique session ID.
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
    session_id = session['session_id']
    # Create a user-specific folder inside the BASE_UPLOAD_FOLDER.
    user_folder = os.path.join(BASE_UPLOAD_FOLDER, session_id)
    if not os.path.exists(user_folder):
        os.makedirs(user_folder)
    return user_folder

def save_uploaded_file(file_storage, default_filename, folder):
    if file_storage:
        unique_filename = f"{uuid.uuid4()}_{default_filename}"
        filepath = os.path.join(folder, unique_filename)
        file_storage.save(filepath)
        return filepath
    return None

def cleanup_upload_folder(user_folder, max_age_seconds=86400):  # For example, 24 hours retention
    now = time.time()
    if os.path.exists(user_folder):
        for filename in os.listdir(user_folder):
            file_path = os.path.join(user_folder, filename)
            if os.path.isfile(file_path):
                file_age = now - os.path.getmtime(file_path)
                if file_age > max_age_seconds:
                    try:
                        os.remove(file_path)
                        print(f"Removed old file: {file_path}")
                    except Exception as e:
                        print(f"Error removing file {file_path}: {e}")

@app.route('/', methods=['GET', 'POST'])
def index():
    user_folder = get_user_upload_folder()  # Get the session-specific folder
    # Clean up old files in the user's folder
    cleanup_upload_folder(user_folder)
    
    if request.method == 'POST':
        print("=== Starting new submission ===")
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

        # Use updated field names (_path)
        atac_bw_path = request.form.get('atac_bw_path')
        print("Absolute ATAC file path:", os.path.abspath(atac_bw_path))
        ctcf_bw_path = (request.form.get('ctcf_bw_path') or "").strip()
        peaks_file = request.form.get('peaks_file_path', "").strip()

        # Handle file uploads; save them in the user_folder.
        uploaded_atac = request.files.get('atac_bw_file')
        uploaded_ctcf = request.files.get('ctcf_bw_file')
        uploaded_peaks = request.files.get('peaks_file')
        if uploaded_atac:
            print("Using uploaded ATAC file.")
            atac_bw_path = save_uploaded_file(uploaded_atac, "atac.bw", user_folder)
        if uploaded_ctcf:
            print("Using uploaded CTCF file.")
            ctcf_bw_path = save_uploaded_file(uploaded_ctcf, "ctcf.bw", user_folder)
        if uploaded_peaks:
            print("Using uploaded Peaks file.")
            peaks_file = save_uploaded_file(uploaded_peaks, "peaks.narrowPeak", user_folder)

        norm_atac = request.form.get('norm_atac')
        norm_ctcf = request.form.get('norm_ctcf')
        print(f"ATAC normalization method: {norm_atac}")
        print(f"CTCF normalization method: {norm_ctcf}")

        if norm_atac != "none":
            print(f"Normalizing ATAC file using {norm_atac} normalization.")
            try:
                bw_in = pyBigWig.open(atac_bw_path, "r")
            except Exception as e:
                return f"Error opening ATAC file: {e}"
            atac_values = bw_in.values(region_chr, region_start, region_end)
            bw_in.close()
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
            normalized_atac = os.path.join(DEFAULT_OUTPUT_DIR, f"normalized_atac_{norm_atac}.bw")
            try:
                bw_in = pyBigWig.open(atac_bw_path, "r")
                chrom_lengths = bw_in.chroms()
                bw_in.close()
            except Exception as e:
                return f"Error reading header from ATAC file: {e}"
            chr_length = chrom_lengths.get(region_chr, region_end)
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

        # Read ds_option to choose the route.
        ds_option = request.form.get('ds_option', 'none')
        print("Selected ds_option:", ds_option)

        if ctcf_bw_path == "none" or ctcf_bw_path == "":
            print(f"No CTCF file provided; generating CTCF file using {norm_ctcf} normalization...")
            roi_file = os.path.join(DEFAULT_OUTPUT_DIR, "temp_roi.bed")
            with open(roi_file, "w") as f:
                f.write(f"{region_chr}\t{region_start}\t{region_end}\n")
            ctcf_generated = os.path.join(DEFAULT_OUTPUT_DIR, "predicted_ctcf.bw")
            generate_cmd = [
                "maxatac", "predict", "--tf", "CTCF",
                "--signal", atac_bw_path,
                "--bed", roi_file,
                "--out", DEFAULT_OUTPUT_DIR,
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
        output_dir = request.form.get('output_dir', DEFAULT_OUTPUT_DIR)

        # Branch based on ds_option
        if ds_option == "screening":
            print("Running screening process for 2Mb window...")
            try:
                perturb_width = int(request.form.get('perturb_width', '1000'))
                step_size = int(request.form.get('step_size', '1000'))
            except ValueError:
                return "Invalid perturb width or step size."
            script = "C.Origami/src/corigami/inference/screening.py"
            print("Running screening script:", script)
            cmd = [
                "python", script, "--no-server",
                "--chr", region_chr,
                "--screen-start", str(region_start),
                "--screen-end", str(region_start + WINDOW_WIDTH),
                "--model", model_path,
                "--seq", "./corigami_data/data/hg38/dna_sequence",
                "--atac", atac_bw_path,
                "--out", output_dir,
                "--perturb-width", str(perturb_width),
                "--step-size", str(step_size),
                "--plot-impact-score",
                "--save-pred", "--save-perturbation", "--save-diff", "--save-bedgraph"
            ]
            cmd.extend(["--ctcf", ctcf_bw_path])
            if peaks_file != "none" and peaks_file != "":
                cmd.extend(["--peaks-file", peaks_file])
            print("Command:", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True, env=env)
                print("Screening script completed.")
            except subprocess.CalledProcessError as e:
                return f"Error running screening script: {e.stderr}"
            hi_c_matrix_path = f"{output_dir}/prediction/npy/{region_chr}_{region_start}.npy"
        elif ds_option == "deletion":
            print("Running deletion process...")
            try:
                del_start = int(request.form.get('del_start', '1500000'))
                del_width = int(request.form.get('del_width', '500000'))
            except ValueError:
                return "Invalid deletion parameters."
            script = "C.Origami/src/corigami/inference/editing.py"
            print("Running editing script:", script)
            cmd = [
                "python", script,
                "--chr", region_chr,
                "--start", str(region_start),
                "--model", model_path,
                "--seq", "./corigami_data/data/hg38/dna_sequence",
                "--atac", atac_bw_path,
                "--del-start", str(del_start),
                "--del-width", str(del_width),
                "--out", output_dir
            ]
            cmd.extend(["--ctcf", ctcf_bw_path])
            print("Command:", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
                print("Editing script completed.")
            except subprocess.CalledProcessError as e:
                return f"Error running editing script: {e.stderr}"
            hi_c_matrix_path = f"{output_dir}/deletion/npy/{region_chr}_{region_start}_del_{del_start}_{del_width}_padding_zero.npy"
        else:
            print("Running standard prediction...")
            script = "C.Origami/src/corigami/inference/prediction.py"
            cmd = [
                "python", script,
                "--chr", region_chr,
                "--start", str(region_start),
                "--model", model_path,
                "--seq", "./corigami_data/data/hg38/dna_sequence",
                "--atac", atac_bw_path,
                "--out", output_dir
            ]
            cmd.extend(["--ctcf", ctcf_bw_path])
            print("Command:", " ".join(cmd))
            try:
                subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
                print("Prediction script completed.")
            except subprocess.CalledProcessError as e:
                return f"Error running prediction script: {e.stderr}"
            hi_c_matrix_path = f"{output_dir}/prediction/npy/{region_chr}_{region_start}.npy"

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
            np.save("debug_hi_c_matrix.npy", hi_c_matrix)
        
        hi_c_matrix = rotate(hi_c_matrix, angle=45, reshape=True)
        num_rows = hi_c_matrix.shape[0]
        hi_c_matrix = hi_c_matrix[:num_rows // 2, :]

        # --- Generate Hi-C map image ---
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
            tick_labels = [f"{region_start/1e6:.2f}", f"{del_start_mb:.2f}-{del_end_mb:.2f}", f"{(region_start+WINDOW_WIDTH)/1e6:.2f}"]
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
        square_img_path = "static/hic_square.png"
        plt.savefig(square_img_path, dpi=300, bbox_inches='tight')
        plt.close(fig_square)
        print("Hi-C map image saved.")

        # --- Generate CTCF plot ---
        print("Generating CTCF plot...")
        fig_ctcf, ax_ctcf = create_aligned_figure(figsize=(10,2), left_margin=0.15, bottom_margin=0.15)
        ctcf_positions, ctcf_values = (get_bigwig_signal(ctcf_bw_path, region_chr, region_start, region_start+WINDOW_WIDTH)
                                        if ctcf_bw_path else ([], []))
        ctcf_positions_mb = [p/1e6 for p in ctcf_positions]
        ax_ctcf.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.1f}"))
        ax_ctcf.plot(ctcf_positions_mb, ctcf_values, color='blue')
        ax_ctcf.set_xlim(region_start/1e6, (region_start+WINDOW_WIDTH)/1e6)
        ax_ctcf.set_ylabel("CTCF Signal")
        ax_ctcf.set_xlabel("")
        ax_ctcf.set_xticks([region_start/1e6, (region_start+WINDOW_WIDTH)/1e6])
        ctcf_img_path = "static/ctcf_signal.png"
        plt.savefig(ctcf_img_path, dpi=300, bbox_inches='tight')
        plt.close(fig_ctcf)
        print("CTCF plot saved.")

        # --- Generate ATAC plot ---
        print("Generating ATAC plot...")
        fig_atac, ax_atac = create_aligned_figure(figsize=(10,2), left_margin=0.15, bottom_margin=0.15)
        atac_positions, atac_values = get_bigwig_signal(atac_bw_path, region_chr, region_start, region_start+WINDOW_WIDTH)
        atac_positions_mb = [p/1e6 for p in atac_positions]
        ax_atac.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.1f}"))
        ax_atac.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{y:.1f}"))
        ax_atac.plot(atac_positions_mb, atac_values, color='green')
        ax_atac.set_xlim(region_start/1e6, (region_start+WINDOW_WIDTH)/1e6)
        ax_atac.set_ylabel("ATAC Signal")
        ax_atac.set_xlabel("")
        ax_atac.set_xticks([region_start/1e6, (region_start+WINDOW_WIDTH)/1e6])
        atac_img_path = "static/atac_signal.png"
        plt.savefig(atac_img_path, dpi=300, bbox_inches='tight')
        plt.close(fig_atac)
        print("ATAC plot saved.")

        # Resize images to ensure consistent width.
        hi_img = Image.open(square_img_path)
        hi_width, _ = hi_img.size
        for img_path in [ctcf_img_path, atac_img_path]:
            sig_img = Image.open(img_path)
            sig_width, sig_height = sig_img.size
            if sig_width != hi_width:
                new_height = int(sig_height * hi_width / sig_width)
                sig_img = sig_img.resize((hi_width, new_height), resample=Image.Resampling.LANCZOS)
                sig_img.save(img_path)

        return render_template("index.html",
                               hic_image=url_for('static', filename="hic_square.png"),
                               ctcf_image=url_for('static', filename="ctcf_signal.png"),
                               atac_image=url_for('static', filename="atac_signal.png"),
                               screening_image="")
    else:
        return render_template("index.html")

@app.route('/run_screening', methods=['GET'])
def run_screening_endpoint():
    try:
        region_chr = request.args.get('region_chr', 'chr2')
        screen_start = int(request.args.get('region_start', '500000'))
        screen_end = screen_start + WINDOW_WIDTH
        perturb_width = int(request.args.get('perturb_width', '1000'))
        step_size = int(request.args.get('step_size', '1000'))
        atac_bw_path = request.args.get('atac_bw_path', "./corigami_data/data/hg38/imr90/genomic_features/atac.bw")
        output_dir = request.args.get('output_dir', "./output_new")
        peaks_file = request.args.get('peaks_file', "").strip()
    except Exception as e:
        return jsonify({"error": "Invalid screening parameters", "details": str(e)}), 400

    model_path = "corigami_data/model_weights/atac_ctcf_model.ckpt"
    ctcf_bw_path = request.args.get('ctcf_bw_path', "./corigami_data/data/hg38/imr90/genomic_features/atac.bw")
    ctcf_bw_param = ["--ctcf", ctcf_bw_path]

    results_file = os.path.join(output_dir, "screening", "screening_results.json")
    
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHON_SRC_PATH
    script = "C.Origami/src/corigami/inference/screening.py"
    cmd = [
        "python", script, "--no-server",
        "--chr", region_chr,
        "--model", model_path,
        "--seq", "./corigami_data/data/hg38/dna_sequence",
        "--atac", atac_bw_path,
        "--screen-start", str(screen_start),
        "--screen-end", str(screen_end),
        "--perturb-width", str(perturb_width),
        "--step-size", str(step_size),
        "--plot-impact-score",
        "--save-pred", "--save-perturbation", "--save-diff", "--save-bedgraph",
        "--out", output_dir
    ]
    cmd.extend(ctcf_bw_param)
    if peaks_file and peaks_file != "none":
        cmd.extend(["--peaks-file", peaks_file])
    
    print("Running screening script with command:")
    print(" ".join(cmd))
    try:
        subprocess.run(cmd, env=env, check=True)
        print("Screening script completed.")
    except subprocess.CalledProcessError as e:
        print("Screening script error:", e, file=sys.stderr)
        return jsonify({"error": "Screening script failed.", "details": str(e)}), 500
    except subprocess.TimeoutExpired:
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
        screen_start_mb = results.get("screen_start_mb", screen_start/1e6)
        screen_end_mb = results.get("screen_end_mb", screen_end/1e6)
        bar_width_mb = perturb_width / 1e6

        fig, ax = create_aligned_figure(figsize=(10, 2), left_margin=0.15, bottom_margin=0.15)
        ax.bar(window_midpoints_mb, impact_scores, width=bar_width_mb, align='center', color='dodgerblue')
        ax.set_xlim(screen_start_mb, screen_end_mb)
        ax.set_xlabel("")
        ax.set_ylabel("Impact Score")
        ax.set_xticks([screen_start_mb, screen_end_mb])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.2f}"))
        
        unique_filename = f"screening_plot_{region_chr}_{screen_start}_{screen_end}_{perturb_width}_{step_size}.png"
        plot_path = os.path.join("static", unique_filename)
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("Screening plot saved.")
    except Exception as e:
        return jsonify({"error": "Failed to generate screening plot from JSON.", "details": str(e)}), 500

    results["screening_image"] = url_for('static', filename=unique_filename)
    return jsonify(results)

if __name__ == '__main__':
    app.run(debug=True)
