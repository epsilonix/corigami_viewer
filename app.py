import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import pyBigWig
from flask import Flask, render_template, request
from PIL import Image

app = Flask(__name__)

# Ensure the static directory exists.
if not os.path.exists("static"):
    os.makedirs("static")

# --- Helper Function to Read bigWig Signal ---
def get_bigwig_signal(bw_path, chrom, start, end, bins=256):
    """
    Reads a bigWig file for a given region and averages the signal into `bins` bins.
    Returns a tuple (positions, values).
    """
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

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # --- Set up PYTHONPATH for subprocesses ---
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        PYTHON_SRC_PATH = os.path.join(BASE_DIR, "C.Origami", "src")
        env = os.environ.copy()
        env["PYTHONPATH"] = PYTHON_SRC_PATH

        # Retrieve required user inputs.
        celltype = request.form.get('celltype', 'imr90')
        region_chr = request.form.get('region_chr', 'chr2')
        try:
            region_start = int(request.form.get('region_start', '500000'))
            region_end   = int(request.form.get('region_end', '2500000'))
        except ValueError:
            return "Invalid start or end position. Please enter integer values."
        
        # Get paths for the bigWig files.
        ctcf_bw_path = request.form.get(
            'ctcf_bw_path', 
            f"./corigami_data/data/hg38/{celltype}/genomic_features/ctcf_log2fc.bw"
        )
        atac_bw_path = request.form.get(
            'atac_bw_path', 
            f"./corigami_data/data/hg38/{celltype}/genomic_features/atac.bw"
        )
        output_dir = request.form.get('output_dir', "./output_baseline")
        
        # Determine which optional mode is enabled.
        show_deletion = request.form.get('show_deletion') == "yes"
        show_screening = request.form.get('show_screening') == "yes"
        
        if show_screening:
            try:
                screen_start = int(request.form.get('screen_start', '1250000'))
                screen_end   = int(request.form.get('screen_end', '2250000'))
                perturb_width = int(request.form.get('perturb_width', '1000'))
                step_size = int(request.form.get('step_size', '1000'))
            except ValueError:
                return "Invalid screening parameter(s). Please enter integer values."
            script = "C.Origami/src/corigami/inference/screening.py"
            cmd = [
                "python", script,
                "--chr", region_chr,
                "--celltype", celltype,
                "--model", "./corigami_data/model_weights/corigami_base.ckpt",
                "--seq", "./corigami_data/data/hg38/dna_sequence",
                "--ctcf", ctcf_bw_path,
                "--atac", atac_bw_path,
                "--screen-start", str(screen_start),
                "--screen-end", str(screen_end),
                "--perturb-width", str(perturb_width),
                "--step-size", str(step_size),
                "--plot-impact-score",
                "--save-pred", "--save-perturbation", "--save-diff", "--save-bedgraph",
                "--out", output_dir
            ]
            try:
                subprocess.run(cmd, check=True, env=env)
            except subprocess.CalledProcessError as e:
                return f"Error running screening script: {e}"
            hi_c_matrix_path = f"{output_dir}/{celltype}/screening/imgs/{region_chr}_screen_{screen_start}_{screen_end}_width_{perturb_width}_step_{step_size}.png"
        elif show_deletion:
            del_start_str = request.form.get('del_start', '')
            del_width_str = request.form.get('del_width', '')
            try:
                del_start = int(del_start_str)
                del_width = int(del_width_str)
            except ValueError:
                return "Invalid deletion start or width. Please enter integer values."
            script = "C.Origami/src/corigami/inference/editing.py"
            cmd = [
                "python", script,
                "--chr", region_chr,
                "--celltype", celltype,
                "--start", str(region_start),
                "--model", "./corigami_data/model_weights/corigami_base.ckpt",
                "--seq", "./corigami_data/data/hg38/dna_sequence",
                "--ctcf", ctcf_bw_path,
                "--atac", atac_bw_path,
                "--del-start", str(del_start),
                "--del-width", str(del_width),
                "--out", output_dir
            ]
            try:
                subprocess.run(cmd, check=True, env=env)
            except subprocess.CalledProcessError as e:
                return f"Error running editing script: {e}"
            hi_c_matrix_path = f"{output_dir}/{celltype}/deletion/npy/{region_chr}_{region_start}_del_{del_start}_{del_width}_padding_zero.npy"
        else:
            script = "C.Origami/src/corigami/inference/prediction.py"
            cmd = [
                "python", script,
                "--chr", region_chr,
                "--celltype", celltype,
                "--start", str(region_start),
                "--model", "./corigami_data/model_weights/corigami_base.ckpt",
                "--seq", "./corigami_data/data/hg38/dna_sequence",
                "--ctcf", ctcf_bw_path,
                "--atac", atac_bw_path,
                "--out", output_dir
            ]
            try:
                subprocess.run(cmd, check=True, env=env)
            except subprocess.CalledProcessError as e:
                return f"Error running prediction script: {e}"
            hi_c_matrix_path = f"{output_dir}/{celltype}/prediction/npy/{region_chr}_{region_start}.npy"
        
        if not os.path.exists(hi_c_matrix_path):
            return f"Error: Hi-C matrix not found at {hi_c_matrix_path}"
        try:
            hi_c_matrix = np.load(hi_c_matrix_path)
        except Exception as e:
            return f"Error loading Hi-C matrix: {e}"
        
        # --- Generate the square Hi-C map image (with axes and gridlines as default) ---
        fig_square, ax_square = plt.subplots(figsize=(10, 10))
        ax_square.imshow(hi_c_matrix, cmap='Reds', interpolation='nearest')
        # Do not turn off the axis or grid; leave them as default.
        square_img_path = "static/hic_square.png"
        plt.savefig(square_img_path, dpi=300)
        plt.close(fig_square)
        
        # --- Generate the signal tracks image (separate image) ---
        ctcf_positions, ctcf_values = get_bigwig_signal(ctcf_bw_path, region_chr, region_start, region_end)
        atac_positions, atac_values = get_bigwig_signal(atac_bw_path, region_chr, region_start, region_end)
        if show_deletion:
            deletion_end = del_start + del_width
            ctcf_values = [np.nan if (p >= del_start and p < deletion_end) else v for p, v in zip(ctcf_positions, ctcf_values)]
            atac_values = [np.nan if (p >= del_start and p < deletion_end) else v for p, v in zip(atac_positions, atac_values)]
        fig_signal, ax_signal = plt.subplots(figsize=(10, 4))
        ax_signal.plot(ctcf_positions, ctcf_values, label='CTCF', color='blue')
        ax_signal.plot(atac_positions, atac_values, label='ATAC-seq', color='green')
        ax_signal.set_xlim(region_start, region_end)
        ax_signal.set_xlabel("Genomic Position (bp)")
        ax_signal.set_ylabel("Signal Intensity")
        ax_signal.legend(loc="upper right")
        ax_signal.grid(True)
        plt.tight_layout()
        signal_img_path = "static/signal_tracks.png"
        plt.savefig(signal_img_path, dpi=300)
        plt.close(fig_signal)
        
        # --- Optionally ensure both images have the same width ---
        try:
            hi_img = Image.open(square_img_path)
            signal_img = Image.open(signal_img_path)
            signal_width, _ = signal_img.size
            hi_width, hi_height = hi_img.size
            if hi_width != signal_width:
                new_hi_height = int(hi_height * signal_width / hi_width)
                hi_img = hi_img.resize((signal_width, new_hi_height))
                hi_img.save(square_img_path)
        except Exception as e:
            return f"Error resizing Hi-C image: {e}"
        
        # Render the template with the square Hi-C image and the signal tracks image.
        return render_template("index.html",
                               hic_image="hic_square.png",
                               signal_image="signal_tracks.png")
    else:
        return render_template("index.html")

if __name__ == '__main__':
    app.run(debug=True)
