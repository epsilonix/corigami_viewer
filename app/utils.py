import os
import time
import uuid
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import pyBigWig
import subprocess
from flask import session
from matplotlib.ticker import FuncFormatter

# Constants
WINDOW_WIDTH = 2097152  # 2 Mb window
SESSION_BASE_FOLDER = os.path.join(os.getcwd(), 'user_data')
os.makedirs(SESSION_BASE_FOLDER, exist_ok=True)
MAX_SESSION_AGE = 86400  # 24 hours

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

def generate_peaks_from_bigwig_macs2(bw_path, chrom, start, end, outdir):
    temp_bedgraph = os.path.join(outdir, "temp_region.bedGraph")
    auto_peaks = os.path.join(outdir, "auto_peaks.narrowPeak")
    convert_cmd = [
        "bigWigToBedGraph",
        bw_path,
        temp_bedgraph,
        f"-chrom={chrom}",
        f"-start={start}",
        f"-end={end}"
    ]
    print("Running:", " ".join(convert_cmd))
    result_convert = subprocess.run(convert_cmd, check=True, capture_output=True, text=True)
    print("bigWigToBedGraph stdout:", result_convert.stdout)
    print("bigWigToBedGraph stderr:", result_convert.stderr)
    macs2_cmd = [
        "macs2", "bdgpeakcall",
        "-i", temp_bedgraph,
        "-o", auto_peaks,
        "--cutoff", "2.0"
    ]
    print("Running:", " ".join(macs2_cmd))
    result_macs2 = subprocess.run(macs2_cmd, check=True, capture_output=True, text=True)
    print("macs2 bdgpeakcall stdout:", result_macs2.stdout)
    print("macs2 bdgpeakcall stderr:", result_macs2.stderr)
    return auto_peaks, temp_bedgraph

def get_bigwig_signal(bw_path, chrom, start, end, bins=256):
    try:
        bw = pyBigWig.open(bw_path)
    except Exception as e:
        print(f"Error opening {bw_path}: {e}")
        return [], []

    # If the chromosome doesn't exist in the BigWig, return empty
    if chrom not in bw.chroms():
        print(f"Chromosome {chrom} not found in BigWig {bw_path}")
        bw.close()
        return [], []

    chrom_length = bw.chroms()[chrom]
    bin_width = max(1, (end - start) // bins)

    values, positions = [], []
    for i in range(bins):
        bin_start = start + i * bin_width
        bin_end   = min(bin_start + bin_width, end)

        # --- CLAMP to valid range ---
        if bin_start < 0:
            bin_start = 0
        if bin_end > chrom_length:
            bin_end = chrom_length

        if bin_start >= chrom_length:
            # Beyond the end of the chromosome; fill with 0
            positions.append((start + end)/2)  # or whatever center you'd like
            values.append(0.0)
            continue

        if bin_end <= bin_start:
            # invalid after clamping => 0
            positions.append((start + end)/2)
            values.append(0.0)
            continue

        # Grab the average
        avg = bw.stats(chrom, bin_start, bin_end, type="mean")[0]
        avg = avg if avg is not None else 0.0

        # The midpoint of this bin
        mid = (bin_start + bin_end)/2
        values.append(avg)
        positions.append(mid)

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
        arr_min, arr_max = np.min(arr), np.max(arr)
        norm_arr = (arr - arr_min) / (arr_max - arr_min if arr_max != arr_min else 1)
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

# New functions added for standardized generation and normalization

def normalize_atac(input_bw, chrom, start, end, method, output_folder):
    """
    Normalize an ATAC-seq bigWig file using the given method.
    method: one of 'log', 'minmax', or 'none'.
    If method is 'none', the original file is returned.
    """
    if method == "none":
        return input_bw
    return normalize_file(input_bw, chrom, start, end, method, output_folder, "normalized_atac")

def normalize_ctcf(input_bw, chrom, start, end, method, output_folder):
    """
    Normalize a CTCF bigWig file using the given method.
    method: one of 'log', 'minmax', or 'none'.
    If method is 'none', the original file is returned.
    """
    if method == "none":
        return input_bw
    return normalize_file(input_bw, chrom, start, end, method, output_folder, "normalized_ctcf")

def predict_ctcf(atac_bw, chrom, start, end, output_folder):
    """
    Generate a synthetic (predicted) CTCF bigWig file from an ATAC bigWig using maxATAC.
    The prediction is normalized using minmax normalization.
    Returns the path to the normalized predicted CTCF file.
    """
    roi_file = os.path.join(output_folder, "temp_roi.bed")
    with open(roi_file, "w") as f:
        f.write(f"{chrom}\t{start}\t{end}\n")
    ctcf_generated = os.path.join(output_folder, "predicted_ctcf.bw")
    generate_cmd = [
        "maxatac", "predict", "--tf", "CTCF",
        "--signal", atac_bw,
        "--bed", roi_file,
        "--out", output_folder,
        "--name", "predicted_ctcf"
    ]
    print("Running CTCF prediction command:", " ".join(generate_cmd))
    subprocess.run(generate_cmd, check=True, env=os.environ, capture_output=True, text=True)
    # Normalize the predicted CTCF file (using minmax as default for predicted files)
    normalized_ctcf = normalize_ctcf(ctcf_generated, chrom, start, end, "minmax", output_folder)
    if os.path.exists(roi_file):
        os.remove(roi_file)
    return normalized_ctcf

def predict_peaks(bw_path, chrom, start, end, outdir):
    """
    Generate a peaks file from an ATAC bigWig file using bigWigToBedGraph and MACS2.
    Returns the path to the generated peaks file.
    """
    peaks_file, temp_bedgraph = generate_peaks_from_bigwig_macs2(bw_path, chrom, start, end, outdir)
    # Optionally, remove the temporary bedgraph file after peaks generation.
    if os.path.exists(temp_bedgraph):
        os.remove(temp_bedgraph)
    return peaks_file

###############################################################################
# Gene Track Configuration Helper
###############################################################################
def prepare_gene_track_config(genome, region_chr, region_start, region_end,
                              ds_option=None, del_start=None, del_width=None):
    """
    Prepare a gene-track config (for gencode annotation etc.), with optional
    'axisBreak' for deletion mode.
    """
    if genome == "hg38":
        annotation_file = "static/genes.gencode.v38.txt"
    elif genome == "mm10":
        annotation_file = "static/genes.gencode.M21.mm10.txt"
    else:
        annotation_file = ""

    gene_track_config = {
        "annotationFile": annotation_file,
        "region": {"chr": region_chr, "start": region_start, "end": region_end},
        "chart": {"height": 50},
        "xAxis": {
            "min": region_start / 1e6,
            "max": region_end / 1e6,
            "title": ""
        }
    }

    # If there's a deletion, configure axisBreak for the gene track
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

    return gene_track_config

def prepare_chimeric_gene_track_config(annotation_file, chr1, start1, end1, chr2, start2, end2):
    """
    Create a synthetic gene track configuration for a chimeric chromosome by merging gene annotations
    from two source regions.
    
    Parameters:
      annotation_file (str): Path to the gene annotation file.
      chr1, start1, end1: Coordinates for the first region.
      chr2, start2, end2: Coordinates for the second region.
    
    Returns:
      dict: Gene track configuration containing a merged gene list and synthetic region info.
    """
    print("prepare_chimeric_gene_track_config being loaded")
    genes_region1 = []
    genes_region2 = []
    
    # Read and process the gene annotation file.
    with open(annotation_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            # Expecting at least 5 columns: gene, ensembl, "chr:start-end", strand, type
            if len(parts) < 5:
                continue
            gene_name = parts[0]
            ensembl = parts[1]
            coord_str = parts[2]
            strand = parts[3]
            gene_type = parts[4]
            # Only consider protein-coding genes that do not start with a digit.
            if gene_type != "protein_coding" or gene_name[0].isdigit():
                continue
            try:
                chrom, coords = coord_str.split(":")
                gene_start, gene_end = map(int, coords.split("-"))
            except Exception as e:
                continue

            # For region 1 (chr1)
            if chrom == chr1 and gene_start >= start1 and gene_end <= end1:
                # Transform: shift so that region1 starts at 0.
                new_gene = {
                    "gene": gene_name,
                    "ensembl": ensembl,
                    "chrom": "chrCHIM",  # synthetic chromosome name
                    "start": gene_start - start1,
                    "end": gene_end - start1,
                    "strand": strand,
                    "type": gene_type
                }
                genes_region1.append(new_gene)

            # For region 2 (chr2)
            if chrom == chr2 and gene_start >= start2 and gene_end <= end2:
                # Compute an offset equal to the length of the first region.
                offset = end1 - start1
                new_gene = {
                    "gene": gene_name,
                    "ensembl": ensembl,
                    "chrom": "chrCHIM",
                    "start": (gene_start - start2) + offset,
                    "end": (gene_end - start2) + offset,
                    "strand": strand,
                    "type": gene_type
                }
                genes_region2.append(new_gene)
    
    # Merge the two lists and sort by the new start coordinate.
    merged_genes = sorted(genes_region1 + genes_region2, key=lambda x: x["start"])
    print(f"merged_genes: {merged_genes}")
    # Calculate total length of the synthetic chromosome.
    total_length = (end1 - start1) + (end2 - start2)
    # Convert lengths to megabases for the xAxis domains.
    offset_mb = (end1 - start1) / 1e6
    total_length_mb = total_length / 1e6

    gene_track_config = {
        # Here we pass the merged gene list directly.
        "genes": merged_genes,
        "region": {"chr": "chrCHIM", "start": 0, "end": total_length},
        "chart": {"height": 50},
        "xAxis": {
            "axisBreak": True,
            "leftDomain": [0, offset_mb],
            "rightDomain": [offset_mb, total_length_mb],
            "title": ""
        }
    }
    return gene_track_config

