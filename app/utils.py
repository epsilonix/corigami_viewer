# utils.py

import os
import shutil
from subprocess import CalledProcessError
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

_CLEANUP_INTERVAL = 24 * 3600      # run at most once every 24 h
_last_cleanup_ts  = 0              # updated after each sweep

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
    """
    Remove session folders older than MAX_SESSION_AGE (24 h by default).
    • Runs at most once per _CLEANUP_INTERVAL (also 24 h here).
    • Ignores race‑conditions (FileNotFoundError).
    """
    global _last_cleanup_ts
    now = time.time()
    if now - _last_cleanup_ts < _CLEANUP_INTERVAL:
        return            # too soon – skip this request entirely

    _last_cleanup_ts = now
    cutoff = now - MAX_SESSION_AGE

    for entry in os.scandir(SESSION_BASE_FOLDER):
        if not entry.is_dir():
            continue
        try:
            if entry.stat().st_mtime > cutoff:
                continue  # still fresh
        except FileNotFoundError:
            continue      # vanished meanwhile

        shutil.rmtree(entry.path, ignore_errors=True)

def generate_peaks_from_bigwig_macs2(bw_path, chrom, start, end, outdir):
    temp_bedgraph = os.path.join(outdir, "temp_region.bedGraph")
    auto_peaks    = os.path.join(outdir, "auto_peaks.narrowPeak")

    # 1) Debug: print where the binary lives and what your PATH is
    print("🔍 bigWigToBedGraph on PATH →", shutil.which("bigWigToBedGraph"))
    print("🔍 $PATH =", os.environ.get("PATH"))

    convert_cmd = [
        "bigWigToBedGraph", bw_path, temp_bedgraph,
        f"-chrom={chrom}", f"-start={start}", f"-end={end}"
    ]
    print("📦 running:", " ".join(convert_cmd))
    try:
        result = subprocess.run(
            convert_cmd,
            check=True,
            capture_output=True,
            text=True
        )
        print("bigWigToBedGraph stdout:", result.stdout)
        print("bigWigToBedGraph stderr:", result.stderr)
    except CalledProcessError as e:
        # 2) Log the exact error output
        print("❌ bigWigToBedGraph FAILED:", e)
        print("    stdout:", e.stdout)
        print("    stderr:", e.stderr)
        raise

    macs2_cmd = ["macs2", "bdgpeakcall", "-i", temp_bedgraph,
                 "-o", auto_peaks, "--cutoff", "0.5"]
    print("📦 running:", " ".join(macs2_cmd))
    try:
        result = subprocess.run(
            macs2_cmd,
            check=True,
            capture_output=True,
            text=True
        )
        print("macs2 stdout:", result.stdout)
        print("macs2 stderr:", result.stderr)
    except CalledProcessError as e:
        print("❌ macs2 FAILED:", e)
        print("    stdout:", e.stdout)
        print("    stderr:", e.stderr)
        raise

    return auto_peaks, temp_bedgraph

def get_bigwig_signal(
    bw_path, 
    chrom, 
    start, 
    end, 
    bins=None,      # if None, we compute from bins_per_mb
    bins_per_mb=100 # each MB => 100 bins
):
    import pyBigWig
    region_len_bp = end - start
    region_len_mb = (region_len_bp / 1e6)

    # If user didn’t specify bins, compute it from bins_per_mb
    if bins is None:
        # e.g. 20 Mb => 20*100=2000 bins
        # 2 Mb  => 2*100=200 bins
        # round or clamp as you like
        bins = max(int(region_len_mb * bins_per_mb), 1)

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
    # integer bin_width in basepairs
    bin_width = max(1, (end - start) // bins)

    values, positions = [], []
    for i in range(bins):
        bin_start = start + i * bin_width
        bin_end   = min(bin_start + bin_width, end)

        # clamp
        if bin_start < 0:
            bin_start = 0
        if bin_end > chrom_length:
            bin_end = chrom_length
        if bin_start >= chrom_length:
            # beyond the end
            positions.append((start + end)//2)
            values.append(0.0)
            continue
        if bin_end <= bin_start:
            # invalid after clamp
            positions.append((start + end)//2)
            values.append(0.0)
            continue

        # Get average coverage
        avg = bw.stats(chrom, bin_start, bin_end, type="mean")[0]
        avg = avg if avg is not None else 0.0

        # midpoint of that bin
        mid = (bin_start + bin_end) / 2
        values.append(avg)
        positions.append(mid)

    print("ATAC open:", bw_path, "chroms:", list(bw.chroms().keys())[:3], "...", "has chrCHIM?", "chrCHIM" in bw.chroms())
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

# New functions added for standardized generation and normalizatio
def normalize_bigwig(
    input_bw: str,
    chrom: str,
    start: int,
    end: int,
    method: str,
    output_folder: str,
    prefix: str = "normalized"
) -> str:
    """
    Applies one of three normalization methods to the specified region of a BigWig file:
      - 'none':    Return the original input path
      - 'log':     norm_arr = ln(value + 1)
      - 'minmax':  norm_arr = (value - min) / (max - min)

    Creates a new BigWig from [chrom, start, end). 
    If the user wants the entire chromosome, pass the full range.
    
    Returns the path to the newly created normalized BigWig (unless method='none', 
    in which case the original input path is returned).
    """
    if method == "none":
        return input_bw
    
    if method == "log2fc":
        # assume the input is already in log₂‑fold‑change units
        return input_bw

    try:
        bw_in = pyBigWig.open(input_bw, "r")
    except Exception as e:
        raise Exception(f"Error opening file {input_bw}: {e}")

    if chrom not in bw_in.chroms():
        raise Exception(f"Chrom {chrom} not found in BigWig {input_bw}")

    chrom_length = bw_in.chroms()[chrom]
    region_len = end - start
    if region_len <= 0:
        raise Exception(f"Invalid region: start={start}, end={end}")

    values = bw_in.values(chrom, start, end)
    bw_in.close()

    arr = np.array(values, dtype=np.float64)
    arr = np.nan_to_num(arr, nan=0.0)

    if method == "log":
        norm_arr = np.log(arr + 1.0)
    elif method == "minmax":
        min_val, max_val = np.min(arr), np.max(arr)
        if max_val > min_val:
            norm_arr = (arr - min_val) / (max_val - min_val)
        else:
            norm_arr = np.zeros_like(arr)
    else:
        raise ValueError(f"Unknown normalization method: {method}")

    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, f"{prefix}_{method}.bw")

    try:
        bw_out = pyBigWig.open(output_path, "w")
        bw_out.addHeader([(chrom, chrom_length)])
        bw_out.addEntries(
            chrom,
            start,
            values=norm_arr.tolist(),
            span=1,
            step=1
        )
        bw_out.close()
    except Exception as e:
        raise Exception(f"Error writing normalized BigWig {output_path}: {e}")

    print(f"Normalization ({method}) done. Output: {output_path}")
    return output_path

def predict_ctcf(atac_bw: str,
                 chrom: str,
                 start: int,
                 end: int,
                 ctcf_out_path: str) -> str:
    """
    Run maxATAC predict after seeding ***all*** reference files it
    hard‑codes internally (sizes, genome 2bit, blacklist) into
    <DATA_PATH>/hg38/.  Works both locally and inside the Docker image.
    """
    import os, subprocess, pathlib, shutil, logging

    # ------------------------------------------------------------------ #
    # 0. baked‑in assets
    # ------------------------------------------------------------------ #
    root   = pathlib.Path(__file__).resolve().parents[1]        # repo root
    static = root / "static"

    model_h5   = static / "CTCF_binary_revcomp99_fullModel_RR0_95.h5"
    sizes_src  = static / "hg38.chrom.sizes"
    twob_src   = static / "hg38.2bit"
    bl_src     = static / "hg38_maxatac_blacklist.bed"   #  ← **add this file**

    # ------------------------------------------------------------------ #
    # 1. ensure <DATA_PATH>/hg38/ has the three files
    # ------------------------------------------------------------------ #
    data_path  = pathlib.Path(os.getenv("MAXATAC_DATA_PATH",
                                        pathlib.Path.home() / "opt" / "maxatac" / "data"))
    hg38_dir   = data_path / "hg38"
    hg38_dir.mkdir(parents=True, exist_ok=True)

    for src, dest in (
        (sizes_src, hg38_dir / "hg38.chrom.sizes"),
        (twob_src,  hg38_dir / "hg38.2bit"),
        (bl_src,    hg38_dir / "hg38_maxatac_blacklist.bed"),
    ):
        if not dest.exists():
            shutil.copy(src, dest)      # copy ≈ tiny; use .link_to for hard‑link

    # ------------------------------------------------------------------ #
    # 2. one‑line ROI BED
    # ------------------------------------------------------------------ #
    roi_bed = f"{ctcf_out_path}.tmp.bed"
    with open(roi_bed, "w") as f:
        f.write(f"{chrom}\t{start}\t{end}\n")

    # ------------------------------------------------------------------ #
    # 3. build & run command
    # ------------------------------------------------------------------ #
    out_dir   = os.path.dirname(ctcf_out_path)
    base_name = pathlib.Path(ctcf_out_path).stem
    cmd = [
        "maxatac", "predict",
        "--signal",     atac_bw,
        "--bed",        roi_bed,
        "--out",        out_dir,
        "--name",       base_name,
        "--model",      str(model_h5),
        "--DATA_PATH",  str(data_path)      # now contains all refs
    ]

    logging.info("maxatac cmd: %s", " ".join(cmd))

    try:
        proc = subprocess.run(cmd, check=True, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.debug("maxatac stdout:\n%s", proc.stdout)
        logging.debug("maxatac stderr:\n%s", proc.stderr)
    except subprocess.CalledProcessError as e:
        logging.error("maxatac FAILED (exit %s)\nstdout:\n%s\nstderr:\n%s",
                      e.returncode, e.stdout, e.stderr)
        raise
    finally:
        try:
            os.remove(roi_bed)
        except FileNotFoundError:
            pass

    return ctcf_out_path




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
        "region1": {"chr": region_chr, "start": region_start, "end": region_end},
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
    from two source regions. Genes that overlap each region *partially* are also included (clipped).
    
    Parameters:
      annotation_file (str): Path to the gene annotation file (TSV).
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

            gene_name  = parts[0]
            ensembl_id = parts[1]
            coord_str  = parts[2]
            strand     = parts[3]
            gene_type  = parts[4]

            # Only consider protein-coding genes that do not start with a digit in the name.
            if gene_type != "protein_coding" or gene_name[0].isdigit():
                continue

            try:
                chrom, coords = coord_str.split(":")
                gene_start, gene_end = map(int, coords.split("-"))
            except Exception:
                continue

            # ------------------------------
            # Handle Region 1 (chr1)
            # Include genes that overlap region1 in any way:
            #    gene_end >= start1 and gene_start <= end1
            if chrom == chr1 and (gene_end >= start1) and (gene_start <= end1):
                # Clip to region1 boundaries
                clipped_start = max(gene_start, start1)
                clipped_end   = min(gene_end,   end1)

                new_gene = {
                    "gene": gene_name,
                    "ensembl": ensembl_id,
                    "chrom": "chrCHIM",  # synthetic chromosome name
                    # shift so that region1 starts at 0
                    "start": clipped_start - start1,
                    "end":   clipped_end - start1,
                    "strand": strand,
                    "type": gene_type
                }
                genes_region1.append(new_gene)

            # ------------------------------
            # Handle Region 2 (chr2)
            # Similarly, include partially overlapping genes:
            if chrom == chr2 and (gene_end >= start2) and (gene_start <= end2):
                clipped_start2 = max(gene_start, start2)
                clipped_end2   = min(gene_end,   end2)
                offset = end1 - start1  # length of region1 in bp

                new_gene = {
                    "gene": gene_name,
                    "ensembl": ensembl_id,
                    "chrom": "chrCHIM",
                    # shift so that region2 starts after region1
                    "start": (clipped_start2 - start2) + offset,
                    "end":   (clipped_end2 - start2)   + offset,
                    "strand": strand,
                    "type": gene_type
                }
                genes_region2.append(new_gene)
    
    # Merge the two lists and sort by the new start coordinate
    merged_genes = sorted(genes_region1 + genes_region2, key=lambda x: x["start"])

    # Calculate total length of the synthetic chromosome
    total_length = (end1 - start1) + (end2 - start2)
    offset_mb = (end1 - start1) / 1e6
    total_length_mb = total_length / 1e6

    gene_track_config = {
        "genes": merged_genes,  # the actual gene entries
        "region1": {"chr": "chrCHIM", "start": 0, "end": total_length},
        "chart": {"height": 50},
        "xAxis": {
            "axisBreak": True,
            "leftDomain": [0, offset_mb],
            "rightDomain": [offset_mb, total_length_mb],
            "title": ""
        }
    }
    return gene_track_config

# In utils.py

import pyBigWig
import numpy as np

def extract_bigwig_region_array(
    bw_path, 
    chrom, 
    start, 
    end, 
    do_reverse=False,
    default_val=0.0
):
    """
    Extract raw bigwig values as a 1D numpy array from [start..end).
    If `do_reverse=True`, reverse the array.
    Raises an Exception if the file isn't recognized by pyBigWig,
    or if it lacks the .bw/.bigwig extension.
    """
    import pyBigWig
    import numpy as np

    # (A) Basic extension check
    bw_path_lower = bw_path.lower()
    if not (bw_path_lower.endswith(".bw") or bw_path_lower.endswith(".bigwig")):
        raise Exception(
            f"[extract_bigwig_region_array] File '{bw_path}' does not have a .bw or .bigwig extension."
        )

    length = end - start
    if length <= 0:
        return np.zeros(0, dtype=np.float32)

    # (B) Attempt to open with pyBigWig
    try:
        bw = pyBigWig.open(bw_path, "r")
        if bw is None or bw.chroms() is None:
            bw.close()
            raise Exception(
                f"[extract_bigwig_region_array] '{bw_path}' cannot be opened as a valid bigWig."
            )
    except Exception as e:
        raise Exception(
            f"[extract_bigwig_region_array] Error opening '{bw_path}' with pyBigWig: {e}. "
            "Likely not a valid bigWig file."
        )

    if chrom not in bw.chroms():
        print(f"[extract_bigwig_region_array] Chrom {chrom} not found in {bw_path}")
        bw.close()
        return np.zeros(length, dtype=np.float32)

    chr_len = bw.chroms()[chrom]
    s = max(0, start)
    e = min(end, chr_len)

    arr = np.full(length, default_val, dtype=np.float32)
    if s < e:
        vals = bw.values(chrom, s, e)
        vals = np.nan_to_num(vals, nan=default_val)
        arr[0:(e - s)] = vals

    bw.close()

    if do_reverse:
        arr = arr[::-1]

    return arr

def write_chimeric_bigwig(
    atac_array,
    ctcf_array,
    outdir,
    chim_len,
    chim_name="chrCHIM"
):
    """
    Writes 2 bigWigs (ATAC + CTCF) each of length = chim_len on chrom=chim_name.
    Returns (atac_bw_path, ctcf_bw_path).
    """
    import pyBigWig
    os.makedirs(outdir, exist_ok=True)

    atac_bw_path = os.path.join(outdir, f"{chim_name}_atac.bw")
    ctcf_bw_path = os.path.join(outdir, f"{chim_name}_ctcf.bw")

    # Write ATAC
    bw_atac = pyBigWig.open(atac_bw_path, "w")
    bw_atac.addHeader([(chim_name, chim_len)])
    chunk_size = 1_000_000
    idx = 0
    while idx < chim_len:
        end_idx = min(idx + chunk_size, chim_len)
        subarr = atac_array[idx:end_idx].tolist()
        starts = list(range(idx, end_idx))
        ends   = list(range(idx+1, end_idx+1))
        bw_atac.addEntries([chim_name]*(end_idx-idx), starts, ends=ends, values=subarr)
        idx = end_idx
    bw_atac.close()

    # Write CTCF
    bw_ctcf = pyBigWig.open(ctcf_bw_path, "w")
    bw_ctcf.addHeader([(chim_name, chim_len)])
    idx = 0
    while idx < chim_len:
        end_idx = min(idx + chunk_size, chim_len)
        subarr = ctcf_array[idx:end_idx].tolist()
        starts = list(range(idx, end_idx))
        ends   = list(range(idx+1, end_idx+1))
        bw_ctcf.addEntries([chim_name]*(end_idx-idx), starts, ends=ends, values=subarr)
        idx = end_idx
    bw_ctcf.close()

    return atac_bw_path, ctcf_bw_path

########################################################################
#  FASTA EXTRACTION + CHIMERIC FASTA
########################################################################

def extract_fasta_region(fasta_path, start, end, do_reverse=False, do_flip=False):
    """
    Extract substring from a FASTA (or .fa.gz). 
      - if do_reverse=True -> reverse the string
      - if do_flip=True -> complement the string (like reverse-complement).
        If you only want reverse complement, set do_reverse=True & do_flip=True.
    Returns a string of ACTG in uppercase. 
    """
    import gzip

    # Open gz or plain
    opener = gzip.open if fasta_path.endswith(".gz") else open

    # 1) Read entire sequence (assuming it's a single-chrom FASTA)
    seq_chunks = []
    with opener(fasta_path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq_chunks.append(line)
    full_seq = "".join(seq_chunks).upper()

    # 2) Substring
    sub_seq = full_seq[start:end]
    if not sub_seq:
        return ""

    # 3) Reverse if needed
    if do_reverse:
        sub_seq = sub_seq[::-1]

    # 4) Flip (complement) if needed
    if do_flip:
        complement = {
            "A": "T", "T": "A",
            "G": "C", "C": "G",
            "N": "N"
        }
        # for safety, map unknown => N
        sub_seq = "".join(complement.get(base, "N") for base in sub_seq)

    return sub_seq

def write_chrom_fasta(out_fa, chrom_name, seq_str):
    """
    Write `seq_str` to out_fa (.fa.gz) with a 70-char line wrap.
    Returns final path to gz file.
    """
    import gzip

    if not out_fa.endswith(".gz"):
        out_fa += ".gz"

    line_len = 70
    with gzip.open(out_fa, "wt") as f:
        f.write(f">{chrom_name}\n")
        for i in range(0, len(seq_str), line_len):
            chunk = seq_str[i:i+line_len]
            f.write(chunk + "\n")

    return out_fa

def assemble_chimeric_fasta(
    fa1_path, chr1_start, chr1_end,
    fa2_path, chr2_start, chr2_end,
    *,                          # keyword-only from here
    chr1_reverse=False,
    chr1_flip=False,
    chr2_reverse=False,
    chr2_flip=False,
    output_folder,
    chim_name="chrCHIM",
):
    os.makedirs(output_folder, exist_ok=True)

    seq1 = extract_fasta_region(
        fa1_path, chr1_start, chr1_end,
        do_reverse=chr1_reverse, do_flip=chr1_flip
    )
    seq2 = extract_fasta_region(
        fa2_path, chr2_start, chr2_end,
        do_reverse=chr2_reverse, do_flip=chr2_flip
    )

    chim_seq = seq1 + seq2
    out_fa_path = os.path.join(output_folder, f"{chim_name}.fa")
    gz_path = write_chrom_fasta(out_fa_path, chim_name, chim_seq)
    return gz_path

########################################################################
#  Putting it all together for chimeric signals
########################################################################

def assemble_chimeric_arrays(
    raw_atac_path,
    raw_ctcf_path,
    chr1, chr1_start, chr1_end,
    chr2, chr2_start, chr2_end,
    chr1_reverse=False, chr1_flip=False,
    chr2_reverse=False, chr2_flip=False
):
    arr_atac_1 = extract_bigwig_region_array(
        raw_atac_path, chr1, chr1_start, chr1_end,
        do_reverse=(chr1_reverse or chr1_flip)
    )
    arr_atac_2 = extract_bigwig_region_array(
        raw_atac_path, chr2, chr2_start, chr2_end,
        do_reverse=(chr2_reverse or chr2_flip)
    )
    chim_atac = np.concatenate([arr_atac_1, arr_atac_2])

    arr_ctcf_1 = extract_bigwig_region_array(
        raw_ctcf_path, chr1, chr1_start, chr1_end,
        do_reverse=(chr1_reverse or chr1_flip)
    )
    arr_ctcf_2 = extract_bigwig_region_array(
        raw_ctcf_path, chr2, chr2_start, chr2_end,
        do_reverse=(chr2_reverse or chr2_flip)
    )
    chim_ctcf = np.concatenate([arr_ctcf_1, arr_ctcf_2])
    return chim_atac, chim_ctcf

def clear_folder(folder):
    """
    Delete all files and subdirectories in the given folder.
    """
    if os.path.exists(folder):
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    import shutil
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")
