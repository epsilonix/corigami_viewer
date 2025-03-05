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

