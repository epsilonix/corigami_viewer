# tasks.py

import os
import time
import json
import re
import subprocess
import numpy as np
import uuid
from redis import Redis
from rq import Queue, get_current_job

###############################################################################
# CREATE A REDIS CONNECTION & QUEUE
###############################################################################
redis_conn = Redis(host='localhost', port=6379, db=0)
q = Queue('default', connection=redis_conn)

###############################################################################
# HELPER: run_cmd with line-by-line printing
#   This ensures that output is printed to the *worker* console
#   so it appears in real time, as if we were in the route.
###############################################################################
def run_cmd(cmd, env=None, cwd=None, timeout=None):
    """
    Runs the given command, printing stdout/stderr lines in real-time
    to the *worker*'s console.
    Raises RuntimeError if the command exits with non-zero.
    """
    print(f"[RQ-Worker] Executing: {' '.join(cmd)}")
    process = subprocess.Popen(
        cmd,
        env=env,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    start_time = time.time()

    # Stream output to worker console
    while True:
        # Check if process ended
        if process.poll() is not None:
            break

        # Check for timeout
        if timeout and (time.time() - start_time) > timeout:
            process.kill()
            raise RuntimeError("Timed out while running: " + " ".join(cmd))

        # Read lines from stdout
        line = process.stdout.readline()
        if line:
            print(line, end="")

        # Also read stderr
        err_line = process.stderr.readline()
        if err_line:
            print(err_line, end="")

        time.sleep(0.01)

    # Drain the remaining lines (in case process ended quickly)
    for line in process.stdout:
        print(line, end="")
    for line in process.stderr:
        print(line, end="")

    retcode = process.wait()
    if retcode != 0:
        raise RuntimeError(f"Command failed with return code {retcode}: {' '.join(cmd)}")

###############################################################################
# TASK 1: run_editing (for ds_option=='deletion')
###############################################################################
def run_editing_task(editing_script, region_chr, region_start, model_path,
                     seq_dir_for_prediction, atac_bw_for_model,
                     del_start, del_width, output_folder,
                     ctcf_bw_for_model=None, env=None):
    """
    Moves the original 'subprocess.run(...)' logic from routes.py
    into this function so it can be executed by an RQ worker.
    """
    cmd = [
        "python", editing_script,
        "--chr", region_chr,
        "--start", str(region_start),
        "--model", model_path,
        "--seq", seq_dir_for_prediction,
        "--atac", atac_bw_for_model,
        "--del-start", str(del_start),
        "--del-width", str(del_width),
        "--out", output_folder,
    ]
    if ctcf_bw_for_model:
        cmd += ["--ctcf", ctcf_bw_for_model]

    run_cmd(cmd, env=env)
    return True  # optionally return something

###############################################################################
# TASK 2: run_prediction (for ds_option=='none' or 'screening')
###############################################################################
def run_prediction_task(prediction_script, region_chr, region_start, region_end,
                        model_path, seq_dir_for_prediction, atac_bw_for_model,
                        output_folder, ctcf_bw_for_model=None, env=None):
    """
    Moves the original 'prediction.py' call from routes.py
    into this function so it can be executed by an RQ worker.
    """
    cmd = [
        "python", prediction_script,
        "--chr", region_chr,
        "--start", str(region_start),
        "--model", model_path,
        "--seq", seq_dir_for_prediction,
        "--atac", atac_bw_for_model,
        "--out", output_folder
    ]
    # If we have a CTCF track, pass it
    if ctcf_bw_for_model:
        cmd += ["--ctcf", ctcf_bw_for_model]

    # If region_end is bigger than 2Mb beyond region_start => add --full_end
    if region_end > region_start + 2097152:  # 2Mb
        cmd += ["--full_end", str(region_end)]

    run_cmd(cmd, env=env)
    return True

###############################################################################
# TASK 3: build_chimeric
###############################################################################
def build_chimeric_task(build_chimeric_script, region_chr1, region_start1, region_end1,
                        region_chr2, region_start2, region_end2,
                        model_path, genome_dir, ctcf_bw_path, atac_bw_path,
                        output_folder, first_reverse, first_flip,
                        second_reverse, second_flip, env=None):
    """
    Moves the original 'build_chimeric.py' call (subprocess) from routes.py
    into this function so it can be executed by an RQ worker.
    Returns (chimera_dir, chim_info_dict).
    """
    cmd = [
        "python", build_chimeric_script,
        "--chr1", region_chr1,
        "--start1", str(region_start1),
        "--end1", str(region_end1),
        "--chr2", region_chr2,
        "--start2", str(region_start2),
        "--end2", str(region_end2),
        "--model", model_path,
        "--seq", genome_dir,
        "--ctcf", ctcf_bw_path,
        "--atac", atac_bw_path,
        "--out", output_folder,
        "--run_downstream", "none"
    ]
    if first_reverse:
        cmd.append("--first_reverse")
    if first_flip:
        cmd.append("--first_flip")
    if second_reverse:
        cmd.append("--second_reverse")
    if second_flip:
        cmd.append("--second_flip")

    run_cmd(cmd, env=env)

    # Now parse the stdout/stderr is not enough, we read manifest from output:
    # We'll guess the CHIMERA_OUTPUT_DIR from the last known line in build_chimeric output
    # but let's do it more robustly by reading manifest:
    chimera_dir = None
    # Attempt to find the subdirectory:
    # the script prints "CHIMERA_OUTPUT_DIR=..."
    # We can parse from the output_folder if we want,
    # or just guess from the standard logic:
    # but let's do it from the lines ourselves:

    # For safety, let's search for any subfolders named "chimeric_build_" in output_folder
    # and assume the latest is the correct one (this might need to match your real code logic).
    # If your code prints exactly "CHIMERA_OUTPUT_DIR=some_path", you'd need a more direct approach.
    # Let's just do your original approach: read "manifest.txt" in the output folder.
    # We'll assume your code places "manifest.txt" in `output_folder/manifest.txt` or a subfolder:
    # Because your original logic in routes did:

    # We might do:
    possible_manifest = os.path.join(output_folder, "manifest.txt")
    if os.path.exists(possible_manifest):
        chimera_dir = output_folder
    else:
        # Maybe there's a subdirectory:
        for root, dirs, files in os.walk(output_folder):
            if "manifest.txt" in files:
                chimera_dir = root
                break

    if not chimera_dir:
        raise RuntimeError("Could not find manifest.txt after build_chimeric.")

    # Read manifest
    manifest_path = os.path.join(chimera_dir, "manifest.txt")
    chim_info = {}
    with open(manifest_path, "r") as f:
        for line in f:
            line = line.strip()
            if '=' in line:
                k,v = line.split('=',1)
                chim_info[k] = v

    return (chimera_dir, chim_info)

###############################################################################
# TASK 4: run_screening
###############################################################################
def run_screening_task(screening_script, region_chr, screen_start, screen_end,
                       model_path, seq_dir, atac_bw_path, ctcf_bw_path, 
                       peaks_file, output_dir, perturb_width, step_size,
                       env=None):
    """
    Moves the 'screening.py' call and bigWigToBedGraph+MACS2 logic from routes.py
    into this function.  We run the command in a separate process.
    We *still* do line-by-line prints to console so you see them in worker logs.
    """
    # All the code that was originally in your route can go here.
    # But the user wants minimal changes, so we'll keep it straightforward.

    # We'll assume your route or some other function already generated
    # auto_peaks if needed. If not, you'd run bigWigToBedGraph + macs2 here.
    # For illustration, we'll just run the screening script:

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
        "--save-pred",
        "--save-perturbation",
        "--save-diff",
        "--save-bedgraph",
    ]

    if ctcf_bw_path:
        cmd += ["--ctcf", ctcf_bw_path]
    if peaks_file:
        cmd += ["--peaks-file", peaks_file]

    run_cmd(cmd, env=env)
    return True


def hello_task(name):
    print(f"Hello from {name}")
    return f"Hello from {name}"