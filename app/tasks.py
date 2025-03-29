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
    print(f"[RQ-Worker] Executing: {' '.join(cmd)}")
    from rq import get_current_job
    
    process = subprocess.Popen(cmd, env=env, cwd=cwd)
    start_time = time.time()

    while True:
        # 1) Check if the job ended
        ret = process.poll()
        if ret is not None:
            # process is done
            break

        # 2) Check if RQ says "canceled"
        job = get_current_job()
        if job and job.get_status() == 'canceled':
            process.kill()
            raise RuntimeError("Job was canceled by user.")

        # 3) Check timeout
        if timeout and (time.time() - start_time) > timeout:
            process.kill()
            raise RuntimeError("Timed out: " + " ".join(cmd))

        # short sleep so we don’t spin CPU
        time.sleep(0.2)

    # After the process ends, get the final exit code
    ret = process.wait()
    if ret != 0:
        raise RuntimeError(f"Command failed with return code {ret}: {' '.join(cmd)}")


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
    Calls prediction.py with explicit start and end parameters.
    """
    print(f"tasks.py: Sending start={region_start}, end={region_end} to prediction.py")

    cmd = [
        "python", prediction_script,
        "--chr", region_chr,
        "--start", str(region_start),
        "--end", str(region_end),  # Always pass the explicit end coordinate
        "--model", model_path,
        "--seq", seq_dir_for_prediction,
        "--atac", atac_bw_for_model,
        "--out", output_folder,
        "--ctcf", ctcf_bw_for_model
    ]
    run_cmd(cmd, env=env)
    return True



###############################################################################
# TASK 4: run_screening
###############################################################################
def run_screening_task(screening_script, region_chr, screen_start, screen_end,
                       model_path, seq_dir, atac_bw_path, ctcf_bw_path, 
                       peaks_file, output_dir, perturb_width,
                       env=None):
    """
    Runs the screening.py script. Now requires:
    - A valid CTCF bigWig file (ctcf_bw_path)
    - A valid peaks file (peaks_file).
    """
    # 1) Validate inputs
    if not ctcf_bw_path or not os.path.isfile(ctcf_bw_path):
        raise ValueError(f"CTCF file is required and must exist: {ctcf_bw_path}")
    if not peaks_file or not os.path.isfile(peaks_file):
        raise ValueError(f"Peaks file is required and must exist: {peaks_file}")

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
        "--plot-impact-score",
        "--save-pred",
        "--save-perturbation",
        "--save-diff",
        "--save-bedgraph",
        "--ctcf", ctcf_bw_path,
        "--peaks-file", peaks_file
    ]

    run_cmd(cmd, env=env)
    return True

