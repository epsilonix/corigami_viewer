"""
app/api.py  – JSON API that farms heavy work out to an RQ worker
─────────────────────────────────────────────────────────────────
The front‑end now:

  1. POSTs its entire form as JSON →  /api/predict
  2. polls                 /api/job/<id>
  3. once finished, GETs   /api/job/<id>/html   (fragment to insert)

This file contains everything needed; no hidden imports from routes.py
besides the heavy `run_prediction_and_render` function.
"""
from __future__ import annotations
from flask import Blueprint, request, jsonify, current_app, session, render_template
from inspect import signature
import os, shutil
import numpy as np

from app.utils import (
    get_user_output_folder, WINDOW_WIDTH,
    normalize_bigwig, predict_ctcf,
    extract_bigwig_region_array, write_chimeric_bigwig,
    assemble_chimeric_fasta,
    prepare_gene_track_config, prepare_chimeric_gene_track_config,
)  

api = Blueprint("api", __name__, url_prefix="/api")

# ═════════════════════════════════ helpers ════════════════════════════════ #

def _allowed_kwargs(func):
    """names of kwargs accepted by the worker function"""
    return {
        p.name
        for p in signature(func).parameters.values()
        if p.kind in (p.KEYWORD_ONLY, p.POSITIONAL_OR_KEYWORD)
    }

_ALLOWED: set[str] = set()

_MODEL_TO_W = {
    "IMR90": "corigami_data/model_weights/v1_jimin.ckpt",
    "BALL":  "corigami_data/model_weights/v4_javier.ckpt",
}

def _genome_to_seq_dir(genome: str) -> str:
    return f"./corigami_data/data/{genome}/dna_sequence"


# ---------------------------------------------------------------- prepare_inputs

def _truthy(v): return str(v).lower() in ("1", "true", "on", "yes")

def prepare_inputs(form: dict, outdir: str) -> tuple[dict, dict, dict]:
    """
    Re‑implements all heavy preprocessing formerly in routes.index():  
      ▸ derive final ATAC / CTCF bigWigs (prediction + normalisation)  
      ▸ handle chimeric splice, deletion bookkeeping  
    Returns:
        run_args ........ kwargs for run_prediction_and_render()
        axis_cfg ........ custom axis dict
        gene_track_cfg .. gene‑track dict
    """
    # ───────── basic selections
    model   = form.get("model_select", "IMR90")
    genome  = form.get("genome_select", "hg38")
    model_w = _MODEL_TO_W.get(model, _MODEL_TO_W["IMR90"])
    seq_dir = _genome_to_seq_dir(genome)

    req_atac_norm = "log"    if model == "IMR90" else "minmax"
    req_ctcf_norm = "log2fc" if model == "IMR90" else "minmax"

    chr1   = form["region_chr1"]
    start1 = int(form["region_start1"])
    end1   = int(form.get("region_end1") or start1 + WINDOW_WIDTH)

    ds_opt = form.get("ds_option", "none")
    chimeric = bool(form.get("region_start2", "").strip())

    # ───────── raw files + flags
    atac_raw   = form["atac_bw_path"].strip()
    ctcf_raw   = form["ctcf_bw_path"].strip()

    apply_atac_norm = form.get("apply_atac_norm") == "true"
    apply_ctcf_norm = form.get("apply_ctcf_norm") == "true"
    predict_ctcf_f  = form.get("predict_ctcf")   == "true"

    # clean target folder
    if os.path.exists(outdir):
        for entry in os.listdir(outdir):
            p = os.path.join(outdir, entry)
            if os.path.isfile(p) or os.path.islink(p):
                os.unlink(p)
            elif os.path.isdir(p):
                shutil.rmtree(p)
    else:
        os.makedirs(outdir, exist_ok=True)

    # ════════════ ATAC & CTCF bigWigs ─────────────────────────────────────── #
    atac_pre_norm = atac_raw

    # Get chimeric params early if needed
    if chimeric:
        chr2        = form["region_chr2"]
        chr2_start  = int(form["region_start2"])
        chr2_end    = int(form.get("region_end2") or chr2_start + WINDOW_WIDTH)
        
        chr1_start  = start1
        chr1_end    = end1
        
        chr1_reverse = _truthy(form.get("first_reverse"))
        chr1_flip    = _truthy(form.get("first_flip"))
        chr2_reverse = _truthy(form.get("second_reverse"))
        chr2_flip    = _truthy(form.get("second_flip"))

    # ----- CTCF prediction if requested/needed
    if predict_ctcf_f or ctcf_raw.lower() == "none":
        if chimeric:
            # Need to predict CTCF for BOTH regions
            ctcf_pred1 = os.path.join(outdir, "predicted_ctcf_region1.bw")
            ctcf_pred2 = os.path.join(outdir, "predicted_ctcf_region2.bw")
            
            # Predict for each region
            predict_ctcf(atac_raw, chr1, chr1_start, chr1_end, ctcf_pred1)
            predict_ctcf(atac_raw, chr2, chr2_start, chr2_end, ctcf_pred2)
            
            # Mark that we have predicted files for chimeric assembly
            ctcf_pre_norm = None  # Signal for chimeric path
            ctcf_predicted_files = (ctcf_pred1, ctcf_pred2)
        else:
            # Single region prediction
            ctcf_pred = os.path.join(outdir, "predicted_ctcf.bw")
            ctcf_pre_norm = predict_ctcf(atac_raw, chr1, start1, end1, ctcf_pred)
            ctcf_predicted_files = None
    else:
        ctcf_pre_norm = ctcf_raw
        ctcf_predicted_files = None

    # ════════════ deletion or chimeric bookkeeping ───────────────────────── #
    del_start = del_width = None
    if ds_opt == "deletion":
        del_start = int(form.get("del_start", 1_500_000))
        del_width = int(form.get("del_width",   500_000))

    if chimeric:
        # Extract ATAC arrays
        a1 = extract_bigwig_region_array(atac_raw, chr1, chr1_start, chr1_end, do_reverse=chr1_reverse)
        a2 = extract_bigwig_region_array(atac_raw, chr2, chr2_start, chr2_end, do_reverse=chr2_reverse)
        
        # Extract CTCF arrays - handle both predicted and provided cases
        if ctcf_predicted_files:
            # Use the separate predicted files
            c1 = extract_bigwig_region_array(ctcf_predicted_files[0], chr1, chr1_start, chr1_end, do_reverse=chr1_reverse)
            c2 = extract_bigwig_region_array(ctcf_predicted_files[1], chr2, chr2_start, chr2_end, do_reverse=chr2_reverse)
        else:
            # Use the provided CTCF file
            c1 = extract_bigwig_region_array(ctcf_pre_norm, chr1, chr1_start, chr1_end, do_reverse=chr1_reverse)
            c2 = extract_bigwig_region_array(ctcf_pre_norm, chr2, chr2_start, chr2_end, do_reverse=chr2_reverse)

        # Write chimeric bigWigs
        atac_final, ctcf_final = write_chimeric_bigwig(
            np.concatenate([a1, a2]),
            np.concatenate([c1, c2]),
            outdir, len(a1) + len(a2), "chrCHIM",
        )

        chim_len = len(a1) + len(a2)
        
        # Apply normalization if needed
        if apply_atac_norm:
            atac_final = normalize_bigwig(
                atac_final, "chrCHIM", 0, chim_len, req_atac_norm, outdir, prefix="chrCHIM_atac"
            )
        if apply_ctcf_norm:
            ctcf_final = normalize_bigwig(
                ctcf_final, "chrCHIM", 0, chim_len, req_ctcf_norm, outdir, prefix="chrCHIM_ctcf"
            )

        # synthetic FASTA
        assemble_chimeric_fasta(
            f"{seq_dir}/{chr1}.fa.gz", chr1_start, chr1_end,
            f"{seq_dir}/{chr2}.fa.gz", chr2_start, chr2_end,
            chr1_reverse=chr1_reverse, chr1_flip=chr1_flip,
            chr2_reverse=chr2_reverse, chr2_flip=chr2_flip,
            output_folder=os.path.join(outdir, "chim_fa"),
            chim_name="chrCHIM",
        )
        seq_dir = os.path.join(outdir, "chim_fa")
        chr_final, start_final, end_final = "chrCHIM", 0, chim_len

        gene_cfg = prepare_chimeric_gene_track_config(
            "static/genes.gencode.v38.txt",
            chr1, chr1_start, chr1_end,
            chr2, chr2_start, chr2_end,
        )
        axis_cfg = {
            "region1": {"chrom": chr1, "startMb": chr1_start/1e6, "endMb": chr1_end/1e6},
            "region2": {"chrom": chr2, "startMb": chr2_start/1e6, "endMb": chr2_end/1e6},
        }

    else:
        # NON-CHIMERIC PATH
        if apply_atac_norm:
            atac_final = normalize_bigwig(
                atac_pre_norm, chr1, start1, end1, req_atac_norm, outdir
            )
        else:
            atac_final = atac_pre_norm

        if apply_ctcf_norm:
            ctcf_final = normalize_bigwig(
                ctcf_pre_norm, chr1, start1, end1, req_ctcf_norm, outdir
            )
        else:
            ctcf_final = ctcf_pre_norm

        chr_final, start_final, end_final = chr1, start1, end1
        gene_cfg = prepare_gene_track_config(
            genome, chr1, start1, end1, ds_opt, del_start, del_width
        )
        axis_cfg = {
            "region1": {"chrom": chr1, "startMb": start1/1e6, "endMb": end1/1e6},
        }
        if ds_opt == "deletion":
            axis_cfg.update({
                "deletionStartMb": del_start/1e6,
                "deletionEndMb":   (del_start+del_width)/1e6,
            })

    # ════════════ build final kwargs for worker ═════════════════════════════ #
    run_args = dict(
        region_chr        = chr_final,
        region_start      = start_final,
        region_end        = end_final,
        ds_option         = ds_opt,
        model_path        = model_w,
        genome            = genome,
        seq_dir           = seq_dir,
        atac_bw_for_model = atac_final,
        ctcf_bw_for_model = ctcf_final,
        output_folder     = outdir,
        del_start         = del_start,
        del_width         = del_width,
        screening_requested = ds_opt == "screening",
    )
    return run_args, axis_cfg, gene_cfg
# ═════════════════════════════════ routes ═════════════════════════════════ #
# ──────────────────────────────────────────────────────────
#  app/api.py
# ──────────────────────────────────────────────────────────
@api.route("/predict", methods=["POST"])
def predict():
    """
    • creates a per‑user output folder
    • enqueues the fast prediction job
    • if ds_option == 'screening' also enqueues the long screening job
      and writes that child‑job id into the parent’s meta
    """
    if not request.is_json:
        return jsonify({"error": "expecting JSON"}), 415

    form_json = request.get_json()
    outdir    = get_user_output_folder()

    q = current_app.q

    # ① fast prediction job
    pred_job = q.enqueue(
        "app.tasks.worker_predict_and_norm",
        form_json,
        outdir,
        job_timeout=3600,
        ttl         = 86400,         # keep queued/started job 24 h
        result_ttl  = 86400,         # keep finished job 24 h
    )

    # ② optional long screening job
    if form_json.get("ds_option") == "screening":
        scr_job = q.enqueue(
            "app.tasks.worker_run_screening",
            pred_job.id,            # pass parent‑id
            depends_on=pred_job,
            job_timeout=7200,
            ttl         = 86400,         # keep queued/started job 24 h
            result_ttl  = 86400,         # keep finished job 24 h
        )
        pred_job.meta["screening_job_id"] = scr_job.id
        pred_job.save_meta()

    session["current_job_id"] = pred_job.id
    return jsonify({"job_id": pred_job.id}), 202


@api.route("/job/<job_id>")
def job_status(job_id):
    """Polling endpoint."""
    job = current_app.q.fetch_job(job_id)
    if not job:
        return jsonify({"status": "unknown"}), 404
    return jsonify({
        "status": (
            "queued"    if job.is_queued
            else "started"  if job.is_started
            else "failed"   if job.is_failed
            else "finished" if job.is_finished
            else job.get_status()
        )
    })


# ────────────────────────────────────────────────────────────────────────────
#  /api/job/<id>/html   – return the fragment that the front‑end injects
# ────────────────────────────────────────────────────────────────────────────
@api.route("/job/<job_id>/html")
def job_html(job_id):
    """
    Returns the HTML fragment once the *prediction* job is finished.
    The optional screening job may still be running, so every
    attribute is fetched with .get().
    """
    job = current_app.q.fetch_job(job_id)
    if not job:
        return "Unknown job‑id", 404
    if not job.is_finished:
        return "Not ready", 425          # 425 Too Early

    m = job.meta
    return render_template(
        "plots_partial.html",
        hi_c_config       = m.get("hi_c_config"),
        ctcf_config       = m.get("ctcf_config"),
        atac_config       = m.get("atac_config"),
        gene_track_config = m.get("gene_track_config"),
        screening_config  = m.get("screening_config"),
        screening_mode    = m.get("screening_mode", False),
        custom_axis_config= m.get("custom_axis_config"),
        current_job_id    = job_id,          # ← this one is read by the JS
    )

# ─────────────────────────────────────────────────────────────────
# /api/run_screening – long-poll endpoint
# ─────────────────────────────────────────────────────────────────
@api.route("/run_screening")
def run_screening():
    parent_id = request.args.get("parent")
    if not parent_id:
        return jsonify(error="missing parent param"), 400

    parent = current_app.q.fetch_job(parent_id)
    if not parent:
        return jsonify(error="unknown parent"), 404

    # ♦ NEW: if the parent already holds a config we’re DONE
    cfg = parent.meta.get("screening_config")
    if cfg is not None:
        return jsonify(status="done", config=cfg), 200

    # screening wasn’t requested at all
    if not parent.meta.get("screening_mode"):
        return jsonify(status="disabled"), 200

    # we haven’t finished yet – fall back to child status
    child_id = parent.meta.get("screening_job_id")
    if not child_id:                          # still queued / meta lost
        return jsonify(status="running"), 202

    child = current_app.q.fetch_job(child_id)
    if child is None or not child.is_finished:
        return jsonify(status="running"), 202
    if child.is_failed:
        return jsonify(status="error",
                       message="screening worker crashed"), 500

    # child finished, but parent.meta didn’t get the config for some reason
    cfg = child.meta.get("screening_config")
    if cfg is None:
        return jsonify(status="error",
                       message="screening_config missing"), 500

    # propagate the config back to the parent for next time
    parent.meta["screening_config"] = cfg
    parent.save_meta()
    return jsonify(status="done", config=cfg), 200