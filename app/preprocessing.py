# app/api.py  – keep everything else exactly as‑is
# ----------------------------------------------------------------------
_ALLOWED = {
    "region_chr", "region_start", "region_end",
    "ds_option", "del_start", "del_width",
    "model_path", "genome", "seq_dir",
    "atac_bw_for_model", "ctcf_bw_for_model",
    "raw_atac_path", "output_folder",
    "screening_requested",
}

def preprocess(form_json):
    """Return kwargs that match *exactly* what run_prediction_and_render expects."""
    f   = form_json          # shortcut
    out = {}

    # -------- direct int fields --------
    for raw, dest in {
        "region_start1": "region_start",
        "region_end1":   "region_end",
        "del_start":     "del_start",
        "del_width":     "del_width",
    }.items():
        if f.get(raw):
            out[dest] = int(f[raw])

    # default end = start + 2 Mb
    if "region_start" in out and "region_end" not in out:
        out["region_end"] = out["region_start"] + WINDOW_WIDTH

    # -------- direct str fields --------
    for raw, dest in {
        "region_chr1":   "region_chr",
        "ds_option":     "ds_option",
        "atac_bw_path":  "raw_atac_path",
        "atac_bw_path":  "atac_bw_for_model",   # same value reused
        "ctcf_bw_path":  "ctcf_bw_for_model",
        "genome_select": "genome",
    }.items():
        if f.get(raw) not in ("", None):
            out[dest] = str(f[raw])

    # -------- model_path & required seq_dir --------
    model = f.get("model_select", "")
    if model   == "IMR90":
        out["model_path"] = "corigami_data/model_weights/v1_jimin.ckpt"
    elif model == "BALL":
        out["model_path"] = "corigami_data/model_weights/v4_javier.ckpt"
    else:
        raise ValueError("unknown model_select")

    genome = out.get("genome", "hg38")
    out["seq_dir"] = f"./corigami_data/data/{genome}/dna_sequence"

    # -------- output folder (per‑user tmp) --------
    out["output_folder"] = get_user_output_folder()

    # -------- screening flag --------
    out["screening_requested"] = (f.get("ds_option") == "screening")

    # finally: drop anything not accepted
    return {k: v for k, v in out.items() if k in _ALLOWED}
