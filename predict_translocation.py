#!/usr/bin/env python3
# assemble_chim_and_predict.py

from pathlib import Path
from datetime import datetime
import argparse
import os
import gzip
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from app.utils import (
    assemble_chimeric_fasta,
    extract_bigwig_region_array,
    write_chimeric_bigwig,
    normalize_bigwig,
)
import subprocess

MODEL_WINDOW = 2_097_152


def pad_fasta_to_length(fasta_path: str, target_len: int) -> int:
    with gzip.open(fasta_path, "rt") as fh:
        lines = [line.rstrip() for line in fh]

    if not lines:
        raise RuntimeError(f"Empty FASTA: {fasta_path}")

    header, seq_lines = lines[0], lines[1:]
    sequence = "".join(seq_lines)

    if len(sequence) < target_len:
        sequence = sequence + ("N" * (target_len - len(sequence)))

    with gzip.open(fasta_path, "wt") as fh:
        fh.write(header + "\n")
        for i in range(0, len(sequence), 70):
            fh.write(sequence[i:i+70] + "\n")

    return len(sequence)

def parse_args():
    script_dir = Path(__file__).resolve().parent
    default_seq_dir = script_dir / "corigami_data/data/hg38/dna_sequence"
    default_atac_bw = script_dir / "corigami_data/data/hg38/ball/MCG023_ATAC_merged_IS_slop20_RP20M_minmax01.bw"
    default_ctcf_bw = script_dir / "corigami_data/data/hg38/ball/MCG023_CTCF_merged_maxatac_predict.bw"
    default_model = script_dir / "corigami_data/model_weights/v1_jimin.ckpt"
    default_output = script_dir / "predict_translocation_output"

    parser = argparse.ArgumentParser(description="Assemble a chimeric window and run C.Origami prediction.")
    parser.add_argument("--chr1", default="chr14", help="Chromosome for the first segment")
    parser.add_argument("--start1", type=int, default=190_000_000, help="Start coordinate for the first segment")
    parser.add_argument("--end1", type=int, default=192_000_000, help="End coordinate for the first segment")
    parser.add_argument("--chr2", default="chr18", help="Chromosome for the second segment")
    parser.add_argument("--start2", type=int, default=60_000_000, help="Start coordinate for the second segment")
    parser.add_argument("--end2", type=int, default=62_000_000, help="End coordinate for the second segment")
    parser.add_argument("--chr1-reverse", action="store_true", help="Reverse the first segment before concatenation")
    parser.add_argument("--chr1-flip", action="store_true", help="Reverse-complement the first segment (sequence only)")
    parser.add_argument("--chr2-reverse", action="store_true", help="Reverse the second segment before concatenation")
    parser.add_argument("--chr2-flip", action="store_true", help="Reverse-complement the second segment (sequence only)")
    parser.add_argument("--seq-dir", default=str(default_seq_dir), help="Directory containing per-chromosome FASTA (.fa.gz) files")
    parser.add_argument("--atac-bw", default=str(default_atac_bw), help="ATAC bigWig path")
    parser.add_argument("--ctcf-bw", default=str(default_ctcf_bw), help="CTCF bigWig path")
    parser.add_argument("--model", default=str(default_model), help="Model checkpoint (.ckpt) path")
    parser.add_argument("--output-base", default=str(default_output), help="Base directory to store per-run outputs")
    parser.add_argument("--atac-norm", default="log", choices=["none", "log", "minmax", "log2fc"], help="Normalization to apply to the stitched ATAC bigWig")
    parser.add_argument("--ctcf-norm", default="none", choices=["none", "log", "minmax", "log2fc"], help="Normalization to apply to the stitched CTCF bigWig")
    return parser.parse_args()


def main():
    args = parse_args()

    script_dir = Path(__file__).resolve().parent
    seq_dir = Path(args.seq_dir).expanduser().resolve()
    atac_bw = Path(args.atac_bw).expanduser().resolve()
    ctcf_bw = Path(args.ctcf_bw).expanduser().resolve()
    model_ckpt = Path(args.model).expanduser().resolve()
    output_base = Path(args.output_base).expanduser().resolve()
    output_base.mkdir(parents=True, exist_ok=True)

    run_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = output_base / run_stamp
    inputs_dir = run_dir / "inputs"
    seq_output_dir = inputs_dir / "seq"
    outputs_dir = run_dir / "outputs"

    inputs_dir.mkdir(parents=True, exist_ok=True)
    seq_output_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    chim_fa = assemble_chimeric_fasta(
        str(seq_dir / f"{args.chr1}.fa.gz"), args.start1, args.end1,
        str(seq_dir / f"{args.chr2}.fa.gz"), args.start2, args.end2,
        output_folder=seq_output_dir,
        chim_name="chrCHIM",
        chr1_reverse=args.chr1_reverse,
        chr1_flip=args.chr1_flip,
        chr2_reverse=args.chr2_reverse,
        chr2_flip=args.chr2_flip,
    )

    reverse1 = args.chr1_reverse or args.chr1_flip
    reverse2 = args.chr2_reverse or args.chr2_flip

    a1 = extract_bigwig_region_array(str(atac_bw), args.chr1, args.start1, args.end1, do_reverse=reverse1)
    a2 = extract_bigwig_region_array(str(atac_bw), args.chr2, args.start2, args.end2, do_reverse=reverse2)
    c1 = extract_bigwig_region_array(str(ctcf_bw), args.chr1, args.start1, args.end1, do_reverse=reverse1)
    c2 = extract_bigwig_region_array(str(ctcf_bw), args.chr2, args.start2, args.end2, do_reverse=reverse2)

    chim_atac = np.concatenate([a1, a2])
    chim_ctcf = np.concatenate([c1, c2])

    if chim_atac.shape[0] < MODEL_WINDOW:
        pad_len = MODEL_WINDOW - chim_atac.shape[0]
        chim_atac = np.pad(chim_atac, (0, pad_len), mode="constant")
    if chim_ctcf.shape[0] < MODEL_WINDOW:
        pad_len = MODEL_WINDOW - chim_ctcf.shape[0]
        chim_ctcf = np.pad(chim_ctcf, (0, pad_len), mode="constant")

    chim_len = max(chim_atac.shape[0], chim_ctcf.shape[0])
    chim_len = pad_fasta_to_length(chim_fa, max(chim_len, MODEL_WINDOW))

    chim_atac_bw, chim_ctcf_bw = write_chimeric_bigwig(
        atac_array=chim_atac,
        ctcf_array=chim_ctcf,
        outdir=str(inputs_dir),
        chim_len=chim_len,
        chim_name="chrCHIM",
    )

    atac_norm_method = args.atac_norm.lower()
    ctcf_norm_method = args.ctcf_norm.lower()

    final_atac_bw = chim_atac_bw
    if atac_norm_method != "none":
        final_atac_bw = normalize_bigwig(
            chim_atac_bw,
            "chrCHIM",
            0,
            chim_len,
            atac_norm_method,
            str(inputs_dir),
            prefix="chrCHIM_atac",
        )

    final_ctcf_bw = chim_ctcf_bw
    if ctcf_norm_method != "none":
        final_ctcf_bw = normalize_bigwig(
            chim_ctcf_bw,
            "chrCHIM",
            0,
            chim_len,
            ctcf_norm_method,
            str(inputs_dir),
            prefix="chrCHIM_ctcf",
        )

    metadata = {
        "run_timestamp": run_stamp,
        "run_dir": str(run_dir),
        "parameters": {
            "chr1": {
                "chrom": args.chr1,
                "start": args.start1,
                "end": args.end1,
                "reverse": args.chr1_reverse,
                "flip": args.chr1_flip,
            },
            "chr2": {
                "chrom": args.chr2,
                "start": args.start2,
                "end": args.end2,
                "reverse": args.chr2_reverse,
                "flip": args.chr2_flip,
            },
        },
        "normalization": {
            "atac": atac_norm_method,
            "ctcf": ctcf_norm_method,
        },
        "source_files": {
            "seq_dir": str(seq_dir),
            "atac_bw": str(atac_bw),
            "ctcf_bw": str(ctcf_bw),
            "model": str(model_ckpt),
        },
        "model_window_bp": MODEL_WINDOW,
        "prediction_span_bp": chim_len,
        "inputs": {
            "folder": str(inputs_dir),
            "fasta": chim_fa,
            "atac_bw": final_atac_bw,
            "ctcf_bw": final_ctcf_bw,
            "pre_normalization": {
                "atac_bw": chim_atac_bw,
                "ctcf_bw": chim_ctcf_bw,
            },
        },
    }

    cmd = [
        "python",
        str(script_dir / "C.Origami/src/corigami/inference/prediction.py"),
        "--chr", "chrCHIM",
        "--start", "0",
        "--end", str(chim_len),
        "--model", str(model_ckpt),
        "--seq", str(seq_output_dir),
        "--ctcf", str(final_ctcf_bw),
        "--atac", str(final_atac_bw),
        "--out", str(outputs_dir),
    ]

    env = dict(os.environ)
    env["PYTHONPATH"] = str(script_dir / "C.Origami/src")

    result_path = outputs_dir / "result.npy"
    png_path = outputs_dir / "result.png"

    error_exc = None
    error_info = None

    try:
        subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError as exc:
        error_exc = exc
        error_info = {
            "returncode": exc.returncode,
            "cmd": exc.cmd,
        }

    png_path_str = None
    if result_path.exists():
        matrix = np.load(result_path)
        fig, ax = plt.subplots(figsize=(6, 6))
        im = ax.imshow(matrix, cmap="Reds", origin="lower")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title("Predicted Hi-C (chrCHIM)")
        fig.tight_layout()
        fig.savefig(png_path, dpi=200)
        plt.close(fig)
        png_path_str = str(png_path)

    metadata["outputs"] = {
        "folder": str(outputs_dir),
        "result_npy": str(result_path) if result_path.exists() else None,
        "result_png": png_path_str,
    }

    if error_info:
        metadata["error"] = error_info

    with open(run_dir / "metadata.json", "w") as fh:
        json.dump(metadata, fh, indent=2)

    if error_exc:
        print(f"Prediction failed. See {run_dir / 'metadata.json'} for details.")
        raise error_exc

    print(f"Translocation run outputs stored in {run_dir}")
    if png_path_str:
        print(f"Saved heatmap to {png_path_str}")


if __name__ == "__main__":
    main()