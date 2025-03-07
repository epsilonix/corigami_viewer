#!/usr/bin/env python
import os
import sys
import argparse
import math
import numpy as np
from skimage.transform import resize  # for resizing consensus matrix
from PIL import Image

# Import your inference utilities
import corigami.inference.utils.inference_utils as infer
from corigami.inference.utils import plot_utils

# Default window size is 2 Mb (2097152 bp)
DEFAULT_WINDOW_SIZE = 2097152

def single_prediction(chr_name, start, model_path, seq_path, ctcf_path, atac_path):
    """
    Run a prediction for a single 2-Mb window.
    """
    seq_region, ctcf_region, atac_region = infer.load_region(chr_name, start, seq_path, ctcf_path, atac_path)
    pred = infer.prediction(seq_region, ctcf_region, atac_region, model_path)
    return pred

def pil_resize(matrix, new_shape):
    """
    Resize a 2D NumPy array (matrix) to new_shape using nearest-neighbor interpolation.
    """
    matrix = matrix.astype(np.float32)
    img = Image.fromarray(matrix, mode='F')
    img_resized = img.resize(new_shape, resample=Image.NEAREST)
    return np.array(img_resized, dtype=np.float32)

def consensus_prediction(chr_name, full_start, full_end, model_path, seq_path, ctcf_path, atac_path,
                         window_size=DEFAULT_WINDOW_SIZE, step_size=1000000, resolution_target=10000):
    """
    Build a consensus Hi-C matrix across a wider region.
    """
    full_length = full_end - full_start
    sample_pred = single_prediction(chr_name, full_start, model_path, seq_path, ctcf_path, atac_path)
    n_bins_model = sample_pred.shape[0]
    effective_resolution = window_size / n_bins_model
    num_bins_consensus = int(math.ceil(full_length / effective_resolution))
    consensus = np.zeros((num_bins_consensus, num_bins_consensus), dtype=np.float32)
    count = np.zeros((num_bins_consensus, num_bins_consensus), dtype=np.float32)
    
    for sub_start in range(full_start, full_end - window_size + 1, step_size):
        pred = single_prediction(chr_name, sub_start, model_path, seq_path, ctcf_path, atac_path)
        start_bin = int((sub_start - full_start) / effective_resolution)
        end_bin = start_bin + n_bins_model
        if end_bin > num_bins_consensus:
            excess = end_bin - num_bins_consensus
            pred = pred[:-excess, :-excess]
            end_bin = num_bins_consensus
        consensus[start_bin:end_bin, start_bin:end_bin] += pred
        count[start_bin:end_bin, start_bin:end_bin] += 1

    valid = count > 0
    consensus[valid] /= count[valid]
    num_bins_target = int(math.ceil(full_length / resolution_target))
    consensus_resized = pil_resize(consensus, (num_bins_target, num_bins_target))
    return consensus_resized

def combine_regions(args):
    """
    Load two regions, apply reverse if specified, and combine them.
    This implementation assumes that the sequence, CTCF, and ATAC signals are 1D arrays.
    When the "flip" flag is selected, only the DNA sequence is flipped (i.e. reversed),
    while the ATAC and CTCF signals remain unchanged.
    """
    # Load first region (all as 1D arrays)
    seq1, ctcf1, atac1 = infer.load_region(args.chr_name, args.start, args.seq_path, args.ctcf_path, args.atac_path)
    if args.first_reverse:
        seq1 = np.flip(seq1)
        ctcf1 = np.flip(ctcf1)
        atac1 = np.flip(atac1)
    if args.first_flip:
        seq1 = seq1[::-1]
    
    # Load second region (all as 1D arrays)
    seq2, ctcf2, atac2 = infer.load_region(args.chr2, args.start2, args.seq_path, args.ctcf_path, args.atac_path)
    if args.second_reverse:
        seq2 = np.flip(seq2)
        ctcf2 = np.flip(ctcf2)
        atac2 = np.flip(atac2)
    if args.second_flip:
        seq2 = seq2[::-1]
    
    # For 1D data, concatenate end-to-end (axis=0)
    combined_seq = np.concatenate([seq1, seq2], axis=0)
    combined_ctcf = np.concatenate([ctcf1, ctcf2], axis=0)
    combined_atac = np.concatenate([atac1, atac2], axis=0)
    
    print(f"First region: {args.chr_name} from {args.start} to {args.start + DEFAULT_WINDOW_SIZE}")
    print(f"Second region: {args.chr2} from {args.start2} to {args.end2}")
    
    # Run prediction on the combined input.
    pred = infer.prediction(combined_seq, combined_ctcf, combined_atac, args.model_path)
    return pred

def main():
    parser = argparse.ArgumentParser(description='C.Origami Prediction Module')
    parser.add_argument('--out', dest='output_path', default='outputs',
                        help='Output path for storing results (default: %(default)s)')
    parser.add_argument('--celltype', dest='celltype', help='Sample cell type for prediction')
    parser.add_argument('--chr', dest='chr_name', help='Chromosome for prediction', required=True)
    parser.add_argument('--start', dest='start', type=int,
                        help='Starting coordinate for prediction (left boundary)', required=True)
    parser.add_argument('--model', dest='model_path', help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', help='Path to the folder with sequence .fa.gz files', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', help='Path to the CTCF .bw file/folder', required=True)
    parser.add_argument('--atac', dest='atac_path', help='Path to the ATAC .bw file/folder', required=True)
    # Optional parameters for consensus prediction
    parser.add_argument('--full_end', dest='full_end', type=int, default=None,
                        help='End coordinate for full region prediction; if > start+window_size, consensus is built.')
    parser.add_argument('--step_size', dest='step_size', type=int, default=253952,
                        help='Step size between windows (default: 253952)')
    parser.add_argument('--resolution_target', dest='resolution_target', type=int, default=10000,
                        help='Target bin size (bp) for consensus matrix (default: 10000)')
    
    # Optional parameters for combining two regions
    parser.add_argument('--chr2', dest='chr2', help='Second chromosome for prediction', default=None)
    parser.add_argument('--start2', dest='start2', type=int, help='Starting coordinate for second region', default=None)
    parser.add_argument('--end2', dest='end2', type=int, help='End coordinate for second region', default=None)
    parser.add_argument('--first_reverse', action='store_true', help='Reverse first chromosome region', default=False)
    parser.add_argument('--first_flip', action='store_true', help='Flip first chromosome region', default=False)
    parser.add_argument('--second_reverse', action='store_true', help='Reverse second chromosome region', default=False)
    parser.add_argument('--second_flip', action='store_true', help='Flip second chromosome region', default=False)
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    celltype_for_path = args.celltype.strip() if args.celltype and args.celltype.strip() != "" else ""
    
    # Determine which branch to run:
    if args.chr2 and args.start2 is not None and args.end2 is not None:
        print("Running combined prediction on two regions.")
        print(f"First region: {args.chr_name} from {args.start} to {args.start + DEFAULT_WINDOW_SIZE}")
        print(f"Second region: {args.chr2} from {args.start2} to {args.end2}")
        pred = combine_regions(args)
        out_fname = f"{args.chr_name}_{args.start}_{args.chr2}_{args.start2}_combined.npy"
    elif args.full_end is not None and args.full_end > args.start + DEFAULT_WINDOW_SIZE:
        print("Running consensus (wider window) prediction.")
        print(f"Region: {args.chr_name} from {args.start} to {args.full_end}")
        consensus_matrix = consensus_prediction(
            chr_name=args.chr_name,
            full_start=args.start,
            full_end=args.full_end,
            model_path=args.model_path,
            seq_path=args.seq_path,
            ctcf_path=args.ctcf_path,
            atac_path=args.atac_path,
            window_size=DEFAULT_WINDOW_SIZE,
            step_size=args.step_size,
            resolution_target=args.resolution_target
        )
        out_fname = f"{args.chr_name}_{args.start}_{args.full_end}_consensus.npy"
        os.makedirs(args.output_path, exist_ok=True)
        consensus_path = os.path.join(args.output_path, out_fname)
        np.save(consensus_path, consensus_matrix)
        print(f"Consensus prediction saved to: {consensus_path}")
        plot = plot_utils.MatrixPlot(args.output_path, consensus_matrix, 'consensus_prediction', celltype_for_path, args.chr_name, args.start)
        plot.plot()
        return
    else:
        print("Running standard 2 Mb prediction.")
        print(f"Region: {args.chr_name} from {args.start} to {args.start + DEFAULT_WINDOW_SIZE}")
        pred = single_prediction(args.chr_name, args.start, args.model_path, args.seq_path, args.ctcf_path, args.atac_path)
        out_fname = f"{args.chr_name}_{args.start}.npy"
    
    os.makedirs(args.output_path, exist_ok=True)
    pred_path = os.path.join(args.output_path, out_fname)
    np.save(pred_path, pred)
    print(f"Prediction saved to: {pred_path}")
    plot = plot_utils.MatrixPlot(args.output_path, pred, 'prediction', celltype_for_path, args.chr_name, args.start)
    plot.plot()

if __name__ == '__main__':
    main()
