#!/usr/bin/env python
import os
import sys
import argparse
import math
import numpy as np
from skimage.transform import resize  # for resizing consensus matrix

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

from PIL import Image
import numpy as np

def pil_resize(matrix, new_shape):
    """
    Resize a 2D NumPy array (matrix) to new_shape using nearest-neighbor interpolation.
    """
    # Ensure matrix is in float32 mode; for 'F' mode in PIL it needs to be a float array.
    matrix = matrix.astype(np.float32)
    # Create an image from the array (using mode 'F' for 32-bit floating point)
    img = Image.fromarray(matrix, mode='F')
    # Resize the image with nearest-neighbor resampling
    img_resized = img.resize(new_shape, resample=Image.NEAREST)
    # Convert the resized image back to a NumPy array
    return np.array(img_resized, dtype=np.float32)


def consensus_prediction(chr_name, full_start, full_end, model_path, seq_path, ctcf_path, atac_path,
                         window_size=2097152, step_size=1000000, resolution_target=10000):
    """
    Build a consensus Hi-C matrix across a wider region.
    This function divides the full region into overlapping windows (of size window_size)
    and for each window runs a prediction. It then maps each predicted matrix into a consensus
    matrix using the effective resolution of the prediction (computed as window_size / n_bins_model)
    and averages overlapping bins.
    """
    full_length = full_end - full_start

    # Run one prediction to determine the model's output shape and effective resolution.
    sample_pred = single_prediction(chr_name, full_start, model_path, seq_path, ctcf_path, atac_path)
    n_bins_model = sample_pred.shape[0]  # e.g. 256
    effective_resolution = window_size / n_bins_model  # effective bp per bin (e.g. ~8192 bp)

    # Determine consensus matrix dimensions using effective resolution.
    num_bins_consensus = int(math.ceil(full_length / effective_resolution))
    consensus = np.zeros((num_bins_consensus, num_bins_consensus), dtype=np.float32)
    count = np.zeros((num_bins_consensus, num_bins_consensus), dtype=np.float32)

    # Loop over subwindows in the full region.
    for sub_start in range(full_start, full_end - window_size + 1, step_size):
        pred = single_prediction(chr_name, sub_start, model_path, seq_path, ctcf_path, atac_path)
        # Use the known model output size
        # Map genomic coordinates to consensus matrix indices:
        start_bin = int((sub_start - full_start) / effective_resolution)
        end_bin = start_bin + n_bins_model
        # If the submatrix goes beyond the consensus dimensions, truncate it.
        if end_bin > num_bins_consensus:
            excess = end_bin - num_bins_consensus
            pred = pred[:-excess, :-excess]
            end_bin = num_bins_consensus
        consensus[start_bin:end_bin, start_bin:end_bin] += pred
        count[start_bin:end_bin, start_bin:end_bin] += 1

    # Average overlapping predictions.
    valid = count > 0
    consensus[valid] /= count[valid]

    # Optionally, resize consensus to a target resolution if desired.
    num_bins_target = int(math.ceil(full_length / resolution_target))
    # Here we use PIL for resizing (or skimage if installed)
    consensus_resized = pil_resize(consensus, (num_bins_target, num_bins_target))
    
    return consensus_resized

def main():
    parser = argparse.ArgumentParser(description='C.Origami Prediction Module with Consensus for Wider Regions')
    parser.add_argument('--out', dest='output_path', default='outputs',
                        help='Output path for storing results (default: %(default)s)')
    parser.add_argument('--celltype', dest='celltype', help='Sample cell type for prediction')
    parser.add_argument('--chr', dest='chr_name', help='Chromosome for prediction', required=True)
    parser.add_argument('--start', dest='start', type=int,
                        help='Starting coordinate for prediction (used as left boundary for window(s))', required=True)
    parser.add_argument('--model', dest='model_path', help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', help='Path to the folder where the sequence .fa.gz files are stored', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', help='Path to the CTCF .bw file/folder', required=True)
    parser.add_argument('--atac', dest='atac_path', help='Path to the ATAC .bw file/folder', required=True)
    # New optional parameters for consensus prediction over a wider region
    parser.add_argument('--full_end', dest='full_end', type=int, default=None,
                        help='End coordinate for full region prediction. If provided and greater than start+window_size, a consensus is built.')
    parser.add_argument('--step_size', dest='step_size', type=int, default=253952,
                        help='Step size (in bp) between consecutive windows for consensus prediction (default: 1000000)')
    parser.add_argument('--resolution_pred', dest='resolution_pred', type=int, default=10000,
                        help='Bin size (bp) used in individual prediction matrices (default: 10000)')
    parser.add_argument('--resolution_target', dest='resolution_target', type=int, default=10000,
                        help='Target bin size (bp) for consensus matrix (default: 10000)')
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    celltype_for_path = args.celltype.strip() if args.celltype and args.celltype.strip() != "" else ""
    
    # Check if we need to run a consensus prediction (wider than 2 Mb)
    if args.full_end is not None and args.full_end > args.start + DEFAULT_WINDOW_SIZE:
        print("Running consensus (wider window) prediction...")
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

        # Save the consensus matrix
        out_fname = f"{args.chr_name}_{args.start}_{args.full_end}_consensus.npy"
        consensus_path = os.path.join(args.output_path, out_fname)
        np.save(consensus_path, consensus_matrix)
        print(f"Consensus prediction saved to: {consensus_path}")
        # Optionally, generate a plot using your plotting utilities:
        plot = plot_utils.MatrixPlot(args.output_path, consensus_matrix, 'consensus_prediction', celltype_for_path, args.chr_name, args.start)
        plot.plot()
    else:
        # Default single-window prediction (2 Mb)
        print("Running standard 2 Mb prediction...")
        pred = single_prediction(args.chr_name, args.start, args.model_path, args.seq_path, args.ctcf_path, args.atac_path)
        out_fname = f"{args.chr_name}_{args.start}.npy"
        pred_path = os.path.join(args.output_path, out_fname)
        np.save(pred_path, pred)
        print(f"Prediction saved to: {pred_path}")
        plot = plot_utils.MatrixPlot(args.output_path, pred, 'prediction', celltype_for_path, args.chr_name, args.start)
        plot.plot()

if __name__ == '__main__':
    main()
