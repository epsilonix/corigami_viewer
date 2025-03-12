# prediction.py
#!/usr/bin/env python
import os
import sys
import argparse
import math
import numpy as np
from skimage.transform import resize  # still imported for pil_resize if needed

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
    matrix = matrix.astype(np.float32)
    img = Image.fromarray(matrix, mode='F')
    img_resized = img.resize(new_shape, resample=Image.NEAREST)
    return np.array(img_resized, dtype=np.float32)

def consensus_prediction(chr_name, full_start, full_end, model_path, seq_path, ctcf_path, atac_path,
                         window_size=DEFAULT_WINDOW_SIZE, consensus_window_step=253952):
    """
    Build a consensus Hi-C matrix across a wider region.
    This function divides the full region into overlapping windows (of size window_size)
    with a fixed step (consensus_window_step) and for each window runs a prediction.
    It then maps each predicted matrix into a consensus matrix using the effective resolution
    of the prediction (computed as window_size / n_bins_model) and averages overlapping bins.
    """
    full_length = full_end - full_start

    # Run one prediction to determine the model's output shape and effective resolution.
    sample_pred = single_prediction(chr_name, full_start, model_path, seq_path, ctcf_path, atac_path)
    n_bins_model = sample_pred.shape[0]  # e.g. 256
    effective_resolution = window_size / n_bins_model  # effective bp per bin

    # Determine consensus matrix dimensions using effective resolution.
    num_bins_consensus = int(math.ceil(full_length / effective_resolution))
    consensus = np.zeros((num_bins_consensus, num_bins_consensus), dtype=np.float32)
    count = np.zeros((num_bins_consensus, num_bins_consensus), dtype=np.float32)

    # Loop over subwindows using the fixed consensus_window_step.
    for sub_start in range(full_start, full_end - window_size + 1, consensus_window_step):
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

    # Removed resizing step based on resolution_target.
    return consensus

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
    # Optional parameter for consensus prediction over a wider region.
    parser.add_argument('--full_end', dest='full_end', type=int, default=None,
                        help='End coordinate for full region prediction. If provided and greater than start+window_size, a consensus is built.')
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    celltype_for_path = args.celltype.strip() if args.celltype and args.celltype.strip() != "" else ""
    
 
    out_path = os.path.join(args.output_path, "result.npy")
    
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
            consensus_window_step=253952  # fixed step size
        )
        np.save(out_path, consensus_matrix)
        print(f"Consensus prediction saved to: {out_path}")
        # Plot generation commented out:
        # plot = plot_utils.MatrixPlot(args.output_path, consensus_matrix, 'consensus_prediction', celltype_for_path, args.chr_name, args.start)
        # plot.plot()
    else:
        print("Running standard 2 Mb prediction...")
        pred = single_prediction(args.chr_name, args.start, args.model_path, args.seq_path, args.ctcf_path, args.atac_path)
        np.save(out_path, pred)
        print(f"Prediction saved to: {out_path}")
        # Plot generation commented out:
        # plot = plot_utils.MatrixPlot(args.output_path, pred, 'prediction', celltype_for_path, args.chr_name, args.start)
        # plot.plot()

if __name__ == '__main__':
    main()
