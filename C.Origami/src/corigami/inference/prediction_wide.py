import corigami.inference.utils.inference_utils as infer
from corigami.inference.utils import plot_utils 
import argparse
import sys
import numpy as np
import os

# New parameters for large window segmentation
WINDOW_SIZE = 2097152  # Prediction window size
STEP_SIZE = 253952     # Overlap/step between windows

def main():
    parser = argparse.ArgumentParser(description='C.Origami Prediction Module with Large Window Support.')
    parser.add_argument('--out', dest='output_path', default='outputs',
                        help='Output path for storing results (default: %(default)s)')
    parser.add_argument('--celltype', dest='celltype', 
                        help='Sample cell type for prediction, used for output separation')
    parser.add_argument('--chr', dest='chr_name', 
                        help='Chromosome for prediction', required=True)
    # Instead of a single start, now require both a start and end to define a large window.
    parser.add_argument('--start', dest='start', type=int,
                        help='Start position of the large window', required=True)
    parser.add_argument('--end', dest='end', type=int,
                        help='End position of the large window', required=True)
    parser.add_argument('--model', dest='model_path', 
                        help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', 
                        help='Path to the folder where the sequence .fa.gz files are stored', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', 
                        help='Path to the folder where the CTCF ChIP-seq .bw files are stored', required=True)
    parser.add_argument('--atac', dest='atac_path', 
                        help='Path to the folder where the ATAC-seq .bw files are stored', required=True)
    parser.add_argument('--norm', dest='norm', default=None,
                        help='Normalization method for ATAC (e.g., "log" or "minmax"). (Currently not used)')

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    # Call the new function that supports a large window.
    large_window_prediction(args.output_path, args.celltype, args.chr_name, 
                              args.start, args.end,
                              args.model_path, args.seq_path, 
                              args.ctcf_path, args.atac_path)

def large_window_prediction(output_path, celltype, chr_name, large_start, large_end, 
                            model_path, seq_path, ctcf_path, atac_path):
    # Create an output directory including cell type if provided
    celltype_for_path = celltype.strip() if celltype and celltype.strip() != "" else ""
    region_output_dir = os.path.join(output_path, celltype_for_path, chr_name)
    os.makedirs(region_output_dir, exist_ok=True)
    
    # Determine segmentation: create list of window start positions from large_start to large_end.
    window_starts = list(range(large_start, large_end, STEP_SIZE))
    
    # To store predictions for each segment.
    segment_predictions = []
    segment_positions = []
    
    for start in window_starts:
        # Ensure the segment does not extend beyond the large window.
        seg_start = start
        seg_end = min(start + WINDOW_SIZE, large_end)
        
        # Load the region-specific data.
        seq_region, ctcf_region, atac_region = infer.load_region(chr_name, seg_start, seq_path, ctcf_path, atac_path)
        
        # Get prediction for the segment.
        pred = infer.prediction(seq_region, ctcf_region, atac_region, model_path)
        
        # Save prediction for later stitching.
        segment_predictions.append(pred)
        segment_positions.append(seg_start)
        
        # Optionally, plot each individual segment:
        seg_plot = plot_utils.MatrixPlot(region_output_dir, pred, 'segment', celltype_for_path, chr_name, seg_start)
        seg_plot.plot()
    
    # Stitch predictions into a consensus Hi-C matrix.
    # Here, we assume each prediction returns a 2D numpy array.
    # You could simply average overlapping regions.
    consensus_matrix = stitch_predictions(segment_predictions, segment_positions, WINDOW_SIZE, large_start, large_end)
    
    # Save and plot the consensus Hi-C map.
    consensus_file = os.path.join(region_output_dir, f"{chr_name}_{large_start}_{large_end}_consensus.npy")
    np.save(consensus_file, consensus_matrix)
    
    consensus_plot = plot_utils.MatrixPlot(region_output_dir, consensus_matrix, 'consensus', celltype_for_path, chr_name, large_start)
    consensus_plot.plot()
    
def stitch_predictions(predictions, positions, window_size, large_start, large_end):
    """
    Stitch the overlapping predictions into a single consensus matrix.
    This example uses a simple averaging in the overlapping areas.
    Assumes predictions are square matrices of fixed size.
    """
    # Determine the number of bins based on the resolution of the prediction.
    # Here we assume each prediction matrix has the same shape.
    pred_shape = predictions[0].shape[0]
    resolution = window_size // pred_shape  # approximate base pairs per pixel
    
    # Compute total number of bins in the consensus matrix.
    total_bins = (large_end - large_start) // resolution
    consensus = np.zeros((total_bins, total_bins), dtype=np.float32)
    counts = np.zeros((total_bins, total_bins), dtype=np.float32)
    
    for pred, pos in zip(predictions, positions):
        # Convert genomic coordinate (pos) to consensus matrix index.
        idx_start = (pos - large_start) // resolution
        idx_end = idx_start + pred_shape
        
        # Define slice in consensus matrix.
        consensus_slice = (slice(idx_start, idx_end), slice(idx_start, idx_end))
        
        # Add prediction and update counts.
        consensus[consensus_slice] += pred
        counts[consensus_slice] += 1
    
    # Avoid division by zero.
    counts[counts == 0] = 1
    consensus = consensus / counts
    
    return consensus

if __name__ == '__main__':
    main()
