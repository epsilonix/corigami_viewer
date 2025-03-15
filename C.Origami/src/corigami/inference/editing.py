#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import sys
import torch
import argparse
import json

import corigami.inference.utils.inference_utils as infer
from corigami.inference import editing 
#from corigami.inference.utils import plot_utils   # Removed plotting dependency

def main():
    parser = argparse.ArgumentParser(description='C.Origami Editing Module.')
    
    # Output location
    parser.add_argument('--out', dest='output_path', 
                        default='outputs',
                        help='Output path for storing results (default: %(default)s)')
    # Location related parameters
    parser.add_argument('--chr', dest='chr_name', 
                        help='Chromosome for prediction', required=True)
    parser.add_argument('--start', dest='start', type=int,
                        help='Starting point for prediction (input window width is 2097152 bp)', required=True)
    parser.add_argument('--model', dest='model_path', 
                        help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', 
                        help='Path to the folder where the sequence .fa.gz files are stored', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', 
                        help='Path to the folder where the CTCF ChIP-seq .bw files are stored', required=True)
    parser.add_argument('--atac', dest='atac_path', 
                        help='Path to the folder where the ATAC-seq .bw files are stored', required=True)
    # Deletion related parameters
    parser.add_argument('--del-start', dest='deletion_start', type=int,
                        help='Starting point for deletion', required=True)
    parser.add_argument('--del-width', dest='deletion_width', type=int,
                        help='Width for deletion', required=True)
    parser.add_argument('--padding', dest='end_padding_type', 
                        default='zero',
                        help='Padding type: "zero" (pad with zeros/N) or "follow" (pad with following region features) (default: %(default)s)')
    parser.add_argument('--hide-line', dest='hide_deletion_line', 
                        action='store_true',
                        help='Remove the line showing the deletion site')
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    single_deletion(args.output_path, args.chr_name, args.start, 
                    args.deletion_start, args.deletion_width, 
                    args.model_path, args.seq_path, args.ctcf_path, args.atac_path, 
                    show_deletion_line = not args.hide_deletion_line,
                    end_padding_type = args.end_padding_type)

def single_deletion(output_path, chr_name, start, deletion_start, deletion_width, 
                    model_path, seq_path, ctcf_path, atac_path, 
                    show_deletion_line=True, end_padding_type='zero'):
    # Define a window that accommodates the deletion.
    window = 2097152 + deletion_width
    seq_region, ctcf_region, atac_region = infer.load_region(chr_name, start, seq_path, ctcf_path, atac_path, window=window)
    # Apply deletion (with padding)
    seq_region, ctcf_region, atac_region = deletion_with_padding(start, deletion_start, deletion_width, 
                                                                  seq_region, ctcf_region, atac_region, end_padding_type)
    # Run prediction on the modified (deleted) input.
    pred = infer.prediction(seq_region, ctcf_region, atac_region, model_path)
    
    # Save the predicted Hi-C matrix as an npy file.
    npy_dir = os.path.join(output_path)
    if not os.path.exists(npy_dir):
        os.makedirs(npy_dir)
    npy_filename = "result.npy"
    npy_filepath = os.path.join(npy_dir, npy_filename)
    np.save(npy_filepath, pred)
    print(f"Saved predicted Hi-C matrix to {npy_filepath}")
    
    # Build a JSON dictionary with summary results.
    results = {
        "chr": chr_name,
        "start": start,
        "deletion_start": deletion_start,
        "deletion_width": deletion_width,
        "window": window,
        "output_matrix_shape": pred.shape,
        "npy_file": npy_filepath  # Include path to saved npy file.
    }
    
    # Save the results JSON to the output folder directly.
    results_file = os.path.join(output_path, "editing_results.json")
    with open(results_file, "w") as f:
        json.dump(results, f)
    print(f"Saved editing results to {results_file}")

def deletion_with_padding(start, deletion_start, deletion_width, seq_region, ctcf_region, atac_region, end_padding_type):
    """ Delete signals from the input region and apply padding.
        end_padding_type: 'zero' or 'follow'
    """
    if end_padding_type == 'zero':
        seq_region, ctcf_region, atac_region = zero_region(seq_region, ctcf_region, atac_region)
    elif end_padding_type == 'follow':
        # Not implemented; fallback to zero padding.
        seq_region, ctcf_region, atac_region = zero_region(seq_region, ctcf_region, atac_region)
    else:
        raise Exception('Unknown padding type: ' + end_padding_type)
    # Perform deletion: remove the region between deletion_start - start and deletion_start - start + deletion_width.
    seq_region, ctcf_region, atac_region = delete(deletion_start - start, deletion_start - start + deletion_width,
                                                    seq_region, ctcf_region, atac_region)
    return seq_region, ctcf_region, atac_region

def zero_region(seq, ctcf, atac, window=2097152):
    """ Replace features after the window with zeros (or N's for sequence). """
    seq[window:] = [0, 0, 0, 0, 1]  # Here 0 represents a replacement; adjust as needed.
    ctcf[window:] = 0
    atac[window:] = 0
    return seq, ctcf, atac

def delete(start, end, seq, ctcf, atac, window=2097152):
    seq = np.delete(seq, np.s_[start:end], axis=0)
    ctcf = np.delete(ctcf, np.s_[start:end])
    atac = np.delete(atac, np.s_[start:end])
    return seq[:window], ctcf[:window], atac[:window]

if __name__ == '__main__':
    main()
