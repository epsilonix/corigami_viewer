#!/usr/bin/env python
import os
import sys
sys.path.insert(0, "/Users/everett/corigami_viewer/C.Origami/src")
import argparse
import json
import numpy as np
import pandas as pd
from tqdm import tqdm

import corigami.inference.utils.inference_utils as infer
from corigami.inference import editing 
from corigami.inference.utils import plot_utils, model_utils

def main():
    parser = argparse.ArgumentParser(description='C.Origami Screening Module.')
    
    # Output location
    parser.add_argument('--out', dest='output_path', 
                        default='outputs',
                        help='Output path for storing results (default: %(default)s)')

    # Location related parameters (celltype is now optional)
    parser.add_argument('--celltype', dest='celltype', 
                        help='Sample cell type for prediction, used for output separation (optional)')
    parser.add_argument('--chr', dest='chr_name', 
                        help='Chromosome for prediction', required=True)
    parser.add_argument('--model', dest='model_path', 
                        help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', 
                        help='Path to the folder where the sequence .fa.gz files are stored', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', 
                        help='Path to the folder where the CTCF ChIP-seq .bw files are stored', required=True)
    parser.add_argument('--atac', dest='atac_path', 
                        help='Path to the folder where the ATAC-seq .bw files are stored', required=True)
    
    # Optional parameters for find_peaks (fallback)
    parser.add_argument('--use-find-peaks', dest='use_find_peaks', action='store_true',
                        help='Use scipy.signal.find_peaks to filter windows based on ATAC signal.')
    parser.add_argument('--peak-height', dest='peak_height', type=float,
                        help='Minimum height for a peak to be considered.')
    parser.add_argument('--peak-distance', dest='peak_distance', type=int, default=50,
                        help='Minimum distance (in bins) between peaks (default: %(default)s).')
    parser.add_argument('--peak-prominence', dest='peak_prominence', type=float,
                        help='Minimum prominence for peaks (optional).')

    # New: Peaks file argument
    parser.add_argument('--peaks-file', dest='peaks_file', type=str, default=None,
                        help='Path to a narrowPeak file containing peaks coordinates.')
    
    # Screening parameters
    parser.add_argument('--screen-start', dest='screen_start', type=int,
                        help='Starting point for screening.', required=True)
    parser.add_argument('--screen-end', dest='screen_end', type=int,
                        help='Ending point for screening.', required=True)
    parser.add_argument('--perturb-width', dest='perturb_width', type=int,
                        help='Width of perturbation used for screening.', required=True)
    parser.add_argument('--step-size', dest='step_size', type=int,
                        help='Step size of perturbations in screening.', required=True)

    # Saving related parameters
    parser.add_argument('--plot-impact-score', dest='plot_impact_score', action='store_true',
                        help='Plot impact score and save png.')
    parser.add_argument('--save-pred', dest='save_pred', action='store_true',
                        help='Save prediction tensor.')
    parser.add_argument('--save-perturbation', dest='save_perturbation', action='store_true',
                        help='Save perturbed tensor.')
    parser.add_argument('--save-diff', dest='save_diff', action='store_true',
                        help='Save difference tensor.')
    parser.add_argument('--save-impact-score', dest='save_impact_score', action='store_true',
                        help='Save impact score array.')
    parser.add_argument('--save-bedgraph', dest='save_bedgraph', action='store_true',
                        help='Save bedgraph file for impact score.')
    parser.add_argument('--save-frames', dest='plot_frames', action='store_true',
                        help='Save each deletion instance with png and npy (not recommended for large scale screening).')

    # Perturbation parameter
    parser.add_argument('--padding', dest='end_padding_type', default='zero',
                        help='Padding type, either zero or follow. (default: %(default)s)')

    # Flag to run in batch mode.
    parser.add_argument('--no-server', dest='no_server', action='store_true',
                        help='Do not run the Flask server; run screening in batch mode.')

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    if args.no_server:
        screening(args.output_path, args.celltype, args.chr_name, args.screen_start, 
                  args.screen_end, args.perturb_width, args.step_size,
                  args.model_path, args.seq_path, args.ctcf_path, args.atac_path,
                  args.use_find_peaks, args.peak_height, args.peak_distance, args.peak_prominence,
                  args.save_pred, args.save_perturbation, args.save_diff,
                  args.save_impact_score, args.save_bedgraph, args.plot_impact_score,
                  args.plot_frames, args.peaks_file)
    else:
        from flask import Flask
        app = Flask(__name__)
        @app.route('/')
        def index():
            return "C.Origami Screening Module"
        app.run(debug=True, use_reloader=False)

def screening(output_path, celltype, chr_name, screen_start, screen_end, perturb_width, step_size,
              model_path, seq_path, ctcf_path, atac_path,
              use_find_peaks, peak_height, peak_distance, peak_prominence,
              save_pred=False, save_perturbation=False, save_diff=True, save_impact_score=True,
              save_bedgraph=True, plot_impact_score=True, plot_frames=False, peaks_file=None):
    # Load data and model.
    print(f'screening: {screen_start}, {screen_end}, {perturb_width}, {step_size}', flush=True)
    print(f'screening peaks_file: {peaks_file}', flush=True)
    seq, ctcf, atac = infer.load_data_default(chr_name, seq_path, ctcf_path, atac_path)
    model = model_utils.load_default(model_path)
    
    # Determine screening windows. (Same logic you already have.)
    if peaks_file:
        # If the peaks file comes from auto-generation, skip one header row.
        skip_rows = 1 if "auto_peaks.narrowPeak" in peaks_file else 0
        peaks_df = pd.read_csv(peaks_file, sep='\t', header=None, skiprows=skip_rows)
        # Filter to the selected chromosome.
        peaks_df = peaks_df[peaks_df[0] == chr_name]
        # Compute the midpoint for each peak (using columns 1 and 2)
        peaks_df['mid'] = peaks_df[[1, 2]].mean(axis=1).astype(int)
        # Filter peaks within [screen_start, screen_end]
        peaks_df = peaks_df[(peaks_df['mid'] >= screen_start) & (peaks_df['mid'] <= screen_end)]
        
        # Determine number of peaks to select.
        num_peaks = int((screen_end - screen_start) // 400000)
        print("Selecting top {} peaks by score from {} available peaks".format(num_peaks, len(peaks_df)))
        
        # Sort peaks by score (column index 4) in descending order and take the top num_peaks.
        peaks_df = peaks_df.sort_values(by=4, ascending=False).head(num_peaks)
        windows = peaks_df['mid'].tolist()
        print("Total windows based on peaks file (after selection): {}".format(len(windows)))
    else:
        windows = [w * step_size + screen_start for w in range(int((screen_end - screen_start) / step_size))]
        print("Total windows (uniform): {}".format(len(windows)))

    # Prepare arrays for storing results.
    preds = np.empty((0, 256, 256))
    preds_deletion = np.empty((0, 256, 256))
    diff_maps = np.empty((0, 256, 256))
    perturb_starts = []
    perturb_ends = []

    processed_windows = 0
    skipped_windows = 0
    print('Screening...')

    # The total width the model needs is 2,097,152 plus the deletion width.
    # (Because `predict_difference` calls `end = start + 2_097_152 + deletion_width`.)
    # We'll call that "MODEL_WINDOW_SIZE".
    MODEL_BASE_SIZE = 2_097_152

    for center in tqdm(windows):
        # Define the perturbation window around 'center'.
        w_start = center - (perturb_width // 2)
        window_end = w_start + perturb_width

        # Adjust boundaries if necessary (optional).
        if w_start < screen_start:
            w_start = screen_startif
            window_end = w_start + perturb_width
        if window_end > screen_end:
            w_start = screen_end - perturb_width
            window_end = screen_end

        # The model uses `pred_start` as the "anchor" for the 2M window:
        pred_start = int(w_start + perturb_width / 2 - MODEL_BASE_SIZE / 2)

        # *** Here we skip partial windows if [pred_start, pred_start + 2M + deletion_width] is out of range. ***

        region_start = pred_start
        region_end = pred_start + MODEL_BASE_SIZE + perturb_width  # same as 'end' in predict_difference

        # If region_start or region_end is outside [screen_start, screen_end], skip
        if region_start < screen_start or region_end > screen_end:
            skipped_windows += 1
            continue  # skip this window entirely

        # Otherwise, we can safely process it
        processed_windows += 1
        pred, pred_deletion, diff_map = predict_difference(
            chr_name, pred_start, int(w_start), perturb_width,
            model, seq, ctcf, atac,
            screen_start, screen_end
        )

        # Optional frame plotting
        if plot_frames:
            plot_combination(
                output_path, celltype, chr_name,
                pred_start, w_start, perturb_width,
                pred, pred_deletion, diff_map, 'screening'
            )

        # Append to arrays
        preds = np.append(preds, np.expand_dims(pred, 0), axis=0)
        preds_deletion = np.append(preds_deletion, np.expand_dims(pred_deletion, 0), axis=0)
        diff_maps = np.append(diff_maps, np.expand_dims(diff_map, 0), axis=0)
        perturb_starts.append(w_start)
        perturb_ends.append(window_end)

    print("Total windows processed: {}".format(processed_windows))
    print("Total windows skipped (partial): {}".format(skipped_windows))

    # ... Your existing code to save JSON, etc. ...
    window_midpoints = (np.array(perturb_starts) + np.array(perturb_ends)) / 2
    window_midpoints_mb = (window_midpoints / 1e6).tolist()
    results = {
        "window_midpoints_mb": window_midpoints_mb,
        "impact_scores": np.abs(diff_maps).mean(axis=(1,2)).tolist(),
        "screen_start_mb": screen_start / 1e6,
        "screen_end_mb": screen_end / 1e6,
        "perturb_width": perturb_width,
        "step_size": step_size
    }
    
    # Build the results directory, etc.
    if celltype and celltype.strip() != "":
        results_dir = os.path.join(output_path, celltype, "screening")
    else:
        results_dir = os.path.join(output_path, "screening")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    results_file = os.path.join(results_dir, "screening_results.json")
    with open(results_file, "w") as f:
        json.dump(results, f)
    print("Saved screening results to {}".format(results_file))

def predict_difference(chr_name, start, deletion_start, deletion_width, model, seq, ctcf, atac, screen_start, screen_end):
    
    if chr_name == "chrCHIM":
        chrom_len = screen_end - screen_start

        # Clip start and end to avoid out-of-bounds on chrCHIM
        if start < 0:
            start = 0
        end = start + 2_097_152 + deletion_width
        if end > chrom_len:
            end = chrom_len
    else:
        # For any other chromosome, do whatever logic you want
        end = start + 2_097_152 + deletion_width
   
   

    seq_region, ctcf_region, atac_region = infer.get_data_at_interval(chr_name, start, end, seq, ctcf, atac)
    inputs = preprocess_prediction(chr_name, start, seq_region, ctcf_region, atac_region)
    pred = model(inputs)[0].detach().cpu().numpy()
    inputs_deletion = preprocess_deletion(chr_name, start, deletion_start, deletion_width, seq_region, ctcf_region, atac_region)
    pred_deletion = model(inputs_deletion)[0].detach().cpu().numpy()
    diff_map = pred_deletion - pred
    return pred, pred_deletion, diff_map

def plot_combination(output_path, celltype, chr_name, start, deletion_start, deletion_width, pred, pred_deletion, diff_map, plot_type='point_screening'):
    plot = plot_utils.MatrixPlot(output_path, pred, plot_type, celltype, chr_name, start)
    plot.plot()
    plot = plot_utils.MatrixPlotDeletion(output_path, pred_deletion, plot_type, celltype, chr_name, start, deletion_start, deletion_width, padding_type='zero', show_deletion_line=True)
    plot.plot()
    plot = plot_utils.MatrixPlotPointScreen(output_path, diff_map, plot_type, celltype, chr_name, start, deletion_start, deletion_width, padding_type='zero', show_deletion_line=False)
    plot.plot()

def preprocess_prediction(chr_name, start, seq_region, ctcf_region, atac_region):
    seq_region, ctcf_region, atac_region = trim(seq_region, ctcf_region, atac_region)
    inputs = infer.preprocess_default(seq_region, ctcf_region, atac_region)
    return inputs

def preprocess_deletion(chr_name, start, deletion_start, deletion_width, seq_region, ctcf_region, atac_region):
    seq_region, ctcf_region, atac_region = editing.deletion_with_padding(start, deletion_start, deletion_width, seq_region, ctcf_region, atac_region, end_padding_type='zero')
    inputs = infer.preprocess_default(seq_region, ctcf_region, atac_region)
    return inputs

def trim(seq, ctcf, atac, window=2097152):
    return seq[:window], ctcf[:window], atac[:window]

if __name__ == '__main__':
    import json
    main()
