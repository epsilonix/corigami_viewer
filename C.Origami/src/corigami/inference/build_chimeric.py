#!/usr/bin/env python
import os
import sys
import argparse
import tempfile
import json
import numpy as np
import pyBigWig
import gzip
from PIL import Image
import subprocess

################################################################################
# 1) HELPER FUNCTIONS FOR BUILDING A CHIMERIC CHROMOSOME
################################################################################

def extract_sequence(fa_path, start, end, do_reverse=False, do_flip=False):
    """
    Extract substring from a gzipped or plain-text FASTA file.
    - If do_reverse/do_flip is True, reverse the extracted sequence.
    Returns a lower-case string.
    """
    seq_str = []
    opener = gzip.open if fa_path.endswith(".gz") else open
    with opener(fa_path, "rt") as f:
        # Skip header lines
        lines = []
        for line in f:
            if line.startswith(">"):
                continue
            lines.append(line.strip().lower())
    full_seq = "".join(lines)
    
    chunk = full_seq[start:end]
    if (do_reverse or do_flip):
        chunk = chunk[::-1]  # simple pythonic string reverse
    return chunk

def extract_bigwig_region(bw_path, chr_name, start, end, default_val=0, do_reverse=False, do_flip=False):
    """
    Extract values from a bigwig file for a specified region.
    Returns an array of values, one per base pair.
    If do_reverse/do_flip is True, reverse the extracted signal array.
    """
    try:
        bw = pyBigWig.open(bw_path)
        
        # Check if chromosome exists in the bigwig file
        if chr_name not in bw.chroms():
            print(f"Warning: Chromosome {chr_name} not found in bigwig file {bw_path}")
            print(f"Available chromosomes: {list(bw.chroms().keys())}")
            # Return an array of default values
            return np.full(end - start, default_val)
        
        # Check if region is within chromosome bounds
        chr_length = bw.chroms()[chr_name]
        if start >= chr_length or end > chr_length:
            print(f"Warning: Region {start}-{end} exceeds chromosome {chr_name} length ({chr_length})")
            # Adjust the region to fit within chromosome bounds
            valid_start = min(start, chr_length-1)
            valid_end = min(end, chr_length)
            
            if valid_end <= valid_start:
                print(f"Error: No valid region remains within chromosome bounds")
                return np.full(end - start, default_val)
            
            # Create array with default values
            values = np.full(end - start, default_val)
            
            # Fill in valid portion with actual values
            valid_values = bw.values(chr_name, valid_start, valid_end)
            valid_length = valid_end - valid_start
            values[:valid_length] = valid_values
            
            return values
        
        # Normal case - region is valid
        vals = bw.values(chr_name, start, end)
        bw.close()
        
        # Replace None/nan with default value
        vals = np.array(vals)
        vals[np.isnan(vals)] = default_val
        
        # Reverse the array if requested
        if do_reverse or do_flip:
            vals = vals[::-1]
            
        return vals
        
    except RuntimeError as e:
        print(f"Error accessing bigwig file {bw_path} for region {chr_name}:{start}-{end}: {e}")
        return np.full(end - start, default_val)
    except Exception as e:
        print(f"Unexpected error processing bigwig file {bw_path}: {e}")
        return np.full(end - start, default_val)

def write_chrom_fasta(out_fa_path, chromosome_name, sequence):
    """
    Write a minimal FASTA for the new chromosome, line-wrapped at 70 chars.
    Creates a gzipped file with .fa.gz extension.
    """
    line_width = 70
    
    # Ensure the output path ends with .fa.gz
    if not out_fa_path.endswith('.gz'):
        out_fa_path = out_fa_path + '.gz'
    
    with gzip.open(out_fa_path, "wt") as out:
        out.write(f">{chromosome_name}\n")
        for i in range(0, len(sequence), line_width):
            out.write(sequence[i:i+line_width] + "\n")
    return out_fa_path

def write_bigwig_from_array(out_bw_path, chromosome_name, chrom_length, signal_array):
    """
    Create a BigWig with a single chromosome 'chromosome_name' of length 'chrom_length'.
    signal_array should be length == chrom_length (1 value per bp).
    """
    header = [(chromosome_name, chrom_length)]
    
    with pyBigWig.open(out_bw_path, "w") as bw:
        bw.addHeader(header)
        
        # For very large arrays, we might need to process in chunks to avoid memory issues
        chunk_size = 1000000  # Process 1M entries at a time
        
        for i in range(0, len(signal_array), chunk_size):
            end_idx = min(i + chunk_size, len(signal_array))
            chunk = signal_array[i:end_idx]
            
            # Create lists for starts and ends
            starts = list(range(i, end_idx))
            ends = list(range(i + 1, end_idx + 1))
            
            # Convert numpy values to a regular list if needed
            values = chunk.tolist() if isinstance(chunk, np.ndarray) else chunk
            
            # Add entries using lists of chroms, starts, ends, values
            chroms = [chromosome_name] * len(starts)
            bw.addEntries(chroms, starts, ends=ends, values=values)

################################################################################
# 2) MAIN: BUILD CHIMERA AND OPTIONALLY RUN DOWNSTREAM ANALYSIS
################################################################################

def main():
    """
    Example usage:
      python build_chimeric.py \
        --chr1 chr1 --start1 35000000 --end1 37097152 \
        --chr2 chr4 --start2 35000000 --end2 37097152 \
        --model corigami_data/model_weights/v1_jimin.ckpt \
        --seq ./corigami_data/data/hg38/dna_sequence \
        --ctcf ./corigami_data/data/hg38/imr90/genomic_features/ctcf_log2fc.bw \
        --atac ./my_atac.bw \
        --out ./output \
        [--first_reverse] [--second_reverse] ... \
        [--run_downstream prediction|screening]
    """
    parser = argparse.ArgumentParser(description="Build a chimeric chromosome from two regions.")
    parser.add_argument('--chr1', required=True, help="First chromosome name (e.g., chr1)")
    parser.add_argument('--start1', required=True, type=int)
    parser.add_argument('--end1', required=True, type=int)

    parser.add_argument('--chr2', required=True, help="Second chromosome name (e.g., chr4)")
    parser.add_argument('--start2', required=True, type=int)
    parser.add_argument('--end2', required=True, type=int)

    parser.add_argument('--model', required=True, help="Path to the model checkpoint")
    parser.add_argument('--seq', required=True, help="Folder with .fa or .fa.gz files")
    parser.add_argument('--ctcf', required=True, help="Path to CTCF bigwig")
    parser.add_argument('--atac', required=True, help="Path to ATAC bigwig")
    parser.add_argument('--out', default="outputs", help="Output folder")

    # Optional flipping/reversal
    parser.add_argument('--first_reverse', action='store_true', default=False)
    parser.add_argument('--first_flip',    action='store_true', default=False)
    parser.add_argument('--second_reverse', action='store_true', default=False)
    parser.add_argument('--second_flip',    action='store_true', default=False)

    # Optional consensus: only full_end is provided
    parser.add_argument('--full_end', type=int, default=None,
                        help="End coordinate for full region prediction (consensus). If not set, do standard 2 Mb.")

    # New argument: whether to run downstream analysis and which one to run.
    # Options: "none" (default, just build the synthetic files), "prediction", or "screening".
    parser.add_argument('--run_downstream', choices=['none','prediction','screening'], default='none',
                        help="Optionally run a downstream analysis on the synthetic chromosome.")

    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # 1) Extract sequences & signals
    fa1 = os.path.join(args.seq, f"{args.chr1}.fa.gz")
    fa2 = os.path.join(args.seq, f"{args.chr2}.fa.gz")
    seq1 = extract_sequence(fa1, args.start1, args.end1, 
                            do_reverse=args.first_reverse, 
                            do_flip=args.first_flip)
    seq2 = extract_sequence(fa2, args.start2, args.end2,
                            do_reverse=args.second_reverse,
                            do_flip=args.second_flip)
    chim_seq = seq1 + seq2
    chim_len = len(chim_seq)

    # CTCF 
    ctcf1 = extract_bigwig_region(args.ctcf, args.chr1, args.start1, args.end1,
                                  do_reverse=args.first_reverse or args.first_flip,
                                  do_flip=False)
    ctcf2 = extract_bigwig_region(args.ctcf, args.chr2, args.start2, args.end2,
                                  do_reverse=args.second_reverse or args.second_flip,
                                  do_flip=False)
    chim_ctcf = np.concatenate([ctcf1, ctcf2])

    # ATAC
    atac1 = extract_bigwig_region(args.atac, args.chr1, args.start1, args.end1,
                                  do_reverse=args.first_reverse or args.first_flip,
                                  do_flip=False)
    atac2 = extract_bigwig_region(args.atac, args.chr2, args.start2, args.end2,
                                  do_reverse=args.second_reverse or args.second_flip,
                                  do_flip=False)
    
    chim_atac = np.concatenate([atac1, atac2])


    # 2) Write out the synthetic chromosome files ("chrCHIM")
    chim_name = "chrCHIM"
    
    # Use the specified output directory instead of a temporary one
    output_dir = args.out
    temp_dir = os.path.join(output_dir, "chimera_files")
    os.makedirs(temp_dir, exist_ok=True)

    # Write FASTA - note that write_chrom_fasta() appends .gz if needed
    chim_fa = os.path.join(temp_dir, "chrCHIM.fa")
    chim_fa_gz = write_chrom_fasta(chim_fa, chim_name, chim_seq)
    
    # Write BigWig files
    chim_ctcf_bw = os.path.join(temp_dir, "chrCHIM_ctcf.bw")
    chim_atac_bw = os.path.join(temp_dir, "chrCHIM_atac.bw")
    write_bigwig_from_array(chim_ctcf_bw, chim_name, chim_len, chim_ctcf)
    write_bigwig_from_array(chim_atac_bw, chim_name, chim_len, chim_atac)

    print(f"[Chimera] Synthetic chromosome built: {chim_name} (length={chim_len:,} bp)")
    print(f"[Chimera] Files written to directory: {temp_dir}")
    print(f"[Chimera] FASTA at: {chim_fa_gz}")
    print(f"[Chimera] CTCF bigwig at: {chim_ctcf_bw}")
    print(f"[Chimera] ATAC bigwig at: {chim_atac_bw}")
    
    # 2.5) Generate peaks from ATAC signal for screening
    print(f"[Chimera] Generating peaks from ATAC signal...")
    chim_peaks_file = os.path.join(temp_dir, "chrCHIM_peaks.bed")
    
    # Method 1: Use MACS2 peak calling if available
    try:
        from corigami.utils.peak_calling import generate_peaks_from_bigwig
        generate_peaks_from_bigwig(chim_atac_bw, chim_peaks_file, chim_name, 0, chim_len)
        print(f"[Chimera] Generated peaks using peak caller at: {chim_peaks_file}")
    except ImportError:
        # Method 2: Simple threshold-based peak calling as fallback
        print(f"[Chimera] No peak caller found, using simple threshold-based approach")
        with open(chim_peaks_file, 'w') as f:
            # Find peaks using a simple threshold approach
            threshold = np.percentile(chim_atac, 95)  # Use top 5% as peaks
            peaks = []
            in_peak = False
            peak_start = 0
            
            for i, val in enumerate(chim_atac):
                if not in_peak and val > threshold:
                    # Start of peak
                    in_peak = True
                    peak_start = i
                elif in_peak and val <= threshold:
                    # End of peak
                    in_peak = False
                    peak_end = i
                    if peak_end - peak_start >= 50:  # Minimum peak width
                        f.write(f"{chim_name}\t{peak_start}\t{peak_end}\t.\t{val}\t.\n")
                        peaks.append((peak_start, peak_end))
            
            # Handle case where we're still in a peak at the end
            if in_peak:
                peak_end = len(chim_atac)
                if peak_end - peak_start >= 50:
                    f.write(f"{chim_name}\t{peak_start}\t{peak_end}\t.\t{chim_atac[peak_start]}\t.\n")
                    peaks.append((peak_start, peak_end))
        
        print(f"[Chimera] Generated {len(peaks)} peaks using threshold method: {chim_peaks_file}")
    
    # Write a simple manifest file with all paths
    manifest_path = os.path.join(temp_dir, "manifest.txt")
    with open(manifest_path, 'w') as f:
        f.write(f"CHIM_DIR={temp_dir}\n")
        f.write(f"CHIM_NAME={chim_name}\n")
        f.write(f"CHIM_LENGTH={chim_len}\n")
        f.write(f"CHIM_FASTA={chim_fa_gz}\n")
        f.write(f"CHIM_CTCF={chim_ctcf_bw}\n")
        f.write(f"CHIM_ATAC={chim_atac_bw}\n")
        f.write(f"CHIM_PEAKS={chim_peaks_file}\n")
    
    # Just print the directory path - simple and unambiguous
    print(f"CHIMERA_OUTPUT_DIR={temp_dir}")

    # 3) Either run downstream analysis or exit
    if args.run_downstream != "none":
        # Determine which downstream script to run.
        if args.run_downstream == "prediction":
            # Assuming prediction.py is in the same directory as this script.
            downstream_script = os.path.join(os.path.dirname(__file__), "prediction.py")
            # Build command for prediction.
            if chim_len > 2000000:
                full_end = chim_len if not args.full_end else args.full_end
                cmd = [
                    "python", downstream_script,
                    "--chr", chim_name,
                    "--start", "0",
                    "--full_end", str(full_end),
                    "--model", args.model,
                    "--seq", os.path.dirname(chim_fa_gz),
                    "--ctcf", chim_ctcf_bw,
                    "--atac", chim_atac_bw,
                    "--out", args.out
                ]
            else:
                cmd = [
                    "python", downstream_script,
                    "--chr", chim_name,
                    "--start", "0",
                    "--model", args.model,
                    "--seq", os.path.dirname(chim_fa_gz),
                    "--ctcf", chim_ctcf_bw,
                    "--atac", chim_atac_bw,
                    "--out", args.out
                ]
        elif args.run_downstream == "screening":
            # Assuming screening.py is in the same directory as this script.
            downstream_script = os.path.join(os.path.dirname(__file__), "screening.py")
            # Build command for screening.
            # Note: screening.py may require additional parameters (e.g., screen-start, screen-end, perturb_width, step_size).
            # You might choose to add extra arguments here or set defaults.
            # For example purposes, we set screen_start=0 and screen_end equal to chim_len.
            cmd = [
                "python", downstream_script,
                "--chr", chim_name,
                "--screen-start", "0",
                "--screen-end", str(chim_len),
                "--model", args.model,
                "--seq", os.path.dirname(chim_fa_gz),
                "--ctcf", chim_ctcf_bw,
                "--atac", chim_atac_bw,
                "--out", args.out,
                "--perturb-width", "1000",
                "--step-size", "1000"
            ]
        cmd_str = " ".join(cmd)
        print(f"[Chimera] Running downstream {args.run_downstream} with command:")
        print(cmd_str)
    
        # Set up the environment (e.g., PYTHONPATH) as needed.
        env = os.environ.copy()
        src_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        if "PYTHONPATH" in env:
            env["PYTHONPATH"] = f"{src_dir}:{env['PYTHONPATH']}"
        else:
            env["PYTHONPATH"] = src_dir
        print(f"[Chimera] Setting PYTHONPATH to: {env['PYTHONPATH']}")
    
        try:
            ret = subprocess.run(cmd_str, shell=True, env=env, check=True, capture_output=True, text=True)
            print(f"[Chimera] Downstream {args.run_downstream} completed successfully with exit code {ret.returncode}")
            if hasattr(ret, 'stdout'):
                print(f"[Chimera] Command stdout: {ret.stdout}")
        except subprocess.CalledProcessError as e:
            print(f"[Chimera] ERROR: Downstream {args.run_downstream} failed with exit code {e.returncode}")
            print(f"[Chimera] STDERR: {e.stderr if hasattr(e, 'stderr') else 'No stderr captured'}")
            print(f"[Chimera] STDOUT: {e.stdout if hasattr(e, 'stdout') else 'No stdout captured'}")
    else:
        # No downstream analysis requested; output the file paths in JSON.
        result_info = {
            "chromosome": chim_name,
            "length": chim_len,
            "fasta": chim_fa_gz,
            "ctcf_bigwig": chim_ctcf_bw,
            "atac_bigwig": chim_atac_bw,
            "temp_dir": temp_dir
        }
        print(json.dumps(result_info, indent=2))

if __name__ == "__main__":
    main()
