#!/usr/bin/env python3
"""
GPU Optimized Step 1: Generate haploblock boundaries from a recombination map
- Deterministic CPU/GPU output
- Vectorized parsing and peak detection
- Optional GPU acceleration
"""

import sys
import logging
import argparse
import pathlib
import numpy as np
import scipy.ndimage

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------
# Core functions
# -------------------------------------------------------------------------

def haploblocks_to_tsv(haploblocks, chrom, out_dir: pathlib.Path):
    """Save haploblock boundaries to a TSV file."""
    out_dir.mkdir(parents=True, exist_ok=True)
    output_file = out_dir / f"haploblock_boundaries_chr{chrom}.tsv"
    with open(output_file, 'w') as f:
        f.write("START\tEND\n")
        f.writelines(f"{start}\t{end}\n" for start, end in haploblocks)
    logger.info(f"Haploblock boundaries written to {output_file}")
    return output_file


def parse_recombination_file(recombination_file: pathlib.Path, chrom: str):
    """
    Deterministic parsing of recombination file into positions & rates.
    Vectorized via NumPy for speed.
    """
    dt = np.dtype([
        ("chr", "U6"),
        ("start", np.int64),
        ("end", np.int64),
        ("rate", np.float64),
        ("cm", np.float64),
    ])
    try:
        data = np.genfromtxt(recombination_file, dtype=dt, comments="#", delimiter="\t")
    except Exception as e:
        logger.error("Failed to read recombination file %s: %s", recombination_file, e)
        raise

    if not str(chrom).startswith("chr"):
        chrom = "chr" + str(chrom)

    mask = data["chr"] == chrom
    positions = np.ascontiguousarray(data["start"][mask], dtype=np.int64)
    rates = np.ascontiguousarray(data["rate"][mask], dtype=np.float64)

    if len(positions) == 0:
        logger.warning("No data for chromosome %s", chrom)
        return np.array([]), np.array([])

    return positions, rates


def smooth_rates(rates: np.ndarray, use_gpu=False, gpu_id=0):
    """Gaussian smoothing (CPU or GPU), deterministic."""
    sigma = 5.0

    if use_gpu:
        try:
            import cupy as cp
            import cupyx.scipy.ndimage as cpx_ndimage

            cp.cuda.Device(gpu_id).use()
            rates_gpu = cp.asarray(rates)
            smoothed_gpu = cpx_ndimage.gaussian_filter1d(rates_gpu, sigma=sigma, mode="reflect")
            smoothed = cp.asnumpy(smoothed_gpu)

            # clean GPU memory
            rates_gpu = None
            smoothed_gpu = None
            cp.get_default_memory_pool().free_all_blocks()

            logger.info("Smoothing completed on GPU")
            return smoothed

        except Exception as e:
            logger.warning("GPU smoothing failed (%s). Falling back to CPU.", e)

    # CPU fallback
    smoothed = scipy.ndimage.gaussian_filter1d(rates, sigma=sigma, mode="reflect")
    logger.info("Smoothing completed on CPU")
    return smoothed


def detect_peaks(smoothed: np.ndarray):
    """Vectorized peak detection."""
    if len(smoothed) < 3:
        return np.array([], dtype=int)
    left = smoothed[1:-1] > smoothed[:-2]
    right = smoothed[1:-1] > smoothed[2:]
    peaks = np.where(left & right)[0] + 1
    return peaks


def build_haploblocks(positions: np.ndarray, peaks: np.ndarray):
    """Generate haploblock boundaries from peak positions."""
    if len(positions) == 0:
        return []

    if len(peaks) == 0:
        return [(1, int(positions[-1]))]

    peak_positions = positions[peaks]
    haploblocks = [(1, int(peak_positions[0]))]

    for i in range(1, len(peak_positions)):
        haploblocks.append((int(peak_positions[i-1]), int(peak_positions[i])))

    haploblocks.append((int(peak_positions[-1]), int(positions[-1])))

    logger.info("Detected %d haploblocks", len(haploblocks))
    return haploblocks


def run_haploblocks(recombination_file, chrom, out_dir, use_gpu=False, gpu_id=0):
    """Main Step 1 pipeline."""
    logger.info("Parsing recombination file: %s", recombination_file)
    positions, rates = parse_recombination_file(recombination_file, chrom)

    if len(positions) == 0:
        haploblocks = []
    else:
        smoothed = smooth_rates(rates, use_gpu=use_gpu, gpu_id=gpu_id)
        peaks = detect_peaks(smoothed)
        haploblocks = build_haploblocks(positions, peaks)

    return haploblocks_to_tsv(haploblocks, chrom, out_dir)


# -------------------------------------------------------------------------
# CLI wrapper
# -------------------------------------------------------------------------
if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=logging.INFO
    )

    parser = argparse.ArgumentParser(
        description="Step 1: Generate haploblock boundaries"
    )
    parser.add_argument('--recombination_file', type=pathlib.Path, required=True)
    parser.add_argument('--chr', type=str, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    parser.add_argument('--gpu', action='store_true', help="Use GPU if available")
    parser.add_argument('--gpu_id', type=int, default=0)

    args = parser.parse_args()

    try:
        run_haploblocks(
            args.recombination_file,
            args.chr,
            args.out,
            use_gpu=args.gpu,
            gpu_id=args.gpu_id
        )
    except Exception as e:
        logger.error("Step 1 failed: %s", e)
        sys.exit(1)


# -------------------------------------------------------------------------
# Pipeline entry point for other scripts
# -------------------------------------------------------------------------
def run(recombination_file, chr=None, chrom=None, out=None, threads=None, gpu=False, gpu_id=0):
    """
    Pipeline-compatible run function.

    Accepts either `chr` or `chrom` from config for backwards compatibility.
    """
    chrom = chrom if chrom is not None else chr
    if chrom is None:
        raise ValueError("chromosome number must be provided via 'chrom' or 'chr'")
    
    return run_haploblocks(
        pathlib.Path(recombination_file),
        str(chrom),
        pathlib.Path(out),
        use_gpu=gpu,
        gpu_id=gpu_id
    )
