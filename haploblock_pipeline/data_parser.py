#!/usr/bin/env python3
import logging
import pathlib
import subprocess
import os

import numpy as np
from scipy.ndimage import correlate1d

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Recombination map parsing (MAIN CPU HOTSPOT – optimized)
# ---------------------------------------------------------------------
def parse_recombination_rates(recombination_file, chromosome):
    """
    Fast CPU-optimized parsing of Halldorsson2019 recombination map.
    Returns list of (start, end) tuples.

    Original scipy-based reference implementation, left untouched so
    it always exists as ground truth for verify_accelerated_matches_cpu()
    below, independent of whatever happens on the GPU/accelerated path.
    """
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"

    data = np.fromiter(
        ((int(line[1]), float(line[3]))
         for line in map(str.split, open(recombination_file))
         if not line[0].startswith("#") and line[0] != "Chr" and line[0] == chromosome),
        dtype=[('start', 'i8'), ('rate', 'f8')]
    )

    positions = data['start']
    rates = data['rate']

    if len(rates) < 3:
        raise ValueError(f"Not enough data points for {chromosome}")

    from scipy.ndimage import gaussian_filter1d
    smoothed = gaussian_filter1d(rates, sigma=5)
    peaks = np.where(
        (smoothed[1:-1] > smoothed[:-2]) &
        (smoothed[1:-1] > smoothed[2:])
    )[0] + 1

    if len(peaks) == 0:
        raise ValueError(f"No recombination peaks detected for {chromosome}")

    peak_positions = positions[peaks]

    haploblocks = [(1, peak_positions[0])]
    haploblocks.extend(zip(peak_positions[:-1], peak_positions[1:]))
    haploblocks.append((peak_positions[-1], positions[-1]))

    logger.info("Found %d haploblocks for %s", len(haploblocks), chromosome)
    return haploblocks


# ---------------------------------------------------------------------
# GPU-accelerated recombination map parsing, with deterministic
# CPU fallback.
# ---------------------------------------------------------------------
try:
    import cupy as cp
    _GPU_AVAILABLE = bool(cp.cuda.runtime.getDeviceCount())
except Exception:
    cp = None
    _GPU_AVAILABLE = False

# Below this many points, GPU transfer/launch overhead outweighs any
# benefit, so "auto" mode stays on CPU regardless of GPU presence.
_GPU_MIN_POINTS = 20_000


def _gaussian_kernel1d(sigma, truncate=4.0):
    """
    Build the same 1D Gaussian kernel scipy.ndimage.gaussian_filter1d
    uses internally (truncate=4.0, order=0), as an explicit weight
    array, so we can hand the identical weights to correlate1d on
    either backend instead of relying on gaussian_filter1d's own
    kernel construction.
    """
    radius = int(truncate * sigma + 0.5)
    sigma2 = sigma * sigma
    x = np.arange(-radius, radius + 1)
    phi = np.exp(-0.5 / sigma2 * x ** 2)
    phi /= phi.sum()
    return radius, phi


def _gaussian_smooth(rates_xp, sigma, xp):
    """
    Gaussian smoothing via a single correlate1d call with our own
    explicit kernel weights (see _gaussian_kernel1d), rather than a
    hand-rolled per-tap Python loop.

    On CPU (xp=numpy): calls scipy.ndimage.correlate1d directly.
    gaussian_filter1d is itself implemented as "build this exact
    Gaussian kernel, then call correlate1d with mode='reflect'" - so
    this CPU path is provably identical to parse_recombination_rates
    above, not just empirically matched to it.

    On GPU (xp=cupy): calls cupyx.scipy.ndimage.correlate1d, a single
    fused CUDA kernel. A hand-rolled per-tap loop would issue one
    kernel launch per tap (~41 launches at this kernel radius);
    this issues one.
    """
    _, weights = _gaussian_kernel1d(sigma)
    weights_xp = xp.asarray(weights)
    if xp is np:
        return correlate1d(rates_xp, weights_xp, mode="reflect")
    else:
        from cupyx.scipy.ndimage import correlate1d as gpu_correlate1d
        return gpu_correlate1d(rates_xp, weights_xp, mode="reflect")


def _select_backend(force_gpu, n_points):
    """
    force_gpu: True -> require GPU, raise if unavailable.
               False -> force CPU.
               None -> auto: GPU if available and n_points >= threshold.
    Returns (xp_module, backend_name).
    """
    if force_gpu is False:
        return np, "cpu"
    if not _GPU_AVAILABLE:
        if force_gpu is True:
            raise RuntimeError(
                "GPU requested (gpu=True / --gpu on) but no CUDA device "
                "or cupy installation was found"
            )
        return np, "cpu"
    if force_gpu is True:
        return cp, "gpu"
    return (cp, "gpu") if n_points >= _GPU_MIN_POINTS else (np, "cpu")


def parse_recombination_rates_accelerated(recombination_file, chromosome, gpu=None):
    """
    GPU-accelerated version of parse_recombination_rates, with automatic
    CPU fallback. gpu=None (default) auto-detects; True forces GPU and
    errors if unavailable; False forces CPU.

    Returns the same list[(start, end)] as parse_recombination_rates.
    """
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"

    data = np.fromiter(
        ((int(line[1]), float(line[3]))
         for line in map(str.split, open(recombination_file))
         if not line[0].startswith("#") and line[0] != "Chr" and line[0] == chromosome),
        dtype=[('start', 'i8'), ('rate', 'f8')]
    )

    positions = data['start']
    rates = data['rate']

    if len(rates) < 3:
        raise ValueError(f"Not enough data points for {chromosome}")

    xp, backend = _select_backend(force_gpu=gpu, n_points=len(rates))
    logger.info("Using %s backend for %s (%d points)", backend, chromosome, len(rates))

    rates_xp = xp.asarray(rates)
    smoothed_xp = _gaussian_smooth(rates_xp, sigma=5, xp=xp)
    smoothed = cp.asnumpy(smoothed_xp) if backend == "gpu" else smoothed_xp

    peaks = np.where(
        (smoothed[1:-1] > smoothed[:-2]) &
        (smoothed[1:-1] > smoothed[2:])
    )[0] + 1

    if len(peaks) == 0:
        raise ValueError(f"No recombination peaks detected for {chromosome}")

    peak_positions = positions[peaks]

    haploblocks = [(1, peak_positions[0])]
    haploblocks.extend(zip(peak_positions[:-1], peak_positions[1:]))
    haploblocks.append((peak_positions[-1], positions[-1]))

    logger.info("Found %d haploblocks for %s (%s)", len(haploblocks), chromosome, backend)
    return haploblocks


def verify_accelerated_matches_cpu(recombination_file, chromosome):
    """
    Cross-checks the original scipy-based CPU implementation against
    the accelerated CPU path (provably identical - see _gaussian_smooth
    docstring), and, if a GPU is present, the GPU path against its own
    CPU fallback (empirically checked - GPU/CPU are separate CUDA/C
    implementations, so this is the practical determinism guarantee,
    not a mathematical one). Raises AssertionError on any mismatch.
    """
    reference = parse_recombination_rates(recombination_file, chromosome)
    cpu_new = parse_recombination_rates_accelerated(recombination_file, chromosome, gpu=False)
    if reference != cpu_new:
        raise AssertionError(
            "Accelerated CPU path does not match the original scipy-based "
            "parse_recombination_rates - do not trust the GPU path until "
            "this is resolved."
        )

    if _GPU_AVAILABLE:
        gpu_result = parse_recombination_rates_accelerated(recombination_file, chromosome, gpu=True)
        if cpu_new != gpu_result:
            raise AssertionError(
                "GPU result does not match its own CPU fallback - check "
                "cupy/CUDA version and consider tolerance-based comparison "
                "instead of exact equality if this is a genuine ULP-level "
                "difference rather than a bug."
            )
        logger.info("Verified: scipy-CPU, accelerated-CPU, and GPU all agree for %s", chromosome)
    else:
        logger.info("Verified: scipy-CPU and accelerated-CPU agree for %s; no GPU present to check", chromosome)


# ---------------------------------------------------------------------
# Simple TSV parsers
# ---------------------------------------------------------------------
def parse_haploblock_boundaries(boundaries_file):
    """
    Parses haploblock boundaries TSV (header: START\tEND)
    Returns integer tuples.
    """
    with open(boundaries_file) as f:
        header = next(f)
        if not header.startswith("START\t"):
            raise ValueError("Boundaries file missing header")
        return [tuple(map(int, line.rstrip().split("\t"))) for line in f]


def parse_samples(samples_file):
    with open(samples_file) as f:
        header = next(f)
        if not header.startswith("Sample name\t"):
            raise ValueError("Samples file missing header")

        samples = []
        for line in f:
            sample = line.split("\t", 1)[0]
            if not (sample.startswith("HG") or sample.startswith("NA")):
                raise ValueError(f"Invalid sample line: {line}")
            samples.append(sample)

    logger.info("Found %d samples", len(samples))
    return samples


def parse_samples_from_vcf(vcf):
    samples = subprocess.run(
        ["bcftools", "query", "-l", vcf],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()

    logger.info("Found %d samples", len(samples))
    return samples


def parse_variants_of_interest(variants_file):
    variants = []
    with open(variants_file) as f:
        for line in f:
            chr_pos = line.rstrip().split(":")
            if len(chr_pos) != 2:
                raise ValueError(f"Bad variant line: {line}")
            variants.append(chr_pos[1])
    return variants


# ---------------------------------------------------------------------
# VCF / FASTA extraction (unchanged, I/O bound)
# ---------------------------------------------------------------------
def extract_region_from_vcf(vcf, chr, chr_map, start, end, out):
    if chr.startswith("chr"):
        chr = chr.replace("chr", "")

    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    temp_vcf = tmp_dir / f"{chr}_region_{start}-{end}.vcf.gz"

    subprocess.run(
        ["bcftools", "view",
         "-r", f"{chr}:{start}-{end}",
         "--min-af", "0.05",
         vcf,
         "-o", temp_vcf],
        check=True,
    )
    subprocess.run(["bcftools", "index", temp_vcf], check=True)

    output_vcf = tmp_dir / f"chr{chr}_region_{start}-{end}.vcf"
    subprocess.run(
        ["bcftools", "annotate", "--rename-chrs", chr_map, temp_vcf],
        stdout=open(output_vcf, "w"),
        check=True,
    )
    subprocess.run(["bgzip", output_vcf], check=True)
    subprocess.run(
        ["bcftools", "index", "-c", output_vcf.with_suffix(".vcf.gz")],
        check=True,
    )

    return output_vcf.with_suffix(".vcf.gz")


def extract_sample_from_vcf(vcf, sample, out):
    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    output_vcf = tmp_dir / f"{sample}_{vcf.stem}.gz"

    subprocess.run(
        ["bcftools", "view",
         "--force-samples",
         "-s", sample,
         "-o", output_vcf,
         str(vcf)],
        check=True,
    )
    subprocess.run(["bcftools", "index", output_vcf], check=True)
    return output_vcf


def extract_region_from_fasta(fasta, chr, start, end, out):
    subprocess.run(["samtools", "faidx", fasta], check=True)

    tmp_dir = pathlib.Path(out) / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    output_fasta = tmp_dir / f"chr{chr}_region_{start}-{end}.fa"
    subprocess.run(
        ["samtools", "faidx", fasta, f"chr{chr}:{start}-{end}"],
        stdout=open(output_fasta, "w"),
        check=True,
    )
    return output_fasta


# ---------------------------------------------------------------------
# Cluster parsing (single-pass optimized)
# ---------------------------------------------------------------------
def parse_clusters(clusters_file):
    representative2cluster = {}
    individual2cluster = {}
    clusters = []
    next_cluster = 0

    with open(clusters_file) as f:
        for line in f:
            rep, indiv = line.rstrip().split("\t")
            if rep not in representative2cluster:
                representative2cluster[rep] = next_cluster
                clusters.append(next_cluster)
                next_cluster += 1
            individual2cluster[indiv] = representative2cluster[rep]

    return individual2cluster, clusters
