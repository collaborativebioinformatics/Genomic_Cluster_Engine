#!/usr/bin/env python3
from __future__ import annotations

import os
import logging
import pathlib
import subprocess
from typing import Dict, List, Optional

import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool, cpu_count

import data_parser

logger = logging.getLogger(__name__)

CLUSTER_HASH_LENGTH = 20
HAPLOBLOCK_HASH_LENGTH = 20
PARALLEL_THRESHOLD = 1000  # parallelize > 1000 haploblocks

# ---------------------------------------------------------------------
# Optional GPU backend (mirrors data_parser.py's guarded import so
# this module is fully importable, and fully usable, on CPU-only nodes)
# ---------------------------------------------------------------------
try:
    import cupy as cp
    _GPU_AVAILABLE = bool(cp.cuda.runtime.getDeviceCount())
except Exception:
    cp = None
    _GPU_AVAILABLE = False


def _resolve_gpu_choice(gpu):
    """
    Normalizes the gpu argument to True (force GPU), False (force CPU),
    or None (auto). Accepts "auto"/"on"/"off" strings as well as legacy
    booleans, so existing pipeline callers using gpu=True/False keep
    working unchanged.
    """
    if isinstance(gpu, bool):
        return gpu
    if gpu in (None, "auto"):
        return None
    if gpu == "on":
        return True
    if gpu == "off":
        return False
    raise ValueError(f"Invalid gpu option: {gpu!r} (expected auto/on/off or bool)")


def _cpu_workers(threads=None):
    """
    os.cpu_count() can return None on some containerized/restricted
    environments, which raised TypeError before `- 1` even ran. Fall
    back to 2 in that case so the `- 1 or 1` floor still works.
    """
    if threads:
        return threads
    return (os.cpu_count() or 2) - 1 or 1


def _check_hash_capacity(n, width, what):
    """
    np.binary_repr(i, width=width) *widens* past `width` characters for
    i >= 2**width rather than truncating, while the GPU bit-shift path
    keeps a fixed width and silently wraps (i % 2**width). Left
    unguarded, that's a real CPU/GPU divergence once n exceeds the
    hash capacity - raise loudly instead of letting it happen quietly.
    """
    capacity = 1 << width
    if n > capacity:
        raise ValueError(
            f"{n:,} {what} exceed the {width}-bit hash capacity "
            f"(2^{width} = {capacity:,}); CPU and GPU hashes would "
            f"diverge past this point - increase the hash length "
            f"constant instead of proceeding."
        )


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------
def _make_hash_cpu(i: int, width: int) -> str:
    return np.binary_repr(i, width=width)


def _make_hash_gpu(indices, width: int):
    """Convert integer indices to binary arrays on GPU (0/1)."""
    bin_array = ((indices[:, None] & (1 << cp.arange(width)[::-1])) > 0).astype(cp.uint8)
    return bin_array


def chromosome_to_int(chrom: str) -> int:
    chrom = chrom.replace("chr", "")  # safety

    if chrom.isdigit():
        return int(chrom)

    mapping = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25
    }

    if chrom in mapping:
        return mapping[chrom]

    raise ValueError(f"Unknown chromosome: {chrom}")


# ---------------------------------------------------------------------
# Backend cross-check
# ---------------------------------------------------------------------
def verify_hash_backends_match(n, width):
    """
    Confirms CPU- and GPU-generated index->binary hashes agree for the
    first n indices at the given bit width. Bitwise integer ops are
    exact on both backends (no floating point involved, unlike the
    Gaussian smoothing in step 1), so a mismatch here should only ever
    come from the capacity-overflow case or an accidental bit-order
    change upstream - this exists to catch exactly that.
    """
    _check_hash_capacity(n, width, "indices")
    cpu_hashes = [_make_hash_cpu(i, width) for i in range(n)]

    if not _GPU_AVAILABLE:
        logger.info("Verified: %d CPU hashes generated; no GPU present to cross-check", n)
        return

    indices = cp.arange(n, dtype=cp.uint32)
    gpu_hashes = ["".join(map(str, row)) for row in cp.asnumpy(_make_hash_gpu(indices, width))]

    if cpu_hashes != gpu_hashes:
        mismatches = [i for i, (a, b) in enumerate(zip(cpu_hashes, gpu_hashes)) if a != b][:5]
        raise AssertionError(f"CPU/GPU hash mismatch at indices {mismatches} (showing up to 5)")

    logger.info("Verified: CPU and GPU hashes agree for all %d indices (width=%d)", n, width)


# ---------------------------------------------------------------------
# Hash generators
# ---------------------------------------------------------------------
def generate_haploblock_hashes(haploblock_boundaries: list[tuple[int, int]], use_gpu=False):
    n = len(haploblock_boundaries)
    logger.info(f"Generating hashes for {n:,} haploblocks (GPU={use_gpu})")
    _check_hash_capacity(n, HAPLOBLOCK_HASH_LENGTH, "haploblocks")

    if use_gpu:
        indices = cp.arange(n, dtype=cp.uint32)
        bin_array = _make_hash_gpu(indices, HAPLOBLOCK_HASH_LENGTH)
        haploblock2hash = {hap: "".join(map(str, cp.asnumpy(bin_array[i]))) for i, hap in enumerate(haploblock_boundaries)}
    else:
        if n > PARALLEL_THRESHOLD:
            with Pool(cpu_count()) as pool:
                hashes = pool.starmap(_make_hash_cpu, [(i, HAPLOBLOCK_HASH_LENGTH) for i in range(n)])
        else:
            hashes = [_make_hash_cpu(i, HAPLOBLOCK_HASH_LENGTH) for i in range(n)]
        haploblock2hash = dict(zip(haploblock_boundaries, hashes))
    return haploblock2hash


def generate_cluster_hashes(clusters: list[int], use_gpu=False):
    n = len(clusters)
    _check_hash_capacity(n, CLUSTER_HASH_LENGTH, "clusters")

    if use_gpu:
        indices = cp.arange(n, dtype=cp.uint32)
        bin_array = _make_hash_gpu(indices, CLUSTER_HASH_LENGTH)
        cluster2hash = {cluster: "".join(map(str, cp.asnumpy(bin_array[i]))) for i, cluster in enumerate(clusters)}
    else:
        cluster2hash = {cluster: np.binary_repr(i, width=CLUSTER_HASH_LENGTH) for i, cluster in enumerate(clusters)}
    return cluster2hash


def generate_variant_hashes(variants: List[str],
                            vcf: pathlib.Path,
                            chrom: str,
                            start: int,
                            end: int,
                            samples: Optional[List[str]]) -> Dict[str, str]:
    """
    Generate binary variant presence hashes for all samples and
    haplotypes, scoped to a single haploblock's own (start, end).

    arguments:
    - variants: list of string variant positions (whole-chromosome list)
    - vcf
    - chrom
    - start, end: this haploblock's own boundaries
    - samples

    returns:
    - variant2hash: dict, key=individual, values=hash. Empty dict if
      no variants of interest fall inside this block.
    """
    if not samples:
        logger.warning("No samples provided for variant hash generation. Returning empty dict.")
        return {}

    block_variants = [v for v in variants if start <= int(v) <= end]
    if not block_variants:
        return {}

    idx_map = {str(v): i for i, v in enumerate(block_variants)}
    variant2hash = {
        f"{sample}_chr{chrom}_region_{start}-{end}_hap{h}": ["0"] * len(block_variants)
        for sample in samples for h in (0, 1)
    }

    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n",
         "-s", ",".join(samples), "--force-samples", "-r", region, str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    )

    for line in result.stdout.splitlines():
        _, pos, *gts = line.split("\t")
        if pos not in idx_map:
            continue

        i = idx_map[pos]
        for sample_idx, gt in enumerate(gts):
            if "|" not in gt:
                continue
            a0, a1 = gt.split("|")
            sample = samples[sample_idx]
            if a0 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap0"][i] = "1"
            if a1 == "1":
                variant2hash[f"{sample}_chr{chrom}_region_{start}-{end}_hap1"][i] = "1"

    variant2hash = {k: "".join(v) for k, v in variant2hash.items()}

    return(variant2hash)


# ---------------------------------------------------------------------
# Individual hash generation (CPU)
# ---------------------------------------------------------------------
def generate_individual_hash(individual,
                             individual2cluster,
                             cluster2hash,
                             haploblock2hash,
                             chr_hash,
                             variant2hash=None):
    strand = individual[-1]
    if strand == "0":
        strand_hash = "0001"
    elif strand == "1":
        strand_hash = "0010"
    else:
        raise ValueError(f"Individual ID {individual!r} does not end in '0' or '1' (strand): {strand!r}")

    individual_split = individual.split("_")
    region_str = individual_split[3].replace(".fa", "").replace(".fasta", "").replace(".vcf", "")
    start, end = map(int, region_str.split("-"))
    haploblock_hash = haploblock2hash[(start, end)]
    cluster_hash = cluster2hash[individual2cluster[individual]]
    hash_val = strand_hash + chr_hash + haploblock_hash + cluster_hash
    if variant2hash:
        hash_val += variant2hash[individual]

    return individual, hash_val


def generate_individual_hashes_parallel(individual2cluster, cluster2hash, haploblock2hash,
                                        chr_hash, variant2hash=None, max_workers=None):
    max_workers = _cpu_workers(max_workers)
    individual2hash = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(generate_individual_hash, ind, individual2cluster, cluster2hash,
                                   haploblock2hash, chr_hash, variant2hash)
                   for ind in individual2cluster]
        for fut in as_completed(futures):
            ind, h = fut.result()
            individual2hash[ind] = h
    return {ind: individual2hash[ind] for ind in individual2cluster}


# ---------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------
def run_hashes(boundaries_file: pathlib.Path,
               clusters_dir: pathlib.Path,
               chrom: str,
               out: pathlib.Path,
               variants_file: Optional[pathlib.Path] = None,
               vcf: Optional[pathlib.Path] = None,
               samples_file: Optional[pathlib.Path] = None,
               threads: Optional[int] = None,
               gpu="auto",
               gpu_id: Optional[int] = None,
               verify: bool = False):
    """
    gpu: "auto" (GPU if available, default), "on" (force GPU, error if
         unavailable), "off" (force CPU) - or a bool, as main.py's YAML
         config passes.
    verify: if True, cross-check CPU/GPU hash generation for the exact
            haploblock and cluster counts encountered in this run
            before writing output.
    """
    gpu_choice = _resolve_gpu_choice(gpu)
    use_gpu = False
    if gpu_choice is not False:
        if _GPU_AVAILABLE:
            try:
                cp.cuda.Device(gpu_id or 0).use()
                use_gpu = True
                logger.info("GPU enabled (ID=%d)", gpu_id or 0)
            except Exception as e:
                if gpu_choice is True:
                    raise RuntimeError(f"GPU requested (--gpu on) but device setup failed: {e}") from e
                logger.warning("GPU device setup failed (%s); falling back to CPU.", e)
        elif gpu_choice is True:
            raise RuntimeError("GPU requested (--gpu on) but no CUDA device / cupy install found")
        else:
            logger.info("No GPU available; using CPU.")

    chr_hash = np.binary_repr(chromosome_to_int(chrom), width=5)

    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)
    if verify:
        verify_hash_backends_match(len(haploblock_boundaries), HAPLOBLOCK_HASH_LENGTH)
    haploblock2hash = generate_haploblock_hashes(haploblock_boundaries, use_gpu=use_gpu)

    out.mkdir(parents=True, exist_ok=True)
    haploblock_hashes_file = out / "haploblock_hashes.tsv"
    with haploblock_hashes_file.open("w") as f:
        f.write("START\tEND\tHASH\n")
        for (start, end), h in haploblock2hash.items():
            f.write(f"{start}\t{end}\t{h}\n")

    variants = None
    samples = None
    if variants_file:
        samples = data_parser.parse_samples(samples_file) if samples_file else data_parser.parse_samples_from_vcf(vcf)
        variants = data_parser.parse_variants_of_interest(variants_file)

    for (start, end) in haploblock_boundaries:
        cluster_file = clusters_dir / f"chr{chrom}_{start}-{end}_cluster.tsv"
        try:
            individual2cluster, clusters = data_parser.parse_clusters(cluster_file)
            if verify:
                verify_hash_backends_match(len(clusters), CLUSTER_HASH_LENGTH)
            cluster2hash = generate_cluster_hashes(clusters, use_gpu=use_gpu)

            cluster_hash_file = out / f"cluster_hashes_{start}-{end}.tsv"
            with cluster_hash_file.open("w") as f:
                f.write("CLUSTER\tHASH\n")
                for cl, h in cluster2hash.items():
                    f.write(f"{cl}\t{h}\n")

            variant2hash = None
            if variants_file:
                variant2hash = generate_variant_hashes(variants, vcf, chrom, start, end, samples)

            max_workers = _cpu_workers(threads)
            individual2hash = generate_individual_hashes_parallel(
                individual2cluster, cluster2hash, haploblock2hash, chr_hash,
                variant2hash=variant2hash, max_workers=max_workers
            )

            out_file = out / f"individual_hashes_{start}-{end}.tsv"
            with out_file.open("w") as f:
                f.write("INDIVIDUAL\tHASH\n")
                for ind, h in individual2hash.items():
                    f.write(f"{ind}\t{h}\n")

        except Exception as e:
            logger.error(f"Failed processing haploblock {start}-{end}: {e}")

    logger.info("All hashes generated successfully")


# ---------------------------------------------------------------------
# Wrapper for pipeline
# ---------------------------------------------------------------------
def run(boundaries_file, clusters_dir, chr, out,
        variants_file=None, vcf=None, samples_file=None,
        threads=None, gpu="auto", gpu_id=None, verify=False):
    run_hashes(
        boundaries_file=pathlib.Path(boundaries_file),
        clusters_dir=pathlib.Path(clusters_dir),
        chrom=str(chr),
        out=pathlib.Path(out),
        variants_file=pathlib.Path(variants_file) if variants_file else None,
        vcf=pathlib.Path(vcf) if vcf else None,
        samples_file=pathlib.Path(samples_file) if samples_file else None,
        threads=threads,
        gpu=gpu,
        gpu_id=gpu_id,
        verify=verify,
    )


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )

    parser = argparse.ArgumentParser(
        prog="step5_variant_hashes",
        description="Generate haploblock, cluster, and individual hashes (Step 5)"
    )
    parser.add_argument('--boundaries_file', type=pathlib.Path, required=True)
    parser.add_argument('--clusters_dir', type=pathlib.Path, required=True)
    parser.add_argument('--chr', type=str, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    parser.add_argument('--variants', type=pathlib.Path, default=None)
    parser.add_argument('--vcf', type=pathlib.Path, default=None)
    parser.add_argument('--samples', type=pathlib.Path, default=None)
    parser.add_argument('--threads', type=int, default=None)
    parser.add_argument(
        '--gpu',
        choices=["auto", "on", "off"],
        default="auto",
        help="auto: use GPU if available (default); on: force GPU, error "
             "if none; off: force CPU",
    )
    parser.add_argument('--gpu_id', type=int, default=0, help="GPU ID (if multiple GPUs)")
    parser.add_argument(
        '--verify',
        action='store_true',
        help="cross-check CPU/GPU hash generation for this run's actual "
             "haploblock and cluster counts before writing output",
    )

    args = parser.parse_args()

    run(
        boundaries_file=args.boundaries_file,
        clusters_dir=args.clusters_dir,
        chr=args.chr,
        out=args.out,
        variants_file=args.variants,
        vcf=args.vcf,
        samples_file=args.samples,
        threads=args.threads,
        gpu=args.gpu,
        gpu_id=args.gpu_id,
        verify=args.verify,
    )
