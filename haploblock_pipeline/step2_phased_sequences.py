#!/usr/bin/env python3
import sys
import os
import logging
import argparse
import pathlib
import subprocess
from multiprocessing import Pool, cpu_count
from functools import partial
from typing import Tuple, Dict, Optional

import numpy as np

import data_parser

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------
# Variant counting
# -------------------------------------------------------------------------
def count_variants(vcf: pathlib.Path) -> Tuple[int, int]:
    """Count variants in CPU mode using bcftools query."""
    result = subprocess.run(
        ["bcftools", "query", "-f", "[ %GT]", str(vcf)],
        capture_output=True,
        text=True,
        check=True
    )

    count_0 = count_1 = 0
    for token in result.stdout.split():
        if token in {".", "./.", ".|."}:
            continue
        if "|" in token:
            a, b = token.split("|")
        elif "/" in token:
            a, b = token.split("/")
        else:
            continue
        count_0 += a == "1"
        count_1 += b == "1"

    return count_0, count_1


def gpu_count_variants(vcf: pathlib.Path, gpu_id=0) -> Tuple[int, int]:
    """Count variants on GPU using CuPy."""
    try:
        import cupy as cp
        cp.cuda.Device(gpu_id).use()
    except Exception as e:
        logger.warning("Cannot select GPU %d: %s. Falling back to CPU.", gpu_id, e)
        return count_variants(vcf)

    result = subprocess.run(
        ["bcftools", "query", "-f", "[ %GT]", str(vcf)],
        capture_output=True, text=True, check=True
    )

    hap0, hap1 = [], []
    for t in result.stdout.split():
        if t in {".", "./.", ".|."}:
            continue
        if "|" in t:
            a,b = t.split("|")
        elif "/" in t:
            a,b = t.split("/")
        else:
            continue
        hap0.append(int(a))
        hap1.append(int(b))

    return int(cp.sum(cp.asarray(hap0))), int(cp.sum(cp.asarray(hap1)))


# -------------------------------------------------------------------------
# VCF processing
# -------------------------------------------------------------------------
def normalize_and_filter_vcf(vcf: pathlib.Path, ref: pathlib.Path, out_dir: pathlib.Path) -> pathlib.Path:
    """Normalize and filter VCF in one pipeline, return filtered BCF path."""
    tmp_dir = out_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    filtered_bcf = tmp_dir / f"{vcf.stem}.norm_flt.bcf"

    with open(filtered_bcf, "wb") as out_f:
        p1 = subprocess.Popen(["bcftools", "norm", "-f", str(ref), "-Ob", str(vcf)], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["bcftools", "filter", "--IndelGap", "5", "-e", "QUAL<40", "-Ob"],
                              stdin=p1.stdout, stdout=out_f)
        p2.communicate()

    # Index the BCF for bcftools consensus
    subprocess.run(["bcftools", "index", "-c", str(filtered_bcf)], check=True)

    return filtered_bcf



# -------------------------------------------------------------------------
# Consensus FASTA generation
# -------------------------------------------------------------------------
def generate_consensus_fasta(ref: pathlib.Path,
                             vcf: pathlib.Path,
                             out_dir: pathlib.Path,
                             sample: str):
    """Generate hap0 and hap1 consensus FASTAs for a sample."""
    tmp_dir = out_dir / "tmp" / "consensus_fasta"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    hap0_fasta = tmp_dir / f"{sample}_{vcf.stem}_hap0.fa"
    hap1_fasta = tmp_dir / f"{sample}_{vcf.stem}_hap1.fa"

    subprocess.run(
        ["bcftools", "consensus", "-H", "1", "-f", str(ref), str(vcf)],
        stdout=open(hap0_fasta, "w"),
        check=True
    )
    subprocess.run(
        ["bcftools", "consensus", "-H", "2", "-f", str(ref), str(vcf)],
        stdout=open(hap1_fasta, "w"),
        check=True
    )

    return hap0_fasta, hap1_fasta


# -------------------------------------------------------------------------
# Worker function for Pool
# -------------------------------------------------------------------------
def _process_sample_for_region(args, gpu=False, gpu_id=0, tmp_dir=None):
    sample, region_vcf, region_fasta, ref, out_dir = args
    logger.info("Processing sample %s on %s", sample, "GPU" if gpu else "CPU")

    try:
        # Extract single sample VCF
        sample_vcf = data_parser.extract_sample_from_vcf(region_vcf, sample, out_dir)

        # Normalize + filter
        filtered_vcf = normalize_and_filter_vcf(sample_vcf, ref, out_dir)

        # Variant counting
        if gpu:
            try:
                count_0, count_1 = gpu_count_variants(filtered_vcf, gpu_id)
            except Exception as e:
                logger.warning("GPU counting failed: %s. Falling back to CPU.", e)
                count_0, count_1 = count_variants(filtered_vcf)
        else:
            count_0, count_1 = count_variants(filtered_vcf)

        # Consensus FASTA
        generate_consensus_fasta(region_fasta, filtered_vcf, out_dir, sample)

        return sample, count_0, count_1

    except Exception as e:
        logger.exception("Error processing sample %s: %s", sample, e)
        return sample, 0, 0


# -------------------------------------------------------------------------
# Write outputs
# -------------------------------------------------------------------------
def variant_counts_to_TSV(haploblock2count: Dict[tuple, Tuple[float, float]], out_dir: pathlib.Path):
    out_file = out_dir / "variant_counts.tsv"
    with open(out_file, "w") as f:
        f.write("START\tEND\tMEAN\tSTDEV\n")
        for (s, e), (mean, stdev) in haploblock2count.items():
            f.write(f"{s}\t{e}\t{mean:.3g}\t{stdev:.3g}\n")


# -------------------------------------------------------------------------
# Main pipeline function
# -------------------------------------------------------------------------
def run_phased_sequences(boundaries_file: pathlib.Path,
                         vcf: pathlib.Path,
                         ref: pathlib.Path,
                         chr_map: pathlib.Path,
                         chrom: str,
                         out_dir: pathlib.Path,
                         samples_file: Optional[pathlib.Path] = None,
                         workers: Optional[int] = None,
                         gpu: bool = False,
                         gpu_id: int = 0):

    logger.info("Parsing haploblock boundaries")
    haploblocks = data_parser.parse_haploblock_boundaries(boundaries_file)
    logger.info("Found %d haploblocks", len(haploblocks))

    samples = data_parser.parse_samples(samples_file) if samples_file else data_parser.parse_samples_from_vcf(vcf)
    if not samples:
        logger.warning("No samples found; proceeding with empty sample list.")

    workers = workers or cpu_count()
    logger.info("Using %d worker(s) (available cores: %d)", workers, cpu_count())
    logger.info("GPU mode: %s (gpu_id=%d)", gpu, gpu_id)

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "tmp"
    (tmp_dir / "consensus_fasta").mkdir(parents=True, exist_ok=True)

    # Prepare work list: all (sample, haploblock)
    work = []
    for start, end in haploblocks:
        region_vcf = data_parser.extract_region_from_vcf(vcf, chrom, chr_map, start, end, out_dir)
        region_fasta = data_parser.extract_region_from_fasta(ref, chrom, start, end, out_dir)
        work.extend([(s, region_vcf, region_fasta, ref, out_dir) for s in samples])

    # Multiprocessing
    if workers == 1:
        results = [_process_sample_for_region(w, gpu=gpu, gpu_id=gpu_id, tmp_dir=tmp_dir) for w in work]
    else:
        with Pool(workers) as pool:
            results = pool.map(partial(_process_sample_for_region, gpu=gpu, gpu_id=gpu_id, tmp_dir=tmp_dir), work)

    # Aggregate variant counts per haploblock
    haploblock2count = {}
    idx = 0
    for start, end in haploblocks:
        counts = []
        for _ in samples:
            _, c0, c1 = results[idx]
            counts.extend([c0, c1])
            idx += 1
        haploblock2count[(start, end)] = (float(np.mean(counts)), float(np.std(counts)) if counts else 0.0)

    variant_counts_to_TSV(haploblock2count, out_dir)
    logger.info(f"Variant counts written to {out_dir}")


# -------------------------------------------------------------------------
# CLI wrapper
# -------------------------------------------------------------------------
def run(boundaries_file, vcf, ref, chr_map, chr, out, samples_file, threads=None, gpu=False, gpu_id=0):
    run_phased_sequences(
        pathlib.Path(boundaries_file),
        pathlib.Path(vcf),
        pathlib.Path(ref),
        pathlib.Path(chr_map),
        str(chr),
        pathlib.Path(out),
        pathlib.Path(samples_file) if samples_file else None,
        threads,
        gpu=gpu,
        gpu_id=gpu_id
    )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )
    logger = logging.getLogger(pathlib.Path(sys.argv[0]).name)

    parser = argparse.ArgumentParser(description="Generate haploblock phased sequences")
    parser.add_argument("--boundaries_file", type=pathlib.Path, required=True)
    parser.add_argument("--vcf", type=pathlib.Path, required=True)
    parser.add_argument("--ref", type=pathlib.Path, required=True)
    parser.add_argument("--chr_map", type=pathlib.Path, required=True)
    parser.add_argument("--chr", type=str, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--samples", type=pathlib.Path, default=None)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu_id", type=int, default=0)

    args = parser.parse_args()

    try:
        run_phased_sequences(
            args.boundaries_file,
            args.vcf,
            args.ref,
            args.chr_map,
            args.chr,
            args.out,
            args.samples,
            args.workers,
            gpu=args.gpu,
            gpu_id=args.gpu_id
        )
    except Exception as e:
        sys.stderr.write(f"ERROR in {pathlib.Path(sys.argv[0]).name}: {repr(e)}\n")
        sys.exit(1)


