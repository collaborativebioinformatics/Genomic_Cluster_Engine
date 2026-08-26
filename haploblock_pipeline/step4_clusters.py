#!/usr/bin/env python3
import os
import time
import logging
import pathlib
import argparse
import subprocess
import signal
from concurrent.futures import ThreadPoolExecutor, as_completed

import data_parser

logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------
# MMseqs2 params
# ----------------------------------------------------------------------
def calculate_mmseq_params(variant_counts_file: pathlib.Path):

    haploblock2min_id = {}
    haploblock2cov_fraction = {}

    with open(variant_counts_file, "r") as f:
        header = f.readline()

        if not header.startswith("START\t"):
            raise ValueError(
                f"Variant counts file missing header: {header.strip()}"
            )

        for line in f:
            start, end, mean, stdev = line.strip().split("\t")

            start = int(start)
            end = int(end)

            hap_len = end - start

            haploblock2min_id[(start, end)] = (
                1 - (float(mean) / hap_len)
            )

            haploblock2cov_fraction[(start, end)] = (
                1 - (682 / hap_len)
            )

    return haploblock2min_id, haploblock2cov_fraction


# ----------------------------------------------------------------------
# Run clustering per FASTA (STABLE v2 FIXED)
# ----------------------------------------------------------------------
def compute_clusters(
    input_fasta: str,
    out: str,
    min_seq_id: float,
    cov_fraction: float,
    cov_mode: int,
    chrom: str,
    start: int,
    end: int,
    mmseq_threads: int
):

    t0 = time.time()

    input_fasta = str(pathlib.Path(input_fasta).resolve())

    output_prefix = (
        pathlib.Path(out)
        / "clusters"
        / f"chr{chrom}_{start}-{end}"
    )
    output_prefix = str(output_prefix.resolve())

    # ------------------------------------------------------------------
    # isolated MMseqs workspace
    # ------------------------------------------------------------------
    tmp_dir = (
        pathlib.Path(out)
        / "mmseqs_tmp"
        / f"chr{chrom}_{start}_{end}"
    )
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = str(tmp_dir.resolve())

    # ------------------------------------------------------------------
    # skip if done
    # ------------------------------------------------------------------
    if pathlib.Path(f"{output_prefix}_cluster.tsv").exists():
        logger.info("Skipping chr%s:%s-%s", chrom, start, end)
        return True

    # ------------------------------------------------------------------
    # force isolated env
    # ------------------------------------------------------------------
    env = os.environ.copy()
    env["TMPDIR"] = tmp_dir
    env["MMSEQS_TMPDIR"] = tmp_dir
    env["APPTAINER_TMPDIR"] = tmp_dir
    env["SINGULARITY_TMPDIR"] = tmp_dir

    cmd = [
        "mmseqs",
        "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(cov_fraction),
        "--cov-mode",
        str(cov_mode),
        "--remove-tmp-files",
        "1",
        "--threads",
        str(mmseq_threads)
    ]

    logger.debug("Running: %s", " ".join(cmd))

    # ================================================================
    # 🔥 STABLE EXECUTION (PROCESS GROUP CONTROL)
    # ================================================================
    try:

        p = subprocess.Popen(
            cmd,
            env=env,
            preexec_fn=os.setsid
        )

        timeout_sec = 3600  # 1h per haploblock safety cap

        while True:

            if p.poll() is not None:
                break

            if time.time() - t0 > timeout_sec:
                logger.error("TIMEOUT chr%s:%s-%s", chrom, start, end)
                os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                return False

            time.sleep(5)

        returncode = p.poll()

        if returncode != 0:
            logger.error(
                "MMseqs failed chr%s:%s-%s | exit=%s",
                chrom, start, end, returncode
            )
            return False

        runtime = time.time() - t0

        logger.info(
            "Finished chr%s:%s-%s | runtime=%.1fs",
            chrom, start, end, runtime
        )

        return True

    except Exception as e:

        runtime = time.time() - t0

        logger.error(
            "Unexpected error chr%s:%s-%s | runtime=%.1fs | error=%s",
            chrom, start, end, runtime, str(e)
        )

        return False


# ----------------------------------------------------------------------
# Main workflow (unchanged logic)
# ----------------------------------------------------------------------
def run_clusters(
    boundaries_file: pathlib.Path,
    merged_consensus_dir: pathlib.Path,
    variant_counts_file: pathlib.Path,
    chrom: str,
    out_dir: pathlib.Path,
    cov_mode: int,
    threads: int | None,
    max_retries: int = 2
):

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "clusters").mkdir(exist_ok=True)

    haploblock_boundaries = data_parser.parse_haploblock_boundaries(boundaries_file)

    haploblock2min_id, haploblock2cov_fraction = calculate_mmseq_params(variant_counts_file)

    logger.info("Found %d haploblocks", len(haploblock_boundaries))

    cpu_budget = threads or (os.cpu_count() or 1)
    mmseq_threads = 16

    max_workers = max(1, cpu_budget // mmseq_threads)
    max_workers = min(max_workers, 12)

    logger.info(
        "CPU=%d mmseq_threads=%d workers=%d",
        cpu_budget, mmseq_threads, max_workers
    )

    remaining_blocks = sorted(haploblock_boundaries, key=lambda x: x[1] - x[0])

    retries = 0

    while remaining_blocks and retries <= max_retries:

        futures = {}
        failed = []

        with ThreadPoolExecutor(max_workers=max_workers) as ex:

            for start, end in remaining_blocks:

                inp = merged_consensus_dir / f"chr{chrom}_region_{start}-{end}.fa"

                fut = ex.submit(
                    compute_clusters,
                    str(inp),
                    str(out_dir),
                    haploblock2min_id[(start, end)],
                    haploblock2cov_fraction[(start, end)],
                    cov_mode,
                    chrom,
                    start,
                    end,
                    mmseq_threads
                )

                futures[fut] = (start, end)

            for f in as_completed(futures):

                start, end = futures[f]

                try:
                    if not f.result():
                        failed.append((start, end))
                except Exception:
                    failed.append((start, end))

        if failed:
            retries += 1
            remaining_blocks = failed
        else:
            break

    if remaining_blocks:
        logger.error("FAILED BLOCKS:")
        for s, e in remaining_blocks:
            logger.error("%s:%s-%s", chrom, s, e)
    else:
        logger.info("ALL DONE")


# ----------------------------------------------------------------------
# Wrapper for pipeline
# ----------------------------------------------------------------------
def run(boundaries_file, merged_consensus_dir, variant_counts, chr, out,
        threads=None, cov_mode=0, gpu=None, gpu_id=None, max_retries=2):
    """
    gpu/gpu_id accepted-and-ignored: MMseqs2 has no GPU-accelerated
    clustering path (only its search/profile-alignment functions do),
    so this step stays CPU-only regardless of pipeline.gpu in the config.

    """
    run_clusters(
        boundaries_file=pathlib.Path(boundaries_file),
        merged_consensus_dir=pathlib.Path(merged_consensus_dir),
        variant_counts_file=pathlib.Path(variant_counts),
        chrom=str(chr),
        out_dir=pathlib.Path(out),
        cov_mode=cov_mode,
        threads=threads,
        max_retries=max_retries,
    )
