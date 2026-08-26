#!/usr/bin/env python3
import sys
import os
import yaml
import argparse
import subprocess
from pathlib import Path
import re
import step1_haploblocks
import step2_phased_sequences
import step3_merge_fasta
import step4_clusters
import step5_variant_hashes
from utils.logging import setup_logger

def load_config(config_path):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def main():
    # ----------------------------
    # Argument parsing
    # ----------------------------
    parser = argparse.ArgumentParser(description="Run the Haploblocks Pipeline")
    parser.add_argument("--config", type=Path, required=True,
                        help="Path to YAML configuration file")
    parser.add_argument("--step", type=str,
                        help="Pipeline step to run ('all' or select one between 1 and 5)")
    parser.add_argument("--threads", type=int, help="Number of CPU threads")
    args = parser.parse_args()

    logger = setup_logger()
    logger.info(f"Loading configuration from {args.config}")
    cfg = load_config(args.config)

    # ----------------------------
    # GPU detection and CuPy install
    # ----------------------------
    gpu_available, detected_gpu_id = False, None
    try:
        subprocess.run(["nvidia-smi"], capture_output=True, check=True)
        gpu_available = True
        detected_gpu_id = 0
    except (subprocess.CalledProcessError, FileNotFoundError):
        gpu_available = False
        detected_gpu_id = None

    use_gpu = cfg["pipeline"].get("gpu", gpu_available)
    gpu_id = cfg["pipeline"].get("gpu_id", detected_gpu_id if gpu_available else None)

    if use_gpu and gpu_id is None:
        logger.warning("GPU requested but not detected, falling back to CPU")
        use_gpu = False

    logger.info(f"GPU enabled: {use_gpu} (gpu_id={gpu_id})")


    # ----------------------------
    # Threads configuration
    # ----------------------------
    step = args.step or cfg["pipeline"]["step"]
    threads = args.threads if args.threads else cfg["pipeline"]["threads"]
    if threads == "auto":
        threads = max(1, (os.cpu_count() or 2) - 1)

    logger.info(f"Starting Haploblocks pipeline (Step: {step}, Threads: {threads})")

    try:
        # ---- STEP 1 ---------------------------------------------------------
        if step in ["1", "all"]:
            step1_haploblocks.run(
                recombination_file=cfg["data"]["recombination_file"],
                chr=cfg["chromosome"]["number"],
                out=cfg["outputs"]["out_dir"],
                threads=threads,
                gpu=use_gpu
            )

        # ---- STEP 2 ---------------------------------------------------------
        if step in ["2", "all"]:
            samples_file = cfg["data"].get("samples_file")
            samples_file = Path(samples_file) if samples_file else None
            step2_phased_sequences.run(
                boundaries_file=cfg["outputs"]["boundaries_file"],
                vcf=cfg["data"]["vcf"],
                ref=cfg["data"]["ref"],
                chr_map=cfg["data"]["chr_map"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["step2_out"]),
                samples_file=samples_file,
                threads=threads,
                gpu=use_gpu,
                gpu_id=gpu_id
            )

        # ---- STEP 3 ---------------------------------------------------------
        if step in ["3", "all"]:
            step3_merge_fasta.run(
                input_dir=Path(cfg["outputs"]["step3_input"]),
                output_dir=Path(cfg["outputs"]["step3_output"]),
                threads=threads,
                gpu=use_gpu
            )

        # ---- STEP 4 ---------------------------------------------------------
        if step in ["4", "all"]:
            step4_clusters.run(
                boundaries_file=cfg["outputs"]["boundaries_file"],
                merged_consensus_dir=cfg["outputs"]["merged_consensus_dir"],
                variant_counts=cfg["outputs"]["variant_counts"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["out_dir"]),
                threads=threads
            )

        # ---- STEP 5 ---------------------------------------------------------
        if step in ["5", "all"]:
            variants_file = cfg["data"].get("variants")
            variants_file = Path(variants_file) if variants_file else None
            vcf = Path(cfg["data"]["vcf"]) if variants_file else None
            samples_file = cfg["data"].get("samples_file")
            samples_file = Path(samples_file) if samples_file else None

            step5_variant_hashes.run(
                boundaries_file=cfg["outputs"]["boundaries_file"],
                clusters_dir=cfg["outputs"]["clusters_dir"],
                chr=cfg["chromosome"]["number"],
                out=Path(cfg["outputs"]["out_dir"]),
                variants_file=variants_file,
                vcf=vcf,
                samples_file=samples_file if samples_file else None,
                threads=threads,
                gpu=use_gpu,
                gpu_id=gpu_id
            )

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)

    logger.info("Pipeline finished successfully!")


if __name__ == "__main__":
    main()
