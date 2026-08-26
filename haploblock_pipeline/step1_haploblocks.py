#!/usr/bin/env python3
import sys
import logging
import argparse
import pathlib
import data_parser
logger = logging.getLogger(__name__)


def _resolve_gpu_choice(gpu):
    """
    Normalizes the gpu argument to True (force GPU), False (force CPU),
    or None (auto). Accepts "auto"/"on"/"off" strings (used by this
    module's own --gpu CLI flag) as well as plain booleans - main.py's
    YAML config passes gpu as a real bool (gpu: true/false), and a bare
    {"auto": None, "on": True, "off": False}[gpu] dict lookup raises
    KeyError the moment a real True/False comes through instead of one
    of those three strings.
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


def haploblocks_to_tsv(haploblocks, chrom, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    output_file = out_dir / f"haploblock_boundaries_chr{chrom}.tsv"
    with output_file.open("w") as f:
        f.write("START\tEND\n")
        for start, end in haploblocks:
            f.write(f"{start}\t{end}\n")


def run_haploblocks(recombination_file, chrom, out_dir, gpu="auto", verify=False):
    """
    gpu: "auto" (GPU if available and file is large enough, else CPU),
         "on" (force GPU, error if unavailable), "off" (force the
         deterministic CPU path) - or a plain bool (True/False), as
         passed by main.py's YAML config.
    verify: if True, cross-check against the original scipy-based CPU
            implementation and abort before writing output on mismatch.
    """
    logger.info("Parsing recombination file %s (chr %s)", recombination_file, chrom)
    if verify:
        data_parser.verify_accelerated_matches_cpu(recombination_file, chrom)
    gpu_flag = _resolve_gpu_choice(gpu)
    haploblocks = data_parser.parse_recombination_rates_accelerated(
        recombination_file, chrom, gpu=gpu_flag
    )
    haploblocks_to_tsv(haploblocks, chrom, out_dir)
    logger.info("Wrote haploblock boundaries to %s", out_dir)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main():
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.INFO,
    )
    parser = argparse.ArgumentParser(
        description="Generate haploblock boundaries from recombination maps"
    )
    parser.add_argument(
        "--recombination_file",
        type=pathlib.Path,
        required=True,
    )
    parser.add_argument("--chr", required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument(
        "--gpu",
        choices=["auto", "on", "off"],
        default="auto",
        help="auto: use GPU if available and worthwhile (default); "
             "on: force GPU, error if none present; off: force CPU",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="cross-check the accelerated result against the original "
             "scipy-based CPU implementation before writing output",
    )
    args = parser.parse_args()
    try:
        run_haploblocks(args.recombination_file, args.chr, args.out, args.gpu, args.verify)
    except Exception as e:
        logger.exception("Haploblock generation failed")
        sys.exit(1)


# Pipeline alias
def run(recombination_file, chr, out, threads=None, gpu="auto"):
    run_haploblocks(
        pathlib.Path(recombination_file),
        str(chr),
        pathlib.Path(out),
        gpu=gpu,
    )


if __name__ == "__main__":
    main()
