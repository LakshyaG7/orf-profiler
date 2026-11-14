"""Command-line interface: ``orf-profiler scan`` and ``orf-profiler plot``."""

from __future__ import annotations

import argparse
import sys
import time

from . import __version__
from .report import make_figures
from .scan import scan_fasta


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="orf-profiler",
        description="Find and characterize canonical and non-canonical ORFs in transcript sequences.",
    )
    parser.add_argument("--version", action="version", version=f"orf-profiler {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    scan = commands.add_parser("scan", help="scan a transcript FASTA file for ORFs")
    scan.add_argument("fasta", help="transcript FASTA file (.fa or .fa.gz), e.g. from GENCODE")
    scan.add_argument("-o", "--out", required=True, help="output folder")
    scan.add_argument(
        "--min-aa", type=int, default=10, help="shortest ORF to report, in amino acids (default: 10)"
    )
    scan.add_argument(
        "--max-aa", type=int, default=None, help="longest ORF to report, in amino acids (default: no limit)"
    )
    scan.add_argument(
        "--near-cognate",
        action="store_true",
        help="also report ORFs starting at CTG, GTG or TTG upstream of the first ATG",
    )
    scan.add_argument("--threads", type=int, default=1, help="worker processes (default: 1)")
    scan.add_argument("--gzip", action="store_true", help="compress orfs.tsv and peptides.fa")
    scan.add_argument("--no-peptides", action="store_true", help="do not write peptides.fa")

    plot = commands.add_parser("plot", help="draw summary figures from scan output")
    plot.add_argument("results", nargs="+", help="one or more scan output folders (their counts are added up)")
    plot.add_argument("-o", "--out", required=True, help="folder for the figures")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "scan":
        if args.min_aa < 1:
            print("error: --min-aa must be at least 1", file=sys.stderr)
            return 2
        if args.max_aa is not None and args.max_aa < args.min_aa:
            print("error: --max-aa must not be smaller than --min-aa", file=sys.stderr)
            return 2
        started = time.time()
        info = scan_fasta(
            args.fasta,
            args.out,
            min_aa=args.min_aa,
            max_aa=args.max_aa,
            near_cognate=args.near_cognate,
            threads=max(1, args.threads),
            compress=args.gzip,
            write_peptides=not args.no_peptides,
        )
        check = info["cds_check"]
        print(f"Scanned {info['transcripts']:,} transcripts in {time.time() - started:.0f} s")
        print(f"Found {info['orfs']:,} ORFs; results written to {args.out}")
        if check["complete_cds"]:
            found = check["recovered_exact"] + check["recovered_as_extension"]
            print(
                f"Annotated CDS recovered for {found:,} of {check['complete_cds']:,} "
                f"complete CDSs ({check['recovered_exact']:,} exactly, "
                f"{check['recovered_as_extension']:,} as an N-terminal extension)"
            )
        return 0

    if args.command == "plot":
        for path in make_figures(args.results, args.out):
            print(f"Wrote {path}")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
