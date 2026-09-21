"""Command line entry point: `proteotrace check`."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .core import normalise, pair_table, similarity, threshold_from_pairs
from .io import read_matrix, read_samples
from .report import run_report, summarise


def _check(args: argparse.Namespace) -> int:
    matrix, fmt = read_matrix(args.matrix, args.format)
    samples = read_samples(args.samples)

    shared = matrix.index.intersection(samples.index)
    if len(shared) < 3:
        unmatched = list(matrix.index[:3])
        raise SystemExit(
            f"only {len(shared)} runs appear in both files. The 'run' column must "
            f"match the row names of the matrix exactly; the matrix starts with "
            f"{unmatched}.")
    dropped = len(matrix.index) - len(shared)
    matrix = matrix.loc[shared]
    samples = samples.loc[shared]

    single_batch = samples["batch"].nunique() == 1
    if single_batch:
        print("note: every run is in one batch, so no batch offset is removed. "
              "Scores will be optimistic relative to a multi-batch study.",
              file=sys.stderr)

    normalised = normalise(matrix, samples["batch"], args.min_coverage)
    scores = similarity(normalised, args.min_overlap)
    pairs = pair_table(scores, samples["donor"], samples["batch"])
    rivals_per_run = args.rival_adjust
    threshold, rival_threshold = threshold_from_pairs(
        pairs, args.specificity, cross_batch_only=not single_batch,
        rivals_per_run=rivals_per_run)
    report = run_report(scores, samples["donor"], samples["batch"], threshold,
                        args.margin, cross_batch_only=not single_batch,
                        rival_threshold=rival_threshold)
    summary = summarise(report, threshold, normalised.n_features, args.specificity,
                        rival_threshold)
    summary.update({"input_format": fmt, "runs_dropped_no_annotation": int(dropped),
                    "features_dropped_low_coverage": normalised.dropped_features,
                    "proteotrace_version": __version__})

    args.out.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.out / "run_report.tsv", sep="\t", index=False)
    pairs.to_csv(args.out / "pair_scores.tsv", sep="\t", index=False)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                           encoding="utf-8")

    counts = summary["verdicts"]
    print(f"\n{summary['runs']} runs · {summary['donors']} donors · "
          f"{summary['features_used']} features")
    print(f"match threshold {threshold:.4f} at {args.specificity:.0%} specificity")
    if rivals_per_run > 1:
        print(f"rival threshold {rival_threshold:.4f}, spread over "
              f"{rivals_per_run} comparisons per run")
    print(f"swap margin {args.margin:.2f}\n")
    for verdict in ("swap_suspected", "review", "no_reference", "ok"):
        print(f"  {verdict:<16} {counts[verdict]:>5}")
    flagged = report.loc[report["verdict"] == "swap_suspected"]
    if len(flagged):
        print("\nRuns whose best match is a different donor:\n")
        print(flagged[["run", "claimed_donor", "score_to_claimed_donor",
                       "best_other_donor", "score_to_best_other"]]
              .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print("\nInvestigate these before trusting the affected donors. Note that "
              "identical twins are not separable by this method.")
    print(f"\nwritten to {args.out}/")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="proteotrace",
        description="Check whether the runs in a cohort match the donors they "
                    "are labelled with, using the quantitative data alone.")
    parser.add_argument("--version", action="version", version=f"proteotrace {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    check = sub.add_parser("check", help="score a cohort and flag inconsistent runs")
    check.add_argument("--matrix", type=Path, required=True,
                       help="quantitative table: generic runs x features, a MaxQuant "
                            "proteinGroups.txt (or zip), or a DIA-NN report")
    check.add_argument("--samples", type=Path, required=True,
                       help="annotation with columns run, donor and batch")
    check.add_argument("--out", type=Path, default=Path("proteotrace_out"),
                       help="output directory (default: proteotrace_out)")
    check.add_argument("--format", choices=["auto", "generic", "maxquant", "diann"],
                       default="auto", help="matrix format (default: auto)")
    check.add_argument("--min-coverage", type=float, default=0.99,
                       help="keep features observed in at least this fraction of "
                            "runs (default: 0.99)")
    check.add_argument("--specificity", type=float, default=0.99,
                       help="false-flag rate to calibrate the threshold at "
                            "(default: 0.99)")
    check.add_argument("--margin", type=float, default=0.20,
                       help="how far another donor must beat the claimed donor "
                            "before a swap is called. This is the main knob. "
                            "Measured on our two cohorts: 0.05 catches more "
                            "swaps but flags 16-22%% of clean runs, 0.20 catches "
                            "81-93%% of swaps and flags 4-13%% (default: 0.20)")
    check.add_argument("--rival-adjust", type=int, default=1,
                       help="spread the error budget over this many rival "
                            "comparisons per run. Raising it cuts false alarms "
                            "but costs far more sensitivity than tightening "
                            "--margin does, so it is off by default (default: 1)")
    check.add_argument("--min-overlap", type=int, default=30,
                       help="minimum shared features for a pair to be scored "
                            "(default: 30)")
    check.set_defaults(func=_check)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except ValueError as error:
        raise SystemExit(f"proteotrace: {error}")


if __name__ == "__main__":
    raise SystemExit(main())
