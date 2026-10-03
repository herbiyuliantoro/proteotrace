#!/usr/bin/env python3
"""Third cohort: PXD062541, run through the shipped tool at its defaults.

PXD062541 (Otton lab, Cruzeiro do Sul University; Orbitrap Exploris 480,
DIA-NN): 28 obese adults, plasma before and after epicatechin
supplementation, two samples each. A different lab, instrument and software
from both validation cohorts, so it is the first cohort the method was not
developed on.

Rules fixed before looking at the results (decided 2026-10-02):

* Primary result = the tool exactly as `proteotrace check` runs it, at its
  defaults, whatever it shows. Other coverage thresholds are reported as a
  sensitivity analysis only.
* No batch information is deposited, so the cohort runs in single-batch mode:
  no batch correction, all pairs are scored. This is a stated limitation.
* The DIA-NN report names one run `WB_1`; the SDRF has no `WB_1` but a
  `W1_B.raw` that is otherwise absent from the report. `WB_1` is therefore
  read as W1_B (Women 1, before). This is an inference from file names, not
  a confirmation from the submitters.
* Planted swaps use the same planter and accounting as validate_tool.py,
  with 6 pairs per repeat (12 of 28 people) because the cohort is small.

    python analysis_pxd062541.py
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_primary_estimand import cluster_bootstrap, operating
from proteotrace.core import normalise, pair_table, similarity, threshold_from_pairs
from proteotrace.io import read_matrix
from proteotrace.report import VERDICTS, run_report
from validate_tool import account, plant_swaps

RENAMED = {"WB_1": "W1_B"}            # report run name -> SDRF sample name
COVERAGES = (0.99, 0.95, 0.90, 0.80)  # 0.99 is the tool default


def samples_table(report_runs, sdrf: Path) -> pd.DataFrame:
    """run (as named in the report), donor and visit, checked against the SDRF."""
    with sdrf.open(encoding="utf-8") as handle:
        deposited = {r["comment[data file]"].removesuffix(".raw")
                     for r in csv.DictReader(handle, delimiter="\t")}
    rows = []
    for run in sorted(report_runs):
        sample = RENAMED.get(run, run)
        if sample not in deposited:
            raise ValueError(f"report run {run!r} has no matching SDRF raw file")
        donor, visit = sample.split("_")
        rows.append({"run": run, "donor": donor, "visit": visit, "batch": "single_batch"})
    table = pd.DataFrame(rows)
    unmatched = deposited - {RENAMED.get(r, r) for r in report_runs}
    if unmatched:
        raise ValueError(f"SDRF samples missing from the report: {sorted(unmatched)}")
    return table


def verdicts(scores, samples, margin=0.20, specificity=0.99):
    """Exactly the CLI's single-batch path."""
    annotation = samples.set_index("run")
    pairs = pair_table(scores, annotation["donor"], annotation["batch"])
    match, rival = threshold_from_pairs(pairs, specificity, cross_batch_only=False)
    report = run_report(scores, annotation["donor"], annotation["batch"], match, margin,
                        cross_batch_only=False, rival_threshold=rival)
    return report.set_index("run").loc[samples["run"]]


def pair_metrics(scores: pd.DataFrame, samples: pd.DataFrame, n_boot: int) -> dict:
    annotation = samples.set_index("run")
    pairs = pair_table(scores, annotation["donor"], annotation["batch"])
    labels = pairs["same_donor"].to_numpy()
    values = pairs["score"].to_numpy()
    out = operating(labels, values)
    out.update({"n_positive": int(labels.sum()), "n_negative": int((~labels).sum())})
    if n_boot:
        out["ci"] = cluster_bootstrap(labels, values,
                                      annotation.loc[pairs["left"], "donor"].to_numpy(),
                                      annotation.loc[pairs["right"], "donor"].to_numpy(), n_boot)
    return out


def raw_scores(matrix: pd.DataFrame, min_coverage: float) -> pd.DataFrame:
    """Ablation baseline: the tool's feature rule, log2 and run-median centring only."""
    positive = matrix.replace(0, np.nan)
    keep = positive.notna().mean(axis=0) >= min_coverage
    keep |= positive.isna().sum(axis=0) <= 1
    values = np.log2(positive.loc[:, keep])
    values = values.subtract(values.median(axis=1), axis=0)
    return values.T.corr(min_periods=30)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", type=Path, default=Path("data/PXD062541_report.tsv"))
    ap.add_argument("--sdrf", type=Path, default=Path("data/PXD062541.sdrf.tsv"))
    ap.add_argument("--repeats", type=int, default=30)
    ap.add_argument("--pairs", type=int, default=6)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--out", type=Path, default=Path("data/PXD062541_primary.json"))
    ap.add_argument("--samples-out", type=Path, default=Path("examples/PXD062541_samples.tsv"))
    ap.add_argument("--pairs-out", type=Path, default=Path("data/PXD062541_pairs.tsv"))
    args = ap.parse_args()

    matrix, fmt = read_matrix(args.report)
    samples = samples_table(matrix.index, args.sdrf)
    matrix = matrix.loc[samples["run"]]
    args.samples_out.parent.mkdir(parents=True, exist_ok=True)
    samples[["run", "donor", "batch"]].to_csv(args.samples_out, sep="\t", index=False)
    annotation = samples.set_index("run")

    result = {"cohort": "PXD062541", "format": fmt, "runs": int(len(samples)),
              "donors": int(samples["donor"].nunique()), "renamed_runs": RENAMED,
              "batch_information": "none deposited; single-batch mode", "coverage": {}}
    for coverage in COVERAGES:
        normalised = normalise(matrix, annotation["batch"], coverage)
        scores = similarity(normalised)
        entry = {"features": normalised.n_features,
                 "default": pair_metrics(scores, samples, args.n_boot if coverage == 0.99 else 0),
                 "raw": pair_metrics(raw_scores(matrix, coverage), samples, 0)}
        if coverage == 0.99:
            # pair scores at the defaults, so the supplementary figure builds
            # without the 106 MB report (which is not tracked)
            pairs = pair_table(scores, annotation["donor"], annotation["batch"])
            raw = raw_scores(matrix, coverage)
            pairs["score_raw"] = [raw.loc[a, b] for a, b in zip(pairs["left"], pairs["right"])]
            pairs.rename(columns={"score": "score_default"})[
                ["left", "right", "same_donor", "score_raw", "score_default"]
            ].to_csv(args.pairs_out, sep="\t", index=False, float_format="%.6f")
            clean = verdicts(scores, samples)
            counts = clean["verdict"].value_counts()
            entry["clean_cohort_verdicts"] = {v: int(counts.get(v, 0)) for v in VERDICTS}
            entry["flagged_runs"] = clean.loc[clean["verdict"] == "swap_suspected"].index.tolist()
            rows = []
            for r in range(args.repeats):
                planted, swapped = plant_swaps(samples, args.pairs, np.random.default_rng(100 + r))
                rows.append(account(verdicts(scores, planted), samples, swapped))
            table = pd.DataFrame(rows)
            entry["planted_swaps"] = {k: {"mean": float(table[k].mean()),
                                          "mc_se": float(table[k].std(ddof=1) / np.sqrt(len(table)))}
                                      for k in table.columns}
        result["coverage"][f"{coverage:.2f}"] = entry
        d, w = entry["default"], entry["raw"]
        print(f"coverage {coverage:.2f}: features={entry['features']:4d}  "
              f"raw AUC {w['auc']:.4f} | default AUC {d['auc']:.4f} sens95 {d['sens_95']:.3f} "
              f"sens99 {d['sens_99']:.3f}  (pos {d['n_positive']}, neg {d['n_negative']})")

    primary = result["coverage"]["0.99"]
    print("primary CI", primary["default"]["ci"])
    print("clean cohort", primary["clean_cohort_verdicts"], "flagged", primary["flagged_runs"])
    print("planted swaps", {k: round(v["mean"], 3) for k, v in primary["planted_swaps"].items()})
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.out, args.samples_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
