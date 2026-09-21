#!/usr/bin/env python3
"""Layer-1 PXD004242 validation using one run per biological sample.

The script deliberately keeps one technical run for each individual x timepoint.
It reports all-pair and cross-plate metrics separately. A fixed technical
replicate is a useful literal "one run" sensitivity analysis, but if same-donor
pairs never cross plates the cross-plate AUC is undefined and must not be treated
as an external validation equivalent to PXD001064.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

TIME_ORDER = {"-8W": 1, "0W": 2, "4W": 3, "13W": 4, "26W": 5, "39W": 6, "52W": 7}


def read_sample_table(path: Path, cohort: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rows = [r for r in rows if r["cohort"] == cohort]
    result = []
    for row in rows:
        raw_files = row["raw_files"].split(";")
        labels = row["maxquant_labels"].split(";")
        techs = row["technical_replicates"].split(",")
        if not (len(raw_files) == len(labels) == len(techs) == int(row["n_runs"])):
            raise ValueError(f"sample row is internally inconsistent: {row['sample_id']}")
        for raw, label, tech in zip(raw_files, labels, techs):
            row_copy = dict(row)
            row_copy.update({"raw_file": raw, "run": label, "technical_replicate": tech})
            plate = label.split("_", 1)[0]
            row_copy["plate"] = plate
            date_match = re.match(r"(\d{8})_", raw)
            if not date_match:
                raise ValueError(f"cannot parse date from {raw}")
            date_text = date_match.group(1)
            day = datetime.strptime(date_text, "%Y%m%d")
            row_copy["date"] = date_text
            row_copy["day"] = str(day.toordinal())
            result.append(row_copy)
    return result


def select_runs(rows: list[dict[str, str]], selection: str, technical_replicate: int | None, seed: int) -> pd.DataFrame:
    by_sample: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_sample.setdefault(row["sample_id"], []).append(row)
    selected = []
    rng = np.random.default_rng(seed)
    for sample_id in sorted(by_sample):
        candidates = by_sample[sample_id]
        if selection == "fixed":
            matches = [r for r in candidates if int(r["technical_replicate"]) == int(technical_replicate)]
            if len(matches) != 1:
                raise ValueError(f"expected one fixed replicate for {sample_id}, got {len(matches)}")
            chosen = matches[0]
        elif selection == "rotate_by_time":
            ordinal = TIME_ORDER[candidates[0]["sampling_time"]]
            target = ((ordinal - 1) % len(candidates)) + 1
            matches = [r for r in candidates if int(r["technical_replicate"]) == target]
            if len(matches) != 1:
                raise ValueError(f"expected one rotated replicate for {sample_id}, got {len(matches)}")
            chosen = matches[0]
        elif selection == "random":
            chosen = candidates[int(rng.integers(0, len(candidates)))]
        else:
            raise ValueError(selection)
        selected.append(chosen)
    meta = pd.DataFrame(selected).set_index("run")
    if meta.index.duplicated().any():
        raise ValueError("selected run labels are duplicated")
    return meta


def pair_table(scores: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    names = list(scores.index)
    records = []
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            score = scores.loc[left, right]
            if not np.isfinite(score):
                continue
            records.append({
                "left": left,
                "right": right,
                "score": float(score),
                "same_donor": bool(metadata.loc[left, "individual"] == metadata.loc[right, "individual"]),
                "same_plate": bool(metadata.loc[left, "plate"] == metadata.loc[right, "plate"]),
                "date_gap_days": abs(int(metadata.loc[left, "day"]) - int(metadata.loc[right, "day"])),
            })
    return pd.DataFrame.from_records(records)


def metric(table: pd.DataFrame) -> dict[str, object]:
    if table.empty:
        return {"n_pairs": 0, "n_same_donor": 0, "n_different_donor": 0, "auc": None}
    same = table["same_donor"].astype(bool)
    scores_same = table.loc[same, "score"]
    scores_diff = table.loc[~same, "score"]
    result: dict[str, object] = {
        "n_pairs": int(len(table)),
        "n_same_donor": int(same.sum()),
        "n_different_donor": int((~same).sum()),
        "same_donor_median": float(scores_same.median()) if len(scores_same) else None,
        "different_donor_median": float(scores_diff.median()) if len(scores_diff) else None,
        "same_donor_same_plate": int(table.loc[same, "same_plate"].sum()) if len(scores_same) else 0,
        "same_donor_cross_plate": int((~table.loc[same, "same_plate"]).sum()) if len(scores_same) else 0,
    }
    if len(scores_same) and len(scores_diff):
        result["auc"] = float(roc_auc_score(same.astype(int), table["score"]))
    else:
        result["auc"] = None
    return result


def analyze(matrix: pd.DataFrame, metadata: pd.DataFrame) -> tuple[dict[str, object], pd.DataFrame]:
    # Same preprocessing contract as PXD001064 pilot.
    positive = matrix.replace(0, np.nan)
    selected = positive.notna().mean(axis=0) >= 0.99
    log_matrix = np.log2(positive.loc[:, selected])
    log_matrix = log_matrix.subtract(log_matrix.median(axis=1), axis=0)
    raw_scores = log_matrix.T.corr(min_periods=30)
    raw_pairs = pair_table(raw_scores, metadata)

    batch_medians = log_matrix.groupby(metadata["plate"]).transform("median")
    residual = log_matrix - batch_medians
    residual = (residual - residual.mean(axis=0)) / residual.std(axis=0).replace(0, np.nan)
    corrected_scores = residual.T.corr(min_periods=30)
    corrected_pairs = pair_table(corrected_scores, metadata)
    corrected_pairs = corrected_pairs.merge(
        raw_pairs[["left", "right", "score"]].rename(columns={"score": "raw_score"}),
        on=["left", "right"], how="left", validate="one_to_one"
    )

    cross = ~corrected_pairs["same_plate"]
    summary = {
        "n_runs": int(len(metadata)),
        "n_donors": int(metadata["individual"].nunique()),
        "n_plates": int(metadata["plate"].nunique()),
        "plates": sorted(metadata["plate"].unique().tolist()),
        "protein_groups_selected": int(selected.sum()),
        "raw_all_pairs": metric(raw_pairs),
        "corrected_all_pairs": metric(corrected_pairs),
        "raw_cross_plate_pairs": metric(raw_pairs.loc[cross]),
        "corrected_cross_plate_pairs": metric(corrected_pairs.loc[cross]),
        "same_donor_cross_plate_available": bool(corrected_pairs.loc[corrected_pairs["same_donor"], "same_plate"].eq(False).any()),
        "same_donor_cross_plate_pairs": int((corrected_pairs["same_donor"] & ~corrected_pairs["same_plate"]).sum()),
    }
    return summary, corrected_pairs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--sample-table", type=Path, required=True)
    ap.add_argument("--selection", choices=["fixed", "rotate_by_time", "random"], default="fixed")
    ap.add_argument("--technical-replicate", type=int, choices=[1, 2, 3, 4])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-json", type=Path, required=True)
    ap.add_argument("--out-pairs", type=Path)
    args = ap.parse_args()
    if args.selection == "fixed" and args.technical_replicate is None:
        ap.error("--technical-replicate is required for --selection fixed")

    matrix = pd.read_csv(args.matrix, sep="\t", index_col=0)
    rows = read_sample_table(args.sample_table, "longitudinal")
    metadata = select_runs(rows, args.selection, args.technical_replicate, args.seed)
    missing = [r for r in metadata.index if r not in matrix.index]
    if missing:
        raise SystemExit(f"selected runs missing from matrix: {missing[:5]} ({len(missing)})")
    matrix = matrix.loc[metadata.index]
    summary, pairs = analyze(matrix, metadata)
    summary.update({
        "selection": args.selection,
        "technical_replicate": args.technical_replicate,
        "seed": args.seed if args.selection == "random" else None,
        "cohort": "longitudinal_only",
        "interpretation": "LOPO/plate transfer is not donor-held-out because donor timepoints share a plate family; cross-plate positives are required for a PXD001064-comparable AUC.",
    })
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.out_pairs:
        args.out_pairs.parent.mkdir(parents=True, exist_ok=True)
        pairs.to_csv(args.out_pairs, sep="\t", index=False)
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
