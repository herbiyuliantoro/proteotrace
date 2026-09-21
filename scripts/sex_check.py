#!/usr/bin/env python3
"""Reproducible, prespecified plasma sex-signal feasibility check.

The metadata labels are read from the deposited sample_metadata file.  The
report is streamed, and only SHBG/PZP are used so the experiment is explicit
and auditable rather than a post-hoc high-dimensional classifier.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SAMPLE_RE = re.compile(r"([MW])(\d+)_([AB])")
MARKERS = ("SHBG", "PZP")


def parse_number(value: str) -> float:
    value = (value or "").strip().strip('"')
    if not value or value.lower() in {"nan", "na", "n/a", "null"}:
        return float("nan")
    # DIA-NN exports decimal commas in this report.  Scientific notation uses
    # a comma only as the decimal separator, not as a list separator here.
    if "," in value and "." not in value:
        value = value.replace(",", ".")
    try:
        number = float(value)
    except ValueError:
        return float("nan")
    return number if math.isfinite(number) and number > 0 else float("nan")


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    # The PRIDE sample metadata is UTF-16 with a BOM.  utf-16 handles both
    # little- and big-endian variants and removes the BOM for DictReader.
    with path.open("r", encoding="utf-16", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        result: dict[str, dict[str, str]] = {}
        for row in rows:
            sample_text = (row.get("Sample_ID") or row.get("Sample ID") or "").strip()
            match = SAMPLE_RE.search(sample_text)
            if not match:
                raise ValueError(f"cannot parse participant/timepoint from {sample_text!r}")
            sex, person_number, timepoint = match.groups()
            sample = f"{sex}{person_number}_{timepoint}"
            result[sample] = {
                "participant": f"{sex}{person_number}",
                "sex": sex,
                "timepoint": timepoint,
                "metadata_label_source": "sample_metadata_description",
                "description": (row.get("Description") or "").strip(),
            }
    return result


def read_marker_values(path: Path, value_column: str) -> dict[tuple[str, str], list[float]]:
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"File.Name", "Genes", value_column}
        missing = required.difference(rows.fieldnames or [])
        if missing:
            raise ValueError(f"report lacks columns: {sorted(missing)}")
        for row in rows:
            gene = (row.get("Genes") or "").strip()
            if gene not in MARKERS:
                continue
            sample = (row.get("File.Name") or "").strip()
            value = parse_number(row.get(value_column, ""))
            if math.isfinite(value):
                values[(sample, gene)].append(value)
    return values


def make_table(metadata: dict[str, dict[str, str]], values: dict[tuple[str, str], list[float]]) -> pd.DataFrame:
    rows = []
    samples = sorted(metadata, key=lambda x: (metadata[x]["participant"], metadata[x]["timepoint"]))
    for sample in samples:
        info = metadata[sample]
        row = {"sample": sample, **info}
        for marker in MARKERS:
            raw = values.get((sample, marker), [])
            median = float(np.median(raw)) if raw else float("nan")
            row[marker] = median
            row[f"log10_{marker}"] = math.log10(median) if median > 0 else float("nan")
            row[f"n_{marker}_rows"] = len(raw)
        rows.append(row)
    return pd.DataFrame(rows)


def evaluate_loso(table: pd.DataFrame, features: list[str]) -> tuple[dict, pd.DataFrame]:
    cols = [f"log10_{x}" for x in features]
    complete = table.dropna(subset=cols).copy()
    predictions = []
    for participant in sorted(complete["participant"].unique()):
        test = complete[complete["participant"] == participant]
        train = complete[complete["participant"] != participant]
        # The baseline table has one row/person; this guard keeps the function
        # correct if a future dataset includes repeated baseline records.
        if train["sex"].nunique() < 2:
            continue
        model = make_pipeline(
            SimpleImputer(strategy="median"),
            StandardScaler(),
            LogisticRegression(solver="liblinear", random_state=1),
        )
        model.fit(train[cols], (train["sex"] == "M").astype(int))
        probability = model.predict_proba(test[cols])[:, 1]
        for (_, row), prob in zip(test.iterrows(), probability):
            predictions.append({
                "sample": row["sample"],
                "participant": row["participant"],
                "timepoint": row["timepoint"],
                "sex": row["sex"],
                "prediction_probability_male": float(prob),
                "prediction": "M" if prob >= 0.5 else "W",
                "mismatch": bool((prob >= 0.5) != (row["sex"] == "M")),
            })
    pred = pd.DataFrame(predictions)
    if pred.empty:
        metrics = {"n": 0, "participants": 0, "accuracy": None, "auc": None, "confusion_matrix_W_M": None}
    else:
        y = (pred["sex"] == "M").astype(int)
        p = pred["prediction_probability_male"].to_numpy()
        metrics = {
            "n": int(len(pred)),
            "participants": int(pred["participant"].nunique()),
            "accuracy": float(accuracy_score(y, (p >= 0.5).astype(int))),
            "auc": float(roc_auc_score(y, p)) if y.nunique() == 2 else None,
            "confusion_matrix_W_M": confusion_matrix(y, (p >= 0.5).astype(int), labels=[0, 1]).tolist(),
        }
    return metrics, pred


def fit_baseline_predict_after(table: pd.DataFrame, features: list[str]) -> tuple[dict, pd.DataFrame]:
    cols = [f"log10_{x}" for x in features]
    baseline = table[(table["timepoint"] == "B")].dropna(subset=cols)
    after = table[(table["timepoint"] == "A")].dropna(subset=cols)
    if baseline.empty or after.empty or baseline["sex"].nunique() < 2:
        return {"n": int(len(after)), "accuracy": None, "confusion_matrix_W_M": None}, pd.DataFrame()
    model = make_pipeline(
        SimpleImputer(strategy="median"),
        StandardScaler(),
        LogisticRegression(solver="liblinear", random_state=1),
    )
    model.fit(baseline[cols], (baseline["sex"] == "M").astype(int))
    probability = model.predict_proba(after[cols])[:, 1]
    pred = after[["sample", "participant", "timepoint", "sex"]].copy()
    pred["prediction_probability_male"] = probability
    pred["prediction"] = np.where(probability >= 0.5, "M", "W")
    pred["mismatch"] = pred["prediction"] != pred["sex"]
    metrics = {
        "n": int(len(pred)),
        "accuracy": float(np.mean(~pred["mismatch"])),
        "confusion_matrix_W_M": confusion_matrix((pred["sex"] == "M").astype(int), (probability >= 0.5).astype(int), labels=[0, 1]).tolist(),
    }
    return metrics, pred


def json_default(value):
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    raise TypeError(type(value).__name__)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--value-column", default="PG.Normalised")
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-tsv", type=Path, required=True)
    args = parser.parse_args()

    metadata = read_metadata(args.metadata)
    values = read_marker_values(args.report, args.value_column)
    table = make_table(metadata, values)
    baseline = table[table["timepoint"] == "B"]
    after = table[table["timepoint"] == "A"]
    results = {
        "dataset": "PXD062541",
        "report": str(args.report),
        "metadata": str(args.metadata),
        "value_column": args.value_column,
        "label_source": "sample_metadata_description (W/M in Sample_ID and Description)",
        "markers": list(MARKERS),
        "samples": int(len(table)),
        "participants": int(table["participant"].nunique()),
        "baseline_samples": int(len(baseline)),
        "after_samples": int(len(after)),
        "marker_coverage": {
            marker: {
                "baseline": int(baseline[marker].notna().sum()),
                "after": int(after[marker].notna().sum()),
            }
            for marker in MARKERS
        },
        "two_marker": {},
    }
    output_rows = []
    for feature_set in (list(MARKERS), ["PZP"]):
        key = "_".join(feature_set)
        baseline_metrics, baseline_pred = evaluate_loso(baseline, feature_set)
        after_metrics, after_pred = fit_baseline_predict_after(table, feature_set)
        results["two_marker" if len(feature_set) == 2 else "pzp_only"] = {
            "features": feature_set,
            "baseline_complete": int(baseline.dropna(subset=[f"log10_{x}" for x in feature_set]).shape[0]),
            "baseline_loso": baseline_metrics,
            "after_complete": int(after.dropna(subset=[f"log10_{x}" for x in feature_set]).shape[0]),
            "after_using_baseline_model": after_metrics,
        }
        if not baseline_pred.empty:
            baseline_pred["evaluation"] = f"baseline_loso_{key}"
            output_rows.extend(baseline_pred.to_dict("records"))
        if not after_pred.empty:
            after_pred["evaluation"] = f"after_using_baseline_model_{key}"
            output_rows.extend(after_pred.to_dict("records"))

    table.to_csv(args.out_tsv, sep="\t", index=False, float_format="%.8g")
    results["prediction_rows"] = output_rows
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(results, indent=2, default=json_default) + "\n")
    print(json.dumps(results, indent=2, default=json_default))


if __name__ == "__main__":
    main()
