#!/usr/bin/env python3
"""Reconcile PXD004242 PRIDE files, SDRF rows, and MaxQuant run columns.

This is a structure check only. It intentionally does not compute donor
similarity or AUC. The unit it emits for downstream analysis is one biological
sample (individual x sampling time), after technical replicates are identified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

PARTICIPANT_RE = re.compile(r"Plate(?P<plate>[^_]+)_(?P<sample>[^_]+)_(?P<rep>[^.]+)\.raw$", re.I)
REFERENCE_RE = re.compile(r"DoubleTop14_(?P<sample>[FM]\d)_(?P<rep>[^.]+)\.raw$", re.I)
DATE_RE = re.compile(r"^(?P<date>\d{8})_")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_manifest(path: Path):
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    required = {"file_name", "file_category", "file_size_bytes", "sha1", "ftp_url", "api_accession"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"manifest must contain {sorted(required)}")
    return rows


def maxquant_intensity_labels(path: Path):
    with zipfile.ZipFile(path) as z:
        with z.open("proteinGroups.txt") as f:
            header = f.readline().decode("utf-8-sig").rstrip("\r\n").split("\t")
    prefix = "Intensity "
    return [col[len(prefix):] for col in header if col.startswith(prefix)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sdrf", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--maxquant-zip", type=Path)
    ap.add_argument("--out-summary", type=Path, required=True)
    ap.add_argument("--out-samples", type=Path, required=True)
    args = ap.parse_args()

    errors = []
    manifest = read_manifest(args.manifest)
    manifest_by_name = {r["file_name"]: r for r in manifest}
    if len(manifest_by_name) != len(manifest):
        errors.append("manifest contains duplicate file_name values")
    raw_manifest = {r["file_name"]: r for r in manifest if r["file_category"] == "RAW"}

    with args.sdrf.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        sdrf_rows = list(reader)
        fields = set(reader.fieldnames or [])
    required_sdrf = {
        "characteristics[individual]", "characteristics[sampling time]", "comment[data file]",
        "comment[technical replicate]", "comment[fraction identifier]",
        "characteristics[biological replicate]",
    }
    missing_fields = sorted(required_sdrf - fields)
    if missing_fields:
        raise ValueError(f"SDRF is missing columns: {missing_fields}")

    sdrf_files = [r["comment[data file]"].strip() for r in sdrf_rows]
    if len(set(sdrf_files)) != len(sdrf_files):
        errors.append("SDRF contains duplicate comment[data file] values")
    sdrf_raw = set(sdrf_files)
    if sdrf_raw != set(raw_manifest):
        missing_in_manifest = sorted(sdrf_raw - set(raw_manifest))
        missing_in_sdrf = sorted(set(raw_manifest) - sdrf_raw)
        errors.append(
            "SDRF/API raw-file mismatch: "
            f"missing_in_manifest={missing_in_manifest[:5]} ({len(missing_in_manifest)}), "
            f"missing_in_sdrf={missing_in_sdrf[:5]} ({len(missing_in_sdrf)})"
        )

    groups = defaultdict(list)
    participant_ids = set()
    reference_ids = set()
    parsed_runs = []
    for row in sdrf_rows:
        individual = row["characteristics[individual]"].strip()
        time = row["characteristics[sampling time]"].strip()
        raw_name = row["comment[data file]"].strip()
        tech = row["comment[technical replicate]"].strip()
        fraction = row["comment[fraction identifier]"].strip()
        biological = row["characteristics[biological replicate]"].strip()
        m = PARTICIPANT_RE.search(raw_name)
        ref = REFERENCE_RE.search(raw_name)
        if m:
            if not individual.isdigit():
                errors.append(f"numeric run pattern has nonnumeric individual: {raw_name} -> {individual}")
            sample_id = f"participant:{individual}|time:{time}"
            participant_ids.add(individual)
            plate = m.group("plate")
            sample_number = m.group("sample")
            file_rep = m.group("rep")
            maxquant_label = f"{plate}_{sample_number}_{file_rep}"
            cohort_hint = "participant"
        elif ref:
            sample_id = f"reference:{individual}"
            reference_ids.add(individual)
            plate = "DoubleTop14"
            sample_number = ref.group("sample")
            file_rep = ref.group("rep")
            maxquant_label = f"Top14_{sample_number}_{file_rep}"
            cohort_hint = "reference_library"
        else:
            errors.append(f"cannot parse raw filename pattern: {raw_name}")
            plate = sample_number = file_rep = maxquant_label = ""
            sample_id = f"unparsed:{raw_name}"
            cohort_hint = "unparsed"

        if raw_name not in raw_manifest:
            errors.append(f"SDRF raw file is absent from API manifest: {raw_name}")
        if m:
            # In participant filenames the final suffix is the timepoint
            # ordinal (1..7), while SDRF technical replicate identifies the
            # repeated plate/run (1..4). The two are deliberately different.
            time_ordinal = {"-8W": "1", "0W": "2", "4W": "3", "13W": "4", "26W": "5", "39W": "6", "52W": "7"}
            expected_ordinal = time_ordinal.get(time)
            if expected_ordinal != file_rep:
                errors.append(f"timepoint suffix mismatch for {raw_name}: filename={file_rep}, sampling_time={time!r}, expected={expected_ordinal}")
        elif ref:
            if file_rep != tech:
                errors.append(f"technical replicate mismatch for reference {raw_name}: filename={file_rep}, SDRF={tech}")
        if fraction != "1":
            errors.append(f"unexpected fraction for {raw_name}: {fraction!r}")
        if biological != "1":
            errors.append(f"unexpected biological replicate for {raw_name}: {biological!r}")
        date_match = DATE_RE.match(raw_name)
        run_date = date_match.group("date") if date_match else ""
        run = {
            "sample_id": sample_id,
            "individual": individual,
            "sampling_time": time,
            "raw_file": raw_name,
            "technical_replicate": tech,
            "plate": plate,
            "run_date": run_date,
            "maxquant_label": maxquant_label,
            "cohort_hint": cohort_hint,
        }
        groups[sample_id].append(run)
        parsed_runs.append(run)

    participant_timepoints = defaultdict(set)
    for r in parsed_runs:
        if r["cohort_hint"] == "participant":
            participant_timepoints[r["individual"]].add(r["sampling_time"])
    longitudinal_ids = sorted(i for i, times in participant_timepoints.items() if len(times) >= 6)
    baseline_only_ids = sorted(i for i, times in participant_timepoints.items() if len(times) == 2)
    other_ids = sorted(i for i, times in participant_timepoints.items() if len(times) not in (2,) and len(times) < 6)
    for i in other_ids:
        errors.append(f"participant has unexpected number of timepoints: {i} -> {sorted(participant_timepoints[i])}")
    cohort_by_individual = {i: "longitudinal" for i in longitudinal_ids}
    cohort_by_individual.update({i: "baseline_only" for i in baseline_only_ids})
    cohort_by_individual.update({i: "other" for i in other_ids})

    sample_rows = []
    for sample_id in sorted(groups):
        runs = groups[sample_id]
        first = runs[0]
        if first["cohort_hint"] == "reference_library":
            cohort = "reference_library"
        else:
            cohort = cohort_by_individual.get(first["individual"], "other")
        techs = sorted({r["technical_replicate"] for r in runs}, key=lambda x: (len(x), x))
        plates = sorted({r["plate"] for r in runs})
        dates = sorted({r["run_date"] for r in runs})
        if len(techs) != len(runs) or techs != [str(i) for i in range(1, len(runs) + 1)]:
            errors.append(f"non-contiguous technical replicates for {sample_id}: {techs}")
        sample_rows.append({
            "sample_id": sample_id,
            "individual": first["individual"],
            "sampling_time": first["sampling_time"],
            "cohort": cohort,
            "n_runs": len(runs),
            "technical_replicates": ",".join(techs),
            "plates": ",".join(plates),
            "run_dates": ",".join(dates),
            "raw_files": ";".join(r["raw_file"] for r in runs),
            "maxquant_labels": ";".join(r["maxquant_label"] for r in runs),
        })

    maxquant_info = {}
    if args.maxquant_zip:
        labels = maxquant_intensity_labels(args.maxquant_zip)
        expected = {r["maxquant_label"] for r in parsed_runs}
        observed = set(labels)
        maxquant_info = {
            "proteinGroups_intensity_columns": len(labels),
            "expected_labels": len(expected),
            "labels_intersection": len(expected & observed),
            "missing_expected_labels": sorted(expected - observed),
            "unexpected_observed_labels": sorted(observed - expected),
            "zip_sha256": sha256(args.maxquant_zip),
        }
        if len(labels) != len(observed):
            errors.append("proteinGroups.txt has duplicate Intensity labels")
        if expected != observed:
            errors.append(
                "MaxQuant/SDRF label mismatch: "
                f"missing={sorted(expected-observed)[:5]} ({len(expected-observed)}), "
                f"unexpected={sorted(observed-expected)[:5]} ({len(observed-expected)})"
            )

    participant_runs = [r for r in parsed_runs if r["cohort_hint"] == "participant"]
    reference_runs = [r for r in parsed_runs if r["cohort_hint"] == "reference_library"]
    longitudinal_runs = [r for r in participant_runs if cohort_by_individual.get(r["individual"]) == "longitudinal"]
    baseline_runs = [r for r in participant_runs if cohort_by_individual.get(r["individual"]) == "baseline_only"]
    def count_by_n(rs):
        return {str(k): v for k, v in sorted(Counter(len(groups[r["sample_id"]]) for r in rs).items())}

    summary = {
        "project_accession": "PXD004242",
        "inputs": {
            "sdrf": str(args.sdrf),
            "sdrf_sha256": sha256(args.sdrf),
            "manifest": str(args.manifest),
            "manifest_sha256": sha256(args.manifest),
        },
        "api_manifest": {
            "files_total": len(manifest),
            "raw_files": len(raw_manifest),
            "non_raw_files": len(manifest) - len(raw_manifest),
            "raw_size_bytes": sum(int(r["file_size_bytes"]) for r in raw_manifest.values()),
        },
        "sdrf": {
            "rows": len(sdrf_rows),
            "unique_data_files": len(set(sdrf_files)),
            "participant_runs": len(participant_runs),
            "reference_runs": len(reference_runs),
        },
        "run_reconciliation": {
            "sdrf_raw_manifest_intersection": len(sdrf_raw & set(raw_manifest)),
            "sdrf_raw_missing_from_manifest": len(sdrf_raw - set(raw_manifest)),
            "manifest_raw_missing_from_sdrf": len(set(raw_manifest) - sdrf_raw),
        },
        "participants": {
            "numeric_participants": len(participant_ids),
            "longitudinal": len(longitudinal_ids),
            "longitudinal_ids": longitudinal_ids,
            "baseline_only": len(baseline_only_ids),
            "baseline_only_ids": baseline_only_ids,
            "other": len(other_ids),
            "other_ids": other_ids,
        },
        "reference_library": {
            "samples": len(reference_ids),
            "ids": sorted(reference_ids),
        },
        "biological_sample_units": {
            "all_groups": len(sample_rows),
            "participant_groups": sum(1 for x in sample_rows if x["cohort"] != "reference_library"),
            "longitudinal_groups": sum(1 for x in sample_rows if x["cohort"] == "longitudinal"),
            "baseline_only_groups": sum(1 for x in sample_rows if x["cohort"] == "baseline_only"),
            "reference_groups": sum(1 for x in sample_rows if x["cohort"] == "reference_library"),
            "participant_technical_replicate_distribution": count_by_n(participant_runs),
            "reference_technical_replicate_distribution": count_by_n(reference_runs),
        },
        "maxquant": maxquant_info,
        "errors": errors,
    }

    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_samples.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fields_out = ["sample_id", "individual", "sampling_time", "cohort", "n_runs", "technical_replicates", "plates", "run_dates", "raw_files", "maxquant_labels"]
    with args.out_samples.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields_out, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(sample_rows)

    print(f"api_files={len(manifest)} raw={len(raw_manifest)} sdrf_rows={len(sdrf_rows)}")
    print(f"participants={len(participant_ids)} longitudinal={len(longitudinal_ids)} baseline_only={len(baseline_only_ids)} references={len(reference_ids)}")
    print(f"biological_units participant={sum(1 for x in sample_rows if x['cohort'] != 'reference_library')} longitudinal={sum(1 for x in sample_rows if x['cohort'] == 'longitudinal')} baseline_only={sum(1 for x in sample_rows if x['cohort'] == 'baseline_only')} references={sum(1 for x in sample_rows if x['cohort'] == 'reference_library')}")
    print(f"participant_runs={len(participant_runs)} reference_runs={len(reference_runs)}")
    if maxquant_info:
        print(f"maxquant_intensity_columns={maxquant_info['proteinGroups_intensity_columns']} label_intersection={maxquant_info['labels_intersection']}")
    if errors:
        print(f"VALIDATION_ERRORS={len(errors)}", file=sys.stderr)
        for e in errors[:20]:
            print(f"- {e}", file=sys.stderr)
        return 1
    print("validation=PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
