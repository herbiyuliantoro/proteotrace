#!/usr/bin/env python3
"""Peptide-to-protein ratio matrices for idea 3 of docs/EXTENSION_PLAN.md.

For every run and peptide: log2(peptide intensity) minus the median log2
intensity of the peptides of the same protein in that run. A common genetic
variant that changes a peptide's sequence lowers the signal of the reference
peptide in carriers, so this ratio behaves like a genotype read-out. Proteins
need at least three quantified peptides; peptides must be observed in at least
80% of the runs of the cohort.

Deviation from the plan, stated here and in the plan's results section: the
protein reference is the median of its peptides in each run, not the protein
group intensity reported by each search engine. This uses the same definition
for MaxQuant, DIA-NN and OpenSWATH, whose protein summaries differ.

Writes data/peptides/<cohort>_peptide_ratios.tsv.gz (runs x peptides) and
<cohort>_peptide_info.tsv (sequence, protein, positions where known).

    python extract_peptide_ratios.py PXD004242 PXD062541 PXD001064
"""
from __future__ import annotations

import io
import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_pxd004242_layer1 import read_sample_table, select_runs
from proteotrace.io import _parse_numbers

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "data/peptides"
MIN_RUN_FRACTION = 0.80
MIN_PEPTIDES_PER_PROTEIN = 3


def ratios(long: pd.DataFrame, runs: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """long: columns run, sequence, protein, intensity (>0). Returns (ratio matrix, info)."""
    long = long[long["run"].isin(runs) & (long["intensity"] > 0)].copy()
    long["log"] = np.log2(long["intensity"])
    wide = long.pivot_table(index="run", columns="sequence", values="log", aggfunc="max").reindex(runs)
    protein_of = long.groupby("sequence")["protein"].first()
    coverage = wide.notna().mean()
    eligible = coverage[coverage >= MIN_RUN_FRACTION].index
    n_per_protein = protein_of.loc[eligible].value_counts()
    proteins = n_per_protein[n_per_protein >= MIN_PEPTIDES_PER_PROTEIN].index
    keep = [s for s in eligible if protein_of[s] in set(proteins)]
    wide = wide[keep]
    ref = pd.DataFrame(index=wide.index)
    for prot, seqs in protein_of.loc[keep].groupby(protein_of.loc[keep]):
        ref[prot] = wide[seqs.index].median(axis=1)
    ratio = wide - ref[protein_of.loc[keep].to_numpy()].to_numpy()
    info = pd.DataFrame({"sequence": keep, "protein": protein_of.loc[keep].to_numpy(),
                         "run_coverage": coverage.loc[keep].to_numpy()})
    return ratio, info


def pxd004242():
    base = ROOT / "data/PXD004242"
    rows = read_sample_table(base / "PXD004242_biological_samples.tsv", "longitudinal")
    runs = list(select_runs(rows, "random", None, 20261090).index)
    with zipfile.ZipFile(base / "MaxQuant_Output.zip") as z, z.open("peptides.txt") as h:
        header = io.TextIOWrapper(h, encoding="utf-8").readline().rstrip("\n").split("\t")
    cols = ["Sequence", "Leading razor protein", "Start position", "End position", "Reverse",
            "Potential contaminant"] + [f"Intensity {r}" for r in runs]
    missing = [c for c in cols if c not in header]
    if missing:
        sys.exit(f"missing columns: {missing[:5]}")
    with zipfile.ZipFile(base / "MaxQuant_Output.zip") as z, z.open("peptides.txt") as h:
        pep = pd.read_csv(h, sep="\t", usecols=cols, low_memory=False)
    pep = pep[(pep["Reverse"] != "+") & (pep["Potential contaminant"] != "+")]
    long = pep.melt(id_vars=["Sequence", "Leading razor protein"], value_vars=[f"Intensity {r}" for r in runs],
                    var_name="run", value_name="intensity")
    long["run"] = long["run"].str[len("Intensity "):]
    long = long.rename(columns={"Sequence": "sequence", "Leading razor protein": "protein"})
    ratio, info = ratios(long, runs)
    pos = pep.set_index("Sequence")[["Start position", "End position"]]
    info["start"] = pos.loc[info["sequence"], "Start position"].to_numpy()
    info["end"] = pos.loc[info["sequence"], "End position"].to_numpy()
    return ratio, info


def pxd062541():
    samples = pd.read_csv(ROOT / "examples/PXD062541_samples.tsv", sep="\t", dtype=str)
    rep = pd.read_csv(ROOT / "data/PXD062541_report.tsv", sep="\t", dtype=str,
                      usecols=["Run", "Stripped.Sequence", "Protein.Group", "Precursor.Normalised"])
    rep["intensity"], _ = _parse_numbers(rep["Precursor.Normalised"])
    rep = rep.groupby(["Run", "Stripped.Sequence"], as_index=False).agg(
        intensity=("intensity", "sum"), protein=("Protein.Group", "first"))
    rep["protein"] = rep["protein"].str.split(";").str[0]
    long = rep.rename(columns={"Run": "run", "Stripped.Sequence": "sequence"})
    return ratios(long, list(samples["run"]))


def pxd001064():
    import analysis_primary_estimand as PE
    _, meta = PE.pxd001064(ROOT)
    path = ROOT / "data/PXD001064/OpenSWATHResult_feature_alignment_requant.tsv"
    head = pd.read_csv(path, sep="\t", nrows=0).columns
    file_col = next(c for c in ("filename", "align_origfilename", "run_id") if c in head)
    cols = [file_col, "Sequence", "ProteinName", "Intensity", "decoy"] + (["m_score"] if "m_score" in head else [])
    parts = []
    for chunk in pd.read_csv(path, sep="\t", usecols=cols, chunksize=500_000, low_memory=False):
        chunk = chunk[chunk["decoy"] == 0]
        if "m_score" in chunk:      # requantified features carry no m_score; keep them
            chunk = chunk[chunk["m_score"].isna() | (chunk["m_score"] <= 0.01)]
        chunk["run"] = chunk[file_col].astype(str).map(
            lambda f: (re.search(r"(yanliu_L\d{6}_\d+_SW)", f) or [None])[0])
        parts.append(chunk[["run", "Sequence", "ProteinName", "Intensity"]])
    df = pd.concat(parts)
    df = df[df["run"].notna()]
    df = df.groupby(["run", "Sequence"], as_index=False).agg(
        intensity=("Intensity", "sum"), protein=("ProteinName", "first"))
    df["protein"] = df["protein"].astype(str).str.replace(r"^\d+/", "", regex=True).str.split("/").str[0]
    long = df.rename(columns={"Sequence": "sequence"})
    runs = list(meta.index)
    found = set(long["run"])
    missing = [r for r in runs if r not in found]
    if missing:
        print(f"warning: {len(missing)} runs absent from the OpenSWATH file, e.g. {missing[:3]}")
    return ratios(long, [r for r in runs if r in found])


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    for name in sys.argv[1:] or ["PXD004242", "PXD062541", "PXD001064"]:
        ratio, info = {"PXD004242": pxd004242, "PXD062541": pxd062541, "PXD001064": pxd001064}[name]()
        ratio.to_csv(OUT / f"{name}_peptide_ratios.tsv.gz", sep="\t", float_format="%.5f")
        info.to_csv(OUT / f"{name}_peptide_info.tsv", sep="\t", index=False)
        print(f"{name}: {ratio.shape[0]} runs x {ratio.shape[1]} peptides from "
              f"{info['protein'].nunique()} proteins; missing values {ratio.isna().mean().mean():.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
