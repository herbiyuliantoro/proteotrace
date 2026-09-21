"""Readers for the quantitative output people actually have.

A generic runs x features table is the contract. The other readers exist so
nobody has to reshape a MaxQuant or DIA-NN file by hand before trying the tool.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import pandas as pd

REQUIRED_SAMPLE_COLUMNS = ("run", "donor")


def read_samples(path: Path) -> pd.DataFrame:
    """Sample annotation: one row per run, with `run`, `donor`, and `batch`.

    `batch` may be omitted, in which case every run is treated as one batch
    and no correction is applied. That is rarely what you want, and the CLI
    says so, but it keeps the tool usable on a single-batch pilot.
    """
    table = pd.read_csv(path, sep=None, engine="python", dtype=str)
    table.columns = [c.strip().lower() for c in table.columns]
    missing = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(
            f"{path.name} is missing the column(s) {missing}. It needs at least "
            "'run' and 'donor', and normally 'batch' as well.")
    if "batch" not in table.columns:
        table["batch"] = "single_batch"
    table = table.dropna(subset=["run", "donor"])
    duplicated = table["run"][table["run"].duplicated()].tolist()
    if duplicated:
        raise ValueError(f"these runs appear more than once: {duplicated[:5]}")
    return table.set_index("run")[["donor", "batch"]]


def read_generic(path: Path) -> pd.DataFrame:
    """A table whose first column names the run and whose others are features."""
    table = pd.read_csv(path, sep=None, engine="python", index_col=0)
    return table.apply(pd.to_numeric, errors="coerce")


def read_maxquant(path: Path) -> pd.DataFrame:
    """MaxQuant `proteinGroups.txt`, or a zip containing one.

    Uses the `Intensity <run>` columns, and drops the reverse, contaminant
    and site-only rows that MaxQuant flags, because scoring on them would be
    scoring on decoys.
    """
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            name = next((n for n in archive.namelist()
                         if n.endswith("proteinGroups.txt")), None)
            if name is None:
                raise ValueError(f"{path.name} contains no proteinGroups.txt")
            with archive.open(name) as handle:
                table = pd.read_csv(handle, sep="\t", low_memory=False)
    else:
        table = pd.read_csv(path, sep="\t", low_memory=False)

    for flag in ("Reverse", "Potential contaminant", "Only identified by site"):
        if flag in table.columns:
            table = table.loc[table[flag].isna()]

    prefix = "Intensity "
    columns = [c for c in table.columns if c.startswith(prefix)]
    if not columns:
        raise ValueError("no 'Intensity <run>' columns found in proteinGroups.txt")
    index = table["Protein IDs"] if "Protein IDs" in table.columns else table.index
    values = table[columns].copy()
    values.columns = [c[len(prefix):] for c in columns]
    values.index = index
    return values.T.apply(pd.to_numeric, errors="coerce")


def read_diann(path: Path) -> pd.DataFrame:
    """A DIA-NN main report, pivoted to runs x protein groups."""
    table = pd.read_csv(path, sep="\t", low_memory=False)
    run_col = next((c for c in ("Run", "File.Name") if c in table.columns), None)
    group_col = next((c for c in ("Protein.Group", "Protein.Ids") if c in table.columns), None)
    value_col = next((c for c in ("PG.MaxLFQ", "PG.Quantity", "PG.Normalised")
                      if c in table.columns), None)
    if not (run_col and group_col and value_col):
        raise ValueError(
            "this does not look like a DIA-NN report: expected a run column, a "
            "protein-group column and a PG quantity column")
    wide = table.pivot_table(index=run_col, columns=group_col, values=value_col,
                             aggfunc="max")
    return wide.apply(pd.to_numeric, errors="coerce")


READERS = {"generic": read_generic, "maxquant": read_maxquant, "diann": read_diann}


def sniff(path: Path) -> str:
    """Guess the format from the file name, falling back to generic."""
    name = path.name.lower()
    if "proteingroups" in name or name.endswith(".zip"):
        return "maxquant"
    if "report" in name and name.endswith((".tsv", ".txt")):
        header = pd.read_csv(path, sep="\t", nrows=0).columns
        if {"Protein.Group", "Protein.Ids"} & set(header):
            return "diann"
    return "generic"


def read_matrix(path: Path, fmt: str = "auto") -> tuple[pd.DataFrame, str]:
    resolved = sniff(path) if fmt == "auto" else fmt
    if resolved not in READERS:
        raise ValueError(f"unknown format {resolved!r}; choose from {sorted(READERS)}")
    return READERS[resolved](path), resolved
