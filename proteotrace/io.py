"""Readers for the quantitative output people actually have.

A generic runs x features table is the contract. The other readers exist so
nobody has to reshape a MaxQuant or DIA-NN file by hand before trying the tool.
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

import pandas as pd

REQUIRED_SAMPLE_COLUMNS = ("run", "donor")


def _parse_numbers(values: pd.Series) -> tuple[pd.Series, bool]:
    """Numbers from a column, accepting a decimal comma (`7,04E+12`).

    Exports made under a European locale write decimal commas, which plain
    `to_numeric` turns into NaN without complaint. The comma reading is used
    only when the column never mixes "." and "," and the comma reading parses
    essentially every value the plain reading could not. Returns the parsed
    column and whether the decimal-comma reading was used.
    """
    if pd.api.types.is_numeric_dtype(values):
        return values, False
    text = values.astype("string").str.strip()
    plain = pd.to_numeric(text, errors="coerce")
    unparsed = text.notna() & (text != "") & plain.isna()
    if not unparsed.any():
        return plain.astype("float64"), False
    mixed = (text.str.contains(".", regex=False) & text.str.contains(",", regex=False)).any()
    comma = pd.to_numeric(text.str.replace(",", ".", regex=False), errors="coerce")
    if not mixed and comma[unparsed].notna().mean() >= 0.99:
        return comma.astype("float64"), True
    return plain.astype("float64"), False


def _parse_frame(frame: pd.DataFrame) -> pd.DataFrame:
    parsed, used = {}, False
    for column in frame.columns:
        parsed[column], hit = _parse_numbers(frame[column])
        used |= hit
    if used:
        _note_decimal_comma()
    return pd.DataFrame(parsed, index=frame.index)


def _note_decimal_comma() -> None:
    print("note: values use a decimal comma (e.g. 7,04E+12); read as decimal numbers.",
          file=sys.stderr)


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
    table = pd.read_csv(path, sep=None, engine="python", index_col=0, dtype=str)
    return _parse_frame(table)


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
    return _parse_frame(values).T


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
    # Parse before pivoting: taking the max of text values would compare strings.
    table[value_col], used_comma = _parse_numbers(table[value_col])
    if used_comma:
        _note_decimal_comma()
    wide = table.pivot_table(index=run_col, columns=group_col, values=value_col,
                             aggfunc="max")
    return wide.astype(float)


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
    matrix = READERS[resolved](path)
    if not (matrix > 0).to_numpy().any():
        raise ValueError(
            f"no positive numbers could be read from {path.name} (format {resolved}). "
            "Check the decimal separator, the column layout, or pass --format.")
    return matrix, resolved
