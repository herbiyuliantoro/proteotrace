#!/usr/bin/env python3
"""Extract a compact run x protein intensity matrix from MaxQuant proteinGroups."""
from __future__ import annotations

import argparse
import csv
import io
import zipfile

csv.field_size_limit(10**9)

from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--zip', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()

    with zipfile.ZipFile(args.zip) as z, z.open('proteinGroups.txt') as raw:
        text = io.TextIOWrapper(raw, encoding='utf-8-sig', newline='')
        reader = csv.reader(text, delimiter='\t')
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        required = ['Majority protein IDs', 'Protein IDs', 'Only identified by site', 'Reverse', 'Potential contaminant']
        missing = [x for x in required if x not in idx]
        if missing:
            raise SystemExit(f'missing proteinGroups columns: {missing}')
        intensity_cols = [(i, name[len('Intensity '):]) for i, name in enumerate(header) if name.startswith('Intensity ')]
        if not intensity_cols:
            raise SystemExit('no Intensity columns found')
        labels = [label for _, label in intensity_cols]
        values_by_run = [[] for _ in labels]
        protein_ids = []
        seen = set()
        rows_seen = 0
        rows_kept = 0
        filtered = 0
        for row in reader:
            if not row:
                continue
            rows_seen += 1
            if any(row[idx[x]] == '+' for x in ('Only identified by site', 'Reverse', 'Potential contaminant')):
                filtered += 1
                continue
            key = row[idx['Majority protein IDs']] or row[idx['Protein IDs']]
            if not key:
                key = f'protein_group_{rows_seen}'
            if key in seen:
                raise SystemExit(f'duplicate clean protein identifier: {key}')
            seen.add(key)
            protein_ids.append(key)
            for j, (col_idx, _label) in enumerate(intensity_cols):
                values_by_run[j].append(row[col_idx])
            rows_kept += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['run'] + protein_ids)
        for label, values in zip(labels, values_by_run):
            writer.writerow([label] + values)
    print(f'proteinGroups_rows={rows_seen} clean_rows={rows_kept} filtered_rows={filtered}')
    print(f'runs={len(labels)} proteins={len(protein_ids)}')
    print(f'out={args.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
