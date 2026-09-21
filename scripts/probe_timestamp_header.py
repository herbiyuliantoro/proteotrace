#!/usr/bin/env python3
"""Inspect a downloaded raw-file header/range for acquisition timestamps."""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


def values(pattern: bytes, data: bytes) -> list[str]:
    return sorted(set(x.decode("utf-8", "replace") for x in re.findall(pattern, data)))


def inspect(path: Path) -> dict[str, str | int]:
    data = path.read_bytes()
    start = values(rb"startTimeStamp\s*=\s*[\"']([^\"']+)", data)
    retention = values(rb"retentionTime\s*=\s*[\"']([^\"']+)", data)
    completion = values(rb"completionTime\s*=\s*[\"']([^\"']+)", data)
    return {
        "file": path.name,
        "bytes_inspected": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "startTimeStamp": ";".join(start),
        "retentionTime_values_seen": len(retention),
        "retentionTime_first": retention[0] if retention else "",
        "completionTime": ";".join(completion),
        "absolute_acquisition_timestamp_found": int(bool(start)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("headers", nargs="+", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    rows = [inspect(path) for path in args.headers]
    fields = list(rows[0]) if rows else []
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    writer = csv.DictWriter(__import__("sys").stdout, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
