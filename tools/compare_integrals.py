#!/usr/bin/env python3
"""Compare BigDFT integral text files by index key and numerical value."""

from __future__ import annotations

import argparse
from pathlib import Path


def load(path: Path) -> dict[tuple[int, ...], float]:
    values = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            key = tuple(int(value) for value in fields[:-1])
            number = float(fields[-1])
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: invalid integral record") from error
        values[key] = number
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("actual", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--atol", type=float, default=1e-8)
    args = parser.parse_args()

    actual = load(args.actual)
    reference = load(args.reference)
    shared = actual.keys() & reference.keys()
    missing = reference.keys() - actual.keys()
    extra = actual.keys() - reference.keys()
    maximum = max((abs(actual[key] - reference[key]) for key in shared), default=float("inf"))
    passed = not missing and not extra and maximum <= args.atol
    print(
        f"{args.actual.name}: actual={len(actual)}, reference={len(reference)}, "
        f"missing={len(missing)}, extra={len(extra)}, max|delta|={maximum:.3e} "
        f"[{'PASS' if passed else 'FAIL'}]"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
