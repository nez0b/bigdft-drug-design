#!/usr/bin/env python3
"""Static compatibility check for the research toy_model API surface.

The integral extractor uses internal BigDFT/liborbs interfaces.  A release number
alone is therefore insufficient: point this script at an unpacked BigDFT suite
before starting a long build.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Check:
    label: str
    relative_path: str
    pattern: str
    explanation: str


CHECKS = (
    Check(
        "public full_local_potential",
        "bigdft/src/modules/rhopotential.f90",
        r"public\s*::[^\n]*full_local_potential[\s\S]*?subroutine\s+full_local_potential",
        "toy_model gathers the distributed local potential through this internal API",
    ),
    Check(
        "text orbital writer with keys",
        "liborbs/src/locregs.f90",
        r"subroutine\s+dump_array_txt\s*\([^)]*withKeys[^)]*\)",
        "the verified text-wavefunction path needs the optional withKeys argument",
    ),
    Check(
        "legacy local-Hamiltonian work-array API",
        "liborbs/src/locreg_operations.f90",
        r"subroutine\s+initialize_work_arrays_locham_nlr\s*\(\s*nlr\s*,\s*lr\s*,\s*nspinor\s*,\s*w\s*\)",
        "the extractor imports this generic from locreg_operations",
    ),
    Check(
        "gfortran-safe gaussian basis interface",
        "bigdft/src/modules/interfaces.f90",
        r"type\s*\(\s*gaussian_basis\s*\)\s*,\s*intent\s*\(\s*inout\s*\)\s*::\s*G\b",
        "intent(inout) preserves caller-initialized state required before Davidson",
    ),
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", type=Path, help="root of an unpacked bigdft-suite tree")
    args = parser.parse_args()
    suite = args.suite.expanduser().resolve()

    failures = 0
    print(f"BigDFT toy_model compatibility: {suite}")
    for check in CHECKS:
        path = suite / check.relative_path
        if not path.is_file():
            passed = False
            detail = "file missing"
        else:
            text = path.read_text(encoding="utf-8", errors="replace")
            passed = re.search(check.pattern, text, re.IGNORECASE) is not None
            detail = check.explanation
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {check.label}: {detail}")
        failures += not passed

    if failures:
        print(
            f"\nINCOMPATIBLE ({failures} failed check(s)). Do not start a long build; "
            "see docs/COMPATIBILITY.md."
        )
        return 1
    print("\nCOMPATIBLE with the verified extractor API surface.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
