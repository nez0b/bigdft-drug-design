"""Report the 109Asp scaling wall without loading any large tensor.

The repository intentionally carries metadata, not the 1.39 GB complete file or
the misleading 3.6 MB truncated file.  This script derives every resource number
from ``data/109asp/manifest.json`` and simple combinatorics, so its memory use is
constant and tiny.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "109asp" / "manifest.json"


def main() -> None:
    metadata = json.loads(MANIFEST.read_text(encoding="utf-8"))
    n_spatial = metadata["system"]["spatial_orbitals"]
    n_spin = 2 * n_spatial
    n_pairs = n_spatial * (n_spatial + 1) // 2
    expected = n_pairs * (n_pairs + 1) // 2
    truncated = metadata["archived_truncated_run"]["line_count"]
    fraction = 100.0 * truncated / expected
    dense_gb = n_spin**4 * 8 / 1e9

    assert expected == metadata["complete_opt_ks_run"]["line_count"]

    print("109Asp (70-atom ligand complex) -- resource/provenance audit")
    print(f"  spatial / spin orbitals     : {n_spatial} / {n_spin}")
    print(f"  symmetry-unique ERIs        : {expected:,}")
    print(f"  archived truncated run      : {truncated:,} ({fraction:.2f}%)")
    print(f"  dense spin-orbital tensor   : {dense_gb:.1f} GB")
    print("  complete opt-KS file        : recovered locally, deliberately not vendored")
    print()
    print("The complete and truncated files are different runs and must not be mixed.")
    print("No 109Asp energy is reported from the truncated tensor.  The complete file")
    print("is retained only for guarded/cluster-scale follow-up analysis.")


if __name__ == "__main__":
    main()
