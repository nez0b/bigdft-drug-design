#!/usr/bin/env bash
# Run toy_model in a scratch case directory containing input.yaml, posinp.xyz,
# data/, hpq.inp, and hpqrs.inp. Never point this at the vendored data directory.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CASE_DIR="${1:-}"
PREFIX="${BIGDFT_PREFIX:-$HERE/../local/install}"
BINARY="${BIGDFT_TOY_MODEL:-$HERE/toy_model_mac.x}"
REFERENCE_DIR=""

test -n "$CASE_DIR" || { echo "usage: $0 CASE_DIR [REFERENCE_DIR]" >&2; exit 2; }
test "${WAVELET_QC_MEMORY_GUARD:-}" = "1" || {
  echo "ERROR: run this script through tools/run_with_memory_guard.py" >&2
  exit 2
}
test -x "$BINARY" || { echo "ERROR: missing executable: $BINARY" >&2; exit 2; }
test -d "$CASE_DIR/data" || { echo "ERROR: missing $CASE_DIR/data" >&2; exit 2; }

CASE_DIR="$(cd "$CASE_DIR" && pwd)"
if test -n "${2:-}"; then
  REFERENCE_DIR="$(cd "$2" && pwd)"
  test -f "$REFERENCE_DIR/hpq.out" && test -f "$REFERENCE_DIR/hpqrs.out" || {
    echo "ERROR: reference directory needs hpq.out and hpqrs.out: $REFERENCE_DIR" >&2
    exit 2
  }
fi
cd "$CASE_DIR"
DYLD_LIBRARY_PATH="$PREFIX/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}" \
  "$BINARY" > toy_model.run.log 2>&1

wc -l hpq.out hpqrs.out
if test -n "$REFERENCE_DIR"; then
  python3 "$ROOT/tools/compare_integrals.py" hpq.out "$REFERENCE_DIR/hpq.out"
  python3 "$ROOT/tools/compare_integrals.py" hpqrs.out "$REFERENCE_DIR/hpqrs.out"
fi
