#!/bin/bash
# Compile the research toy_model against a compatible Apple-Silicon BigDFT build.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BigDFT="${BIGDFT_PREFIX:-$HERE/../local/install}"

export OMPI_FC=gfortran-13 OMPI_CC=gcc-13   # keep the SAME gfortran that built the suite

SRC="${1:-$HERE/../postprocessing/toy_model_cpu.f90}"
OUT="${2:-$HERE/toy_model_mac.x}"

test "${WAVELET_QC_MEMORY_GUARD:-}" = "1" || {
  echo "ERROR: run this script through tools/run_with_memory_guard.py" >&2
  exit 2
}
test -d "$BigDFT/include" || {
  echo "ERROR: BigDFT headers not found under $BigDFT (set BIGDFT_PREFIX)" >&2
  exit 2
}

mpifort -o "$OUT" -O2 -fopenmp -fallow-argument-mismatch -ffree-line-length-none -fbacktrace \
  -I"$BigDFT/include" -I"$BigDFT/include/futile" \
  "$SRC" \
  -L"$BigDFT/lib" -lbigdft-1 -labinit -lxcf90 -lxc -lGaIn -lPSolver-1 -lCheSS-1 -lfutile-1 \
  -latlab-1 -lyaml -ldicts -lNTPoly -lorbs -framework Accelerate -lstdc++ \
  -Wl,-rpath,"$BigDFT/lib"

old_yaml="$(otool -L "$OUT" | awk '/libyaml-[0-9.]*\.dylib/{print $1; exit}')"
case "$old_yaml" in
  /*) install_name_tool -change "$old_yaml" @rpath/libyaml-0.2.dylib "$OUT" ;;
esac

echo "built: $OUT"
file "$OUT"
otool -L "$OUT"
