#!/usr/bin/env bash
# Build the compatible BigDFT suite into bigdft/local under a memory guard.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="${BIGDFT_SUITE_SOURCE:-}"
WORK="${BIGDFT_BUILD_DIR:-$ROOT/bigdft/local}"
INSTALLER_PYTHON="${BIGDFT_INSTALLER_PYTHON:-}"

test -n "$SOURCE" || {
  echo "ERROR: set BIGDFT_SUITE_SOURCE to the compatible patched suite tree" >&2
  exit 2
}
test "${WAVELET_QC_MEMORY_GUARD:-}" = "1" || {
  echo "ERROR: run this script through tools/run_with_memory_guard.py" >&2
  exit 2
}
for required_command in /opt/homebrew/bin/glibtoolize gfortran-13 gcc-13 g++-13 mpifort; do
  command -v "$required_command" >/dev/null 2>&1 || {
    echo "ERROR: required build command not found: $required_command" >&2
    exit 2
  }
done

# The bundled 2024 JHBuild imports stdlib distutils and old automake helpers use
# imp. Homebrew Python 3.12+ removed both. macOS still ships a suitable Python
# 3.9; callers on other systems can set BIGDFT_INSTALLER_PYTHON explicitly.
if test -z "$INSTALLER_PYTHON" && test -x /usr/bin/python3; then
  if /usr/bin/python3 -c 'import distutils.sysconfig, imp' >/dev/null 2>&1; then
    INSTALLER_PYTHON=/usr/bin/python3
  fi
fi
if test -z "$INSTALLER_PYTHON"; then
  INSTALLER_PYTHON="$(command -v python3 || true)"
fi
test -n "$INSTALLER_PYTHON" && "$INSTALLER_PYTHON" -c \
  'import distutils.sysconfig, imp' >/dev/null 2>&1 || {
  echo "ERROR: BigDFT's bundled JHBuild needs Python with distutils and imp." >&2
  echo "       On macOS use /usr/bin/python3, or set BIGDFT_INSTALLER_PYTHON." >&2
  exit 2
}

python3 "$ROOT/tools/check_bigdft_compatibility.py" "$SOURCE"
mkdir -p "$WORK/shims"
ln -sf /opt/homebrew/bin/glibtoolize "$WORK/shims/libtoolize"

export PATH="$WORK/shims:$PATH"
export OMPI_FC=gfortran-13
export OMPI_CC=gcc-13
export OMPI_CXX=g++-13
export PYTHON="$INSTALLER_PYTHON"

echo "BigDFT installer Python: $($INSTALLER_PYTHON --version 2>&1) ($INSTALLER_PYTHON)"

cd "$WORK"
"$INSTALLER_PYTHON" "$SOURCE/Installer.py" -y -f "$HERE/apple_silicon.rc" build
test -x "$WORK/install/bin/bigdft" && test -d "$WORK/install/include" || {
  echo "ERROR: Installer returned without a usable BigDFT prefix at $WORK/install" >&2
  exit 1
}
echo "verified BigDFT prefix: $WORK/install"
