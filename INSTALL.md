# Installation and execution

This project has two independent layers:

1. the compact Python analysis, which starts from committed `hpq`/`hpqrs` files;
2. the optional native BigDFT layer, which produces fresh wavefunctions and
   integrals on Apple Silicon.

Start with the Python layer. Building BigDFT is only necessary when changing a
geometry or wavelet parameter.

## Memory-safe command wrapper

`tools/run_with_memory_guard.py` watches the RSS of a command and all descendants,
checks system memory pressure, and caps common numerical thread pools. A limit must
be exceeded for three consecutive samples before the process group is terminated.

```bash
python3 tools/run_with_memory_guard.py \
  --max-rss-gb 8 \
  --min-free-percent 15 \
  --threads 2 \
  --log .artifacts/memory/example.log \
  -- COMMAND ARGS...
```

Exit 137 means the guard stopped the job to protect the machine. Logs under
`.artifacts/` are intentionally ignored by Git.

## Python environment

### Prerequisites

- macOS or Linux;
- [`uv`](https://docs.astral.sh/uv/) 0.9 or newer;
- Python 3.12 (installed automatically by uv if necessary).

On macOS:

```bash
brew install uv
```

Create the locked environment, including plotting and test dependencies:

```bash
python3 tools/run_with_memory_guard.py \
  --max-rss-gb 8 --min-free-percent 15 --threads 2 \
  --log .artifacts/memory/uv-sync.log -- \
  uv sync --extra analysis --extra dev
```

Run the tests:

```bash
python3 tools/run_with_memory_guard.py \
  --max-rss-gb 8 --min-free-percent 15 --threads 2 \
  --log .artifacts/memory/pytest.log -- \
  uv run pytest -q
```

Use `uv run` for scripts; activation is unnecessary.

## Reproduce the committed N₂ analysis

The first six stages read only compact committed files:

```bash
uv run python scripts/01_validate_n2.py
uv run python scripts/02_active_space.py
uv run python scripts/03_map_and_diag.py
uv run python scripts/04_vqe.py
uv run python scripts/05_sqd_qsci.py
uv run python scripts/06_dissociation.py
uv run python scripts/gen_headline_plot.py
uv run python scripts/08_scan_analysis.py
uv run python scripts/09_109asp_wall.py
uv run python scripts/10_wavelet_figures.py
```

Stage 6 runs eight PySCF calculations for each of two reference methods. Keep the
guard around it when changing basis sets or scan density:

```bash
python3 tools/run_with_memory_guard.py \
  --max-rss-gb 8 --min-free-percent 15 --threads 2 \
  --log .artifacts/memory/n2-dissociation.log -- \
  uv run python scripts/06_dissociation.py
```

## Apple-Silicon BigDFT layer

### Important source constraint

The custom integral extractor calls private BigDFT/liborbs APIs. It is verified
against the local patched 2024 research tree; it does not compile against stock
BigDFT 1.9.6 or current `devel`. A tag or commit identifier alone does not recreate
the patched tree.

Before spending time on a build:

```bash
python3 tools/check_bigdft_compatibility.py /path/to/bigdft-suite
```

All four checks must pass. See [docs/COMPATIBILITY.md](docs/COMPATIBILITY.md).

### Toolchain

The verified toolchain is macOS arm64 with gfortran 13.4.0. Do not silently
substitute Homebrew’s newest GCC: Fortran module files are compiler-version
specific, and later compilers apply stricter checks to this source era.

```bash
brew install gcc@13 open-mpi openblas autoconf automake libtool cmake pkg-config

export OMPI_FC=gfortran-13
export OMPI_CC=gcc-13
export OMPI_CXX=g++-13

gfortran-13 --version
OMPI_FC=gfortran-13 mpifort --showme:command
```

Both commands should resolve to gfortran 13.

The 2024 BigDFT tree also bundles an old JHBuild that imports `distutils` and
automake helpers that import `imp`. Those modules are absent from Homebrew Python
3.12+. The wrapper therefore uses macOS `/usr/bin/python3` (verified: 3.9.6) for
the native installer while leaving the analysis environment on Python 3.12. To
use another compatible interpreter, set `BIGDFT_INSTALLER_PYTHON` explicitly:

```bash
export BIGDFT_INSTALLER_PYTHON=/path/to/python3
"$BIGDFT_INSTALLER_PYTHON" -c 'import distutils.sysconfig, imp'
```

### Build the suite

Point `BIGDFT_SUITE_SOURCE` at the compatible patched source tree. The wrapper
checks the API surface, creates a `glibtoolize` shim, and runs the upstream
installer into ignored `bigdft/local/`.

```bash
export BIGDFT_SUITE_SOURCE=/path/to/patched/bigdft-suite

python3 tools/run_with_memory_guard.py \
  --max-rss-gb 8 --min-free-percent 15 --threads 2 \
  --interval 2 --status-every 30 \
  --log .artifacts/memory/bigdft-build.log -- \
  bash bigdft/macos/run_build.sh
```

Expected install prefix: `bigdft/local/install`. A clean build historically took
about 10–15 minutes and 2–4 GB RAM on the test machine. The memory guard, not this
estimate, is the safety boundary.

The supplied `apple_silicon.rc` includes:

- `-fallow-argument-mismatch`, `-ffree-line-length-none`, and `-std=legacy` for
  the older Fortran/MPI sources;
- `-finit-local-zero` as a mitigation for the nondeterministic Davidson
  uninitialized-memory failure;
- Apple Accelerate for BLAS/LAPACK;
- `CMAKE_POLICY_VERSION_MINIMUM=3.5` for the older NTPoly CMake project.

If a newer Python is forced instead, the source tree needs modern automake
`py-compile` helpers and a replacement for JHBuild's `distutils` import. The exact
prior-build notes are preserved in `docs/COMPATIBILITY.md`; this repository does
not mutate a source tree automatically beyond the explicit interface patch.

### Compile the extractor

```bash
export BIGDFT_PREFIX="$PWD/bigdft/local/install"

python3 tools/run_with_memory_guard.py \
  --max-rss-gb 4 --min-free-percent 15 --threads 2 \
  --log .artifacts/memory/toy-model-compile.log -- \
  bash bigdft/macos/build_toy_mac.sh
```

The output `bigdft/macos/toy_model_mac.x` is ignored by Git. The build script adds
an rpath and rewrites a stale absolute `libyaml` dependency to
`@rpath/libyaml-0.2.dylib`; verify with:

```bash
otool -L bigdft/macos/toy_model_mac.x
```

### Run a prepared case

Prepare a scratch directory containing:

```text
input.yaml
posinp.xyz
data/wavefunction-*
data/virtuals-*
hpq.inp
hpqrs.inp
```

Then run singleton MPI (plain executable, no `mpirun`):

```bash
export BIGDFT_PREFIX="$PWD/bigdft/local/install"

python3 tools/run_with_memory_guard.py \
  --max-rss-gb 8 --min-free-percent 15 --threads 2 \
  --log .artifacts/memory/toy-model-run.log -- \
  bash bigdft/macos/run_toy_model.sh /absolute/path/to/scratch-case
```

An optional second argument is a directory containing reference `hpq.out` and
`hpqrs.out`; the wrapper then performs a keyed numerical comparison.

The 75 MB archived wavefunction fixture is not committed. The 2026-08-07
verification used that local fixture, produced 240 + 7,260 records, and peaked at
1.62 GiB RSS.

### Fresh parameter scan

`scripts/07_param_scan.py` runs BigDFT, generates index lists in Python, extracts
integrals, evaluates HF/CASCI, and appends to `results/data/scan_results.json`.
It refuses an unguarded launch.

```bash
export BIGDFT_PREFIX="$PWD/bigdft/local/install"

python3 tools/run_with_memory_guard.py \
  --max-rss-gb 8 --min-free-percent 15 --threads 2 \
  --log .artifacts/memory/parameter-scan.log -- \
  uv run python scripts/07_param_scan.py --threads 2
```

For N₂, stay within approximately `crmult=4..7` and `frmult=6..8`. The extractor’s
projector construction failed deterministically at `crmult<=3.5` or `frmult>=9`
in the archived study; those are implementation limits, not evidence of new
physics.

## 109Asp safety rule

Do not run the complete 109Asp tensor through the dense loaders on a laptop. The
spin-orbital tensor alone is 38.9 GB before intermediates. The committed workflow
uses only `data/109asp/manifest.json` and constant-memory combinatorics. Any future
large-system analysis should be streamed, screened, and guarded on a cluster.
