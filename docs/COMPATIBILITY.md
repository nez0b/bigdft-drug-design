# BigDFT compatibility and macOS port

## Bottom line

`bigdft/postprocessing/toy_model_cpu.f90` is verified on Apple Silicon against a
patched 2024 research tree. It is **not** compatible with an unmodified current
BigDFT release. The phrase “ported to the latest BigDFT” in earlier working notes
was too broad: the code was ported along a locally patched internal-API line.

This distinction matters because `toy_model` is not a consumer of a stable public
library interface. It reaches into BigDFT and liborbs implementation details to
reconstruct local/nonlocal Hamiltonian terms and two-electron Coulomb integrals.

## Source revisions audited

| source | revision | date | extractor status |
|---|---|---:|---|
| local research tree (`.gitrepo` metadata) | BigDFT `66dbb1b…`, liborbs `232316e…`, plus local patches | 2024-era | compatible; compiled and run |
| [BigDFT suite 1.9.6](https://gitlab.com/l_sim/bigdft-suite/-/tree/1.9.6) | `1e7d6b3d…` | 2025-06-15 | incompatible |
| [BigDFT suite devel](https://gitlab.com/l_sim/bigdft-suite) | `aee86bd5…` | 2026-03-09 | incompatible |

The public checks were performed on 2026-08-07. Use
`tools/check_bigdft_compatibility.py` again for any later revision.

## Required internal APIs

The static checker tests four load-bearing assumptions:

1. `rhopotential` publicly exposes `full_local_potential`. Current public source
   retains stale tutorial call sites but no callable public implementation.
2. liborbs `dump_array_txt` accepts the optional `withKeys` argument used by the
   verified text-wavefunction path.
3. `initialize_work_arrays_locham` is exported from `locreg_operations` with the
   legacy overload used by this extractor. Newer liborbs reorganizes this API.
4. `gaussian_pswf_basis` declares its `G` argument as `intent(inout)`, matching
   the implementation and protecting a required initialization from gfortran
   optimization.

Run:

```bash
python3 tools/check_bigdft_compatibility.py /path/to/bigdft-suite
```

Do not interpret a passing release/tag check as full numerical validation; it is
only an early fail-fast check. A build plus N₂ invariant test is still required.

## macOS-specific changes

### Toolchain and linking

- Pinned Homebrew `gcc@13`/gfortran 13.4.0 and used the same compiler through
  Open MPI (`OMPI_FC=gfortran-13`).
- Replaced Intel MKL, `iomp5`, Linux `librt`, and Linux linker flags with Apple
  Accelerate and GNU OpenMP.
- Added compatibility flags for old Fortran call sites:
  `-fallow-argument-mismatch`, `-ffree-line-length-none`, and `-std=legacy`.
- Added an executable rpath and rewrote `libyaml` to
  `@rpath/libyaml-0.2.dylib`. The prior library carried an absolute install name
  pointing at a deleted cluster-build path; `DYLD_LIBRARY_PATH` had masked it.
- Run `toy_model` as singleton MPI. `mpirun` was unnecessary and triggered a
  local PMIx/PRRTE version conflict.

### Autotools and CMake

- The bundled JHBuild imports `distutils`, and older automake helpers import
  `imp`; both are absent from current Homebrew Python. The build wrapper selects
  macOS `/usr/bin/python3` (3.9.6 on the verified host), or accepts an explicit
  compatible `BIGDFT_INSTALLER_PYTHON`.
- Homebrew provides `glibtoolize`, so the build wrapper creates a local
  `libtoolize` shim.
- CMake 4 needs `CMAKE_POLICY_VERSION_MINIMUM=3.5` for the older NTPoly project.
- NTPoly needs the same modern-gfortran argument/line-length compatibility flags.
- If a modern Python must be used, replace old `py-compile` helpers with the
  current Homebrew automake helper (or regenerate autotools) and port JHBuild away
  from `distutils` first.

### Fortran interface correctness

The implementation of `gaussian_pswf_basis` reads `G%rxyz` as an already-built
guard and declares `G` `intent(inout)`. The interface declared `intent(out)`.
At `-O2`, gfortran was therefore allowed to discard the caller’s preceding
`nullify(G%rxyz)`, causing the Davidson virtual-orbital stage to read uninitialized
state and fail in `gaudim_check`—sometimes with a Fortran `stop` that returned
zero. `bigdft/macos/gfortran_interface_fix.patch` records the two intent
alignments.

### Remaining Davidson issue

A separate nondeterministic uninitialized-memory path can make the Davidson
gradient NaN from iteration one, leaving all virtual eigenvalues zero. The current
mitigation is:

- compile BigDFT with `-finit-local-zero`;
- detect the all-zero virtual spectrum;
- discard that attempt’s wavefunctions and retry up to three times.

This is a mitigation, not a root-cause fix. Successful attempts were checked
against cluster orbital-energy and HF/CAS invariants. The behavior must not be
hidden in a scientific result record; retry counts are preserved.

## Verified 2026-08-07 checks

The checked source is exactly the committed
`bigdft/postprocessing/toy_model_cpu.f90` (SHA-256
`c490ed61f84bf2ad6bb82d27c4646e811856e06f6d8dfa6be9bcda9ccae202ce`).

1. Static compatibility: all four checks passed on the local research tree; three
   failed on public `devel`.
2. Installer wrapper: an incremental audit of the existing 18-module build passed
   with macOS system Python 3.9.6 after rejecting Homebrew Python 3.14's missing
   `distutils`; a clean suite rebuild was not repeated.
3. Fresh extractor compile: passed under a 4 GiB guard, at most 0.10 GiB observed
   peak RSS; Mach-O arm64 binary; `@rpath` yaml linkage confirmed.
4. N₂ runtime: all 240 one-electron and 7,260 two-electron records were written.
5. Keyed numerical comparison: no missing or extra records; maximum deviations
   `5e-12 Ha` (`hpq`) and `1e-12 Ha` (`hpqrs`).
6. Runtime memory: 1.62 GiB observed process-tree RSS, with system free memory
   remaining above 90% during the run.

Machine-readable details are in `bigdft/macos/VERIFIED_SOURCE.json`.

## What is not claimed

- The public 1.9.6 or current `devel` trees were not forward-ported.
- The local source cannot be reconstructed solely from the two `.gitrepo`
  revisions; the research patches are part of the required source state.
- The entire BigDFT suite was not vendored here. Doing so would add hundreds of
  megabytes and obscure upstream licensing/history.
- Fresh SCF wavefunctions need not reproduce individual virtual-orbital integrals
  byte-for-byte: degenerate π orbitals may rotate or change sign. Compare
  gauge-invariant energies, eigenvalues, and curve properties for fresh SCF runs.

Forward-porting to current BigDFT should replace the private potential-gathering
and work-array calls with supported APIs, then re-establish the N₂ invariants. It
should not be treated as a compiler-flag-only exercise.
