# Provenance and contribution boundaries

This repository spans several stages of a collaboration. The purpose of this note
is to preserve who did what and to distinguish historical artifacts from the later
completion work.

## Timeline

### Original NTU project

The project explored whether BigDFT’s Daubechies-wavelet orbitals could reduce the
orbital/qubit cost of molecular quantum simulations and eventually support a drug-
design workflow.

- Shih-Kai Chou developed the original `toy_model.f90` workflow and the Poisson-
  solver strategy used to generate electron-repulsion integrals.
- Po-Jen Wang and Jing-Yan Huang worked on the electron-integral engine and its
  BigDFT integration.
- The QEE encoding work described in the historical report was Tung’s and Hsiao’s.
- The linear-scaling variant was joint work.
- Po-Jen Wang’s work included the BigDFT API port, three correctness fixes, CUDA
  Poisson offload, N₂ and 109Asp runs, and downstream Qiskit/SQD integration.

The project did not produce a binding affinity, and it did not land a matched
BigDFT-vs-Gaussian N₂ benchmark. Historical presentations say so explicitly.

### Later completion

The later work collected the scattered files, rebuilt the pipeline on modern
Qiskit/PySCF, ported the local build to Apple Silicon, recovered the complete
109Asp two-electron archive, and finished the N₂ comparison using curve-shape
quantities that remain meaningful across a pseudopotential boundary.

The QSCI/SQD code in `src/wavelet_qc/sqd.py` uses the maintained
`qiskit-addon-sqd` API. It replaces older notebooks whose package pins and sampling
provenance were not suitable for a reproducible modern result.

## Code lineage

| component | lineage |
|---|---|
| `bigdft/postprocessing/toy_model_cpu.f90` | BigDFT GPL source plus NTU research modifications; retained notices |
| `src/wavelet_qc/integrals.py` | cleaned port of the original `BigQiskit.py` load/symmetry logic |
| `src/wavelet_qc/reference.py` | cleaned HF reconstruction plus new PySCF checks |
| `src/wavelet_qc/active_space.py` | cleaned port of the original MP2/IEPA1-style selection code |
| `src/wavelet_qc/mapping.py` | modern explicit Jordan–Wigner implementation |
| `src/wavelet_qc/vqe.py` | modern Qiskit statevector demonstration |
| `src/wavelet_qc/sqd.py` | modern QSCI/SQD implementation using qiskit-addon-sqd |
| `src/wavelet_qc/dissociation.py` | later N₂ comparison and spectroscopy analysis |
| `src/wavelet_qc/bigdft_run.py` | later native-run harness and Python `genindex` reconstruction |

The Python `gen_index` routine was checked byte-for-byte against the surviving N₂
index inputs. The original `genindex` source was not found; only a Linux binary
survived, so the implementation is described as a reconstructed port, not the
original source.

## Data lineage

### Equilibrium N₂ fixture

`data/n2/` contains a complete 15-spatial-orbital valence-pseudopotential fixture:

- 240 one-electron records;
- 7,260 symmetry-unique two-electron records;
- the BigDFT log and geometry used to interpret them.

The local Apple-Silicon runtime fixture contains the larger wavefunction files but
is not committed. On 2026-08-07 it regenerated the integral keys within
`5e-12/1e-12 Ha` of the archived reference.

### N₂ scan

`data/n2_scan/` contains eight geometries from 0.2 to 1.6976 Å, with 10 spatial
orbitals per geometry (five occupied plus five virtual). Each complete geometry
contains 110 one-body and 1,540 symmetry-unique two-body records. Two archived
directories lack `posinp.xyz`; their bond lengths are encoded by the historical
directory convention and are decoded explicitly in `dissociation.py`.

### 109Asp

`data/109asp/manifest.json` records two distinct datasets:

- a recovered complete opt-KS file: 38,531,031 records, 1,387,117,116 bytes,
  SHA-256 `6717effd…`;
- an older incomplete file: 100,000 records, 3,600,000 bytes, SHA-256
  `6e37d086…`.

They are different runs and cannot be combined. Neither file is vendored. The
machine-readable `results/data/asp109_mp2_archived.json` came from a prior guarded
analysis of the complete opt-KS dataset; it is kept as an archived result and is
not part of the laptop reproduction path.

## Reports and figures

- `docs/reports/zapata_wavelet_project_deck.pdf` presents the original work and
  its then-unresolved limitations.
- `docs/reports/wavelet_qc_primer.pdf` is the first-principles teaching draft. Its
  comparison section incorrectly calls GTH-PBE the same pseudopotential as the
  BigDFT HGH-K/PBE run; `docs/reports/README.md` records the erratum and the
  curated methodology corrects it.
- `results/figures/` contains generated outputs, not manually redrawn values.
  Their source JSON and plotting scripts are committed alongside them.

## Verification dates are part of provenance

Software stacks and private BigDFT APIs change. Statements such as “current
upstream is incompatible” refer to the revisions and 2026-08-07 audit recorded in
`docs/COMPATIBILITY.md`; they are not timeless claims. Re-run the checker and
record a new manifest when moving to another source revision or platform.
