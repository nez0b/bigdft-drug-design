# Results and interpretation

## Completed N₂ comparison

The later work closes the main gap in the original project: a complete N₂
dissociation pipeline from BigDFT Daubechies-wavelet orbitals to electron
integrals, active-space Hamiltonians, correlated energies, and matched reference
curves.

All correlated curves use ten valence electrons in ten spatial orbitals. The
all-electron cc-pVDZ calculation represents all 14 electrons but freezes two
doubly occupied N 1s-derived molecular orbitals, leaving the same ten correlated
valence electrons as the pseudopotential calculations.

| method | represented electrons | correlated space | rₑ (Å) | Δrₑ vs exp. | ωₑ (cm⁻¹) |
|---|---:|---:|---:|---:|---:|
| BigDFT HGH-PBE wavelet | 10 valence | CAS(10e,10o) | 1.0875 | -0.0102 | 2548.5 |
| PySCF GTH-PBE / gth-dzvp | 10 valence | CAS(10e,10o) | 1.1066 | +0.0089 | 2452.5 |
| PySCF all-electron / cc-pVDZ | 14 total; 4 frozen-core | CAS(10e,10o) | 1.1050 | +0.0074 | 2358.5 |
| experiment | — | — | 1.0977 | — | 2358.6 |

![Minimum-shifted N2 curves](results/figures/n2_dissociation_headline.png)

The wavelet equilibrium distance differs from experiment by about 0.93%. Its
harmonic frequency is about 8% high and should be treated as preliminary: the
archived curve has 0.2 Å spacing near equilibrium and only five virtual orbitals.
A denser near-minimum scan is the clearest next scientific improvement.

### Why absolute energies are not compared

The BigDFT/GTH calculations remove core electrons with pseudopotentials and have
total energies near `-20 Ha`; the all-electron cc-pVDZ values are near `-109 Ha`.
Those numbers have different cores and energy zeros. Subtracting them would not
measure a basis error. Each curve is instead shifted to its own minimum, and the
comparison uses zero-independent shape quantities: `rₑ`, curvature/`ωₑ`, and
relative energy along the bond coordinate.

The GTH-PBE/gth-dzvp comparator aligns the explicit valence-electron count and the
broad pseudopotential treatment, but not the exact pseudopotential or energy zero:
BigDFT uses HGH-K/PBE. It also does not isolate a pure “wavelet versus Gaussian”
basis effect: the BigDFT orbitals begin from PBE Kohn–Sham orbitals, whereas the
PySCF CASCI helper begins from RHF orbitals, and the implementations and finite
active spaces differ.

## Integral and qubit validation at equilibrium

The complete equilibrium fixture contains 15 spatial orbitals, 10 valence
electrons, 240 symmetry-unique one-electron records, and 7,260 symmetry-unique
two-electron records.

| check | result |
|---|---:|
| reconstructed one-electron occupied contribution | `-54.945421 Ha` |
| reconstructed HF electronic energy | `-31.526183 Ha` |
| ion-ion term | `+12.053052 Ha` |
| reconstructed HF total | `-19.473130 Ha` |
| BigDFT PBE Kohn–Sham total | `-19.910391 Ha` |
| HF − PBE gap on shared orbitals | `+0.4373 Ha` |

The HF/PBE gap is a functional difference, not an integral error. For the
MP2-selected CAS(4e,4o), explicit Jordan–Wigner diagonalization and independent
PySCF FCI agree to `6.75e-14 Ha` at `-31.60439547 Ha` electronic energy. Adding
the ion-ion term gives `-19.551343 Ha`, with a physical negative active-space
correlation contribution of `-0.078213 Ha` relative to reconstructed HF.

The corresponding eight-qubit Hamiltonian has 185 Pauli terms. A 64-parameter
hardware-efficient VQE reaches `-31.600342 Ha`, 4.05 mHa above exact. This is a
useful runnable NISQ demonstration, but it is outside the conventional 1.6 mHa
chemical-accuracy threshold.

## Wavelet-parameter sensitivity and cost

At the equilibrium geometry, decreasing grid spacing from `hgrid=0.45` to `0.20`
changes the BigDFT energy from `-19.905076` to `-19.910426 Ha` and the fixed
CAS(10e,10o) energy from `-19.568101` to `-19.572608 Ha`. The `hgrid=0.35` point
is already within 1.6 mHa of the finest tested CAS value, while the extractor cost
grows from roughly 57 seconds/0.5 GB to 406 seconds/2.9 GB in the archived run.

`frmult=6,7,8` produces essentially unchanged energies for this case. Increasing
`crmult=4→7` converges the underlying Kohn–Sham energy but does not make the fixed
CAS energy monotonically lower (`-19.577734→-19.565733 Ha`): changing the spatial
box changes which diffuse virtual orbitals occupy the fixed-size subspace, so
these CAS calculations are not variational with respect to `crmult`.

Three complete curves give stable fitted geometry as resolution tightens:

| hgrid (bohr) | rₑ (Å) | ωₑ (cm⁻¹) |
|---:|---:|---:|
| 0.35 | 1.0876 | 2547 |
| 0.25 | 1.0875 | 2549 |
| 0.20 | 1.0873 | 2552 |

The extractor fails deterministically outside part of the tested support range
(`crmult<=3.5` and `frmult>=9`). These are implementation limits, not physical
conclusions, and the failed records remain in the JSON for auditability.

## QSCI/SQD completion

On the wavelet CAS(6e,6o) Hamiltonian, samples drawn from the exact PySCF FCI
coefficient distribution reduce the selected-CI residual from 5.561 mHa at 5,000
oracle shots to 0.184 mHa at 500,000 shots. The 100,000-shot result is 1.536 mHa
from FCI. This validates the maintained QSCI/SQD classical stages and energy
convention; it is explicitly a classical sampler-oracle study, not QPU data.

The older QURI-QSCI attempt remains an unresolved negative diagnostic at
3,978.770 mHa from its CASCI reference. It is not counted as a result. See
[docs/SQD_QSCI.md](docs/SQD_QSCI.md) for the full interpretation.

## 109Asp scaling wall

The recovered complete 109Asp two-electron archive has 38,531,031 text records
and is not vendored. A separate older file contains only 100,000 records (0.26%)
from a different run and cannot be appended to it. With 132 spatial orbitals, a
dense 264-spin-orbital two-electron tensor alone would require about 38.9 GB
before intermediates. The laptop-safe workflow therefore reports metadata and
combinatorial resource estimates only; it does not manufacture a molecular energy
from the truncated file.

## Scope of the claims

- The native extractor is verified against the local patched 2024 BigDFT API
  line, not stock BigDFT 1.9.6 or current `devel`.
- The N₂ curve is a completed, reproducible small-molecule demonstration, not a
  drug-binding calculation or proof of quantum advantage.
- Reported experimental agreement does not remove discretization, active-space,
  pseudopotential, orbital-choice, or finite-scan uncertainties.
- JSON files under `results/data/` are authoritative; figures are generated views.
