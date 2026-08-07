# QSCI and sample-based quantum diagonalization

## What this stage tests

Quantum-selected configuration interaction (QSCI) uses measured bitstrings to
choose a determinant subspace, then diagonalizes the electronic Hamiltonian in
that subspace on a classical computer. Sample-based quantum diagonalization (SQD)
adds particle-number postselection, batches, and iterative configuration recovery.
Both methods move the large eigensolve away from the full qubit Hilbert space, but
their accuracy still depends on whether sampling reveals the important
determinants.

`scripts/05_sqd_qsci.py` isolates and validates this classical selected-CI stage
on the real BigDFT-wavelet N₂ Hamiltonian. It uses a CAS(6e,6o) problem: 12
Jordan–Wigner qubits and 400 determinants in the fixed `(Nα,Nβ)=(3,3)` sector.
The independent PySCF FCI reference electronic energy, including the frozen-core
constant, is `-31.617439026 Ha`.

## Sampling provenance

The samples are independent draws from the squared coefficients of the exact
PySCF FCI vector. This is a **classical FCI sampler oracle**, chosen so that errors
in state preparation and hardware noise do not obscure validation of the QSCI/SQD
bookkeeping. It is not a quantum-device experiment, a noisy-device simulation, or
evidence of quantum advantage.

Configurations use the `qiskit-addon-sqd` layout
`[beta_0 ... beta_m-1 | alpha_0 ... alpha_m-1]`. QSCI postselects the requested
particle numbers, ranks unique configurations by empirical count, and calls the
same restricted fermionic eigensolver used by the SQD implementation.

## Shot-coverage result

The locked run used sampler seed 2 and SQD seed 3:

| oracle shots | unique sampled determinants | QSCI error (mHa) | SQD error (mHa) |
|---:|---:|---:|---:|
| 5,000 | 25 | 5.561 | 5.561 |
| 20,000 | 36 | 5.558 | 5.558 |
| 100,000 | 40 | 1.536 | 1.536 |
| 500,000 | 43 | 0.184 | 0.184 |

The 100,000-shot selected space falls just inside the conventional 1.6 mHa
chemical-accuracy band; 500,000 shots reduce the residual to 0.184 mHa. This is
coverage convergence: rare but energetically important determinants enter the
sampled subspace as the number of oracle draws grows.

QSCI and SQD agree exactly in this controlled run because every oracle sample
already has the correct particle numbers and both paths ultimately diagonalize
the same full set of unique sampled configurations. Configuration recovery cannot
reconstruct a determinant that is absent from these clean samples when every
batch already spans the available set. Their equality is therefore an internal
consistency result, not a claim that recovery is ineffective under noisy sampling.

The machine-readable record is
[`results/data/sqd_qsci.json`](../results/data/sqd_qsci.json).

## Archived QURI-QSCI comparison

An older cross-stack experiment used N₂ at 1.0 Å in 6-31g with a CAS(6e,8o)
Hamiltonian. The archived `qiskit-addon-sqd` path was 2.355 mHa from PySCF CASCI,
but the QURI-QSCI path stayed 3,978.770 mHa away and did not change with requested
subspace size. That QURI result remains unresolved and is deliberately excluded
from the successful pipeline.

The old report also mislabeled an approximately `-97.04 Ha` core/frozen-orbital
constant as “nuclear repulsion”; its own molecular nuclear-repulsion value was
approximately `+25.93 Ha`. The curated record corrects the label without claiming
that the old run was repaired. See
[`results/data/quri_qsci_diagnostic.json`](../results/data/quri_qsci_diagnostic.json).

## What a hardware follow-up must add

A hardware or noisy-simulator study should record the state-preparation circuit,
backend, transpilation settings, mitigation/postselection rules, shot allocation,
seeds, rejected-sample fraction, and confidence intervals. It should compare
subspace coverage and energy against the oracle curve above; replacing the
sampler without recording those quantities would make the two experiments
scientifically incomparable.
