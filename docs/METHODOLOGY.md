# Methodology

## Why a wavelet basis?

Electronic-structure calculations expand each molecular orbital in a finite basis.
In a Gaussian calculation, basis functions are centered on atoms and accuracy is
improved by choosing a larger hierarchy such as cc-pVDZ → cc-pVTZ. In BigDFT, the
orbitals are represented by compactly supported Daubechies scaling functions and
wavelets on an adaptive real-space grid.

The multiresolution idea is:

- scaling functions represent the smooth, coarse part of an orbital;
- wavelets add localized detail only where the chosen resolution requires it;
- decreasing `hgrids` improves spatial resolution systematically;
- `crmult` controls the coarse support around atoms;
- `frmult` controls the fine-resolution support close to the nuclei.

In three dimensions a fine-grid point has one tensor-product scaling function and
seven tensor-product wavelets: all combinations of φ/ψ except φφφ’s seven
complements. Compact support makes differential operators sparse and prevents
vacuum from carrying the same basis cost as chemically active regions.

For a quantum algorithm, the orbital basis is also a resource decision. After
orbital selection, `m` spatial orbitals normally become `2m` Jordan–Wigner qubits.
The basis therefore connects classical discretization error directly to quantum
width.

## From BigDFT orbitals to electron integrals

For orthonormal orbitals φₚ, the second-quantized electronic Hamiltonian is

```text
H = Σpq hpq a†p aq + 1/2 Σpqrs (pq|rs) a†p a†r as aq + Ecore.
```

The custom `toy_model_cpu.f90` evaluates:

```text
hpq    = <φp | T + Vlocal + Vnonlocal | φq>
(pq|rs)= ∬ φp(r1)φq(r1) 1/|r1-r2| φr(r2)φs(r2) dr1 dr2.
```

The expensive four-index Coulomb term is reduced to repeated Poisson solves. For
one orbital pair, form the real-space product density ρpq(r)=φp(r)φq(r), solve
∇²Vpq(r)=−4πρpq(r) with BigDFT’s Poisson solver, then integrate
`∫ ρrs(r) Vpq(r) dr`. This reuses BigDFT’s wavelet/grid machinery instead of
building an unrelated numerical quadrature.

The emitted text convention is:

```text
hpq.out    p q       value
hpqrs.out  p q r s   value
```

Indices are zero-based interleaved spin orbitals (even α, odd β). The RHF file
stores symmetry-unique α/spatial records; `src/wavelet_qc/integrals.py` restores
matrix symmetry, eightfold real electron-repulsion symmetry, and the allowed spin
blocks.

## Validation ladder

The pipeline does not trust a visually plausible energy. Each representation is
checked before moving downstream:

1. **Completeness:** the equilibrium 15-spatial-orbital fixture has 240 one-body
   and 7,260 two-body records, matching the index combinatorics.
2. **Permutation symmetry:** randomly sampled two-electron permutations agree.
3. **Mean-field reconstruction:** the occupied-orbital HF expression rebuilt from
   `hpq`/`hpqrs` is finite and lies about 0.44 Ha above the BigDFT PBE energy on
   the shared Kohn–Sham orbitals. The gap is expected because HF exact exchange
   and PBE exchange-correlation are different functionals.
4. **Mapping check:** exact diagonalization of the explicit Jordan–Wigner
   `SparsePauliOp` agrees with PySCF FCI on the identical reduced integrals.
5. **Selected-CI limits:** QSCI in the complete fixed-particle configuration
   sector and SQD on the small test problem recover the FCI reference.
6. **Physical curve check:** the N₂ minimum lies near the experimental bond
   length and the fitted curvature is positive.

## MP2-guided active space

The full orbital Hamiltonian is too large for a statevector/VQE demonstration.
The original NTU code attributes the second-order correlation energy to occupied
and virtual orbitals. Orbitals whose percentage contribution falls below a chosen
threshold are frozen (occupied) or removed (virtual).

Frozen occupied orbitals are not simply deleted. Their constant energy and their
mean field on each active electron are folded into `Ecore` and the effective
one-body tensor. The resulting active Hamiltonian is then passed unchanged to
both the qubit mapper and PySCF, so their equality is a meaningful independent
check.

This selection is a pragmatic demonstration, not a generally optimal active-space
method. In particular, changing `crmult` can change which diffuse virtuals fall
among the lowest orbitals; a fixed-size CAS is not variational with respect to that
changing orbital subspace.

## N₂ comparison protocol

### Three curves

The result uses the same eight bond lengths for:

1. **Wavelet:** BigDFT HGH-PBE valence orbitals and integrals, with 10 valence
   electrons in 10 spatial orbitals.
2. **Valence-matched pseudopotential comparator:** PySCF GTH-PBE pseudopotential
   with `gth-dzvp`, followed by RHF and CASCI(10,10).
3. **All electron:** PySCF cc-pVDZ followed by CASCI(10,10). PySCF therefore
   freezes two doubly occupied molecular core orbitals (the N 1s-derived core) and
   correlates the same 10 valence electrons.

The pseudopotential comparator aligns the number of explicitly represented
electrons and removes the all-electron/core mismatch. It does **not** make the
pseudopotentials or energy zeros identical: BigDFT uses HGH-K/PBE, whereas PySCF
uses GTH-PBE here. Nor does it isolate only a mathematical basis-function family:
BigDFT begins from PBE Kohn–Sham orbitals, while the PySCF CASCI helper begins from
RHF orbitals, and implementation/truncation details also differ.

### What can be compared

The valence-pseudopotential and all-electron total energies have different energy
zeros. Therefore each curve is shifted by its own minimum before plotting. We
compare only quantities invariant to that shift:

- equilibrium position `rₑ`;
- local curvature and derived harmonic frequency `ωₑ`;
- qualitative shape near the minimum.

Absolute cross-boundary energy differences are never reported.

### Spectroscopic fit

The committed analysis fits a fourth-degree polynomial in bohr to points between
0.75 and 1.55 Å, finds the interior stationary point with positive second
derivative, and evaluates

```text
ωe = sqrt(k / μ) × (Hartree-to-wavenumber conversion),
```

using the reduced mass of ¹⁴N₂. The quartic is an interpolation device for this
small archived dataset, not a global potential model. `ωₑ` is substantially more
sensitive than `rₑ` to the coarse 0.2 Å spacing, so it is labeled preliminary.

## QSCI and SQD

Both methods use bit rows in qiskit-addon-sqd’s
`[beta_0…beta_m-1 | alpha_0…alpha_m-1]` convention.

- **QSCI:** particle-number postselect the samples, count unique configurations,
  keep the most frequent `K`, and diagonalize the fermionic Hamiltonian in that
  selected subspace.
- **SQD:** recover configurations using current orbital occupations, postselect,
  draw several batches, diagonalize each, and feed refined occupations into the
  next recovery iteration.

The committed result samples the exact PySCF FCI coefficient distribution. This
avoids a dense `2^(2m)` qubit matrix and provides a controlled oracle for testing
the classical selected-CI stages. It does not demonstrate state preparation,
noise, mitigation, or quantum hardware. Those would replace the sampler while
leaving the postselection/diagonalization interface intact.

Selected primary literature for the wavelet basis, BigDFT representation, and
prior molecular benchmarks is listed in [REFERENCES.md](REFERENCES.md).

## Scaling boundary

For `n` spatial orbitals, the number of symmetry-unique real two-electron
integrals is

```text
M = n(n+1)/2
Nunique = M(M+1)/2.
```

For 132 spatial orbitals, `Nunique = 38,531,031`; a dense 264-spin-orbital float64
tensor is 38.9 GB. The old 109Asp archive stopped after 100,000 records (0.26%) and
belongs to a different run than the later recovered complete file. It is used only
to explain the resource wall. No physical 109Asp energy is derived from it.
