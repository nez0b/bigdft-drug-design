"""QSCI and sample-based quantum diagonalization (SQD) on wavelet integrals.

SQD samples electron configurations from a state prepared on the quantum computer,
then classically diagonalizes the Hamiltonian in the subspace those configurations
span -- refining the subspace over a few self-consistent "configuration recovery"
iterations. Here the sampling source is the exact active-space ground state (its
statevector), so the samples are genuine amplitude-weighted draws rather than the
uniform proxy the original notebook used; the classical core is the maintained
``qiskit-addon-sqd`` package.

Convention (matches ``qiskit-addon-sqd`` and the QRE pipeline's reference impl):
each configuration is a bit-row ``[beta_0..beta_{m-1} | alpha_0..alpha_{m-1}]``;
``solve_fermion(batch, h1, h2)`` returns the electronic energy for spatial
integrals ``h1`` (hcore) and ``h2`` (chemist ERI), to which the constant
(frozen-core + ion-ion) energy is added.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import SparsePauliOp


def ground_state_vector(qubit_op: SparsePauliOp) -> np.ndarray:
    """Exact ground-state statevector of the qubit Hamiltonian (small systems)."""
    from .mapping import _guard_dense
    _guard_dense(qubit_op.num_qubits)
    _eigenvalues, v = np.linalg.eigh(qubit_op.to_matrix())
    return v[:, 0]


def sample_configurations(statevector: np.ndarray, n_orb: int, n_shots: int,
                          seed: int = 0) -> np.ndarray:
    """Draw bit-rows from ``|amplitude|^2`` in SQD ``[beta | alpha]`` layout.

    The JW qubit order is blocked alpha (qubits 0..m-1) then beta (m..2m-1); qiskit
    is little-endian, so bit ``j`` of a basis index is spin-orbital ``j``. We split
    into alpha/beta halves and reorder to ``[beta | alpha]``.
    """
    rng = np.random.default_rng(seed)
    probs = np.abs(statevector) ** 2
    probs = probs / probs.sum()
    idx = rng.choice(len(probs), size=n_shots, p=probs)
    n = 2 * n_orb
    bits = ((idx[:, None] >> np.arange(n)[None, :]) & 1).astype(bool)  # bit j = qubit j
    alpha = bits[:, :n_orb]
    beta = bits[:, n_orb:]
    return np.concatenate([beta, alpha], axis=1)  # [beta | alpha]


def sample_configurations_from_fci(
    ci_vector: np.ndarray,
    n_orb: int,
    n_alpha: int,
    n_beta: int,
    n_shots: int,
    seed: int = 0,
) -> np.ndarray:
    """Sample SQD bit rows from a PySCF FCI coefficient array.

    This avoids constructing a dense ``2**(2*n_orb)`` qubit Hamiltonian.  It is
    deliberately described as a classical oracle in all result metadata: useful
    for checking the QSCI/SQD classical stages, but not a hardware experiment.
    """
    from pyscf.fci import cistring

    strings_alpha = np.asarray(
        cistring.make_strings(range(n_orb), n_alpha), dtype=np.int64
    )
    strings_beta = np.asarray(
        cistring.make_strings(range(n_orb), n_beta), dtype=np.int64
    )
    coefficients = np.asarray(ci_vector).reshape(
        len(strings_alpha), len(strings_beta)
    )
    probabilities = np.abs(coefficients.ravel()) ** 2
    probabilities /= probabilities.sum()

    rng = np.random.default_rng(seed)
    flat = rng.choice(probabilities.size, size=n_shots, p=probabilities)
    alpha_address, beta_address = np.unravel_index(flat, coefficients.shape)
    alpha_strings = strings_alpha[alpha_address]
    beta_strings = strings_beta[beta_address]
    bit_positions = np.arange(n_orb, dtype=np.int64)
    alpha = ((alpha_strings[:, None] >> bit_positions) & 1).astype(bool)
    beta = ((beta_strings[:, None] >> bit_positions) & 1).astype(bool)
    return np.concatenate([beta, alpha], axis=1)


@dataclass
class SQDResult:
    energy: float                 # best (lowest) total energy found
    energy_per_iter: list[float]  # best energy after each recovery iteration
    subspace_dims: list[int]
    n_shots: int


@dataclass
class QSCIResult:
    energy: float
    n_configurations: int
    n_postselected_shots: int


def run_qsci(
    h1_act,
    h2_act,
    samples,
    n_alpha,
    n_beta,
    *,
    e_const=0.0,
    max_configurations=None,
) -> QSCIResult:
    """Diagonalize in the most frequently sampled configuration subspace.

    This is the direct QSCI step: postselect the target particle-number sector,
    rank unique bitstrings by empirical frequency, retain at most
    ``max_configurations``, then solve the restricted fermionic Hamiltonian.
    """
    from qiskit_addon_sqd.fermion import solve_fermion

    samples = np.asarray(samples, dtype=bool)
    n_orb = h1_act.shape[0]
    if samples.ndim != 2 or samples.shape[1] != 2 * n_orb:
        raise ValueError(f"samples must have shape (shots, {2 * n_orb})")

    particle_mask = (
        (samples[:, :n_orb].sum(axis=1) == n_beta)
        & (samples[:, n_orb:].sum(axis=1) == n_alpha)
    )
    postselected = samples[particle_mask]
    if len(postselected) == 0:
        raise ValueError("no samples remain after particle-number postselection")

    configurations, counts = np.unique(postselected, axis=0, return_counts=True)
    order = np.argsort(-counts, kind="stable")
    configurations = configurations[order]
    if max_configurations is not None:
        if max_configurations <= 0:
            raise ValueError("max_configurations must be positive")
        configurations = configurations[:max_configurations]

    energy, _state, _occupancies, _spin = solve_fermion(
        configurations, h1_act, h2_act
    )
    return QSCIResult(
        energy=float(energy) + e_const,
        n_configurations=len(configurations),
        n_postselected_shots=len(postselected),
    )


def run_sqd(h1_act, h2_act, samples, n_alpha, n_beta, *, e_const=0.0,
            iterations=5, num_batches=5, samples_per_batch=100,
            seed=0) -> SQDResult:
    """Self-consistent SQD on a batch of sampled configurations.

    ``h1_act``/``h2_act`` are the active-space spatial integrals (hcore / chemist
    ERI); ``e_const`` is the frozen-core + ion-ion constant added to the electronic
    energy. Returns the convergence of the best energy and the subspace sizes used.
    """
    from qiskit_addon_sqd.configuration_recovery import recover_configurations
    from qiskit_addon_sqd.fermion import solve_fermion
    from qiskit_addon_sqd.subsampling import (
        postselect_by_hamming_right_and_left,
        subsample,
    )

    probs = np.ones(len(samples)) / len(samples)
    # prior occupancies: right half = alpha, left half = beta
    n_orb = h1_act.shape[0]
    occ_alpha = samples[:, n_orb:].mean(axis=0)
    occ_beta = samples[:, :n_orb].mean(axis=0)

    best_per_iter: list[float] = []
    dims: list[int] = []
    best = np.inf

    for it in range(iterations):
        recovered, probs_r = recover_configurations(
            samples, probs, (occ_alpha, occ_beta), n_alpha, n_beta,
            rand_seed=seed + it,
        )
        recovered, probs_r = postselect_by_hamming_right_and_left(
            recovered, probs_r, hamming_right=n_alpha, hamming_left=n_beta,
        )
        if len(recovered) == 0:
            best_per_iter.append(best if np.isfinite(best) else float("nan"))
            dims.append(0)
            continue

        spb = min(samples_per_batch, len(recovered))
        batches = subsample(recovered, probs_r, spb, num_batches, rand_seed=seed + it)

        iter_occ_a = np.zeros(n_orb)
        iter_occ_b = np.zeros(n_orb)
        for batch in batches:
            e_sci, _state, (occ_a, occ_b), _s2 = solve_fermion(batch, h1_act, h2_act)
            best = min(best, float(e_sci) + e_const)
            iter_occ_a += occ_a / len(batches)
            iter_occ_b += occ_b / len(batches)
        # feed refined occupancies back into the next recovery iteration
        occ_alpha, occ_beta = iter_occ_a, iter_occ_b
        best_per_iter.append(best)
        dims.append(max(len(b) for b in batches))

    return SQDResult(
        energy=float(best),
        energy_per_iter=best_per_iter,
        subspace_dims=dims,
        n_shots=len(samples),
    )
