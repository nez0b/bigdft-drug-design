"""Fermions to qubits: active-space reduction + Jordan-Wigner -> SparsePauliOp.

The original NTU code built the qubit Hamiltonian with qiskit-nature 0.7
(``ElectronicEnergy`` + ``JordanWignerMapper``), a stack pinned to qiskit 1.x that
no longer installs. This module does the same job directly on modern qiskit
(``qiskit.quantum_info.SparsePauliOp``), so it (a) runs on the current stack and
(b) shows the mapping explicitly rather than hiding it behind a library call.

Two stages:
  1. ``reduce_to_active_space`` -- freeze the correlation-inert orbitals into a
     mean-field core, leaving a small active-space Hamiltonian (chemist notation).
  2. ``jordan_wigner`` -- second-quantize on 2*m spin-orbitals (blocked
     alpha-then-beta) and map each ladder operator to Pauli strings.

Correctness is checked, not asserted: the exact diagonalization of the resulting
SparsePauliOp is compared against an independent PySCF FCI on the same reduced
integrals (``scripts/03_map_and_diag.py``). A sign error in the JW string would
show up there immediately.
"""

from __future__ import annotations

import numpy as np
from qiskit.quantum_info import SparsePauliOp


# --------------------------------------------------------------------------- #
# Stage 1: active-space reduction (frozen-core)
# --------------------------------------------------------------------------- #
def active_indices(n_occ, n_virt, freeze_list, remove_list):
    """Translate interleaved-spin freeze/remove lists into spatial-orbital sets.

    Returns ``(core, active)`` spatial-index lists: ``core`` = doubly-occupied
    orbitals frozen into the mean field; ``active`` = orbitals kept in the CAS.
    """
    core = sorted({i // 2 for i in freeze_list})
    removed = {i // 2 for i in remove_list}
    active = [p for p in range(n_occ + n_virt) if p not in core and p not in removed]
    return core, active


def reduce_to_active_space(h1, h2, core, active, e_const=0.0):
    """Frozen-core reduction of spatial integrals to the active space.

    ``h1`` ``(M,M)`` and chemist-notation ``h2`` ``(M,M,M,M)`` over all M spatial
    orbitals; ``core``/``active`` are spatial-index lists. Returns
    ``(h1_act, h2_act, e_core)`` where ``e_core`` folds in ``e_const`` (e.g. the
    ion-ion energy) plus the frozen electrons' contribution.
    """
    e_core = float(e_const)
    for i in core:
        e_core += 2.0 * h1[i, i]
        for j in core:
            e_core += 2.0 * h2[i, i, j, j] - h2[i, j, j, i]

    # Mean field the core exerts on the active electrons (added to the one-body).
    vhf = np.zeros_like(h1)
    for i in core:
        vhf += 2.0 * h2[:, :, i, i] - np.transpose(h2[:, i, i, :], (0, 1))
    h1_eff = h1 + vhf

    idx = np.ix_(active, active)
    idx4 = np.ix_(active, active, active, active)
    return h1_eff[idx], h2[idx4], e_core


# --------------------------------------------------------------------------- #
# Stage 2: Jordan-Wigner
# --------------------------------------------------------------------------- #
def _annihilation(j: int, n: int) -> SparsePauliOp:
    """JW image of the fermionic annihilation operator a_j on n qubits.

    a_j = (Z_0 ... Z_{j-1}) (X_j + i Y_j)/2, in qiskit's little-endian qubit order
    (qubit 0 is the rightmost character of a label).
    """
    x = ["I"] * n
    y = ["I"] * n
    for k in range(j):
        x[k] = "Z"
        y[k] = "Z"
    x[j] = "X"
    y[j] = "Y"
    lbl_x = "".join(reversed(x))
    lbl_y = "".join(reversed(y))
    return SparsePauliOp.from_list([(lbl_x, 0.5), (lbl_y, 0.5j)])


def jordan_wigner(h1_act, h2_act, e_core=0.0) -> SparsePauliOp:
    """Second-quantized active-space Hamiltonian mapped to a qubit operator.

    Spin-orbitals are blocked: qubits ``0..m-1`` are alpha spatial orbitals,
    ``m..2m-1`` beta. The Hamiltonian (chemist notation ``[pq|rs]``) is

        H = sum_{s,pq} h1[p,q] a+_{p,s} a_{q,s}
          + 1/2 sum_{s,t,pqrs} h2[p,q,r,s] a+_{p,s} a+_{r,t} a_{s',t} a_{q,s}
          + e_core

    with the usual chemist-notation pairing (p,q) and (r,s) on the same electron.
    """
    m = h1_act.shape[0]
    n = 2 * m

    a = [_annihilation(j, n) for j in range(n)]
    adag = [op.adjoint() for op in a]

    def so(p, spin):  # spatial p, spin 0=alpha/1=beta -> blocked spin-orbital
        return p + spin * m

    terms: list[SparsePauliOp] = []

    # one-body
    for spin in (0, 1):
        for p in range(m):
            for q in range(m):
                c = h1_act[p, q]
                if abs(c) < 1e-14:
                    continue
                terms.append(c * (adag[so(p, spin)] @ a[so(q, spin)]))

    # two-body (chemist notation)
    for s1 in (0, 1):
        for s2 in (0, 1):
            for p in range(m):
                for q in range(m):
                    for r in range(m):
                        for ss in range(m):
                            c = h2_act[p, q, r, ss]
                            if abs(c) < 1e-14:
                                continue
                            op = (adag[so(p, s1)] @ adag[so(r, s2)]
                                  @ a[so(ss, s2)] @ a[so(q, s1)])
                            terms.append(0.5 * c * op)

    H = SparsePauliOp.from_list([("I" * n, e_core)])
    for t in terms:
        H = H + t
    return H.simplify()


_MAX_DENSE_QUBITS = 14  # 2^14 dense complex matrix ~ 4 GB; refuse larger


def _guard_dense(n_qubits: int) -> None:
    if n_qubits > _MAX_DENSE_QUBITS:
        raise ValueError(
            f"dense diagonalization of {n_qubits} qubits needs "
            f"{(2 ** n_qubits) ** 2 * 16 / 1e9:.0f} GB; pick a smaller active "
            f"space (raise the MP2 threshold) or use SQD instead.")


def exact_ground_energy(qubit_op: SparsePauliOp) -> float:
    """Lowest eigenvalue of the qubit Hamiltonian (dense, small active spaces)."""
    _guard_dense(qubit_op.num_qubits)
    return float(np.linalg.eigvalsh(qubit_op.to_matrix())[0])
