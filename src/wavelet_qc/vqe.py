"""Variational quantum eigensolver on modern qiskit.

The original ``BigQiskit.HeuristicVQE`` used ``qiskit_algorithms.VQE`` with a V1
``Estimator`` -- both removed from qiskit 2.x. This is the same idea rebuilt on the
qiskit 2 core primitives only (no qiskit-algorithms): a hardware-efficient ansatz,
the exact statevector estimator for expectation values, and a classical optimizer
from scipy. For the small active spaces the wavelet pipeline produces, statevector
VQE reproduces the exact ground state to well within chemical accuracy, which is
the point of the demonstration -- the deck's "VQE matched exact diagonalization"
claim, made runnable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from qiskit.circuit.library import EfficientSU2
from qiskit.quantum_info import SparsePauliOp, Statevector
from scipy.optimize import minimize


@dataclass
class VQEResult:
    energy: float
    optimal_parameters: np.ndarray
    n_iters: int
    history: list[float] = field(default_factory=list)


def run_vqe(qubit_op: SparsePauliOp, *, reps: int = 3, optimizer: str = "COBYLA",
            maxiter: int = 2000, seed: int = 0) -> VQEResult:
    """Minimize <psi(theta)| H |psi(theta)> over a hardware-efficient ansatz.

    Uses ``EfficientSU2`` (RY+RZ rotations + linear CX entanglement). Expectation
    values are exact statevector evaluations (``Statevector.expectation_value``) --
    equivalent to the ideal ``StatevectorEstimator`` primitive but without its
    per-call job overhead, so the optimization loop is fast. Returns the best
    energy and the optimization history.
    """
    n = qubit_op.num_qubits
    ansatz = EfficientSU2(n, reps=reps, entanglement="linear").decompose()

    history: list[float] = []

    def energy(theta):
        bound = ansatz.assign_parameters(theta)
        val = float(Statevector(bound).expectation_value(qubit_op).real)
        history.append(val)
        return val

    rng = np.random.default_rng(seed)
    x0 = rng.uniform(-0.1, 0.1, ansatz.num_parameters)
    res = minimize(energy, x0, method=optimizer, options={"maxiter": maxiter})

    return VQEResult(
        energy=float(res.fun),
        optimal_parameters=np.asarray(res.x),
        n_iters=len(history),
        history=history,
    )
