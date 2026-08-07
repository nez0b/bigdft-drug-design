"""Step 4 -- VQE vs exact diagonalization (qiskit core, no qiskit-algorithms).

A hardware-efficient ansatz optimized with the exact statevector estimator + scipy.
For the small active space, VQE reaches the exact ground state within a few mHa,
making the deck's qualitative comparison runnable without overstating chemical
accuracy (1.6 mHa).
(Modern rebuild of ``BigQiskit.HeuristicVQE``, whose V1-primitive stack is gone.)
"""
import _bootstrap  # noqa: F401

from wavelet_qc import mapping, vqe
from wavelet_qc.pipeline import build_active_problem, load_molecule

mol = load_molecule("n2")
prob = build_active_problem(mol, threshold_pct=40.0)

qop = mapping.jordan_wigner(prob.h1_act, prob.h2_act, e_core=prob.e_core_elec)
e_ed = mapping.exact_ground_energy(qop)

print(f"N2 -- {prob.cas.summary()}  ({qop.num_qubits} qubits)")
print("  running VQE (EfficientSU2, 3x COBYLA restarts)...")
best = min((vqe.run_vqe(qop, reps=3, optimizer="COBYLA", maxiter=2500, seed=s)
            for s in range(3)), key=lambda r: r.energy)

print()
print(f"  exact diagonalization    : {e_ed:.6f} Ha (electronic)")
print(f"  VQE (hardware-efficient) : {best.energy:.6f} Ha")
gap_mha = abs(best.energy - e_ed) * 1e3
print(f"  |VQE - exact|            : {gap_mha:.2f} mHa")
print(f"  ansatz parameters        : {len(best.optimal_parameters)}")
print()
print("  The hardware-efficient ansatz plateaus a few mHa above exact -- the")
print("  expected NISQ behaviour (cf. the deck's 11.6 mHa on 109Asp). On this same")
print("  active space SQD reaches full-CI exactly (test_sqd_converges_to_fci),")
print("  where this ansatz plateaus; step 5 shows SQD on a larger space.")
