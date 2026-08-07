"""Step 3 -- map to qubits and diagonalize (qiskit + pyscf).

Second-quantize the active-space Hamiltonian, Jordan-Wigner map it to a qubit
operator (``SparsePauliOp``), and diagonalize. The exact ground energy is
cross-checked against an independent PySCF full-CI on the *same* reduced integrals
-- if the JW string had a sign error, the two would disagree.
"""
import _bootstrap  # noqa: F401

from wavelet_qc import mapping
from wavelet_qc.pipeline import build_active_problem, load_molecule
from wavelet_qc.reference import pyscf_fci

mol = load_molecule("n2")
prob = build_active_problem(mol, threshold_pct=40.0)
cas = prob.cas

qop = mapping.jordan_wigner(prob.h1_act, prob.h2_act, e_core=prob.e_core_elec)
e_ed = mapping.exact_ground_energy(qop)
e_fci = pyscf_fci(prob.h1_act, prob.h2_act, cas.n_active_orb,
                  (prob.n_alpha, prob.n_beta), e_core=prob.e_core_elec)

print(f"N2 -- {cas.summary()}")
print(f"  qubit Hamiltonian        : {qop.num_qubits} qubits, {len(qop)} Pauli terms")
print()
print(f"  Jordan-Wigner exact-diag : {e_ed:.8f} Ha (electronic)")
print(f"  PySCF full-CI            : {e_fci:.8f} Ha (electronic)")
diff = abs(e_ed - e_fci)
print(f"  |difference|             : {diff:.2e} Ha")
print(f"  [{'PASS' if diff < 1e-8 else 'FAIL'}] two independent methods agree "
      f"-> the mapping is correct")
print()
e_ed_total = e_ed + mol.log.e_ion
print(f"  + E_ion -> total          : {e_ed_total:.6f} Ha")
print(f"  E_HF total               : {mol.e_hf_total:.6f} Ha")
print(f"  correlation (active)     : {e_ed_total - mol.e_hf_total:+.6f} Ha (< 0)")
