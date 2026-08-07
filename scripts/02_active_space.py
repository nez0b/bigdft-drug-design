"""Step 2 -- MP2 active-space selection (numpy only, no qiskit).

15 spatial orbitals is far too many for a near-term quantum computer. Rank each
orbital by its share of the MP2 correlation energy and keep only the important
ones. For N2 the triple bond dominates: the two degenerate pi orbitals and their
antibonding partners carry ~65% each, so a small active space captures the physics.
"""
import _bootstrap  # noqa: F401

from wavelet_qc.pipeline import build_active_problem, load_molecule

mol = load_molecule("n2")

print("N2 -- MP2 correlation share per orbital")
print()
for pct in (40.0, 10.0):
    prob = build_active_problem(mol, threshold_pct=pct)
    cas = prob.cas
    if pct == 40.0:  # print the per-orbital breakdown once
        print("  occupied orbital  MP2 %:",
              "  ".join(f"{p:5.1f}" for p in cas.occ_pct))
        print("  virtual  orbital  MP2 %:",
              "  ".join(f"{p:5.1f}" for p in cas.virt_pct))
        print()
    n_qubits = 2 * cas.n_active_orb
    print(f"  threshold {pct:4.0f}%  ->  {cas.summary()}"
          f"  ->  {n_qubits} qubits (Jordan-Wigner)")
print()
print("  The pi/pi* orbitals (occ 2,3 and virt 0,1) dominate -- the classic N2")
print("  active space -- selected automatically from the wavelet integrals.")
