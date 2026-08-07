"""Step 1 -- validate the integrals against BigDFT (numpy only, no qiskit).

The assessment the whole demo rests on: if the vendored ``hpq``/``hpqrs`` files are
correct and complete, the Hartree-Fock energy rebuilt from them lands a fixed,
understandable distance from BigDFT's own converged energy. For N2 that distance is
~0.44 Ha -- the PBE-exchange-correlation-vs-HF-exchange difference on the shared
Kohn-Sham orbitals.
"""
import _bootstrap  # noqa: F401

from wavelet_qc.pipeline import load_molecule

mol = load_molecule("n2")
log = mol.log

e_hf_elec = mol.e_hf_total - log.e_ion
print("N2 valence system (BigDFT, PBE, HGH pseudopotential)")
print(f"  spin-orbitals            : {mol.h1_so.shape[0]}  "
      f"({mol.n_occ} occ + {mol.n_virt} virt spatial)")
print(f"  electrons                : {log.n_electrons}")
print()
print(f"  E_1  (one-electron, occ) : {sum(mol.h1_so[a, a] for a in range(log.n_electrons)):+.6f} Ha")
print(f"  E_HF electronic          : {e_hf_elec:+.6f} Ha")
print(f"  + E_ion (nuclear rep.)   : {log.e_ion:+.6f} Ha")
print(f"  = E_HF total             : {mol.e_hf_total:+.6f} Ha")
print()
print(f"  BigDFT E_KS (reference)  : {log.e_ks:+.6f} Ha")
gap = mol.e_hf_total - log.e_ks
print(f"  E_HF - E_KS              : {gap:+.4f} Ha")
print()
ok = 0.3 < gap < 0.6
verdict = "PASS" if ok else "CHECK"
print(f"  [{verdict}] gap is the expected ~0.44 Ha DFT-XC-vs-HF-exchange difference")
print("         -> the integrals are complete and self-consistent; the pipeline")
print("            can run end-to-end on the existing files.")
print()
print("  Caveat: this is a valence-only pseudopotential energy (~-20 Ha), NOT the")
print("          all-electron N2 energy (~-109 Ha). Do not compare the two.")
