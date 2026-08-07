"""Reference energies: reconstruct Hartree-Fock from the integrals, read the
BigDFT metadata from ``log.yaml``, and (optionally) cross-check with PySCF.

The point of this module is *validation without re-running BigDFT*. If the
integral files are correct and complete, the Hartree-Fock energy rebuilt from them
should land a fixed, understandable distance from BigDFT's own converged energy.
For the N2 valence system that distance is ~0.44 Ha -- the difference between the
PBE exchange-correlation functional BigDFT minimised and the exact-exchange HF
energy evaluated on those same Kohn-Sham orbitals.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def hartree_fock_energy(n_occ_spin: int, h1: np.ndarray, h2: np.ndarray) -> float:
    """Closed-/open-shell HF electronic energy in the spin-orbital basis.

    E0 = sum_a h_aa + 1/2 sum_ab ( [aa|bb] - [ab|ba] )   (chemist notation)

    ``n_occ_spin`` is the number of occupied spin-orbitals (= number of electrons).
    Lifted from ``BigQiskit.HartreeFockE``.
    """
    one = sum(h1[a, a] for a in range(n_occ_spin))
    two = 0.0
    for a in range(n_occ_spin):
        for b in range(n_occ_spin):
            j = h2[a, a, b, b]      # Coulomb  [aa|bb]
            k = h2[a, b, b, a]      # exchange [ab|ba]
            two += j - k
    return float(one + two / 2.0)


@dataclass
class BigDFTLog:
    """The handful of numbers a validation needs out of a BigDFT ``log.yaml``."""

    e_ks: float                 # final converged Kohn-Sham (DFT) total energy
    e_ion: float                # ion-ion (nuclear) repulsion
    n_electrons: int
    n_occ: int                  # occupied spatial orbitals
    occ_energies: list[float]   # converged occupied orbital energies
    virt_energies: list[float]  # converged virtual orbital energies

    @property
    def orbital_energies(self) -> list[float]:
        """Occupied + virtual orbital energies, in spatial-orbital order -- the
        input MP2 active-space selection expects (replaces PyBigDFT's
        ``Logfile.evals``)."""
        return self.occ_energies + self.virt_energies


def read_bigdft_log(path: str) -> BigDFTLog:
    """Parse the fields we need from a BigDFT ``log.yaml`` with plain regex.

    Uses the ``Complete list of energy eigenvalues`` section (the final converged
    values) rather than the input-guess block, so no BigDFT/PyBigDFT import is
    needed.
    """
    text = Path(path).read_text(encoding="utf-8")

    def last_float(pattern: str) -> float:
        matches = re.findall(pattern, text)
        if not matches:
            raise ValueError(f"pattern not found in {path}: {pattern!r}")
        return float(matches[-1])

    e_ion = last_float(r"Ion-Ion interaction energy\s*:\s*([-\d.E+]+)")
    e_ks = last_float(r"Energy \(Hartree\)\s*:\s*([-\d.E+]+)")
    n_electrons = int(re.findall(r"Total Number of Electrons\s*:\s*(\d+)", text)[-1])
    n_occ = int(re.findall(r"Total Number of Orbitals\s*:\s*(\d+)", text)[-1])

    occ = [float(x) for x in re.findall(r"e_occ\s*:\s*([-\d.E+]+)", text)]
    virt = [float(x) for x in re.findall(r"e_vrt\s*:\s*([-\d.E+]+)", text)]

    return BigDFTLog(
        e_ks=e_ks, e_ion=e_ion, n_electrons=n_electrons, n_occ=n_occ,
        occ_energies=occ, virt_energies=virt,
    )


def pyscf_n2_casci(bond_length_ang, basis, pseudo=None, ncas=10, nelecas=10) -> float:
    """CASCI(ncas, nelecas) total energy of N2 at a given bond length, in PySCF.

    Two references are built from one helper by toggling ``pseudo``:
      * ``pseudo='gth-pbe'`` with a ``gth-*`` basis -> a valence-matched
        pseudopotential comparator: only the 10 valence electrons are treated.
        BigDFT uses an HGH-K/PBE pseudopotential, not this GTH-PBE parameterization,
        so even these two curves retain separate energy zeros. Other differences
        include basis representation, orbital generation (RHF versus BigDFT PBE),
        implementation, and virtual-space truncation.
      * ``pseudo=None`` with an all-electron basis (e.g. ``cc-pVDZ``) -> the
        standard all-electron reference (CASCI freezes the two 1s cores to keep the
        valence CAS(10,10) consistent with the wavelet space).

    Returns the CASCI total energy (includes nuclear repulsion). Import is local so
    the pure-numpy scripts don't require PySCF.
    """
    from pyscf import gto, mcscf, scf

    half = bond_length_ang / 2.0
    mol = gto.M(atom=f"N 0 0 {half}; N 0 0 {-half}", basis=basis, pseudo=pseudo,
                unit="Angstrom", spin=0, verbose=0)
    mf = scf.RHF(mol).run()
    if not mf.converged:
        raise RuntimeError(f"RHF did not converge at R={bond_length_ang:.4f} A")
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.verbose = 0
    return float(mc.kernel()[0])


def pyscf_fci_state(h1_act, h2_act, n_act_orb, n_act_elec, e_core=0.0):
    """Independent full-CI energy and CI vector for a reduced Hamiltonian.

    The CI vector is useful as a deterministic *classical sampler oracle* when
    validating QSCI/SQD.  Results produced from it are algorithm checks, not
    quantum-hardware measurements.
    """
    from pyscf import fci

    energy, ci = fci.direct_spin1.FCI().kernel(
        h1_act, h2_act, n_act_orb, n_act_elec, ecore=e_core
    )
    return float(energy), ci


def pyscf_fci(h1_act, h2_act, n_act_orb, n_act_elec, e_core=0.0) -> float:
    """Independent full-CI electronic energy from PySCF on an already-reduced
    active-space Hamiltonian -- the external cross-check for the qubit-mapping
    exact diagonalization.

    ``h1_act`` is ``(m, m)``, ``h2_act`` the chemist-notation ``(m, m, m, m)``
    electron-repulsion tensor on the *active* orbitals, ``e_core`` the frozen-core +
    (optionally) nuclear constant. ``n_act_elec`` may be an int or ``(na, nb)``.
    Both solvers must be fed the identical reduced integrals for the comparison to
    mean anything. Import is local so the pure-numpy scripts don't require PySCF.
    """
    energy, _ = pyscf_fci_state(
        h1_act, h2_act, n_act_orb, n_act_elec, e_core=e_core
    )
    return energy
