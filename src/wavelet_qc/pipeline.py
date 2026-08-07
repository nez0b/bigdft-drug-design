"""End-to-end convenience wiring, so the scripts stay short.

Loads a vendored molecule (its BigDFT integral files + ``log.yaml``), builds the
spatial-orbital integrals, selects an MP2 active space, reduces to it, and returns
everything the qubit-mapping / VQE / SQD stages need.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import mapping
from .active_space import ActiveSpace, select_active_space
from .integrals import load_one_body, load_two_body, spin_orbital_count
from .reference import BigDFTLog, hartree_fock_energy, read_bigdft_log

DATA = Path(__file__).resolve().parents[2] / "data"

# Vendored molecules: (n_occ spatial, n_virt spatial).
MOLECULES = {
    "n2": {"n_occ": 5, "n_virt": 10, "dir": "n2"},
}


@dataclass
class Loaded:
    name: str
    n_occ: int
    n_virt: int
    h1_so: np.ndarray          # spin-orbital one-body (interleaved)
    h2_so: np.ndarray          # spin-orbital two-body (interleaved, chemist)
    h1_spatial: np.ndarray     # (m, m)
    h2_spatial: np.ndarray     # (m, m, m, m) chemist
    log: BigDFTLog

    @property
    def e_hf_total(self) -> float:
        return hartree_fock_energy(self.log.n_electrons, self.h1_so, self.h2_so) + self.log.e_ion


def load_molecule(name: str, threshold: float = 0.0) -> Loaded:
    spec = MOLECULES[name]
    n_occ, n_virt = spec["n_occ"], spec["n_virt"]
    d = DATA / spec["dir"]
    nso = spin_orbital_count(n_occ, n_virt)
    h1 = load_one_body(str(d / "hpq.out"), nso, threshold)
    h2 = load_two_body(str(d / "hpqrs.out"), nso, threshold)
    log = read_bigdft_log(str(d / "log.yaml"))
    ev = list(range(0, nso, 2))
    h1_sp = h1[np.ix_(ev, ev)]
    h2_sp = h2[np.ix_(ev, ev, ev, ev)]
    return Loaded(name, n_occ, n_virt, h1, h2, h1_sp, h2_sp, log)


@dataclass
class ActiveProblem:
    cas: ActiveSpace
    h1_act: np.ndarray
    h2_act: np.ndarray
    e_core_elec: float         # frozen-core electronic constant (no ion-ion)
    n_alpha: int
    n_beta: int

    @property
    def e_const_total(self) -> float:
        """Constant to add to an electronic active-space energy for a total energy
        (frozen core + ion-ion) -- filled in by ``build_active_problem``."""
        return self._e_const_total

    def with_ion(self, e_ion: float) -> ActiveProblem:
        self._e_const_total = self.e_core_elec + e_ion
        return self


def build_active_problem(mol: Loaded, threshold_pct: float = 40.0) -> ActiveProblem:
    cas = select_active_space(mol.n_occ, mol.n_virt, mol.h2_so,
                              mol.log.orbital_energies, threshold_pct=threshold_pct)
    core, active = mapping.active_indices(mol.n_occ, mol.n_virt,
                                          cas.freeze_list, cas.remove_list)
    h1a, h2a, ecore = mapping.reduce_to_active_space(
        mol.h1_spatial, mol.h2_spatial, core, active, e_const=0.0)
    na = nb = cas.n_active_elec // 2
    prob = ActiveProblem(cas=cas, h1_act=h1a, h2_act=h2a,
                         e_core_elec=ecore, n_alpha=na, n_beta=nb)
    return prob.with_ion(mol.log.e_ion)
