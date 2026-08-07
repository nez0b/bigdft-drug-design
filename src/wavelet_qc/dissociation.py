"""N2 dissociation scan: build the wavelet potential-energy surface from the
BigDFT integral files and extract the zero-independent spectroscopic constants.

This is how the valence-only wavelet method is made comparable to an all-electron
calculation. Absolute total energies cannot be compared across a pseudopotential
boundary (different energy zeros), but the *shape* of the potential-energy curve --
the equilibrium bond length r_e and the harmonic frequency omega_e -- is
independent of the energy zero and can be compared directly, to other methods and to
experiment. (This is exactly why Chou et al. 2023 benchmarked vibrational
frequencies rather than energies; see ``docs/REFERENCES.md``.)

The wavelet curve here is CAS(10,10) full-CI in the 10-orbital valence space the
BigDFT run produced (5 occupied + 5 virtual), evaluated from the vendored
`data/n2_scan/data_*/hpq.out`,`hpqrs.out` at eight bond lengths, plus the analytic
valence ion-ion repulsion. The valence-matched pseudopotential comparator and
all-electron reference curves are built in `reference.pyscf_n2_casci`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .integrals import load_one_body, load_two_body, spin_orbital_count
from .reference import pyscf_fci

DATA = Path(__file__).resolve().parents[2] / "data" / "n2_scan"

# --- physical constants (CODATA) ---
BOHR_PER_ANG = 1.0 / 0.529177210903
HA_TO_CM1 = 219474.6313632
AMU_TO_ME = 1822.888486209
N14_MASS_AMU = 14.003074           # 14N
MU_N2_AU = (N14_MASS_AMU / 2.0) * AMU_TO_ME   # reduced mass of 14N2 in electron masses
Z_VAL_N = 5                         # N valence charge under the HGH pseudopotential

# Experimental N2 constants (Huber & Herzberg).
EXP_RE_ANG = 1.09768
EXP_WE_CM1 = 2358.57

# BigDFT scan: dir name is (half bond length in Angstrom) x 1e4, so R = 2 * z.
SCAN_DIRS = ["data_1000", "data_3000", "data_4000", "data_4488",
             "data_5488", "data_6488", "data_7488", "data_8488"]
N_OCC, N_VIRT = 5, 5                # the scan used norbv=5 -> CAS(10,10)


def bond_length_ang(scan_dir: str) -> float:
    """R (Angstrom) from the directory name, = 2 x (name/1e4)."""
    z = int(scan_dir.split("_")[1]) / 1e4
    return 2.0 * z


def ion_ion(R_ang: float, z_val: int = Z_VAL_N) -> float:
    """Valence ion-ion repulsion Z_val^2 / R (Ha). Matches BigDFT's Ion-Ion energy
    exactly (verified: 25/R at R=1.0976 A gives 12.0530 Ha)."""
    return z_val ** 2 / (R_ang * BOHR_PER_ANG)


def wavelet_energy(scan_dir: str) -> float:
    """CAS(10,10) full-CI total energy (Ha) from the wavelet integrals at one
    geometry: electronic FCI on the 10 spatial orbitals + valence ion-ion."""
    d = DATA / scan_dir
    nso = spin_orbital_count(N_OCC, N_VIRT)          # 20 spin-orbitals
    h1 = load_one_body(str(d / "hpq.out"), nso)
    h2 = load_two_body(str(d / "hpqrs.out"), nso)
    ev = list(range(0, nso, 2))
    h1s = h1[np.ix_(ev, ev)]                          # (10,10) spatial
    h2s = h2[np.ix_(ev, ev, ev, ev)]
    e_elec = pyscf_fci(h1s, h2s, N_OCC + N_VIRT, (N_OCC, N_OCC), e_core=0.0)
    return e_elec + ion_ion(bond_length_ang(scan_dir))


def wavelet_pes() -> tuple[np.ndarray, np.ndarray]:
    """The full wavelet dissociation curve: (R in Angstrom, total energy in Ha)."""
    R = np.array([bond_length_ang(d) for d in SCAN_DIRS])
    E = np.array([wavelet_energy(d) for d in SCAN_DIRS])
    order = np.argsort(R)
    return R[order], E[order]


@dataclass
class Spectroscopy:
    r_e_ang: float
    omega_e_cm1: float
    e_min_ha: float


def fit_spectroscopy(R_ang, E_ha, window=(0.75, 1.55), degree=4) -> Spectroscopy:
    """Fit a polynomial to the near-equilibrium points and read off r_e and omega_e.

    r_e is the position of the minimum; omega_e = sqrt(k / mu) with k the curvature
    at the minimum (in Ha/bohr^2), converted to cm^-1. Only points inside ``window``
    (Angstrom) are used, since the far-repulsive and far-dissociated points would
    bias a low-order polynomial.
    """
    R_ang = np.asarray(R_ang)
    E_ha = np.asarray(E_ha)
    mask = (R_ang >= window[0]) & (R_ang <= window[1])
    Rb = R_ang[mask] * BOHR_PER_ANG                  # fit in bohr -> k in Ha/bohr^2
    E = E_ha[mask]
    coeffs = np.polyfit(Rb, E, degree)
    p = np.poly1d(coeffs)
    dp = p.deriv()
    # minimum: the root of dp inside the fitted range, real, with p'' > 0
    roots = dp.roots
    cands = [r.real for r in roots
             if abs(r.imag) < 1e-9 and Rb.min() < r.real < Rb.max()
             and p.deriv(2)(r.real) > 0]
    if not cands:
        raise RuntimeError("no interior minimum found in the fitted window")
    r_min_bohr = min(cands, key=lambda r: p(r))
    k = float(p.deriv(2)(r_min_bohr))                # Ha/bohr^2
    omega = np.sqrt(k / MU_N2_AU) * HA_TO_CM1        # cm^-1
    return Spectroscopy(
        r_e_ang=r_min_bohr / BOHR_PER_ANG,
        omega_e_cm1=float(omega),
        e_min_ha=float(p(r_min_bohr)),
    )
