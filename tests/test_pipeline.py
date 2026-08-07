"""End-to-end checks for the compact N2 pipeline."""
from itertools import combinations

import numpy as np
import pytest

from wavelet_qc import mapping, sqd
from wavelet_qc.integrals import load_two_body, spin_orbital_count
from wavelet_qc.pipeline import build_active_problem, load_molecule
from wavelet_qc.reference import pyscf_fci, pyscf_fci_state


def test_two_body_symmetry_roundtrip():
    """Loaded two-electron tensor obeys the 8-fold permutation symmetry."""
    nso = spin_orbital_count(5, 10)
    h2 = load_two_body("data/n2/hpqrs.out", nso)
    rng = np.random.default_rng(0)
    for _ in range(200):
        p, q, r, s = rng.integers(0, nso, 4)
        v = h2[p, q, r, s]
        assert h2[q, p, r, s] == pytest.approx(v)
        assert h2[p, q, s, r] == pytest.approx(v)
        assert h2[r, s, p, q] == pytest.approx(v)


def test_hf_reconstruction_matches_bigdft():
    """HF rebuilt from the integrals sits the expected ~0.44 Ha above E_KS."""
    mol = load_molecule("n2")
    gap = mol.e_hf_total - mol.log.e_ks
    assert 0.3 < gap < 0.6


def test_jordan_wigner_matches_pyscf_fci():
    """Exact-diag of the JW qubit operator == independent PySCF full-CI."""
    mol = load_molecule("n2")
    prob = build_active_problem(mol, threshold_pct=40.0)
    qop = mapping.jordan_wigner(prob.h1_act, prob.h2_act, e_core=prob.e_core_elec)
    e_ed = mapping.exact_ground_energy(qop)
    e_fci = pyscf_fci(prob.h1_act, prob.h2_act, prob.cas.n_active_orb,
                      (prob.n_alpha, prob.n_beta), e_core=prob.e_core_elec)
    assert e_ed == pytest.approx(e_fci, abs=1e-8)


def test_correlation_energy_is_negative():
    mol = load_molecule("n2")
    prob = build_active_problem(mol, threshold_pct=40.0)
    qop = mapping.jordan_wigner(prob.h1_act, prob.h2_act, e_core=prob.e_core_elec)
    e_corr = (mapping.exact_ground_energy(qop) + mol.log.e_ion) - mol.e_hf_total
    assert e_corr < 0


def test_wavelet_dissociation_bond_length():
    """The wavelet PES (from the integral scan) gives a bond length within ~1% of
    the experimental N2 value -- the zero-independent property is comparable."""
    from wavelet_qc import dissociation as dis
    R, E = dis.wavelet_pes()
    s = dis.fit_spectroscopy(R, E)
    assert abs(s.r_e_ang - dis.EXP_RE_ANG) < 0.03      # within ~0.03 A of 1.098
    assert 2000 < s.omega_e_cm1 < 3000                 # a sane N2 frequency


def test_sqd_converges_to_fci():
    """SQD on ground-state samples reaches full-CI on the same integrals."""
    mol = load_molecule("n2")
    prob = build_active_problem(mol, threshold_pct=40.0)
    e_fci, ci = pyscf_fci_state(prob.h1_act, prob.h2_act, prob.cas.n_active_orb,
                                (prob.n_alpha, prob.n_beta), e_core=prob.e_core_elec)
    samples = sqd.sample_configurations_from_fci(
        ci, prob.cas.n_active_orb, prob.n_alpha, prob.n_beta, 2000, seed=1)
    res = sqd.run_sqd(prob.h1_act, prob.h2_act, samples, prob.n_alpha, prob.n_beta,
                      e_const=prob.e_core_elec, iterations=3, num_batches=5,
                      samples_per_batch=50, seed=2)
    assert res.energy == pytest.approx(e_fci, abs=1e-6)


def test_qsci_full_configuration_space_matches_fci():
    """QSCI becomes FCI when its selected subspace is the full number sector."""
    mol = load_molecule("n2")
    prob = build_active_problem(mol, threshold_pct=40.0)
    n_orb = prob.cas.n_active_orb
    samples = []
    for occupied_beta in combinations(range(n_orb), prob.n_beta):
        for occupied_alpha in combinations(range(n_orb), prob.n_alpha):
            row = np.zeros(2 * n_orb, dtype=bool)
            row[list(occupied_beta)] = True
            row[n_orb + np.array(occupied_alpha)] = True
            samples.append(row)
    result = sqd.run_qsci(
        prob.h1_act,
        prob.h2_act,
        np.asarray(samples),
        prob.n_alpha,
        prob.n_beta,
        e_const=prob.e_core_elec,
    )
    reference = pyscf_fci(prob.h1_act, prob.h2_act, n_orb,
                          (prob.n_alpha, prob.n_beta), e_core=prob.e_core_elec)
    assert result.energy == pytest.approx(reference, abs=1e-9)
