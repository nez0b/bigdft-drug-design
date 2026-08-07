"""Load BigDFT electron-integral files into dense numpy tensors.

`toy_model.f90` (the BigDFT post-processor) writes the one- and two-electron
integrals over Daubechies-wavelet Kohn-Sham orbitals as sparse text:

    hpq.out    :  p q            h_pq          (one line per surviving index pair)
    hpqrs.out  :  p q r s        h_pqrs        (chemist notation [pq|rs])

Indices are 0-based *interleaved spin-orbital* indices (even = alpha, odd = beta),
and only the symmetry-unique entries are stored -- for a closed-shell RHF system
the file carries only the even/alpha block, from which the beta blocks are
reconstructed here.

The two functions below are the load + symmetry-restoration stage of the original
`BigQiskit.py` (Shih-Kai Chou, NTU), lifted essentially verbatim -- they are pure
numpy. The only changes: they take a path instead of an open handle, they do not
write `.npy` side-effect files, and `restore_two_body` documents the eight
permutation images it fills.
"""

from __future__ import annotations

import numpy as np


def load_one_body(path: str, n_spin_orbitals: int, threshold: float = 0.0) -> np.ndarray:
    """Read ``hpq.out`` into a dense ``(n, n)`` one-electron integral matrix.

    Restores the single permutation symmetry ``h_pq = h_qp``.
    """
    h1 = np.zeros((n_spin_orbitals, n_spin_orbitals))
    with open(path) as fh:
        for line in fh:
            tok = line.split()
            if len(tok) < 3:
                continue
            p, q, val = int(tok[0]), int(tok[1]), float(tok[2])
            if abs(val) <= threshold:
                continue
            h1[p, q] = val
            h1[q, p] = val  # permutation symmetry h_pq = h_qp
    return h1


def load_two_body(path: str, n_spin_orbitals: int, threshold: float = 0.0) -> np.ndarray:
    """Read ``hpqrs.out`` into a dense ``(n, n, n, n)`` two-electron tensor.

    Restores the full eight-fold permutation symmetry of a real chemist-notation
    integral [pq|rs] *and* expands the stored alpha (even-index) block into all
    four spin blocks (the ``p, p+1`` pairs are the alpha/beta partners of one
    spatial orbital). This is the RHF->UHF expansion the downstream qubit
    Hamiltonian expects.

    Memory note: the dense tensor is ``n**4`` float64. This is fine for the N2
    valence system (n=30 -> 6.5 MB) but is why the 70-atom ligand (n=264 ->
    38.9 GB) cannot be materialised -- see ``scripts/09_109asp_wall.py``.
    """
    h2 = np.zeros((n_spin_orbitals,) * 4)
    with open(path) as fh:
        for line in fh:
            tok = line.split()
            if len(tok) < 5:
                continue
            p, q, r, s = (int(tok[i]) for i in range(4))
            val = float(tok[4])
            if abs(val) <= threshold:
                continue
            # Eight-fold permutation symmetry x four spin blocks (even->even/odd).
            for a, b, c, d in (
                (p, q, r, s), (p, q, s, r), (q, p, r, s), (q, p, s, r),
                (r, s, p, q), (r, s, q, p), (s, r, p, q), (s, r, q, p),
            ):
                h2[a, b, c, d] = val
                h2[a, b, c + 1, d + 1] = val
                h2[a + 1, b + 1, c, d] = val
                h2[a + 1, b + 1, c + 1, d + 1] = val
    return h2


def spin_orbital_count(n_occ: int, n_virt: int) -> int:
    """Total interleaved spin-orbitals for ``n_occ`` occupied + ``n_virt`` virtual
    spatial orbitals (each spatial orbital -> 2 spin orbitals)."""
    return 2 * (n_occ + n_virt)


def deinterleave_to_spatial(h1: np.ndarray, h2: np.ndarray):
    """Split interleaved spin-orbital integrals into spatial-orbital blocks.

    Returns ``(h1_a, h2_aa, h1_b, h2_bb, h2_ba)`` where the one-body arrays are
    ``(m, m)`` and the two-body arrays are ``(m, m, m, m)`` for ``m = n/2`` spatial
    orbitals. The two-body arrays are converted from chemist ``[pq|rs]`` order to
    the physicist ``++--`` order (``ijkl -> iklj``) used by qubit-mapping code.

    This is the de-interleave from ``BigQiskit.QiskitElectronicIntegrals``.
    """
    even = [range(0, i, 2) for i in h1.shape]
    odd = [range(1, i, 2) for i in h1.shape]
    even4 = [range(0, i, 2) for i in h2.shape]
    odd4 = [range(1, i, 2) for i in h2.shape]

    h1_a = h1[np.ix_(*even)]
    h1_b = h1[np.ix_(*odd)]
    h2_aa = np.einsum("ijkl->iklj", h2[np.ix_(*even4)])
    h2_bb = np.einsum("ijkl->iklj", h2[np.ix_(*odd4)])
    # beta-alpha block: first index pair beta, second alpha
    h2_ba = np.einsum(
        "ijkl->iklj",
        h2[np.ix_(*([range(1, i, 2) for i in h1.shape] + [range(0, i, 2) for i in h1.shape]))],
    )
    return h1_a, h2_aa, h1_b, h2_bb, h2_ba
