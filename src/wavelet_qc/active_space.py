"""MP2-based active-space selection.

An accurate wavelet Hamiltonian over 15+ orbitals is still far too big for a
near-term quantum computer, so the pipeline keeps only the orbitals that carry the
electron correlation. The criterion here is each orbital's percentage share of the
total second-order (MP2) correlation energy: occupied orbitals below the threshold
are frozen, virtual orbitals below it are removed. What survives is a compact
active space (n electrons, m orbitals).

Ported from ``BigQiskit.MP2Delta`` / ``MP2ASMO`` (pure numpy). Indices are 0-based
interleaved spin-orbitals throughout. (The report calls its selection criterion
IEPA1; the first-order pair energy of IEPA1 coincides with MP2, which is what the
released code computes.)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def mp2_delta(n_occ, n_virt, h2_chem, orbital_energies,
              remove_list=(), orbital_list=()) -> float:
    """Second-order Moller-Plesset correlation energy.

    ``h2_chem`` is the interleaved-spin two-electron tensor in chemist notation
    (as loaded); it is converted to physicist notation internally.
    ``orbital_list`` restricts the sum to pairs touching those orbitals (used to
    attribute correlation to a single orbital); ``remove_list`` excludes orbitals.
    """
    n_so = n_occ * 2
    n_sv = n_virt * 2
    n_tot = n_so + n_sv

    o_e = np.zeros(n_tot)
    for i, e in enumerate(orbital_energies):
        o_e[2 * i] = e
        o_e[2 * i + 1] = e

    tbi = np.einsum("ijkm->ikmj", h2_chem)  # chemist -> physicist

    occ_sel = [o for o in orbital_list if o < n_so]
    vir_sel = [o for o in orbital_list if o >= n_so]

    delta = 0.0
    for i in range(n_so):
        for j in range(i + 1, n_so):
            if i in remove_list or j in remove_list:
                continue
            if occ_sel and not (i in occ_sel or j in occ_sel):
                continue
            for a in range(n_so, n_so + n_sv):
                for b in range(a + 1, n_so + n_sv):
                    if a in remove_list or b in remove_list:
                        continue
                    if vir_sel and not (a in vir_sel or b in vir_sel):
                        continue
                    num = (tbi[i, j, a, b] - tbi[i, j, b, a]) ** 2
                    denom = o_e[i] + o_e[j] - o_e[a] - o_e[b]
                    delta += num / denom
    return float(delta)


@dataclass
class ActiveSpace:
    freeze_list: list[int]      # occupied spin-orbitals frozen into the core
    remove_list: list[int]      # virtual spin-orbitals discarded
    occ_pct: list[float]        # per-occupied-orbital MP2 %
    virt_pct: list[float]       # per-virtual-orbital MP2 %
    n_active_elec: int
    n_active_orb: int

    def summary(self) -> str:
        kept_occ = [i for i in range(len(self.occ_pct)) if 2 * i not in self.freeze_list]
        kept_virt = [i for i in range(len(self.virt_pct))
                     if 2 * (len(self.occ_pct) + i) not in self.remove_list]
        return (f"CAS[{self.n_active_elec},{self.n_active_orb}] "
                f"(kept occ {kept_occ}, virt {kept_virt})")


def select_active_space(n_occ, n_virt, h2_chem, orbital_energies,
                        threshold_pct=40.0) -> ActiveSpace:
    """Keep orbitals whose MP2 correlation share exceeds ``threshold_pct``.

    Returns the freeze/remove spin-orbital lists plus the resulting active-space
    size. Mirrors ``BigQiskit.MP2ASMO``.
    """
    total = mp2_delta(n_occ, n_virt, h2_chem, orbital_energies)

    occ_pct = []
    for i in range(n_occ):
        d = mp2_delta(n_occ, n_virt, h2_chem, orbital_energies,
                      orbital_list=[2 * i, 2 * i + 1])
        occ_pct.append(100.0 * d / total)

    virt_pct = []
    for i in range(n_virt):
        idx = (n_occ + i)
        d = mp2_delta(n_occ, n_virt, h2_chem, orbital_energies,
                      orbital_list=[2 * idx, 2 * idx + 1])
        virt_pct.append(100.0 * d / total)

    freeze_list, remove_list = [], []
    for i in range(n_occ):
        if occ_pct[i] < threshold_pct:
            freeze_list += [2 * i, 2 * i + 1]
    for i in range(n_virt):
        if virt_pct[i] < threshold_pct:
            remove_list += [2 * (n_occ + i), 2 * (n_occ + i) + 1]

    n_kept_occ = n_occ - len(freeze_list) // 2
    n_kept_virt = n_virt - len(remove_list) // 2
    return ActiveSpace(
        freeze_list=freeze_list, remove_list=remove_list,
        occ_pct=occ_pct, virt_pct=virt_pct,
        n_active_elec=2 * n_kept_occ,
        n_active_orb=n_kept_occ + n_kept_virt,
    )
