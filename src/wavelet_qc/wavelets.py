"""Daubechies scaling-function and wavelet samples used by the explainers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# BigDFT's 16-tap least-asymmetric filter agrees with PyWavelets ``sym8``.
FAMILY = "sym8"


@dataclass(frozen=True)
class Basis1D:
    """One-dimensional scaling function and wavelet on a common bohr axis."""

    x: np.ndarray
    phi: np.ndarray
    psi: np.ndarray
    hgrid: float

    @property
    def support(self) -> float:
        return float(self.x[-1] - self.x[0])


def basis_1d(
    hgrid: float = 0.25,
    level: int = 8,
    family: str = FAMILY,
) -> Basis1D:
    """Sample ``phi`` and ``psi`` with the cascade algorithm, in bohr."""
    import pywt

    wavelet = pywt.Wavelet(family)
    phi, psi, x = wavelet.wavefun(level=level)
    x = (x - 0.5 * (x[0] + x[-1])) * hgrid
    return Basis1D(x=x, phi=phi, psi=psi, hgrid=hgrid)
