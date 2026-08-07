"""Step 10 -- generate the pedagogical basis and adaptive-region figures."""

from pathlib import Path

import _bootstrap  # noqa: F401
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wavelet_qc.wavelets import FAMILY, basis_1d

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "results" / "figures"
HGRID = 0.25
BOHR_TO_ANG = 0.529177210903
R_COARSE_N = 1.37026
R_FINE_N = 0.24549
C = {
    "coarse": "#2878B5",
    "fine": "#A33A32",
    "ink": "#1F2933",
    "muted": "#6B7280",
}


def plot_basis_1d() -> None:
    """Draw the compactly supported 1D generators at the reference spacing."""
    basis = basis_1d(hgrid=HGRID)
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    ax.plot(basis.x, basis.phi, color=C["coarse"], lw=2.2,
            label=r"scaling function $\phi(x)$")
    ax.plot(basis.x, basis.psi, color=C["fine"], lw=2.2,
            label=r"wavelet $\psi(x)$")
    ax.axhline(0.0, color=C["muted"], lw=0.7, ls=":")

    bond = 1.09768 / BOHR_TO_ANG
    ax.set_ylim(-2.05, 1.35)
    ax.annotate("", xy=(bond / 2, -1.62), xytext=(-bond / 2, -1.62),
                arrowprops={"arrowstyle": "<->", "color": C["ink"], "lw": 1.2})
    ax.text(0.0, -1.90, f"N$_2$ bond ({bond:.2f} $a_0$)", ha="center",
            fontsize=9, color=C["ink"])
    ax.set(xlabel="$x$ ($a_0$)", ylabel="amplitude")
    ax.set_title(
        f"Daubechies-16 basis ({FAMILY}) at hgrid = {HGRID} $a_0$: "
        f"support {basis.support:.2f} $a_0$"
    )
    ax.legend(frameon=False)
    fig.tight_layout()
    for extension in ("png", "pdf"):
        fig.savefig(FIGURES / f"wavelet_basis_1d.{extension}", dpi=300)
    plt.close(fig)


def plot_regions() -> None:
    """Draw coarse/fine support masks for reference and coarse settings."""
    panels = ((0.25, 5.0, 8.0), (0.45, 3.5, 6.0))
    z_atom = 0.5488 / BOHR_TO_ANG
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.1))
    for ax, (hgrid, crmult, frmult) in zip(axes, panels, strict=True):
        r_coarse = crmult * R_COARSE_N
        r_fine = frmult * R_FINE_N
        z_max, y_max = z_atom + r_coarse, r_coarse
        yy = np.arange(-y_max - hgrid, y_max + hgrid, hgrid)
        zz = np.arange(-z_max - hgrid, z_max + hgrid, hgrid)
        y_grid, z_grid = np.meshgrid(yy, zz, indexing="ij")
        d1 = np.hypot(y_grid, z_grid - z_atom)
        d2 = np.hypot(y_grid, z_grid + z_atom)
        fine = (d1 < r_fine) | (d2 < r_fine)
        coarse = ((d1 < r_coarse) | (d2 < r_coarse)) & ~fine
        ax.plot(z_grid[coarse], y_grid[coarse], "s", color=C["coarse"],
                ms=2.2, alpha=0.55, mew=0, label="coarse: scaling functions")
        ax.plot(z_grid[fine], y_grid[fine], "s", color=C["fine"], ms=3.2,
                mew=0, label="fine: +7 wavelets per point")
        for z_position in (+z_atom, -z_atom):
            ax.add_patch(plt.Circle((z_position, 0), r_coarse, fill=False,
                                    lw=1.5, ec=C["coarse"]))
            ax.add_patch(plt.Circle((z_position, 0), r_fine, fill=False,
                                    lw=1.5, ec=C["fine"]))
            ax.plot(z_position, 0, "o", color=C["ink"], ms=8, zorder=5)
        ax.set_title(
            f"hgrid = {hgrid} $a_0$, rmult = [{crmult:g}, {frmult:g}]"
        )
        ax.set(xlabel="z ($a_0$, bond axis)", ylabel="y ($a_0$)")
        ax.set_aspect("equal")
        ax.set_xlim(-z_max - 1, z_max + 1)
        ax.set_ylim(-y_max - 1, y_max + 1)
    axes[0].legend(loc="upper left", fontsize=9, framealpha=0.9)
    fig.suptitle(
        "N$_2$ wavelet-grid controls: spacing (hgrid), coarse region (crmult), "
        "fine region (frmult)"
    )
    fig.text(0.5, 0.015, "Schematic 2D support masks; not a molecular-orbital plot.",
             ha="center", color=C["muted"], fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    for extension in ("png", "pdf"):
        fig.savefig(FIGURES / f"wavelet_regions.{extension}", dpi=300)
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plot_basis_1d()
    plot_regions()
    print(f"wrote {FIGURES / 'wavelet_basis_1d.png'} and wavelet_regions.png (+PDF)")


if __name__ == "__main__":
    main()
