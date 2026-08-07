"""Headline figure: the N2 dissociation curves that make the wavelet method
comparable. Reads results/data/dissociation.json (written by 06_dissociation.py) and draws
the three CAS(10,10) potential-energy curves, each shifted to its own minimum so
their shapes overlay, with r_e / omega_e against experiment.

Saves PNG + PDF into results/figures/.
"""
import json
from pathlib import Path

import _bootstrap  # noqa: F401
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wavelet_qc import dissociation as dis

HA_TO_EV = 27.211386245988
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"

# deck palette
C = {"wavelet": "#2F7D4F", "matched": "#2878B5", "allelectron": "#7652C7",
     "exp": "#A33A32", "ink": "#1F2933", "muted": "#6B7280"}

data = json.loads((DATA / "dissociation.json").read_text())
spec = data["spectroscopy"]
exp = data["experiment"]

fig, (ax, axt) = plt.subplots(
    1, 2, figsize=(10.5, 4.6), gridspec_kw={"width_ratios": [1.7, 1.0]})

order = ["wavelet", "matched", "allelectron"]
for key in order:
    c = data["curves"][key]
    R = np.array(c["R"]); E = np.array(c["E"])
    win = (R >= 0.75) & (R <= 1.55)
    Rw, Ew = R[win], E[win]
    # smooth quartic through the near-equilibrium window
    Rb = Rw * dis.BOHR_PER_ANG
    p = np.poly1d(np.polyfit(Rb, Ew, 4))
    rr = np.linspace(Rw.min(), Rw.max(), 200)
    ee = p(rr * dis.BOHR_PER_ANG)
    e0 = ee.min()
    ax.plot(rr, (ee - e0) * HA_TO_EV, color=C[key], lw=2.0, zorder=3,
            label=c["label"])
    ax.plot(Rw, (Ew - e0) * HA_TO_EV, "o", color=C[key], ms=5,
            mec="white", mew=0.8, zorder=4)
    # mark the fitted r_e
    ax.axvline(spec[key]["r_e"], color=C[key], lw=0.8, ls=":", alpha=0.5, zorder=1)

# experimental r_e
ax.axvline(exp["r_e"], color=C["exp"], lw=1.4, ls="--", zorder=2)
ax.text(exp["r_e"] + 0.006, 1.78, "exp. $r_e$", color=C["exp"], fontsize=9,
        rotation=90, va="top")

ax.set_xlim(0.83, 1.52)
ax.set_ylim(-0.03, 1.9)
ax.set_xlabel("N–N bond length  (Å)", fontsize=11)
ax.set_ylabel("energy above minimum  (eV)", fontsize=11)
ax.set_title("N$_2$ dissociation: three CAS(10,10) curves, each shifted to its "
             "own minimum", fontsize=10.5, color=C["ink"])
ax.legend(fontsize=9, loc="upper center", frameon=True, framealpha=0.95)
ax.grid(True, color="0.9", lw=0.5)
for s in ax.spines.values():
    s.set_color("0.5")

# ---- side table: r_e and omega_e vs experiment ----
axt.axis("off")
rows = [("method", "$r_e$ (Å)", "$\\omega_e$ (cm$^{-1}$)")]
labels = {"wavelet": "Wavelet (BigDFT)", "matched": "GTH-PBE valence",
          "allelectron": "All-electron cc-pVDZ"}
for key in order:
    rows.append((labels[key], f"{spec[key]['r_e']:.3f}",
                 f"{spec[key]['omega_e']:.0f}"))
rows.append(("Experiment", f"{exp['r_e']:.3f}", f"{exp['omega_e']:.0f}"))

tab = axt.table(cellText=rows[1:], colLabels=rows[0], cellLoc="center",
                colWidths=[0.52, 0.22, 0.30], loc="center", bbox=[0.0, 0.30, 1.0, 0.60])
tab.auto_set_font_size(False)
tab.set_fontsize(9.2)
tab.scale(1, 1.5)
for (r, col), cell in tab.get_celld().items():
    cell.set_edgecolor("0.8")
    if r == 0:
        cell.set_text_props(weight="bold", color="white")
        cell.set_facecolor(C["ink"])
    elif r <= 3:
        key = order[r - 1]
        cell.set_text_props(color=C[key] if col == 0 else C["ink"])
    else:
        cell.set_text_props(weight="bold", color=C["exp"])

axt.text(0.0, 0.24, "Zero-independent constants — comparable across the\n"
                    "pseudopotential boundary (see note). $r_e$ agrees to ~1%;\n"
                    "$\\omega_e$ is preliminary (coarse 0.2 Å scan, 5 virtuals).",
         fontsize=8.2, color=C["muted"], va="top", transform=axt.transAxes)
axt.text(0.0, 0.02, "Wavelet from BigDFT integrals via the demo pipeline;\n"
                    "references: PySCF CASCI(10,10); cc-pVDZ freezes two\n"
                    "doubly occupied N 1s-derived molecular orbitals.",
         fontsize=7.6, color=C["muted"], va="bottom", transform=axt.transAxes)

fig.tight_layout()
FIGURES.mkdir(parents=True, exist_ok=True)
fig.savefig(FIGURES / "n2_dissociation_headline.png", dpi=300, bbox_inches="tight")
fig.savefig(FIGURES / "n2_dissociation_headline.pdf", bbox_inches="tight")
print(f"wrote {FIGURES / 'n2_dissociation_headline.png'} (+.pdf)")
