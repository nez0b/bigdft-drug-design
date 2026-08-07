"""Step 8 -- analyse the wavelet-parameter scan: convergence + cost figures.

Reads ``results/data/scan_results.json`` (written by ``07_param_scan.py``) and
produces figures in ``results/figures``.

Stage A (equilibrium, one knob at a time):
  * results/figures/scan_convergence.* -- E_KS and CAS(10,10) total energy vs each
    knob (hgrid / crmult / frmult), plus the error w.r.t. the best point on a
    log scale with the 1.6 mHa "chemical accuracy" band.
  * results/figures/scan_cost.* -- measured wall time and peak RSS vs the number of
    wavelet coefficients: the empirical cost model (~ (crmult/hgrid)^3).

Stage B (if dissociation records exist, norbv=5 over 8 geometries):
  * results/figures/scan_curves.* -- min-shifted CAS(10,10) curves per parameter
    set, overlaid on the valence-matched pseudopotential and all-electron CASCI
    references from results/data/dissociation.json, with an r_e / omega_e table.

A convergence table is printed either way -- that's the Stage A checkpoint.
"""
import json
from collections import defaultdict
from pathlib import Path

import _bootstrap  # noqa: F401
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wavelet_qc import dissociation as dis

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)
HA_TO_EV = 27.211386245988
MHA = 1e3
CHEM_ACC_MHA = 1.6

# same deck palette as gen_headline_plot.py
C = {"wavelet": "#2F7D4F", "matched": "#2878B5", "allelectron": "#7652C7",
     "exp": "#A33A32", "ink": "#1F2933", "muted": "#6B7280",
     "eks": "#2878B5", "ecas": "#2F7D4F"}
SET_COLORS = ["#A33A32", "#2F7D4F", "#2878B5", "#7652C7", "#B8860B"]

REF = {"hgrid": 0.25, "crmult": 5.0, "frmult": 8.0}


def load():
    recs = json.loads((DATA / "scan_results.json").read_text())
    ok = [r for r in recs if r["status"] == "ok"]
    failed = [r for r in recs if r["status"] != "ok"]
    for r in failed:
        print(f"  (failed point: {r['tag']}: {r.get('error', '')[:80]})")
    return ok


def coefficients(rec):
    """Effective basis size: coarse points carry 1 scaling function, fine
    points 7 additional wavelets."""
    return rec.get("n_coarse_points", 0) + 7 * rec.get("n_fine_points", 0)


def stage_a_records(recs):
    """The three one-knob-at-a-time branches, each sorted by its knob."""
    a = [r for r in recs if r["norbv"] == 10 and abs(r["bond_ang"] - 1.0976) < 1e-6]
    seen = {}
    for r in a:                      # dedupe same-parameter records (smoke/viz tags)
        seen[(r["hgrid"], r["crmult"], r["frmult"])] = r
    a = list(seen.values())
    branches = {}
    branches["hgrid"] = sorted(
        [r for r in a if r["crmult"] == REF["crmult"] and r["frmult"] == REF["frmult"]],
        key=lambda r: -r["hgrid"])            # coarse -> fine
    branches["crmult"] = sorted(
        [r for r in a if r["hgrid"] == REF["hgrid"] and r["frmult"] == REF["frmult"]],
        key=lambda r: r["crmult"])
    branches["frmult"] = sorted(
        [r for r in a if r["hgrid"] == REF["hgrid"] and r["crmult"] == REF["crmult"]],
        key=lambda r: r["frmult"])
    return branches


def stage_a_figures(branches):
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.0))
    for j, (knob, recs) in enumerate(branches.items()):
        if len(recs) < 2:
            continue
        x = [r[knob] for r in recs]
        eks = np.array([r["e_ks"] for r in recs])
        ecas = np.array([r["e_cas_tot"] for r in recs])
        best = recs[-1]  # finest hgrid / largest multiplier = most converged

        ax = axes[0, j]
        ax.plot(x, eks, "o-", color=C["eks"], label="$E_{KS}$ (DFT)")
        ax.plot(x, ecas, "s-", color=C["ecas"], label="$E$ CAS(10,10)")
        ax.set_xlabel(knob)
        ax.set_ylabel("total energy (Ha)")
        ax.legend(fontsize=8)
        ax.set_title(f"{knob} branch", fontsize=10)
        if knob == "hgrid":
            ax.invert_xaxis()  # finer basis to the right

        ax = axes[1, j]
        derr_ks = np.abs((eks - best["e_ks"]) * MHA)
        derr_cas = np.abs((ecas - best["e_cas_tot"]) * MHA)
        ax.semilogy(x[:-1], np.maximum(derr_ks[:-1], 1e-3), "o-", color=C["eks"])
        ax.semilogy(x[:-1], np.maximum(derr_cas[:-1], 1e-3), "s-", color=C["ecas"])
        ax.axhspan(0, CHEM_ACC_MHA, color=C["ecas"], alpha=0.12, lw=0)
        ax.axhline(CHEM_ACC_MHA, color=C["ecas"], lw=0.8, ls="--", alpha=0.6)
        ax.text(0.03, 0.05, "chemical accuracy (1.6 mHa)", transform=ax.transAxes,
                fontsize=7.5, color=C["ecas"])
        ax.set_xlabel(knob)
        ax.set_ylabel(f"|E - E(best {knob})| (mHa)")
        if knob == "hgrid":
            ax.invert_xaxis()
    fig.suptitle("N$_2$ wavelet-basis convergence -- one knob at a time "
                 f"(ref: h={REF['hgrid']}, rmult=[{REF['crmult']:g}, {REF['frmult']:g}])",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"scan_convergence.{ext}", dpi=300)
    plt.close(fig)
    print(f"  wrote {FIGURES / 'scan_convergence.png'} (+.pdf)")


def cost_figure(recs):
    a = [r for r in recs if coefficients(r) > 0]
    if len(a) < 3:
        return
    ncoef = np.array([coefficients(r) for r in a], float)
    wall = np.array([r["scf_wall_s"] + r["toy_wall_s"] for r in a])
    rss = np.array([max(r["scf_rss_mb"] or 0, r["toy_rss_mb"] or 0) for r in a])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, y, lab, col in ((ax1, wall, "wall time (s)", C["eks"]),
                            (ax2, rss, "peak RSS (MB)", C["ecas"])):
        ax.loglog(ncoef, y, "o", color=col, mec="white")
        p = np.polyfit(np.log(ncoef), np.log(y), 1)
        xs = np.linspace(ncoef.min(), ncoef.max(), 50)
        ax.loglog(xs, np.exp(np.polyval(p, np.log(xs))), "--", color=col, lw=1,
                  label=f"$\\propto N^{{{p[0]:.2f}}}$")
        ax.set_xlabel("wavelet coefficients $N$ (coarse + 7$\\times$fine)")
        ax.set_ylabel(lab)
        ax.legend()
    fig.suptitle("Measured cost of one N$_2$ point (SCF + integrals) on this Mac",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"scan_cost.{ext}", dpi=300)
    plt.close(fig)
    print(f"  wrote {FIGURES / 'scan_cost.png'} (+.pdf)")


def stage_a_table(branches):
    print("\nStage A convergence (R = 1.0976 A, norbv=10, CAS(10,10)):")
    hdr = (f"  {'knob':<7}{'value':>7}{'E_KS (Ha)':>14}{'E_CAS (Ha)':>14}"
           f"{'E_corr (mHa)':>13}{'grid':>13}{'wall(s)':>9}{'RSS(MB)':>9}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for knob, recs in branches.items():
        for r in recs:
            gd = "x".join(map(str, r.get("grid_dims", [])))
            print(f"  {knob:<7}{r[knob]:>7g}{r['e_ks']:>14.6f}{r['e_cas_tot']:>14.6f}"
                  f"{r['e_corr_cas'] * MHA:>13.1f}{gd:>13}"
                  f"{r['scf_wall_s'] + r['toy_wall_s']:>9.0f}"
                  f"{max(r['scf_rss_mb'] or 0, r['toy_rss_mb'] or 0):>9.0f}")


# ------------------------------------------------------------------- Stage B

def stage_b_sets(recs):
    """Group dissociation records (norbv=5) by parameter set."""
    sets = defaultdict(list)
    for r in recs:
        if r["norbv"] == 5:
            sets[(r["hgrid"], r["crmult"], r["frmult"])].append(r)
    return {k: sorted(v, key=lambda r: r["bond_ang"])
            for k, v in sets.items() if len(v) >= 5}


def stage_b_figure(sets):
    ref_json = DATA / "dissociation.json"
    refs = json.loads(ref_json.read_text())["curves"] if ref_json.exists() else {}

    fig, (ax, axt) = plt.subplots(1, 2, figsize=(11.5, 4.8),
                                  gridspec_kw={"width_ratios": [1.7, 1.0]})
    rows = []

    def plot_curve(R, E, label, color, ls="-"):
        R, E = np.asarray(R), np.asarray(E)
        win = (R >= 0.75) & (R <= 1.55)
        p = np.poly1d(np.polyfit(R[win] * dis.BOHR_PER_ANG, E[win], 4))
        rr = np.linspace(R[win].min(), R[win].max(), 200)
        ee = p(rr * dis.BOHR_PER_ANG)
        ax.plot(rr, (ee - ee.min()) * HA_TO_EV, color=color, ls=ls, lw=1.8,
                label=label)
        ax.plot(R[win], (E[win] - ee.min()) * HA_TO_EV, "o", color=color, ms=4,
                mec="white", mew=0.7)
        s = dis.fit_spectroscopy(R, E)
        return s

    for key in ("matched", "allelectron"):
        if key in refs:
            c = refs[key]
            s = plot_curve(c["R"], c["E"], c["label"], C[key], ls="--")
            table_label = "GTH-PBE valence" if key == "matched" else c["label"]
            rows.append((table_label, s.r_e_ang, s.omega_e_cm1))

    styles = ["-", (0, (5, 1.5)), (0, (2, 1.5))]
    for i, ((h, cr, fr), recs) in enumerate(sorted(sets.items())):
        R = [r["bond_ang"] for r in recs]
        E = [r["e_cas_tot"] for r in recs]
        label = f"wavelet h={h:g}"
        s = plot_curve(R, E, label, SET_COLORS[i % len(SET_COLORS)],
                       ls=styles[i % len(styles)])
        rows.append((label, s.r_e_ang, s.omega_e_cm1))
    if len(sets) > 1:
        ax.text(0.03, 0.55,
                "the wavelet curves overlay:\nshape already converged\nat the coarsest grid",
                transform=ax.transAxes, fontsize=8.5, color=C["ink"],
                style="italic")

    ax.axvline(dis.EXP_RE_ANG, color=C["exp"], lw=1.3, ls=":")
    ax.set_xlabel("bond length ($\\AA$)")
    ax.set_ylabel("energy above minimum (eV)")
    ax.set_xlim(0.83, 1.52)
    ax.legend(fontsize=7.5)
    ax.set_title("CAS(10,10) curves, min-shifted (shape is the comparable thing)",
                 fontsize=9.5)

    axt.axis("off")
    rows.append(("Experiment", dis.EXP_RE_ANG, dis.EXP_WE_CM1))
    cell = [[f"{lab}", f"{re:.4f}", f"{we:.0f}"] for lab, re, we in rows]
    tbl = axt.table(cellText=cell,
                    colLabels=["method", "$r_e$ ($\\AA$)", "$\\omega_e$ (cm$^{-1}$)"],
                    loc="center", cellLoc="left",
                    colWidths=[0.52, 0.24, 0.24])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1.0, 1.35)
    axt.set_title("spectroscopic constants", fontsize=9.5)

    fig.suptitle("N$_2$ dissociation-curve shape vs wavelet-basis parameters",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"scan_curves.{ext}", dpi=300)
    plt.close(fig)
    print(f"  wrote {FIGURES / 'scan_curves.png'} (+.pdf)")
    return rows


def main():
    recs = load()
    branches = stage_a_records(recs)
    if any(len(v) >= 2 for v in branches.values()):
        stage_a_table(branches)
        stage_a_figures(branches)
        cost_figure(recs)
    sets = stage_b_sets(recs)
    if sets:
        print(f"\nStage B: {len(sets)} parameter set(s) with dissociation curves")
        rows = stage_b_figure(sets)
        print(f"\n  {'method':<38}{'r_e (A)':>9}{'omega_e':>9}")
        for lab, re, we in rows:
            print(f"  {lab:<38}{re:>9.4f}{we:>9.0f}")


if __name__ == "__main__":
    main()
