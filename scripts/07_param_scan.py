"""Step 7 -- the wavelet-parameter scan: re-run BigDFT itself, vary the basis.

Every earlier step consumed vendored cluster integrals. This one drives the
compatible native BigDFT build to produce *new* integrals while varying the
three wavelet-basis knobs:

    hgrids  -- the grid spacing h (bohr). The wavelet analogue of a plane-wave
               cutoff: cost ~ h^-3, accuracy improves systematically as h -> 0.
    crmult  -- coarse-region multiplier: scaling functions cover a ball of
               crmult * r_loc around each atom (where the exponential tails live).
    frmult  -- fine-region multiplier: the second-resolution shell near the
               nuclei (7 extra wavelets/point where the orbitals oscillate).

Stage A (default): one knob at a time around the reference (h=0.25, [5,8]) at
the equilibrium geometry, cheapest points first. ~9 runs, roughly an hour.

Stage B (--curves "H,C,F" ...): the 8-geometry dissociation scan (same bond
lengths as the vendored ``data/n2_scan``, norbv=5 -> CAS(10,10)) for each
requested parameter set, to see how basis quality moves r_e / omega_e.

Results accumulate in ``results/data/scan_results.json`` (resumable: finished
tags are skipped). Analysis and figures: ``08_scan_analysis.py``. This script
refuses to start unless it is inside ``tools/run_with_memory_guard.py``.
"""
import argparse

import _bootstrap  # noqa: F401

from wavelet_qc.bigdft_run import run_point

HALF_EQ = 0.5488          # equilibrium half bond length (Angstrom), R = 1.0976
SCAN_HALves = [0.1000, 0.3000, 0.4000, 0.4488,
               0.5488, 0.6488, 0.7488, 0.8488]   # matches data/n2_scan dirs

# (hgrid, crmult, frmult) -- one knob at a time around the reference.
STAGE_A = [
    (0.45, 5.0, 8.0), (0.35, 5.0, 8.0), (0.30, 5.0, 8.0),
    (0.25, 5.0, 8.0),                                   # reference point
    (0.20, 5.0, 8.0),
    (0.25, 4.0, 8.0), (0.25, 6.0, 8.0), (0.25, 7.0, 8.0),  # crmult branch
    (0.25, 5.0, 6.0), (0.25, 5.0, 7.0),                 # frmult branch
]


def tag_for(h, c, f, half, norbv):
    return f"h{h:g}_c{c:g}_f{f:g}_r{round(half * 1e4)}_v{norbv}"


def show(rec):
    if rec["status"] == "ok":
        gd = "x".join(map(str, rec.get("grid_dims", []))) or "?"
        print(f"  {rec['tag']:<34} E_KS={rec['e_ks']:.6f}  "
              f"E_CAS={rec['e_cas_tot']:.6f}  grid {gd:<11} "
              f"wall {rec['scf_wall_s'] + rec['toy_wall_s']:6.0f}s  "
              f"rss {max(rec['scf_rss_mb'] or 0, rec['toy_rss_mb'] or 0):6.0f}MB")
    else:
        print(f"  {rec['tag']:<34} FAILED: {rec['error'][:90]}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--curves", nargs="*", metavar="H,C,F",
                    help="run the 8-geometry dissociation scan for each "
                         "hgrid,crmult,frmult triple (Stage B)")
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()

    if args.curves:
        sets = [tuple(float(x) for x in s.split(",")) for s in args.curves]
        print(f"Stage B: dissociation scan for {len(sets)} parameter set(s), "
              f"{len(SCAN_HALves)} geometries each (norbv=5)")
        for h, c, f in sets:
            print(f"\n-- hgrid {h:g}, rmult [{c:g}, {f:g}] --")
            for half in SCAN_HALves:
                rec = run_point(tag_for(h, c, f, half, 5), half, h, c, f,
                                norbv=5, omp_threads=args.threads)
                show(rec)
    else:
        # Cheapest first: cost ~ (crmult/hgrid)^3.
        points = sorted(STAGE_A, key=lambda t: (t[1] / t[0]) ** 3)
        print(f"Stage A: {len(points)} points at R = {2 * HALF_EQ} A (norbv=10)")
        for h, c, f in points:
            rec = run_point(tag_for(h, c, f, HALF_EQ, 10), HALF_EQ, h, c, f,
                            norbv=10, omp_threads=args.threads)
            show(rec)

    print("\nDone. Analyse with: python scripts/08_scan_analysis.py")


if __name__ == "__main__":
    main()
