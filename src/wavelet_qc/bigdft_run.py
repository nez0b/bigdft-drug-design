"""Run fresh BigDFT SCF + toy_model integral-extraction points on this machine.

Everything upstream of the integral files in this repo was originally vendored
from cluster runs. This module drives a compatible native BigDFT install to
produce new points on demand, which is what a wavelet-
parameter study needs: change ``hgrids`` (grid spacing) or ``rmult = [crmult,
frmult]`` (coarse/fine region multipliers), re-run the SCF, re-extract
``hpq.out`` / ``hpqrs.out``, and harvest the energies.

One scan point = one self-contained run directory under ``scan_runs/<tag>/``:

    input.yaml + posinp.xyz          (written from a template, knobs substituted)
    hpq.inp / hpqrs.inp              (index lists from `gen_index`, replacing the
                                      lost-source Linux `genindex` binary)
    bigdft  ->  data/wavefunction-*  data/virtuals-*  log.yaml     (SCF stage)
    toy_model_mac.x  ->  hpq.out  hpqrs.out                        (integral stage)

The wavefunction files are grid-dependent (their size changes with hgrid/rmult),
so run directories are never shared between parameter sets, and the ~60 MB
``data/`` is deleted after harvesting unless ``keep_data=True`` (needed only for
orbital visualisation). Results append to ``results/data/scan_results.json``; a tag that
is already recorded is skipped, so scans are resumable.

Resource accounting: each stage runs under ``/usr/bin/time -l`` so the record
carries measured wall time and peak RSS -- the scan validates its own cost model
(points scale ~ (crmult/hgrid)^3).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from .dissociation import ion_ion
from .integrals import load_one_body, load_two_body
from .reference import hartree_fock_energy, pyscf_fci, read_bigdft_log

PROJECT_ROOT = Path(__file__).resolve().parents[2]
BIGDFT_PREFIX = Path(
    os.environ.get("BIGDFT_PREFIX", PROJECT_ROOT / "bigdft" / "local" / "install")
).expanduser()
BIGDFT_BIN = Path(
    os.environ.get("BIGDFT_BIN", BIGDFT_PREFIX / "bin" / "bigdft")
).expanduser()
TOY_MODEL_BIN = Path(
    os.environ.get(
        "BIGDFT_TOY_MODEL",
        PROJECT_ROOT / "bigdft" / "macos" / "toy_model_mac.x",
    )
).expanduser()
SCAN_RUNS = Path(
    os.environ.get("WAVELET_QC_SCAN_DIR", PROJECT_ROOT / "scratch" / "scan_runs")
).expanduser()
RESULTS_JSON = PROJECT_ROOT / "results" / "data" / "scan_results.json"

N_OCC = 5  # N2 valence occupied spatial orbitals (10 electrons, HGH pseudo)

# Template mirrors the verified macOS N2 reference run;
# only the scanned knobs are substituted. `output_wf: 21` + `orbitals: text` are
# what make BigDFT write the plain-text wavefunction files toy_model reads.
INPUT_TEMPLATE = """\
dft:
  hgrids: {hgrid}
  rmult: [{crmult}, {frmult}]
  ixc: 11
  qcharge: 0
  nspin: 1
  norbv: {norbv}
  nvirt: {norbv}
  nplot: {nplot}
  nrepmax: 10
  gnrm_cv: 1.e-4
  itermax: 50
  gnrm_cv_virt: 1.e-4
  itermax_virt: 50
  disablesym: Yes
lin_general:
  output_wf: 21
output:
  verbosity: 2
  orbitals: {orbitals_fmt}
perf:
  exctxpar: BC
"""


# ---------------------------------------------------------------- genindex port

def gen_index(n_occ: int, n_virt: int, rhf: bool = True):
    """Rebuild the ``hpq.inp`` / ``hpqrs.inp`` index lists.

    Port of the lost-source ``genindex`` C program (only a Linux x86-64 binary
    survives). Verified byte-identical to the shipped reference files for
    (n_occ=5, n_virt=10, RHF): 240 hpq lines, 7260 hpqrs lines.

    Convention (0-based interleaved spin-orbitals, even=alpha odd=beta):
      * hpq   -- same-spin pairs p <= q (same parity), both spins listed.
      * hpqrs -- RHF: alpha block only (all-even indices); pairs (p<=q) ordered
        by p then q, and the (r,s) pair runs from the (p,q) pair onward
        (triangular in the pair index), giving C(npair+1, 2) lines.
    """
    if not rhf:
        raise NotImplementedError("only the RHF index pattern is known/used")
    n_spin = 2 * (n_occ + n_virt)
    hpq = [(p, q) for p in range(n_spin) for q in range(p, n_spin, 2)]
    pairs = [(p, q) for p in range(0, n_spin, 2) for q in range(p, n_spin, 2)]
    hpqrs = [pq + rs for i, pq in enumerate(pairs) for rs in pairs[i:]]
    return hpq, hpqrs


def write_index_files(run_dir: Path, n_occ: int, n_virt: int) -> None:
    hpq, hpqrs = gen_index(n_occ, n_virt)
    (run_dir / "hpq.inp").write_text(
        "".join(f"{p}\t{q}\n" for p, q in hpq))
    (run_dir / "hpqrs.inp").write_text(
        "".join(f"{p}\t{q}\t{r}\t{s}\n" for p, q, r, s in hpqrs))


# ------------------------------------------------------------------ energetics

def energies_from_integrals(run_dir: Path, n_virt: int, bond_ang: float,
                            n_cas_spatial: int = 10) -> dict:
    """Rebuild HF and CAS full-CI energies from a run's integral files.

    Returns electronic *and* total (+ valence ion-ion) energies. The CAS space
    is the first ``n_cas_spatial`` spatial orbitals (5 occupied + lowest
    virtuals), i.e. CAS(10,10) by default -- matching the vendored scan data.
    """
    nso = 2 * (N_OCC + n_virt)
    h1 = load_one_body(str(run_dir / "hpq.out"), nso)
    h2 = load_two_body(str(run_dir / "hpqrs.out"), nso)
    e_hf = hartree_fock_energy(2 * N_OCC, h1, h2)
    ev = list(range(0, 2 * n_cas_spatial, 2))
    e_cas = pyscf_fci(h1[np.ix_(ev, ev)], h2[np.ix_(ev, ev, ev, ev)],
                      n_cas_spatial, (N_OCC, N_OCC), e_core=0.0)
    e_nn = ion_ion(bond_ang)
    return {
        "e_hf_elec": e_hf, "e_cas_elec": e_cas,
        "e_hf_tot": e_hf + e_nn, "e_cas_tot": e_cas + e_nn,
        "e_corr_cas": e_cas - e_hf,
    }


# ------------------------------------------------------------------ subprocess

def _run_timed(cmd, cwd: Path, env: dict, timeout_s: float, logname: str):
    """Run cmd, return (wall_s, peak_rss_mb). Peak RSS is sampled by polling
    ``ps`` -- NOT ``/usr/bin/time -l``, because that is a SIP-protected binary
    and macOS strips DYLD_LIBRARY_PATH from its environment, so the BigDFT
    binaries (which need it for libyaml) would die on launch. Raises on
    non-zero exit (with the log tail in the message) or timeout."""
    logfile = cwd / logname
    command_label = Path(cmd[0]).name
    cwd_label = f"{cwd.parent.name}/{cwd.name}"
    t0 = time.time()
    peak_rss_kb = 0
    with open(logfile, "w") as fh:
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=fh,
                                stderr=subprocess.STDOUT)
        while True:
            try:
                proc.wait(timeout=0.5)
                break
            except subprocess.TimeoutExpired:
                pass
            ps = subprocess.run(
                ["ps", "-o", "rss=", "-p", str(proc.pid)],
                check=False,
                capture_output=True,
                text=True,
            )
            try:
                peak_rss_kb = max(peak_rss_kb, int(ps.stdout.strip()))
            except ValueError:
                pass
            if time.time() - t0 > timeout_s:
                proc.kill()
                proc.wait()
                raise RuntimeError(f"{command_label} timed out after {timeout_s}s")
    wall = time.time() - t0
    if proc.returncode != 0:
        tail = "".join(logfile.read_text(encoding="utf-8").splitlines(keepends=True)[-20:])
        raise RuntimeError(
            f"{command_label} exited {proc.returncode} in {cwd_label}"
            f"\n--- log tail ---\n{tail}"
        )
    return wall, round(peak_rss_kb / 1024) or None


def _parse_grid(log_path: Path) -> dict:
    """Grid-size / cost metadata out of log.yaml (regex, no YAML parse -- the
    log contains fortran-style literals PyYAML chokes on)."""
    text = log_path.read_text(encoding="utf-8")
    out = {}
    m = re.search(r"Grid Spacing Units\s*:\s*\[\s*(\d+),\s*(\d+),\s*(\d+)\s*\]", text)
    if m:
        out["grid_dims"] = [int(g) for g in m.groups()]
    for key, pat in [
        ("n_coarse_points", r"Coarse resolution grid:\s*\n\s*No\. of segments\s*:\s*\d+\s*\n\s*No\. of points\s*:\s*(\d+)"),
        ("n_fine_points", r"Fine resolution grid:\s*\n\s*No\. of segments\s*:\s*\d+\s*\n\s*No\. of points\s*:\s*(\d+)"),
        ("est_mem_peak_mb", r"Estimated Memory Peak \(MB\)\s*:\s*(\d+)"),
    ]:
        m = re.search(pat, text)
        if m:
            out[key] = int(m.group(1))
    gnrms = re.findall(r"gnrm\s*:\s*([\d.Ee+-]+)", text)
    if gnrms:
        out["final_gnrm"] = float(gnrms[-1])
    return out


# -------------------------------------------------------------------- the point

def load_results() -> list[dict]:
    if RESULTS_JSON.exists():
        return json.loads(RESULTS_JSON.read_text())
    return []


def _save_result(rec: dict) -> None:
    RESULTS_JSON.parent.mkdir(parents=True, exist_ok=True)
    results = [r for r in load_results() if r["tag"] != rec["tag"]]
    results.append(rec)
    RESULTS_JSON.write_text(json.dumps(results, indent=2))


def run_point(tag: str, half_bond_ang: float, hgrid: float, crmult: float,
              frmult: float, norbv: int = 10, keep_data: bool = False,
              timeout_s: float = 3600, omp_threads: int = 8,
              force: bool = False, orbitals_fmt: str = "text",
              nplot: int = 0) -> dict:
    """Run one (geometry, wavelet-parameter) point end to end and record it.

    Returns the result record (cached from results/data/scan_results.json when the tag
    already succeeded, unless ``force``). On failure the record carries
    ``status: failed`` + the error, so a scan continues past bad points.
    """
    if os.environ.get("WAVELET_QC_MEMORY_GUARD") != "1":
        raise RuntimeError(
            "fresh BigDFT runs must use tools/run_with_memory_guard.py; "
            "see INSTALL.md"
        )
    if not BIGDFT_BIN.is_file():
        raise FileNotFoundError(
            f"BigDFT executable not found: {BIGDFT_BIN}; set BIGDFT_PREFIX or BIGDFT_BIN"
        )
    if not TOY_MODEL_BIN.is_file():
        raise FileNotFoundError(
            f"toy_model executable not found: {TOY_MODEL_BIN}; run "
            "bigdft/macos/build_toy_mac.sh or set BIGDFT_TOY_MODEL"
        )

    guard_threads = int(os.environ.get("OMP_NUM_THREADS", omp_threads))
    omp_threads = min(omp_threads, guard_threads)

    if not force:
        for r in load_results():
            if r["tag"] == tag and r.get("status") == "ok":
                return r

    run_dir = SCAN_RUNS / tag
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)

    bond_ang = 2.0 * half_bond_ang
    rec = {
        "tag": tag, "bond_ang": bond_ang, "hgrid": hgrid, "crmult": crmult,
        "frmult": frmult, "norbv": norbv, "omp_threads": omp_threads,
    }

    (run_dir / "input.yaml").write_text(INPUT_TEMPLATE.format(
        hgrid=hgrid, crmult=crmult, frmult=frmult, norbv=norbv,
        orbitals_fmt=orbitals_fmt, nplot=nplot))
    (run_dir / "posinp.xyz").write_text(
        f"2 angstroemd0\nfree\n"
        f"N 0.00000 0.00000 {half_bond_ang:.5f}\n"
        f"N 0.00000 0.00000 {-half_bond_ang:.5f}\n")
    write_index_files(run_dir, N_OCC, norbv)

    env = os.environ.copy()
    env["DYLD_LIBRARY_PATH"] = str(BIGDFT_PREFIX / "lib")
    env["OMP_NUM_THREADS"] = str(omp_threads)

    try:
        # Stage 1: SCF + Davidson virtuals (singleton MPI -- mpirun hits a PMIx
        # clash on this Mac). The 2024 BigDFT sources have a NONDETERMINISTIC
        # uninitialised-heap read that makes the Davidson gradient NaN from
        # iteration 1 in a fraction of runs (all virtual eigenvalues stay 0);
        # converged runs are correct (verified against cluster invariants), so
        # detect the loud failure and retry the SCF stage a few times.
        for attempt in range(1, 4):
            wall, rss = _run_timed([str(BIGDFT_BIN)], run_dir, env, timeout_s,
                                   "bigdft.stdout")
            rec["scf_wall_s"], rec["scf_rss_mb"] = round(wall, 1), rss

            log = read_bigdft_log(run_dir / "log.yaml")
            rec["virt_energies"] = (log.virt_energies[-norbv:]
                                    if log.virt_energies else [])
            davidson_dead = rec["virt_energies"] and all(
                abs(e) < 1e-14 for e in rec["virt_energies"])
            if not davidson_dead:
                break
            rec["davidson_retries"] = attempt
            shutil.rmtree(run_dir / "data", ignore_errors=True)
        else:
            raise RuntimeError(
                "Davidson produced all-zero virtual eigenvalues in 3 attempts "
                "(nondeterministic NaN-gradient bug, see INSTALL)")

        data = run_dir / "data"
        wfs = sorted(data.glob("wavefunction-*"))
        virts = sorted(data.glob("virtuals-*"))
        if len(wfs) != N_OCC or len(virts) != norbv:
            raise RuntimeError(
                f"expected {N_OCC} wavefunction + {norbv} virtual files, "
                f"got {len(wfs)} + {len(virts)}")
        rec["wf_file_mb"] = round(wfs[0].stat().st_size / 2**20, 2)

        rec["e_ks"], rec["e_ion_log"] = log.e_ks, log.e_ion
        rec["occ_energies"] = log.occ_energies[-N_OCC:] if log.occ_energies else []
        rec.update(_parse_grid(run_dir / "log.yaml"))

        # Stage 2: integral extraction.
        wall, rss = _run_timed([str(TOY_MODEL_BIN)], run_dir, env, timeout_s,
                               "toy_model.stdout")
        rec["toy_wall_s"], rec["toy_rss_mb"] = round(wall, 1), rss

        rec.update(energies_from_integrals(run_dir, norbv, bond_ang))
        rec["status"] = "ok"
    except Exception as exc:  # noqa: BLE001 - each failed scan point is a result record
        rec["status"] = "failed"
        rec["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        if not keep_data and (run_dir / "data").exists():
            shutil.rmtree(run_dir / "data")

    _save_result(rec)
    return rec
