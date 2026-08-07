# Result artifacts

`data/` contains machine-readable values; `figures/` contains generated plots.
The source of each primary artifact is:

| artifact | producer |
|---|---|
| `data/dissociation.json` | `scripts/06_dissociation.py` |
| `data/sqd_qsci.json` | `scripts/05_sqd_qsci.py` |
| `data/quri_qsci_diagnostic.json` | archived QURI-QSCI cross-stack diagnostic; unresolved negative result |
| `data/scan_results.json` | `scripts/07_param_scan.py` (archived completed runs) |
| `figures/n2_dissociation_headline.*` | `scripts/gen_headline_plot.py` |
| `figures/scan_convergence.*`, `scan_cost.*`, `scan_curves.*` | `scripts/08_scan_analysis.py` |
| `figures/wavelet_basis_1d.*`, `wavelet_regions.*` | `scripts/10_wavelet_figures.py` |

`data/qee_resources.json` and `data/asp109_mp2_archived.json` are archived
resource-analysis outputs. The latter required the non-vendored complete 109Asp
file and is not part of the default laptop reproduction path.

Do not treat PNG text as the authoritative value; use the JSON and regenerate the
figure. `RESULTS.md` states which quantities are comparable and which are only
diagnostic.
