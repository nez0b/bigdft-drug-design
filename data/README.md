# Committed data

Only compact inputs needed for the laptop reproduction are committed.

| path | contents | status |
|---|---|---|
| `n2/` | equilibrium BigDFT log, geometry, 240 `hpq`, 7,260 `hpqrs` records | complete |
| `n2_scan/` | eight bond lengths, 110 `hpq` + 1,540 `hpqrs` records per point | complete for the archived CAS(10,10) scan |
| `109asp/` | provenance/size manifest only | large tensors intentionally omitted |

BigDFT wavefunction files, cube files, build trees, executables, scratch scans,
virtual environments, and the 1.39 GB 109Asp integral file are excluded.

The N₂ text files use zero-based, interleaved spin-orbital indices. See
`docs/METHODOLOGY.md` before consuming them directly.
