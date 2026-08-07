# Report status and errata

The PDFs in this directory are preserved presentation artifacts. They are useful
for the pedagogical story and project history, but the maintained Markdown and
machine-readable results in this repository take precedence.

## Wavelet-QC primer

The primer draft correctly explains the wavelet basis, electron-integral pipeline,
active-space mapping, and why pseudopotential/all-electron absolute energies must
not be subtracted. Its “matched-pseudopotential reference” subsection overstates
the comparison: it says the PySCF GTH-PBE and BigDFT calculations use the same
pseudopotential and that the residual is purely a basis-set effect.

The archived BigDFT log identifies an **HGH-K/PBE** pseudopotential; the PySCF
comparator uses **GTH-PBE** with `gth-dzvp`, and it starts from RHF rather than
BigDFT PBE Kohn–Sham orbitals. The comparator aligns the ten explicitly treated
valence electrons and broad pseudopotential treatment, but neither the exact
pseudopotential nor the energy zero. See
[`docs/METHODOLOGY.md`](../METHODOLOGY.md) for the corrected protocol.

## Zapata interview deck

The deck deliberately captures the project at interview time. At that point the
matched N₂ dissociation comparison had not been completed. The later results are
in [`RESULTS.md`](../../RESULTS.md), including the final curve comparison and the
modern QSCI/SQD oracle study.
