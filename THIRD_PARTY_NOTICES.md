# Third-party notices

## BigDFT-derived Fortran source

`bigdft/postprocessing/toy_model_cpu.f90` begins with the upstream BigDFT notice:
copyright © 2007–2013 BigDFT group and distribution under the GNU General Public
License. The corresponding GNU GPL version 2 text from the compatible BigDFT tree
is included at `licenses/BigDFT-GPL-2.0.txt`.

The file also contains subsequent NTU research modifications and notes. Retaining
the upstream notice and license text does not assign a license to unrelated files
in this repository.

## Python dependencies

NumPy, SciPy, PySCF, Qiskit, qiskit-addon-sqd, Matplotlib, PyWavelets, pytest,
ruff, and uv are third-party projects installed from their own distributions.
They are not vendored here and remain subject to their respective licenses.

## Project-level license status

The repository currently has no single top-level license grant. Public visibility
is not itself permission to reuse code. Anyone preparing a release should first
obtain contributor agreement and select a project license that is compatible with
the GPL-covered Fortran component, or separate that component cleanly.
