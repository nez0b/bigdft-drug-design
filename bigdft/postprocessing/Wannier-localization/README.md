### Important Configuration: Executable Paths

Please edit the following two lines in `KS2wannier.py` to point to the actual paths of the `wannier90.x` and `BigDFT2Wannier` executables in your system:
```python
wann_executable = "/data/wannier90-3.1.0/wannier90.x"
b2w_executable = "/data/bigdft/bigdft-suite/build/bigdft/src/BigDFT2Wannier"
```

---

### **Installing `wannier90.x`:**
To install `wannier90.x`, run the following commands in the terminal:
```bash
wget https://github.com/wannier-developers/wannier90/archive/v3.1.0.tar.gz
tar -xvzf v3.1.0.tar.gz
cd wannier90-3.1.0
cp ./config/make.inc.gfortran ./make.inc
make
```

---

### **Installing `BigDFT2Wannier`:**

1. **Modify the `Makefile` for Compilation:**
   - Navigate to the `bigdft-suite/build/bigdft/src/Makefile` file.
   - Find the `LDFLAGS` line and add `-L/opt/anaconda3/lib` to it:
   
   **Original:**
   ```makefile
   LDFLAGS = -L/home2/r09222059/BigDFT_workspace/bigdft-suite/build/install/lib
   ```

   **Modified:**
   ```makefile
   LDFLAGS = -L/opt/anaconda3/lib -L/home2/r09222059/BigDFT_workspace/bigdft-suite/build/install/lib 
   ```

2. **Replace the `BigDFT2Wannier.f90` File:**
   - Download the correct version of `BigDFT2Wannier.f90` from **this folder** where the README is located.
   - Replace the existing file at `bigdft-suite/bigdft/src/BigDFT2Wannier.f90` with this version.

3. **Compile `BigDFT2Wannier`:**
   - Navigate to the directory containing the `Makefile`:
   ```bash
   cd bigdft-suite/build/bigdft/src
   ```
   - Compile the program by running:
   ```bash
   make BigDFT2Wannier
   ```

---

#### Usage:
To run the script directly from the command line:
```bash
python KS2wannier.py
```

---

#### Purpose:
`KS2wannier.py` prepares the necessary input files for both `BigDFT2Wannier` and `wannier90.x`, and also executes both `BigDFT2Wannier` and `wannier90.x`.

---

#### Outputs:
**`*_u.mat`**: The unitary matrix that transforms the KS orbitals to localized Wannier orbitals. It is written in row-major order, where each line contains the real and imaginary parts of a single matrix element.

---
#### Generating Cube Files:

After compiling `plot_wann_fun.f90` with the same `Makefile` as `toy_model.f90`, use the following command in the terminal:
```bash
./plot_wann_fun.x name n_plot
```

Where:
- `name` is the prefix of your `_u.mat` file (without the `_u.mat` extension).
- `n_plot` is the number of cube files you want to generate.

**Example:**
To generate the cube files for the first 10 Wannier functions from the file `1b_2ASP_u.mat`, run:
```bash
./plot_wann_fun.x 1b_2ASP 10
```

---

#### Required Files and Directory Structure:
Make sure you are in the directory containing the following files and folders:
1. **XYZ file**
2. **`input.yaml` file**
3. **`data/` folder** that contains wavefunction files.
4. **`b2w.yaml` file**

**Additional Requirement (Conditional):**
- If `n_occ: "all"` is set in `b2w.yaml`, you must also include the following file:
   - **`log.yaml`** that generated during the execution of BigDFT.

---

#### Dependencies:
If the `yaml` module is not installed, run the following command:
```bash
python -m pip install pyyaml
```

---

#### Requirements:
- **wannier90 Version**: v3.1.0
- **BigDFT Version**: v1.9.4

---

#### More Information:
For more details about `wannier90`, check the following links:
1. [Wannier90 Support](https://wannier.org/support/)
2. [Wannier90 User Guide (PDF)](https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf)
