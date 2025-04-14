### Purpose

`KS2wannier.py` prepares the necessary input files for both `BigDFT2Wannier` and `wannier90.x`, executes these tools, and produces the unitary matrix file (`*_u.mat`) containing the Wannier localization transformation.

`toy_model_wannier_u_mat.f90` reads the `*_u.mat` file, constructs the Wannier-localized orbitals, and calculates the one- and two-electron integrals for these orbitals. The output files are:
- `hpq_wannier.out` (one-electron integrals)
- `hpqrs_wannier.out` (two-electron integrals)

---

### How to use

#### One-Time Setup (Do This Once)
1. Install `wannier90.x`.
2. Install `BigDFT2Wannier`.
3. Complile `toy_model_wannier_u_mat.f90` and `plot_wann_fun.f90` using the same Makefile of the regular `toy_model.f90`.
4. Run `python -m pip install pyyaml` in the correct environment.
5. In `KS2wannier.py`, modify the following lines to point to the actual paths of the `wannier90.x` and `BigDFT2Wannier` executables on your system:
   ```python
   wann_executable = "/data/wannier90-3.1.0/wannier90.x"
   b2w_executable = "/data/bigdft/bigdft-suite/build/bigdft/src/BigDFT2Wannier"
   ```

#### Run-Time Preparation (Do for Every Run)
1. Make sure the following files and folders are in your current working folder:
   - `*.xyz`
   - `input.yaml`
   - `log.yaml`
   - `data/` folder containing the wavefunction files.
   - `b2w.yaml`
   - `KS2wannier.py`
   - `toy_model_wannier_u_mat.x`
   - `plot_wann_fun.x`
2. Modify `b2w.yaml`
3. Modify `input.yaml`

#### Execution
1. `python KS2wannier.py`
2. `./toy_model_wannier_u_mat.x`

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

### Modifying `b2w.yaml`
   - `n_wann`: The number of Wannier-localized orbitals you want to generate. This should not be greater than `n_occ+n_virt`.
   - `name`: Used to name the output files.
   - `n_occ`: **Currently only supports "all"**. The number of occupied orbitals used to construct Wannier-localized orbitals.
   - `n_virt`: **Currently only supports 0**. The number of virtual orbitals used to construct Wannier-localized orbitals.
   - `xyz`: .xyz file name
   - `input_yaml_name`: input.yaml file name.
   - `projections`: The specified projections. Note that the number of projections should not be less than n_wann, and if it is greater than n_wann, it needs additional settings. Please check the comments in b2w.yaml to learn more about the format.

---

### Modifying `input.yaml`
   - Set `nvirt` to 0.
   - Set `norbv` to the number of virtual orbitals you want to use to calculate the one- and two- electron integrals.

---

#### Generating Cube Files:

After compiling `plot_wann_fun.f90` with the same `Makefile` as for `toy_model.f90`, use the following command in the terminal:
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

#### Requirements:
- **wannier90 Version**: v3.1.0
- **BigDFT Version**: v1.9.4

---

#### More Information:
For more details about `wannier90`, check the following links:
1. [Wannier90 Support](https://wannier.org/support/)
2. [Wannier90 User Guide (PDF)](https://github.com/wannier-developers/wannier90/raw/v3.1.0/doc/compiled_docs/user_guide.pdf)
