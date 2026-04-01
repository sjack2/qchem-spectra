# Q-Chem Computational Spectra Lab

A modular Bash/Python workflow for computing UV-Vis absorption, electronic circular dichroism (ECD), and vibrational circular dichroism (VCD) spectra using Q-Chem density functional theory. Designed for teaching and research, and runs identically on a local workstation or a SLURM-managed HPC cluster.

> **Accompanying article:**  _A Computational Spectroscopy Module for Teaching DFT-Based Spectral Workflows in Graduate Chemistry_ (in preparation for the _Journal of Chemical Education_).

> **License requirement:** Q-Chem is commercial software. You must hold a valid Q-Chem license before using this workflow. Academic site licenses are available from [Q-Chem, Inc.](https://www.q-chem.com/)

---

## Features

- **Six-stage pipeline** from raw XYZ geometry to publication-quality broadened spectra, with parallel Stage 6 branches for electronic (UV-Vis/ECD) and vibrational (IR/VCD) spectroscopy.
- **Dual execution mode:** every Q-Chem-calling script auto-detects whether SLURM is available. Pass `--local` to force direct execution; omit it on a cluster and jobs are submitted via `sbatch`.
- **Two conformer search paths:** Open Babel Confab (fast, force-field-based) or CREST (GFN2-xTB metadynamics, more thorough). Both branches produce the same per-conformer XYZ files consumed by Stage 4.
- **Boltzmann-weighted spectral averaging** with configurable temperature and population threshold.
- **Physically correct broadening:** Gaussian convolution in energy space (eV) for electronic spectra and in wavenumber space (cm-1) for vibrational spectra, with Jacobian correction for ECD. Produces both PNG/PDF plots and CSV data tables.
- **`--dry-run` on every script** to inspect generated Q-Chem inputs without running any calculations.

---

## Pipeline Overview

```
Stage 1                    Gas-phase geometry optimization
  |                        1-qchem-init-opt.sh
  v
Stage 2                    Conformer enumeration
  |                        2-qchem-conf-search.sh  (Confab)
  |                   or   2b-qchem-crest-conf-search.sh (CREST)
  v
Stage 3                    Split conformers into individual XYZ files
  |                        3-qchem-conf-split.sh   (Confab output)
  |                   or   3b-qchem-crest-conf-split.sh  (CREST output)
  v
Stage 4                    Solvent-phase re-optimization (SMD)
  |                        4-qchem-solvent-opt.sh
  v
Stage 5                    Boltzmann weighting & filtering
  |                        5-qchem-boltzmann-weight.sh
  |
  |------------------------------------------+
  v                                          v
Stage 6-tddft                            Stage 6-vcd
  TD-DFT excited-state calculations        Frequency calculations (IR + VCD)
  6-qchem-tddft.sh                         6-qchem-vcd.sh
  |                                          |
  v                                          v
Plot                                     Plot
  Broadened UV-Vis & ECD spectra           Broadened IR & VCD spectra
  qc_ecd_uvvis_tools.py                    qc_vcd_ir_tools.py
```

Stages 1, 4, and 5 are shared. Stages 2/2b and 3/3b are interchangeable: both produce per-conformer XYZ files in `<TAG>/02_conf_search/split_xyz/` consumed by Stage 4. After Boltzmann filtering, the pipeline branches: run `6-qchem-tddft.sh` for electronic spectra, `6-qchem-vcd.sh` for vibrational spectra, or both.

---

## Quick Start

### 1. Install prerequisites

| Software | Version | Purpose |
|----------|---------|---------|
| Q-Chem | 6.0+ | Quantum chemistry engine (commercial license required) |
| Open Babel | >= 3.0 | File conversion & Confab conformer search |
| Python 3 | >= 3.8 | Plotting tools |
| CREST | >= 2.12 | _(optional)_ GFN2-xTB conformer search |

Q-Chem uses shared-memory parallelism (OpenMP threads), so no external MPI installation is required.

### 2. Clone the repository

```bash
git clone https://github.com/sjack2/qchem-spectra.git
cd qchem-spectra
```

### 3. Check Python and install dependencies

First verify that Python 3.8 or later is available:

```bash
python3 --version
```

If the command is not found or the version is below 3.8, install Python before continuing. On an HPC cluster, the easiest approach is [Miniconda](https://docs.conda.io/en/latest/miniconda.html):

```bash
# Download and install Miniconda (Linux, x86-64)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow the prompts, then open a new shell or run: source ~/.bashrc
```

Once Python is available, install the plotting dependencies:

```bash
pip install -r requirements.txt
```

### 4. Make scripts executable

```bash
chmod +x *.sh
```

### 5. Configure for your cluster (HPC users)

Copy the provided template and fill in the paths for your site. This file is read automatically by every Q-Chem-calling script, so you will not need to pass `--qchem-setup` or `--partition` on every invocation.

```bash
cp cluster.cfg.example cluster.cfg
```

Then edit `cluster.cfg` with your site's values, for example:

```bash
QCHEM_SETUP=/path/to/qchem/setqc
CLUSTER_PARTITION=general
# CLUSTER_WALL=04:00:00
```

`cluster.cfg` is gitignored, so site-specific paths are never committed to the repository. The template `cluster.cfg.example` documents all supported variables, including `CLUSTER_LD_LIBRARY_PATH` for sites where Q-Chem depends on Intel or GCC libraries not present on compute nodes by default.

### 6. Set up your molecule

Place a starting XYZ file in `pre_xyz/`. Line 2 must encode the charge and multiplicity:

```
14
charge=0 mult=1
C    0.000000    0.000000    0.000000
...
```

Both `charge=0 mult=1` (key=value) and `0 1` (bare integers) are accepted. If neither is found, defaults are charge=0, mult=1.

### 7. Run the pipeline (local workstation example)

```bash
# Stage 1: Optimize geometry in vacuum
./1-qchem-init-opt.sh --local --cpus 4 pna

# Stage 2: Generate conformers (Confab)
./2-qchem-conf-search.sh --ecut 5 --conf 500 pna
# Stage 3: Split conformers into individual XYZ files
./3-qchem-conf-split.sh pna

# Stage 4: Re-optimize each conformer in solvent
./4-qchem-solvent-opt.sh --local --cpus 4 --solvent water pna

# Stage 5: Boltzmann filter
./5-qchem-boltzmann-weight.sh pna

# --- Electronic spectra branch ---
# Stage 6-tddft: TD-DFT on populated conformers
./6-qchem-tddft.sh --local --cpus 4 --method CAM-B3LYP --roots 30 pna

# Plot UV-Vis & ECD spectra
python3 qc_ecd_uvvis_tools.py pna/05_tddft \
    --bw pna/04_boltzmann/pna_energies.dat \
    --outdir pna --xlim 200 500

# --- Vibrational spectra branch ---
# Stage 6-vcd: Frequency calculations on populated conformers
./6-qchem-vcd.sh --local --cpus 4 --method B3LYP --solvent water pna

# Plot IR & VCD spectra
python3 qc_vcd_ir_tools.py pna/05_vcd \
    --bw pna/04_boltzmann/pna_energies.dat \
    --outdir pna --xlim 800 3500
```

### 8. Run the pipeline (HPC/SLURM example)

On a cluster with SLURM, simply omit `--local`. The scripts detect `sbatch` automatically and submit jobs. With `cluster.cfg` configured, no additional flags are needed:

```bash
./1-qchem-init-opt.sh --cpus 8 pna
# submits SLURM job to the partition set in cluster.cfg; monitor with squeue
```

SLURM jobs use `--cpus-per-task` (shared-memory parallelism). Q-Chem does not require MPI; all parallelism is handled via OpenMP threads within a single task.

---

## Directory Structure

Every script follows the same directory convention. For a molecule tagged `aspirin`:

```
aspirin/
|-- 01_gas_opt/                Stage 1 - gas-phase optimization
|   |-- aspirin.inp            Q-Chem input ($rem/$molecule/$end blocks)
|   |-- aspirin.out            Q-Chem output
|   `-- aspirin.slurm          SLURM script (HPC mode only)
|-- 02_conf_search/            Stage 2/2b - conformer enumeration
|   |-- aspirin.xyz            extracted optimized geometry
|   |-- aspirin.sdf            SDF conversion (Confab path)
|   |-- aspirin_combined.sdf   all conformers (Confab output)
|   |-- crest_conformers.xyz   all conformers (CREST output, Stage 2b)
|   |-- crest.log              CREST log (Stage 2b)
|   |-- split_sdf/             individual SDF files (Stage 3 writes these)
|   `-- split_xyz/             individual conformer XYZ files (Stage 3 or 3b)
|       |-- aspirin_001.xyz
|       |-- aspirin_002.xyz
|       `-- ...
|-- 03_solvent_opt/            Stage 4 - solvent-phase re-optimization
|   |-- aspirin_conf_list.txt  conformer working-dir list (SLURM array input)
|   |-- aspirin_array.slurm    SLURM array job script (HPC mode)
|   |-- aspirin_001/
|   |   |-- aspirin_001.inp
|   |   `-- aspirin_001.out
|   `-- aspirin_002/
|       `-- ...
|-- 04_boltzmann/              Stage 5 - Boltzmann weighting
|   |-- aspirin_energies.dat   full table (CID, E, dE, p)
|   `-- aspirin_bw_labels.dat  conformer IDs above threshold
|-- 05_tddft/                  Stage 6-tddft - TD-DFT excited states
|   |-- aspirin_conf_list.txt  conformer working-dir list (SLURM array input)
|   |-- aspirin_array.slurm    SLURM array job script (HPC mode)
|   |-- aspirin_001/
|   |   |-- aspirin_001.inp
|   |   `-- aspirin_001.out
|   `-- aspirin_002/
|       `-- ...
|-- 05_vcd/                    Stage 6-vcd - frequency calculations
|   |-- aspirin_conf_list.txt  conformer working-dir list (SLURM array input)
|   |-- aspirin_array.slurm    SLURM array job script (HPC mode)
|   |-- aspirin_001/
|   |   |-- aspirin_001.inp
|   |   `-- aspirin_001.out
|   `-- aspirin_002/
|       `-- ...
`-- 06_spectra/                Plotting outputs
    |-- aspirin_uvvis.png
    |-- aspirin_ecd.png
    |-- aspirin_ir.png
    |-- aspirin_vcd.png
    `-- *.csv / *.pdf
```

**Note:** Stage 2 produces a combined SDF ensemble; Stage 3 splits it. If you used the CREST path (Stage 2b), run Stage 3b instead to split `crest_conformers.xyz`. Both splitters write to the same `split_xyz/` directory consumed by Stage 4.

Starting geometries live in `pre_xyz/`:

```
pre_xyz/
|-- pna.xyz
|-- ephedrine.xyz
`-- methyloxirane.xyz
```

---

## CLI Reference

All scripts accept `--help` for full usage. Flags shown with `[default]`.

### 1-qchem-init-opt.sh -- Gas-Phase Geometry Optimization

```
1-qchem-init-opt.sh [OPTIONS] TAG [TAG ...]
1-qchem-init-opt.sh [OPTIONS] --list FILE
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | DFT functional | `B3LYP` |
| `--basis NAME` | `-b` | Basis set | `def2-SVP` |
| `--disp KW` | | Dispersion: `auto`, `none`, `D3BJ` | `auto` |
| `--max-iter N` | | SCF iteration limit | `150` |
| `--cpus N` | `-c` | CPU cores (OpenMP threads) | `12` |
| `--grid KW` | `-g` | Integration grid: `SG1`, `SG2`, `SG3` | `SG3` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--qchem-setup PATH` | | Path to Q-Chem setup script | _auto-detected_ |
| `--list FILE` | | Text file of TAGs (one per line) | |
| `--local` | | Run Q-Chem directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `02:00:00` |

**Output:** `<TAG>/01_gas_opt/<TAG>.out` (Q-Chem output with optimized geometry)

### 2-qchem-conf-search.sh -- Confab Conformer Search

```
2-qchem-conf-search.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--conf N` | Maximum conformers | `1000` |
| `--ecut E` | Energy cutoff (kcal/mol) | `5` |
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo commands without running | |

**Output:** `<TAG>/02_conf_search/<TAG>_combined.sdf`

### 3-qchem-conf-split.sh -- Split Confab Output

```
3-qchem-conf-split.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo actions without running | |

**Input:** `<TAG>/02_conf_search/<TAG>_combined.sdf`
**Output:** `<TAG>/02_conf_search/split_xyz/<TAG>_NNN.xyz`

### 2b-qchem-crest-conf-search.sh -- CREST Conformer Search (alternative to Stage 2)

```
2b-qchem-crest-conf-search.sh [OPTIONS] TAG [TAG ...]
2b-qchem-crest-conf-search.sh [OPTIONS] --list FILE
```

Uses CREST with GFN2-xTB metadynamics to explore the conformational landscape. More thorough than Confab for flexible molecules; requires the standalone `crest` binary (not via conda). The optimized geometry is extracted automatically from the Stage 1 Q-Chem output; alternatively, supply an XYZ path directly or use `--pre-xyz` to fall back to `pre_xyz/`.

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--cpus N` | `-c` | CPU threads for CREST | `4` |
| `--mem MB` | | Memory for SLURM job (MB) | `4096` |
| `--list FILE` | | Text file of TAGs or XYZ paths | |
| `--local` | | Run CREST directly (no SLURM) | |
| `--dry-run` | | Echo actions without running | |
| `--pre-xyz` | | Also search `pre_xyz/` for geometries | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall time _(HPC only)_ | `06:00:00` |

**Input:** `<TAG>/01_gas_opt/<TAG>.out` (geometry extracted automatically) or an explicit `.xyz` path
**Output:** `<TAG>/02_conf_search/crest_conformers.xyz`, `<TAG>/02_conf_search/crest.log`

> **Environment note:** The standalone `crest` binary requires GCC/GFortran runtime libraries. If `crest` was compiled against a specific toolchain, activate the matching conda environment (e.g., `conda activate xtb_rt`) before submitting SLURM jobs so that compute nodes inherit the correct `LD_LIBRARY_PATH`.

### 3b-qchem-crest-conf-split.sh -- Split CREST Output (follows Stage 2b)

```
3b-qchem-crest-conf-split.sh [OPTIONS] TAG [TAG ...]
3b-qchem-crest-conf-split.sh [OPTIONS] --list FILE
```

Reads `crest_conformers.xyz` produced by Stage 2b and splits it into numbered per-conformer XYZ files in the same `split_xyz/` directory consumed by Stage 4.

| Flag | Description | Default |
|------|-------------|---------|
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Echo actions without creating files | |

**Input:** `<TAG>/02_conf_search/crest_conformers.xyz`
**Output:** `<TAG>/02_conf_search/split_xyz/<TAG>_NNN.xyz` (same format as Stage 3)

### 4-qchem-solvent-opt.sh -- Solvent-Phase Re-optimization

```
4-qchem-solvent-opt.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | DFT functional | `B3LYP` |
| `--basis NAME` | `-b` | Basis set | `6-31+G(d)` |
| `--disp KW` | | Dispersion: `auto`, `none`, `D3BJ` | `none` |
| `--solvent NAME` | | SMD solvent keyword | `water` |
| `--max-scf N` | | Max SCF cycles | `150` |
| `--cpus N` | `-c` | CPU cores (OpenMP threads) | `12` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--max-running N` | | Max concurrent array tasks _(HPC only)_ | `10` |
| `--qchem-setup PATH` | | Path to Q-Chem setup script | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run Q-Chem directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `02:00:00` |

**Input:** `<TAG>/02_conf_search/split_xyz/` (from Stage 3 or 3b)
**Output:** `<TAG>/03_solvent_opt/<CID>/<CID>.out` (Q-Chem output per conformer)

On HPC, all conformers are submitted as a single SLURM job array (`<TAG>_array.slurm`) throttled to `--max-running` concurrent tasks. Use `--max-running` to tune cluster load for large ensembles (e.g., CREST output). Set `CLUSTER_MAX_RUNNING` in `cluster.cfg` to apply a site-wide default.

### 5-qchem-boltzmann-weight.sh -- Boltzmann Weighting & Filtering

```
5-qchem-boltzmann-weight.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Description | Default |
|------|-------------|---------|
| `--temp K` | Temperature in Kelvin | `298.15` |
| `--p-cut VAL` | Probability cutoff (e.g., 0.01 = 1%) | `0.01` |
| `--list FILE` | Text file of TAGs | |
| `--dry-run` | Show what would be computed | |

**Input:** `<TAG>/03_solvent_opt/<CID>/<CID>.out` (Q-Chem outputs from Stage 4)
**Output:**
- `<TAG>/04_boltzmann/<TAG>_energies.dat` -- full table: conformer ID, energy (Hartree), relative energy (kcal/mol), Boltzmann probability
- `<TAG>/04_boltzmann/<TAG>_bw_labels.dat` -- conformer IDs above the probability cutoff

### 6-qchem-tddft.sh -- TD-DFT Excited-State Calculations

```
6-qchem-tddft.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | TD-DFT functional | `CAM-B3LYP` |
| `--basis NAME` | `-b` | Basis set | `def2-SVP` |
| `--roots N` | | Number of excited states | `30` |
| `--solvent NAME` | | SMD solvent keyword | `water` |
| `--max-scf N` | | Max SCF cycles | `150` |
| `--cpus N` | `-c` | CPU cores (OpenMP threads) | `12` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--max-running N` | | Max concurrent array tasks _(HPC only)_ | `10` |
| `--qchem-setup PATH` | | Path to Q-Chem setup script | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run Q-Chem directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `04:00:00` |

**Input:** `<TAG>/04_boltzmann/<TAG>_bw_labels.dat` (from Stage 5) + `<TAG>/03_solvent_opt/<CID>/<CID>.out` (from Stage 4)
**Output:** `<TAG>/05_tddft/<CID>/<CID>.out` (Q-Chem TD-DFT output with oscillator strengths and rotatory strengths)

On HPC, all conformers are submitted as a single SLURM job array (`<TAG>_array.slurm`) throttled to `--max-running` concurrent tasks.

### 6-qchem-vcd.sh -- Frequency / IR + VCD Calculations

```
6-qchem-vcd.sh [OPTIONS] TAG [TAG ...]
```

| Flag | Short | Description | Default |
|------|-------|-------------|---------|
| `--method NAME` | `-m` | DFT functional | `B3LYP` |
| `--basis NAME` | `-b` | Basis set | `6-31+G(d)` |
| `--solvent NAME` | | SMD solvent keyword | `water` |
| `--max-scf N` | | Max SCF cycles | `150` |
| `--cpus N` | `-c` | CPU cores (OpenMP threads) | `12` |
| `--mem-per-cpu MB` | | Memory per core (MB) | `2048` |
| `--max-running N` | | Max concurrent array tasks _(HPC only)_ | `10` |
| `--qchem-setup PATH` | | Path to Q-Chem setup script | _auto-detected_ |
| `--list FILE` | | Text file of TAGs | |
| `--local` | | Run Q-Chem directly (no SLURM) | |
| `--dry-run` | | Write inputs without running | |
| `--partition NAME` | | SLURM partition _(HPC only)_ | `general` |
| `--time HH:MM:SS` | | SLURM wall-clock limit _(HPC only)_ | `06:00:00` |

**Input:** `<TAG>/04_boltzmann/<TAG>_bw_labels.dat` (from Stage 5) + `<TAG>/03_solvent_opt/<CID>/<CID>.out` (from Stage 4)
**Output:** `<TAG>/05_vcd/<CID>/<CID>.out` (Q-Chem frequency output with IR intensities and VCD rotatory strengths)

On HPC, all conformers are submitted as a single SLURM job array (`<TAG>_array.slurm`) throttled to `--max-running` concurrent tasks.

### qc_ecd_uvvis_tools.py -- Electronic Spectral Broadening & Plotting

```
qc_ecd_uvvis_tools.py [OPTIONS] LOGS [LOGS ...]
```

`LOGS` can be individual `.out` files, glob patterns, or directories (searched recursively).

| Flag | Description | Default |
|------|-------------|---------|
| `--bw PATH` | Boltzmann weight file (`_energies.dat` from Stage 5) | _equal weights_ |
| `--outdir TAG` | Molecule directory; outputs go to `<TAG>/06_spectra/<TAG>_*` | |
| `--prefix STR` | Explicit output prefix (overrides `--outdir`) | `spectra` |
| `--uv_fwhm EV` | Gaussian FWHM for UV-Vis (eV) | `0.35` |
| `--ecd_fwhm EV` | Gaussian FWHM for ECD (eV) | `0.25` |
| `--xlim MIN MAX` | Wavelength range (nm) | _auto_ |
| `--uv_ylim YMIN YMAX` | UV-Vis y-axis limits | _auto_ |
| `--ecd_ylim YMIN YMAX` | ECD y-axis limits | _auto_ |
| `--stick` | Overlay stick spectrum | _off_ |
| `--flip_x` | Plot long-to-short wavelength | _off_ |
| `--scale FAC` | Multiply intensities after weighting | `1.0` |
| `--no_title` | Suppress figure titles | _off_ |
| `--title STR` | Custom figure title (overrides auto-generated title from prefix) | _auto_ |

**Output:** `<TAG>/06_spectra/<TAG>_uvvis.png/.pdf/.csv`, `<TAG>/06_spectra/<TAG>_ecd.png/.pdf/.csv`

### qc_vcd_ir_tools.py -- Vibrational Spectral Broadening & Plotting

```
qc_vcd_ir_tools.py [OPTIONS] LOGS [LOGS ...]
```

`LOGS` can be individual `.out` files, glob patterns, or directories (searched recursively).

| Flag | Description | Default |
|------|-------------|---------|
| `--bw PATH` | Boltzmann weight file (`_energies.dat` from Stage 5) | _equal weights_ |
| `--outdir TAG` | Molecule directory; outputs go to `<TAG>/06_spectra/<TAG>_*` | |
| `--prefix STR` | Explicit output prefix (overrides `--outdir`) | `vib` |
| `--ir_fwhm CM` | Gaussian FWHM for IR (cm-1) | `10` |
| `--vcd_fwhm CM` | Gaussian FWHM for VCD (cm-1) | `6` |
| `--xlim MIN MAX` | Wavenumber range (cm-1) | _auto_ |
| `--ir_ylim YMIN YMAX` | IR y-axis limits | _auto_ |
| `--vcd_ylim YMIN YMAX` | VCD y-axis limits | _auto_ |
| `--stick` | Overlay stick spectrum | _off_ |
| `--invert_ir` | Plot IR absorption peaks downward | _off_ |
| `--no_title` | Suppress figure titles | _off_ |
| `--title STR` | Custom figure title base (overrides auto-generated title from prefix) | _auto_ |

**Output:** `<TAG>/06_spectra/<TAG>_ir.png/.pdf/.csv`, `<TAG>/06_spectra/<TAG>_vcd.png/.pdf/.csv`

**Examples:**

```bash
# Boltzmann-weighted IR/VCD for methyloxirane
python3 qc_vcd_ir_tools.py methyloxirane/05_vcd \
    --bw methyloxirane/04_boltzmann/methyloxirane_energies.dat \
    --outdir methyloxirane --stick --xlim 800 3200

# Compare computed VCD to experiment with adjusted broadening
python3 qc_vcd_ir_tools.py ephedrine/05_vcd \
    --bw ephedrine/04_boltzmann/ephedrine_energies.dat \
    --outdir ephedrine --vcd_fwhm 8 --ir_fwhm 12
```

---

## Supported Q-Chem Options

### Density Functionals

The `--method` flag accepts any Q-Chem-recognized functional keyword. Common choices for spectroscopy:

| Functional | Type | HF Exchange | Typical Use |
|------------|------|-------------|-------------|
| B3LYP | Hybrid GGA | 20% | Geometry optimization, IR/VCD frequencies |
| PBE0 | Hybrid GGA | 25% | Slightly better than B3LYP for many properties |
| CAM-B3LYP | Range-separated | 19-65% | Charge-transfer excitations, ECD |
| wB97X-D3 | Range-separated | 22-100% | TD-DFT benchmark standard |
| M06-2X | Hybrid meta-GGA | 54% | Main-group thermochemistry |

**Note on VCD:** For vibrational spectra, B3LYP is the most extensively benchmarked functional and is generally recommended as a starting point. Range-separated hybrids (CAM-B3LYP, wB97X-D3) offer no systematic advantage for harmonic frequencies and are significantly more expensive for analytic Hessian calculations.

### Basis Sets

The `--basis` flag accepts any Q-Chem basis set keyword. Q-Chem does not require a separate auxiliary basis for RI-J -- density fitting is handled internally when applicable.

| Basis | Quality | Typical Use |
|-------|---------|-------------|
| `def2-SVP` | Double-zeta | Gas-phase geometry optimization |
| `6-31+G(d)` | Double-zeta + diffuse | Solvent-phase optimization, VCD frequencies |
| `def2-TZVP` | Triple-zeta | TD-DFT production runs |
| `6-311+G(d,p)` | Triple-zeta + diffuse | Basis set convergence checks |

### Solvents

The `--solvent` flag accepts Q-Chem solvent keywords for the SMD implicit solvation model. Solvent parameters are passed via the `$smx` input block. Common options:

`water`, `methanol`, `ethanol`, `acetone`, `acetonitrile`, `dichloromethane`, `trichloromethane`, `benzene`, `toluene`, `hexane`, `cyclohexane`, `dimethylsulfoxide`, `dimethylformamide`, `tetrahydrofuran`, `pyridine`, `diethylether`

**Note:** Q-Chem solvent keywords use full names with no spaces or hyphens (e.g., `1hexanol` not `1-hexanol`, `propanoicacid` not `propanoic acid`). See the SMx table in the Q-Chem manual for the complete list of 179 built-in solvents.

### Dispersion Corrections

The `--disp` flag controls how dispersion is applied. In Q-Chem, dispersion is set via `DFT_D` and related keywords in the `$rem` block:

| Value | Behavior |
|-------|-----------|
| `auto` _(default for Stage 1)_ | Detects whether the functional already includes dispersion (e.g., wB97X-**D3**). If not, adds D3(BJ). |
| `none` _(default for Stage 4)_ | No dispersion correction. |
| `D3BJ` | Grimme's D3 with Becke-Johnson damping (`DFT_D D3_BJ` in `$rem`). |

---

## Q-Chem Input Format

Q-Chem inputs use `$rem`/`$molecule`/`$end` block syntax rather than the keyword-line format used by some other programs. A typical input generated by the scripts looks like:

```
$rem
   JOBTYPE        opt
   METHOD         b3lyp
   BASIS          def2-svp
   DFT_D          D3_BJ
   XC_GRID        3
   MAX_SCF_CYCLES 150
   MEM_TOTAL      8192
$end

$molecule
   0 1
   C    0.000000    0.000000    0.000000
   H    0.000000    0.000000    1.089000
   ...
$end
```

For implicit solvation, an additional `$smx` block is written:

```
$smx
   solvent  water
$end
```

TD-DFT calculations include `RPA TRUE` and `CIS_N_ROOTS` in the `$rem` block. The scripts handle all input generation automatically -- you only need to supply CLI flags.

---

## Charge and Multiplicity

All scripts parse charge and multiplicity from line 2 of the XYZ file. Two formats are supported:

```
14                           # atom count
charge=0 mult=1              # key=value (any order)
C  0.000  0.000  0.000
...
```

```
14                           # atom count
0 1                          # bare integers: charge multiplicity
C  0.000  0.000  0.000
...
```

If neither format is detected, the defaults `charge=0 mult=1` are used.

---

## Molecule List Files

Instead of passing molecule tags as positional arguments, you can supply a text file with `--list`:

```
# molecules.txt
# Lines starting with # are ignored
pna
ephedrine
aspirin
```

```bash
./1-qchem-init-opt.sh --local --cpus 4 --list molecules.txt
```

---

## Example Molecules

The `pre_xyz/` directory includes starting geometries for several test molecules:

**Session 1 -- Rigid chromophores (UV-Vis benchmarking):**
- `pna.xyz` -- para-nitroaniline (intramolecular charge transfer)
- `dmabn.xyz` -- 4-(dimethylamino)benzonitrile (dual fluorescence)

**Session 2 -- Flexible chiral molecules (ECD + conformational averaging):**
- `ephedrine.xyz` -- (1R,2S)-ephedrine (amino alcohol, 4-10 conformers)
- `norephedrine.xyz` -- (1S,2R)-norephedrine (amino alcohol)
- `methyloxirane.xyz` -- (R)-methyloxirane (propylene oxide, small chiral reference)
- `phenylethan1ol.xyz` -- (S)-1-phenylethan-1-ol (alpha-methylbenzyl alcohol)
- `sparteine.xyz` -- (-)-sparteine (tetracyclic alkaloid, large conformer space)

**Session 3 -- Vibrational chiroptical spectroscopy (VCD):**
- `methyloxirane.xyz` -- (R)-methyloxirane (the canonical VCD benchmark, rigid)
- `phenylethan1ol.xyz` -- (S)-1-phenylethan-1-ol (moderate flexibility, VCD literature reference)

---

## Lab Sessions

This workflow supports three independent teaching sessions. Each lab is self-contained; instructors can offer them in any order or combination.

| Lab | Duration |
|-----|----------|
| UV-Vis Benchmarking (PNA) | 3 hours |
| ECD + Conformational Averaging (Ephedrine) | 3 hours |
| VCD/IR Spectroscopy (Methyloxirane) | 3 hours |

---

## Troubleshooting

**Q-Chem not found or wrong version used:** The recommended setup is to configure `cluster.cfg` (see Quick Start step 5). This is the most reliable approach because it is explicit, per-repo, and never conflicts with system `PATH` or shell rc files. Alternatives in order of priority:
```bash
# Option 1 -- cluster.cfg (recommended, set once per clone)
echo 'QCHEM_SETUP=/path/to/qchem/setqc' >> cluster.cfg

# Option 2 -- per-invocation flag
./1-qchem-init-opt.sh --qchem-setup /path/to/qchem/setqc pna
```
The Q-Chem setup script (typically called `setqc` or `qcenv.sh`) sets `QC`, `QCAUX`, and updates `PATH`. The workflow scripts source this file before launching Q-Chem. If the wrong version is being picked up, check that `QCHEM_SETUP` in `cluster.cfg` points to the correct installation.

**Library errors on compute nodes (libiomp5.so, libstdc++.so):** Intel-compiled Q-Chem installations require `libiomp5.so` (Intel OpenMP runtime) and a sufficiently new `libstdc++.so.6`. Set `CLUSTER_LD_LIBRARY_PATH` in `cluster.cfg` to include the Intel compiler lib directory and a compatible GCC lib64 directory. See `cluster.cfg.example` for detailed guidance.

**Convergence failure:** Increase `--max-iter` (Stage 1) or `--max-scf` (Stages 4/6) or try a different initial geometry. Check the `.out` file for SCF convergence messages.

**Empty Boltzmann output:** All conformers fell below `--p-cut`. Lower the threshold (e.g., `--p-cut 0.001`) or check that Stage 4 completed successfully for all conformers.

**Plotting tool import errors:** Install Python dependencies: `pip install -r requirements.txt`

**CREST fails with a shared library error on HPC** (e.g., `version 'GLIBCXX_3.4.XX' not found`): The CREST binary was compiled against a newer libgcc/libgomp than the cluster's system provides. Create a minimal conda environment to supply the missing libraries, then activate it before submitting the Stage 2b job:

```bash
conda create -n xtb_rt -c conda-forge libgcc libgomp libgfortran
conda activate xtb_rt
bash 2b-qchem-crest-conf-search.sh TAG
```

The SLURM job inherits the activated environment's `LD_LIBRARY_PATH`, so no changes to the batch script are needed.

**VCD frequency job runs slowly:** Analytic Hessian calculations scale more steeply than single-point energies. Use `--dry-run` first to check conformer count. B3LYP is recommended over range-separated hybrids for frequency calculations.

**Q-Chem output files:** Q-Chem writes `.out` files (not `.log`). The scripts and Python tools expect this extension. If you have outputs with a different extension, rename them or create symlinks.

---

## Citation

If you use this workflow in published research or teaching, please cite:

> [Author], [Author]. A Computational Spectroscopy Module for Teaching DFT-Based Spectral Workflows in Graduate Chemistry. _J. Chem. Educ._ **2026**, _XX_, XXXX-XXXX.

---

## License

[To be determined -- recommend MIT or BSD-3-Clause for JCE educational software.]
