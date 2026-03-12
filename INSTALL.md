# Installation Guide

This guide covers installing all software required to run the Q-Chem Electronic Spectra Lab workflow. Three environments are supported:

- **Linux workstation** (native Ubuntu, Fedora, or similar)
- **Windows via WSL** (Windows Subsystem for Linux)
- **HPC cluster** (SLURM-managed, modules-based)

---

## Prerequisites at a Glance

| Software | Version | Required? | Purpose |
|----------|---------|-----------|---------|
| Q-Chem | >= 5.4 | Yes | Quantum chemistry calculations |
| Open Babel | >= 3.0 | Yes | XYZ/SDF conversion, Confab conformer search |
| Python 3 | >= 3.8 | Yes | Spectral plotting tools |
| NumPy | >= 1.20 | Yes | Numerical operations |
| Matplotlib | >= 3.4 | Yes | Plot generation |
| Pandas | >= 1.3 | Yes | Data handling |
| CREST | >= 2.12 | Optional | GFN2-xTB conformer search (Stage 2b) |
| Bash | >= 4.0 | Yes | Workflow scripts (ships with all Linux distributions) |
| getopt (GNU) | -- | Yes | CLI flag parsing (ships with `util-linux`) |

> **Note:** Q-Chem is commercial software and requires a paid license. Contact [Q-Chem, Inc.](https://www.q-chem.com/) for academic or commercial pricing. Unlike ORCA, Q-Chem is not free for academic use.

---

## 1. Check Existing Installations

Before installing anything, check what you already have:

```bash
which qchem        # Q-Chem
qchem -h 2>&1 | head -3   # Q-Chem help / version info
which obabel       # Open Babel
obabel -V          # Open Babel version
python3 --version  # Python 3
which crest        # CREST (optional)
```

If a command prints a valid path and version, that tool is already installed and you can skip its section below.

---

## 2. Python 3

Most Linux distributions ship Python 3. Verify:

```bash
python3 --version
```

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
```

### Fedora / CentOS / RHEL

```bash
sudo dnf install python3 python3-pip
```

### Python Dependencies

Install the plotting tools' requirements:

```bash
pip install numpy matplotlib pandas
```

> **Tip:** If you prefer isolated environments, use a virtual environment:
> ```bash
> python3 -m venv .venv
> source .venv/bin/activate
> pip install numpy matplotlib pandas
> ```

---

## 3. Open Babel

Open Babel provides the `obabel` command for file format conversion and the Confab conformer search engine used in Stage 2.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install openbabel
```

### Fedora / CentOS

```bash
sudo dnf install openbabel
```

### Verify

```bash
obabel -V
```

This should print the version (e.g., `Open Babel 3.1.1`). Test Confab:

```bash
obabel -L conformer
```

This should list the Confab conformer search plugin.

---

## 4. Q-Chem

Q-Chem is a commercial quantum chemistry package. You must obtain a license from [Q-Chem, Inc.](https://www.q-chem.com/) before proceeding.

### a) Obtain a License

1. Visit [q-chem.com](https://www.q-chem.com/) and contact sales for academic or commercial licensing.
2. Download the installer or receive installation media from Q-Chem after purchase.

### b) Install

Follow the instructions provided with your Q-Chem distribution. A typical installation places Q-Chem in a directory such as `/opt/qchem` or `/usr/local/qchem`. The installer creates a setup script called `setqc` inside the Q-Chem root directory.

### c) Configure Environment

Q-Chem uses a setup script (`setqc`) that configures all required environment variables (`QC`, `QCAUX`, `QCSCRATCH`, `PATH`, etc.). Source it in your shell:

```bash
source /opt/qchem/setqc
```

To make this permanent, add it to your `~/.bashrc` (or `~/.bash_profile`, `~/.zshrc`, etc.):

```bash
# Q-Chem
source /opt/qchem/setqc
```

Then reload:

```bash
source ~/.bashrc
```

Alternatively, the workflow scripts accept `--qchem-setup /path/to/setqc` or respect the `QCHEM_SETUP` environment variable:

```bash
export QCHEM_SETUP=/opt/qchem/setqc
```

### d) Verify

```bash
which qchem

# Run a quick test (single-point HF on water)
mkdir -p /tmp/qchem_test && cd /tmp/qchem_test
cat > test.in << 'EOF'
$molecule
0 1
O   0.0   0.0   0.0
H   0.0   0.757  0.587
H   0.0  -0.757  0.587
$end

$rem
jobtype      sp
method       hf
basis        sto-3g
$end
EOF
qchem test.in test.out
grep "Total energy" test.out
# Should print a total energy around -74.96 Hartree
cd - && rm -rf /tmp/qchem_test
```

### Parallelism

Q-Chem handles parallelism internally via shared-memory threads -- no external MPI installation is needed. The number of threads is controlled by the `-nt` flag:

```bash
qchem -nt 4 input.in output.out
```

The workflow scripts pass the CPU count to Q-Chem automatically based on the `--cpus` flag.

---

## 5. CREST (Optional)

CREST performs conformer searches using GFN2-xTB metadynamics. It is used in the alternative Stage 2b path and is more thorough than Confab for flexible molecules.

### Download

CREST binaries are available from the [Grimme group GitHub](https://github.com/crest-lab/crest/releases).

```bash
# Download the latest static binary
wget https://github.com/crest-lab/crest/releases/download/v3.0.2/crest-gnu-12-ubuntu-latest.tar.xz
tar -xf crest-gnu-12-ubuntu-latest.tar.xz
sudo mv crest /usr/local/bin/
```

Or install via conda:

```bash
conda install -c conda-forge crest
```

### Verify

```bash
crest --version
```

---

## 6. Windows Users: Installing WSL

If you are on Windows 10 (version 2004+) or Windows 11, you can run the entire workflow inside the Windows Subsystem for Linux.

### a) Check Windows Version

Press `Win+R`, type `winver`, and press Enter. Confirm you have Windows 10 version 2004 or later, or Windows 11.

### b) Install WSL

Open **PowerShell** or **Command Prompt** as Administrator and run:

```powershell
wsl --install
```

This enables WSL, installs the Virtual Machine Platform, and installs Ubuntu by default. Restart your computer when prompted.

If the automatic method does not work, enable the features manually:

```powershell
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```

Then restart, and install Ubuntu from the Microsoft Store.

### c) First Launch

Launch Ubuntu from the Start Menu. You will be prompted to create a Linux username and password.

### d) Verify

In PowerShell:

```powershell
wsl --list --verbose
```

This should show your Ubuntu installation with State = Running.

### e) Install Software Inside WSL

Once inside the WSL terminal, follow Sections 2--5 above as if you were on a native Linux machine. All `apt` and `pip` commands work identically.

> **Performance note:** WSL2 provides near-native performance for CPU-bound tasks like Q-Chem. Store your working files inside the Linux filesystem (`/home/username/`) rather than on the Windows drive (`/mnt/c/`) for best I/O performance.

---

## 7. HPC Cluster Setup

On a SLURM-managed cluster, software is typically loaded via the `module` system rather than installed locally.

### Typical Module Setup

```bash
module load apps/qchem/6.1
module load apps/python/3.11
```

The exact module names vary by institution. Check with `module avail qchem`.

If Q-Chem is not available as a module but is installed in a shared directory, source the setup script directly:

```bash
source /shares/chem/qchem/setqc
```

### Python Dependencies

On clusters, you may not have root access. Use `--user` or a virtual environment:

```bash
pip install --user numpy matplotlib pandas
# or
python3 -m venv ~/.venvs/spectra
source ~/.venvs/spectra/bin/activate
pip install numpy matplotlib pandas
```

### Cluster Configuration File

The workflow reads site-specific paths from `cluster.cfg`. Copy the example and edit it for your cluster:

```bash
cp cluster.cfg.example cluster.cfg
```

The key variable is `QCHEM_SETUP` -- the full path to your Q-Chem `setqc` script:

```bash
# cluster.cfg
QCHEM_SETUP=/path/to/qchem/setqc
CLUSTER_PARTITION=general
```

If Q-Chem was compiled with Intel compilers, you may also need `CLUSTER_LD_LIBRARY_PATH` to provide the Intel OpenMP runtime (`libiomp5.so`) and a compatible `libstdc++.so.6`:

```bash
CLUSTER_LD_LIBRARY_PATH=/path/to/intel/compiler/lib:/path/to/gcc-X.Y.Z/lib64
```

See `cluster.cfg.example` for detailed notes on when this is needed.

### SLURM Configuration

The workflow scripts accept SLURM options:

```bash
./1-qchem-init-opt.sh --partition main --time 02:00:00 --cpus 8 --mem-per-cpu 4096 aspirin
```

Q-Chem uses shared-memory threads, so the scripts request `--cpus-per-task` (not `--ntasks-per-node`) in the generated SLURM batch files. The default partition is `general`. Override with `--partition` to match your institution.

### CREST on HPC

CREST is often available as a module:

```bash
module load apps/crest/3.0
```

If not, download the static binary and place it in your `$HOME/bin`:

```bash
mkdir -p $HOME/bin
wget -O $HOME/bin/crest https://github.com/crest-lab/crest/releases/download/v3.0.2/crest
chmod +x $HOME/bin/crest
export PATH=$HOME/bin:$PATH
```

---

## 8. Final Verification

Run these commands to confirm everything is ready:

```bash
echo "=== Python ==="
python3 --version
python3 -c "import numpy, matplotlib, pandas; print('NumPy', numpy.__version__); print('Matplotlib', matplotlib.__version__); print('Pandas', pandas.__version__)"

echo "=== Open Babel ==="
obabel -V

echo "=== Q-Chem ==="
which qchem
echo "(run a test job to confirm -- see Section 4d)"

echo "=== CREST (optional) ==="
which crest 2>/dev/null && crest --version || echo "CREST not installed (optional)"

echo "=== Workflow scripts ==="
./1-qchem-init-opt.sh --help | head -3
```

If all commands succeed, you are ready to run the pipeline. See [README.md](README.md) for usage instructions.
