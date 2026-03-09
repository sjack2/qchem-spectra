#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 1-qchem-init-opt.sh – Generate (and optionally submit) Q-Chem 6.2 gas-phase
# geometry optimisations.
# -----------------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

# --------------------------- defaults ----------------------------------------
DEFAULT_METHOD="B3LYP"
DEFAULT_BASIS="def2-SVP"
DEFAULT_MAX_SCF=150
DEFAULT_PARTITION="circe"
DEFAULT_CPUS=12
DEFAULT_GRID="SG3"          # SG1–3  → XC_GRID 1–3
DEFAULT_WALL="02:00:00"
DEFAULT_MEM_PER_CPU=2048    # MB
XYZ_DIR="pre_xyz"

# ----------------------------- help ------------------------------------------
show_help() {
cat <<EOF
Usage:
  $(basename "$0") [OPTIONS] tag-or-xyz
  $(basename "$0") --list file.txt [OPTIONS]

Options:
  -m, --method METHOD       Electronic-structure method      [$DEFAULT_METHOD]
  -b, --basis  BASIS        Basis set                         [$DEFAULT_BASIS]
  --disp {auto,none,D3BJ}   Extra dispersion keyword          [auto]
  --grid {SG1,SG2,SG3}      DFT grid (maps to XC_GRID 1–3)    [$DEFAULT_GRID]
  --max-scf N               Max SCF cycles                    [$DEFAULT_MAX_SCF]
  -c, --cpus N              CPUs for SLURM / qchem            [$DEFAULT_CPUS]
  --mem-per-cpu MB          Memory per CPU (MB)               [$DEFAULT_MEM_PER_CPU]
  --partition NAME          SLURM partition                   [$DEFAULT_PARTITION]
  --list FILE               File of tags or XYZ paths
  --dry-run                 Emit files but do not submit
  -h, --help                Show this help and exit
EOF
}

usage_err(){ echo "Error: $1" >&2; show_help; exit 1; }

# ----------------------------- CLI -------------------------------------------
parse_cli() {
  METHOD=$DEFAULT_METHOD
  BASIS=$DEFAULT_BASIS
  DISP_MODE="auto"
  GRID=$DEFAULT_GRID
  MAX_SCF=$DEFAULT_MAX_SCF
  PARTITION=$DEFAULT_PARTITION
  CPUS=$DEFAULT_CPUS
  MEM_MB=$DEFAULT_MEM_PER_CPU
  DRY_RUN=false
  LIST_FILE=""

  local opts
  opts=$(getopt -o hb:m:c: --long help,basis:,method:,disp:,grid:,max-scf:,partition:,list:,cpus:,mem-per-cpu:,dry-run -- "$@")
  eval set -- "$opts"
  while true; do
    case "$1" in
      -m|--method)       METHOD=$2; shift 2 ;;
      -b|--basis)        BASIS=$2; shift 2 ;;
      --disp)            DISP_MODE=$2; shift 2 ;;
      --grid)            GRID=$2; shift 2 ;;
      --max-scf)         MAX_SCF=$2; shift 2 ;;
      --partition)       PARTITION=$2; shift 2 ;;
      -c|--cpus)         CPUS=$2; shift 2 ;;
      --mem-per-cpu)     MEM_MB=$2; shift 2 ;;
      --list)            LIST_FILE=$2; shift 2 ;;
      --dry-run)         DRY_RUN=true; shift ;;
      -h|--help)         show_help; exit 0 ;;
      --) shift; break ;;
      *) usage_err "Unknown option $1" ;;
    esac
  done

  if [[ -n $LIST_FILE ]]; then
    [[ $# -eq 0 ]] || usage_err "XYZ not allowed with --list"
    [[ -f $LIST_FILE ]] || usage_err "List file '$LIST_FILE' not found"
  else
    [[ $# -eq 1 ]] || usage_err "Exactly one tag or XYZ path must be provided"
    if [[ $1 == *.xyz ]]; then
      XYZ_ARG=$1
    else
      XYZ_ARG="${XYZ_DIR}/$1.xyz"
    fi
    [[ -f $XYZ_ARG ]] || usage_err "XYZ '$XYZ_ARG' not found"
  fi
}

# --------------------- charge / multiplicity parser --------------------------
read_xyz_header() {
  read -r _atom_count <"$XYZ_FILE"
  read -r header_line < <(tail -n +2 "$XYZ_FILE" | head -n 1)
  if [[ $header_line =~ ([+-]?[0-9]+)[[:space:]]+([0-9]+) ]]; then
    CHARGE=${BASH_REMATCH[1]}
    MULT=${BASH_REMATCH[2]}
  else
    CHARGE=0
    MULT=1
  fi
}

# --------------------------- dispersion logic --------------------------------
disp_line() {
  case $DISP_MODE in
    none)  echo "" ;;
    D3BJ)  echo "  DFT_D              D3_BJ" ;;
    auto)
      local f=${METHOD^^}            # upper-case
      [[ $f =~ (-D|D3|D4|-V)$ ]] && { echo ""; return; }
      echo "  DFT_D              D3_BJ" ;;
    *) usage_err "--disp must be auto, none, or D3BJ" ;;
  esac
}

grid_number(){ echo "${GRID#SG}"; }     # SG3 → 3, SG1 → 1 …

# --------------------------- Q-Chem input writer -----------------------------
write_qchem_input() {
  local DISP=$(disp_line)
  cat >"$INPUT_FILE" <<EOF
\$comment
Gas-phase optimisation of $MOL
\$end

\$molecule
$CHARGE $MULT
$(tail -n +3 "$XYZ_FILE")
\$end

\$rem
  JOB_TYPE            OPT
  METHOD              $METHOD
  BASIS               $BASIS
  XC_GRID             $(grid_number)
  XC_SMART_GRID       TRUE
  GEN_SCFMAN          TRUE
  SCF_ALGORITHM       DIIS_GDM
  MAX_SCF_CYCLES      $MAX_SCF
  SCF_CONVERGENCE     8
  THRESH              12
$DISP
\$end
EOF
}

# --------------------------- SLURM writer ------------------------------------
write_slurm() {
  cat >"$SLURM_FILE" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=opt_${MOL}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem-per-cpu=${MEM_MB}
#SBATCH --time=${DEFAULT_WALL}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module purge
source /shares/chem_hlw/qchem/7.0.0-qmmm-nt-intel-pre/setqc
export QCSCRATCH=/tmp/\$SLURM_JOBID

qchem -nt ${CPUS} ${MOL}.inp ${MOL}.out
EOF
  chmod +x "$SLURM_FILE"
}

# ------------------------ process one molecule -------------------------------
process_molecule() {
  XYZ_FILE=$1
  MOL=$2

  ROOTDIR=$MOL
  WORKDIR="${ROOTDIR}/${MOL}_qchem_opt"
  mkdir -p "$WORKDIR"

  INPUT_FILE="${WORKDIR}/${MOL}.inp"
  SLURM_FILE="${WORKDIR}/${MOL}.slurm"

  read_xyz_header
  write_qchem_input
  write_slurm

  if $DRY_RUN; then
    printf '\n--- %s (dry run) ---\n' "$MOL"
    cat "$INPUT_FILE"
    printf '\n--- SLURM ---\n'
    cat "$SLURM_FILE"
  else
    sbatch --chdir="$WORKDIR" "$(basename "$SLURM_FILE")"
  fi
}

# -------------------------------- main ---------------------------------------
main() {
  parse_cli "$@"

  if [[ -n $LIST_FILE ]]; then
    while IFS= read -r tag || [[ -n $tag ]]; do
      [[ -z $tag || $tag == \#* || $tag == \;* ]] && continue
      if [[ $tag == *.xyz ]]; then
        XYZ_PATH=$tag; MOL_TAG=$(basename "$tag" .xyz)
      else
        XYZ_PATH="${XYZ_DIR}/${tag}.xyz"; MOL_TAG=$tag
      fi
      [[ -f $XYZ_PATH ]] || { echo "Warning: '$XYZ_PATH' not found, skipping" >&2; continue; }
      process_molecule "$XYZ_PATH" "$MOL_TAG"
    done <"$LIST_FILE"
  else
    XYZ_PATH=$XYZ_ARG
    MOL_TAG=$(basename "$XYZ_ARG" .xyz)
    process_molecule "$XYZ_PATH" "$MOL_TAG"
  fi
}

main "$@"
