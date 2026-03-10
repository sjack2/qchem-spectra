#!/usr/bin/env bash
# ============================================================================
# 6-qchem-tddft.sh — TD-DFT UV-Vis and ECD for Boltzmann-kept conformers
# ============================================================================
#
# OVERVIEW
#   For each conformer retained by Stage 5, this script:
#     1. Extracts the optimised geometry from the Stage 4 output,
#     2. Writes a TD-DFT single-point Q-Chem input with implicit solvent (SMD),
#     3. Submits a SLURM job.
#
#   Each output file contains both UV-Vis (oscillator strengths) and ECD
#   (rotatory strengths, length and velocity gauge) — Q-Chem prints both
#   automatically for every TD-DFT calculation.
#
#   For high-accuracy ECD requiring EOM-CCSD, use 6-qchem-ecd.sh instead
#   (note: significantly more expensive).
#
# Usage:
#   6-qchem-tddft.sh TAG
#   6-qchem-tddft.sh --list mols.txt --method CAM-B3LYP --roots 40 --dry-run
#
# Flags:
#   -m | --method NAME       TD-DFT functional            [CAM-B3LYP]
#   -b | --basis NAME        Basis set                    [def2-SVP]
#        --roots N           Number of excited states     [30]
#        --solvent NAME      SMD solvent                  [water]
#        --max-scf N         Max SCF cycles               [150]
#   -c | --cpus N            CPU cores                    [12]
#        --mem-per-cpu MB    Memory per core (MB)         [2048]
#        --partition NAME    SLURM partition              [general]
#        --time HH:MM:SS     Wall-clock limit             [04:00:00]
#        --list FILE         File of molecule TAGs
#        --dry-run           Write inputs but do not submit
#   -h | --help              Show this help and exit
#
# Cluster configuration (cluster.cfg):
#   See cluster.cfg.example for supported variables.
#
# Directory layout:
#   <TAG>/
#   ├── 03_solvent_opt/<CID>/<CID>.out    ← Stage 4 geometry source
#   ├── 04_boltzmann/<TAG>_bw_labels.dat  ← Stage 5 conformer list
#   └── 05_tddft/
#       └── <CID>/
#           ├── <CID>.inp
#           ├── <CID>.out    (after job completes)
#           └── <CID>.slurm
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_METHOD="CAM-B3LYP"
DEFAULT_BASIS="def2-SVP"
DEFAULT_ROOTS=30
DEFAULT_SOLVENT="water"
DEFAULT_MAX_SCF=150
DEFAULT_CPUS=12
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="general"
DEFAULT_WALL="04:00:00"

XYZ_DIR="pre_xyz"
SOLV_OUT_DIR="03_solvent_opt"
BW_SUBDIR="04_boltzmann"
OUT_SUBDIR="05_tddft"

# ============================================================================
# CLUSTER CONFIG
# ============================================================================
source_cluster_cfg() {
    local script_dir cfg
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cfg="${script_dir}/cluster.cfg"
    if [[ -f $cfg ]]; then
        log "Loaded cluster config: ${cfg}"
        # shellcheck source=/dev/null
        source "$cfg"
    fi
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
die()          { printf 'Error: %s\n' "$*" >&2; exit 1; }
log()          { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
warn()         { printf '[%s] Warning: %s\n' "$(date '+%F %T')" "$*" >&2; }
require_file() { [[ -f $1 ]] || die "File '$1' not found"; }

show_help() {
    sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# \{0,1\}//'
    exit 0
}

# ============================================================================
# Q-CHEM OUTPUT → XYZ EXTRACTOR (last geometry in file)
# ============================================================================
qcout2xyz() {
    python3 - "$1" <<'PY'
import pathlib, re, sys
lines = pathlib.Path(sys.argv[1]).read_text(errors="ignore").splitlines()
hdr = re.compile(
    r'(Standard Nuclear Orientation|Input orientation|Coordinates \(Angstroms\))', re.I
)
start = None
for i in range(len(lines) - 1, -1, -1):
    if hdr.search(lines[i]):
        start = i + 3; break
if start is None:
    sys.exit(f"qcout2xyz: no coordinate block found in {sys.argv[1]}")
geom = []
for l in lines[start:]:
    if not l.strip() or l.lstrip().startswith('-'): break
    p = l.split()
    if len(p) >= 5:   sym, x, y, z = p[1], p[2], p[3], p[4]
    elif len(p) == 4: sym, x, y, z = p
    else: continue
    geom.append(f"{sym} {x} {y} {z}")
if not geom:
    sys.exit(f"qcout2xyz: coordinate parsing failed in {sys.argv[1]}")
print(len(geom)); print(); print("\n".join(geom))
PY
}

# ============================================================================
# CHARGE / MULTIPLICITY FROM pre_xyz
# ============================================================================
read_charge_mult() {
    local tag=$1
    local pxyz="${XYZ_DIR}/${tag}.xyz"
    charge=0; mult=1
    [[ -f $pxyz ]] || { warn "pre_xyz/${tag}.xyz not found — using charge=0 mult=1"; return; }
    local h
    h=$(sed -n '2p' "$pxyz")
    [[ $h =~ charge[[:space:]]*=[[:space:]]*([+-]?[0-9]+) ]] && charge=${BASH_REMATCH[1]}
    if [[ $h =~ mult[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
        mult=${BASH_REMATCH[1]}
    elif [[ $h =~ ^[[:space:]]*([+-]?[0-9]+)[[:space:]]+([0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}; mult=${BASH_REMATCH[2]}
    fi
}

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    method=$DEFAULT_METHOD
    basis=$DEFAULT_BASIS
    roots=$DEFAULT_ROOTS
    solvent=$DEFAULT_SOLVENT
    max_scf=$DEFAULT_MAX_SCF
    cpus=$DEFAULT_CPUS
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=${CLUSTER_PARTITION:-$DEFAULT_PARTITION}
    wall=${CLUSTER_WALL:-$DEFAULT_WALL}
    qchem_setup=${QCHEM_SETUP:-/shares/chem_hlw/qchem/7.0.0-qmmm-nt-intel-pre/setqc}
    dry_run=false
    list_file=""
    single=""

    local opts
    opts=$(getopt -o hb:m:c: \
        --long help,method:,basis:,roots:,solvent:,max-scf:,cpus:,mem-per-cpu:,partition:,time:,list:,dry-run \
        -- "$@") || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case "$1" in
            -m|--method)    method=$2;      shift 2 ;;
            -b|--basis)     basis=$2;       shift 2 ;;
            --roots)        roots=$2;       shift 2 ;;
            --solvent)      solvent=$2;     shift 2 ;;
            --max-scf)      max_scf=$2;     shift 2 ;;
            -c|--cpus)      cpus=$2;        shift 2 ;;
            --mem-per-cpu)  mem_mb=$2;      shift 2 ;;
            --partition)    partition=$2;   shift 2 ;;
            --time)         wall=$2;        shift 2 ;;
            --list)         list_file=$2;   shift 2 ;;
            --dry-run)      dry_run=true;   shift ;;
            -h|--help)      show_help ;;
            --)             shift; break ;;
            *)              die "Unknown option '$1'" ;;
        esac
    done

    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional arguments not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide a molecule TAG"
        single=$1
    fi

    command -v sbatch >/dev/null 2>&1 || die "sbatch not in PATH"
}

# ============================================================================
# Q-CHEM INPUT WRITER
# ============================================================================
write_input() {
    local tag=$1 cid=$2 xyz=$3 inp=$4

    cat >"$inp" <<EOF
\$comment
TD-DFT UV-Vis + ECD (SMD ${solvent}) — generated by 6-qchem-tddft.sh
\$end

\$molecule
${charge} ${mult}
$(tail -n +3 "$xyz")
\$end

\$rem
  JOB_TYPE            SP
  METHOD              ${method}
  BASIS               ${basis}
  CIS_N_ROOTS         ${roots}
  CIS_SINGLETS        TRUE
  CIS_TRIPLETS        FALSE
  RPA                 TRUE
  SOLVENT_METHOD      SMD
  SCF_CONVERGENCE     8
  MAX_SCF_CYCLES      ${max_scf}
  SYM_IGNORE          TRUE
  XC_GRID             3
\$end

\$smx
  solvent ${solvent}
\$end
EOF
}

# ============================================================================
# SLURM SCRIPT WRITER
# ============================================================================
write_slurm() {
    local slurm_file=$1 cid=$2 workdir=$3
    local abs_workdir
    abs_workdir=$(cd "$workdir" && pwd)

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${cid}
#SBATCH --partition=${partition}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem-per-cpu=${mem_mb}
#SBATCH --time=${wall}
#SBATCH --chdir=${abs_workdir}
#SBATCH --output=${abs_workdir}/slurm-%j.out
#SBATCH --error=${abs_workdir}/slurm-%j.err

source ${qchem_setup}
${CLUSTER_LD_LIBRARY_PATH:+export LD_LIBRARY_PATH=${CLUSTER_LD_LIBRARY_PATH}:\$LD_LIBRARY_PATH}
export QCSCRATCH=/tmp/\$SLURM_JOB_ID

qchem -nt ${cpus} ${abs_workdir}/${cid}.inp ${abs_workdir}/${cid}.out
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 6: TD-DFT UV-Vis + ECD
-------------------------------------------------------------
 Q-Chem setup : ${qchem_setup}
 Method       : ${method}
 Basis        : ${basis}
 Roots        : ${roots}
 Solvent      : ${solvent} (SMD)
 Max SCF      : ${max_scf}
 Cores        : ${cpus}
 Mem/core     : ${mem_mb} MB
 Partition    : ${partition}
 Wall time    : ${wall}
 Dry run      : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# SUBMIT ONE CONFORMER
# ============================================================================
submit_one() {
    local tag=$1 cid=$2 outfile=$3
    local workdir="${tag}/${OUT_SUBDIR}/${cid}"
    mkdir -p "$workdir"

    local xyz_tmp
    xyz_tmp=$(mktemp --suffix=.xyz)
    qcout2xyz "$outfile" >"$xyz_tmp" || {
        warn "[${tag}] cannot extract geometry from ${cid}.out"
        rm -f "$xyz_tmp"; return
    }

    local inp="${workdir}/${cid}.inp"
    local slurm_file="${workdir}/${cid}.slurm"

    read_charge_mult "$tag"
    write_input "$tag" "$cid" "$xyz_tmp" "$inp"
    write_slurm "$slurm_file" "$cid" "$workdir"
    rm -f "$xyz_tmp"

    if $dry_run; then
        log "[${tag}] dry run — ${cid}.inp written"
    else
        log "[${tag}] submitting ${cid} (partition=${partition})"
        sbatch "$slurm_file"
    fi
}

# ============================================================================
# PROCESS ONE TAG
# ============================================================================
process_tag() {
    local tag=$1
    local lab="${tag}/${BW_SUBDIR}/${tag}_bw_labels.dat"
    if [[ ! -f $lab ]]; then
        warn "[${tag}] label file not found: ${lab}"
        return
    fi

    mapfile -t cids < <(grep -v '^[[:space:]]*$' "$lab")
    if [[ ${#cids[@]} -eq 0 ]]; then
        warn "[${tag}] label file is empty"
        return
    fi

    log "[${tag}] found ${#cids[@]} conformers in label file"
    for cid in "${cids[@]}"; do
        local outpath="${tag}/${SOLV_OUT_DIR}/${cid}/${cid}.out"
        if [[ ! -f $outpath ]]; then
            warn "[${tag}] missing ${cid}.out — skipping"
            continue
        fi
        submit_one "$tag" "$cid" "$outpath"
    done
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"
    print_banner

    if [[ -n $list_file ]]; then
        while IFS= read -r entry || [[ -n $entry ]]; do
            [[ -z $entry || $entry == \#* || $entry == \;* ]] && continue
            process_tag "$entry"
        done < "$list_file"
    else
        process_tag "$single"
    fi

    log "Stage 6-tddft complete."
}

main "$@"
