#!/usr/bin/env bash
# ============================================================================
# 4-qchem-solvent-opt.sh — Implicit-solvent geometry optimisation (Stage 4)
# ============================================================================
#
# OVERVIEW
#   For every conformer produced by Stage 3, this script:
#     1. Reads charge/multiplicity from the original XYZ in pre_xyz/,
#     2. Writes a Q-Chem input with implicit solvent (SMD),
#     3. Submits a SLURM job.
#
# Usage:
#   4-qchem-solvent-opt.sh TAG
#   4-qchem-solvent-opt.sh --list mols.txt --solvent acetonitrile --dry-run
#
# Flags:
#   -m | --method NAME          DFT functional                  [B3LYP]
#   -b | --basis NAME           Basis set                       [6-31+G(d)]
#        --disp KW              Dispersion: auto|none|D3BJ       [none]
#        --solvent NAME         SMD solvent keyword              [water]
#        --max-scf N            Max SCF cycles                  [150]
#   -c | --cpus N               CPU cores                       [12]
#        --mem-per-cpu MB       Memory per core (MB)            [2048]
#        --partition NAME       SLURM partition                  [general]
#        --time HH:MM:SS        Wall-clock limit                 [02:00:00]
#        --list FILE            File of molecule TAGs
#        --dry-run              Write inputs but do not submit
#   -h | --help                 Show this help and exit
#
# Cluster configuration (cluster.cfg):
#   See cluster.cfg.example for supported variables.
#
# Directory layout:
#   <TAG>/
#   ├── 02_conf_search/split_xyz/       ← Stage 3 output (input here)
#   ├── 03_solvent_opt/
#   │   ├── <TAG>_1/
#   │   │   ├── <TAG>_1.inp
#   │   │   ├── <TAG>_1.out             (after job completes)
#   │   │   └── <TAG>_1.slurm
#   │   └── <TAG>_2/ ...
#   └── pre_xyz/<TAG>.xyz               charge/multiplicity source
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_METHOD="B3LYP"
DEFAULT_BASIS="6-31+G(d)"
DEFAULT_DISP="none"
DEFAULT_SOLVENT="water"
DEFAULT_MAX_SCF=150
DEFAULT_CPUS=12
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="general"
DEFAULT_WALL="02:00:00"

XYZ_DIR="pre_xyz"
CONF_SUBDIR="02_conf_search/split_xyz"
SOLV_OPT_SUBDIR="03_solvent_opt"

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
# CHARGE / MULTIPLICITY FROM pre_xyz
# ============================================================================
read_charge_mult() {
    local tag=$1
    local pxyz="${XYZ_DIR}/${tag}.xyz"
    charge=0; mult=1
    [[ -f $pxyz ]] || { warn "pre_xyz/${tag}.xyz not found — using charge=0 mult=1"; return; }

    local header_line
    header_line=$(sed -n '2p' "$pxyz")

    if [[ $header_line =~ charge[[:space:]]*=[[:space:]]*([+-]?[0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
    fi
    if [[ $header_line =~ mult[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
        mult=${BASH_REMATCH[1]}
    elif [[ $header_line =~ ^[[:space:]]*([+-]?[0-9]+)[[:space:]]+([0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
        mult=${BASH_REMATCH[2]}
    fi
}

# ============================================================================
# DISPERSION LOGIC
# ============================================================================
disp_line() {
    local f=${method^^}
    case $disp_mode in
        none|NONE) printf '' ;;
        D3BJ|d3bj) printf '  DFT_D              D3_BJ' ;;
        auto|AUTO)
            [[ $method =~ (-D|D3|D4|-V)$ ]] && { printf ''; return; }
            case $f in
                WB97X-D|WB97X-D3|WB97X-D4|WB97X-V|WB97XD|WB97M-V|B97-D|B97-D3)
                    printf ''; return ;;
            esac
            printf '  DFT_D              D3_BJ'
            ;;
        *) die "--disp must be auto, none, or D3BJ (got '$disp_mode')" ;;
    esac
}

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    method=$DEFAULT_METHOD
    basis=$DEFAULT_BASIS
    disp_mode=$DEFAULT_DISP
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
        --long help,basis:,method:,disp:,solvent:,max-scf:,cpus:,mem-per-cpu:,partition:,time:,list:,dry-run \
        -- "$@") || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)    method=$2;      shift 2 ;;
            -b|--basis)     basis=$2;       shift 2 ;;
            --disp)         disp_mode=$2;   shift 2 ;;
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
write_qchem_input() {
    local tag=$1 cid=$2 xyz=$3 inp=$4
    local disp
    disp=$(disp_line)

    cat >"$inp" <<EOF
\$comment
Implicit-solvent optimisation of ${cid} (SMD ${solvent})
\$end

\$molecule
${charge} ${mult}
$(tail -n +3 "$xyz")
\$end

\$rem
  JOB_TYPE            OPT
  METHOD              ${method}
  BASIS               ${basis}
  SOLVENT_METHOD      SMD
  MAX_SCF_CYCLES      ${max_scf}
  SCF_CONVERGENCE     8
  SYM_IGNORE          TRUE
  THRESH              11
  XC_GRID             3
${disp}
\$end

\$smx
  solvent             ${solvent}
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
 Stage 4: Implicit-Solvent Geometry Optimisation
-------------------------------------------------------------
 Q-Chem setup : ${qchem_setup}
 Method       : ${method}
 Basis        : ${basis}
 Dispersion   : ${disp_mode}
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
submit_conf() {
    local tag=$1 xyz=$2
    local cid
    cid=$(basename "$xyz" .xyz)
    local workdir="${tag}/${SOLV_OPT_SUBDIR}/${cid}"
    mkdir -p "$workdir"

    local inp="${workdir}/${cid}.inp"
    local slurm_file="${workdir}/${cid}.slurm"

    read_charge_mult "$tag"
    write_qchem_input "$tag" "$cid" "$xyz" "$inp"
    write_slurm "$slurm_file" "$cid" "$workdir"

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
    local split_xyz_dir="${tag}/${CONF_SUBDIR}"
    if [[ ! -d $split_xyz_dir ]]; then
        warn "[${tag}] split_xyz directory not found: ${split_xyz_dir}"
        return
    fi

    shopt -s nullglob
    local xyzs=("${split_xyz_dir}/${tag}_"*.xyz)
    shopt -u nullglob

    if [[ ${#xyzs[@]} -eq 0 ]]; then
        warn "[${tag}] no conformer XYZ files found in ${split_xyz_dir}"
        return
    fi

    log "[${tag}] found ${#xyzs[@]} conformers"
    for xyz in "${xyzs[@]}"; do
        submit_conf "$tag" "$xyz"
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

    log "Stage 4 complete."
}

main "$@"
