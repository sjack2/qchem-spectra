#!/usr/bin/env bash
# ============================================================================
# 1-qchem-init-opt.sh — Gas-phase geometry optimisation (Stage 1)
# ============================================================================
#
# OVERVIEW
#   For each supplied molecule this script:
#     1. Reads charge and multiplicity from the XYZ comment line,
#     2. Writes a Q-Chem geometry-optimisation input file,
#     3. Submits a SLURM job.
#
# Usage:
#   1-qchem-init-opt.sh TAG
#   1-qchem-init-opt.sh path/to/foo.xyz --method PBE0
#   1-qchem-init-opt.sh --list molecules.txt --dry-run
#
# Flags:
#   -m | --method METHOD          DFT functional               [B3LYP]
#   -b | --basis  BASIS           Basis set                    [def2-SVP]
#        --disp {auto,none,D3BJ}  Dispersion correction        [auto]
#        --grid {SG1,SG2,SG3}     XC grid (XC_GRID 1-3)        [SG3]
#        --max-scf N              Max SCF cycles               [150]
#   -c | --cpus N                 CPU cores                    [12]
#        --mem-per-cpu MB         Memory per core (MB)         [2048]
#        --partition NAME         SLURM partition              [general]
#        --time HH:MM:SS          Wall-clock limit             [02:00:00]
#        --list FILE              File of TAGs or XYZ paths
#        --dry-run                Write inputs but do not submit
#   -h | --help                   Show this help and exit
#
# Cluster configuration (cluster.cfg):
#   Create a file named cluster.cfg in the same directory as this script
#   to set site-specific defaults (QCHEM_SETUP, CLUSTER_PARTITION,
#   CLUSTER_WALL). See cluster.cfg.example for details.
#
# XYZ charge/multiplicity convention:
#   Line 2 may specify charge and multiplicity in either format:
#     charge=0 mult=1       (key=value tokens, any order)
#     0 1                   (bare integers: charge multiplicity)
#   If neither is found, defaults charge=0 mult=1 are used.
#
# Directory layout produced:
#   <TAG>/
#   └── 01_gas_opt/
#       ├── <TAG>.inp          Q-Chem input
#       ├── <TAG>.out          Q-Chem output  (after job completes)
#       └── <TAG>.slurm        SLURM script
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_METHOD="B3LYP"
DEFAULT_BASIS="def2-SVP"
DEFAULT_DISP="auto"
DEFAULT_GRID="SG3"
DEFAULT_MAX_SCF=150
DEFAULT_CPUS=12
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="general"
DEFAULT_WALL="02:00:00"

XYZ_DIR="pre_xyz"

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
# CHARGE / MULTIPLICITY PARSER
# ============================================================================
read_xyz_header() {
    local xyz_file=$1
    local header_line
    header_line=$(sed -n '2p' "$xyz_file")

    if [[ $header_line =~ charge[[:space:]]*=[[:space:]]*([+-]?[0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
    else
        charge=0
    fi

    if [[ $header_line =~ mult[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
        mult=${BASH_REMATCH[1]}
    elif [[ $header_line =~ ^[[:space:]]*([+-]?[0-9]+)[[:space:]]+([0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
        mult=${BASH_REMATCH[2]}
    else
        mult=1
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

grid_number() { echo "${grid#SG}"; }    # SG3 → 3

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    method=$DEFAULT_METHOD
    basis=$DEFAULT_BASIS
    disp_mode=$DEFAULT_DISP
    grid=$DEFAULT_GRID
    max_scf=$DEFAULT_MAX_SCF
    cpus=$DEFAULT_CPUS
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=${CLUSTER_PARTITION:-$DEFAULT_PARTITION}
    wall=${CLUSTER_WALL:-$DEFAULT_WALL}
    qchem_setup=${QCHEM_SETUP:-/shares/chem_hlw/qchem/7.0.0-qmmm-nt-intel-pre/setqc}
    dry_run=false
    list_file=""
    xyz_arg=""

    local opts
    opts=$(getopt -o hb:m:c: \
        --long help,basis:,method:,disp:,grid:,max-scf:,cpus:,mem-per-cpu:,partition:,time:,list:,dry-run \
        -- "$@") || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)    method=$2;      shift 2 ;;
            -b|--basis)     basis=$2;       shift 2 ;;
            --disp)         disp_mode=$2;   shift 2 ;;
            --grid)         grid=$2;        shift 2 ;;
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
        [[ $# -eq 1 ]] || die "Provide exactly one TAG or XYZ path (got $#)"
        if [[ $1 == *.xyz ]]; then
            xyz_arg=$1
        else
            xyz_arg="${XYZ_DIR}/$1.xyz"
        fi
        require_file "$xyz_arg"
    fi
}

# ============================================================================
# Q-CHEM INPUT WRITER
# ============================================================================
write_qchem_input() {
    local inp_file=$1 xyz_file=$2
    local disp
    disp=$(disp_line)

    cat >"$inp_file" <<EOF
\$comment
Gas-phase optimisation of $(basename "$xyz_file" .xyz)
\$end

\$molecule
${charge} ${mult}
$(tail -n +3 "$xyz_file")
\$end

\$rem
  JOB_TYPE            OPT
  METHOD              ${method}
  BASIS               ${basis}
  XC_GRID             $(grid_number)
  XC_SMART_GRID       TRUE
  GEN_SCFMAN          TRUE
  SCF_ALGORITHM       DIIS_GDM
  MAX_SCF_CYCLES      ${max_scf}
  SCF_CONVERGENCE     8
  THRESH              12
${disp}
\$end
EOF
}

# ============================================================================
# SLURM SCRIPT WRITER
# ============================================================================
write_slurm() {
    local slurm_file=$1 tag=$2 workdir=$3
    local abs_workdir
    abs_workdir=$(cd "$workdir" && pwd)

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=opt_${tag}
#SBATCH --partition=${partition}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem-per-cpu=${mem_mb}
#SBATCH --time=${wall}
#SBATCH --chdir=${abs_workdir}
#SBATCH --output=${abs_workdir}/slurm-%j.out
#SBATCH --error=${abs_workdir}/slurm-%j.err

module purge
source ${qchem_setup}
export QCSCRATCH=/tmp/\$SLURM_JOB_ID

qchem -nt ${cpus} ${abs_workdir}/${tag}.inp ${abs_workdir}/${tag}.out
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 1: Gas-Phase Geometry Optimisation
-------------------------------------------------------------
 Q-Chem setup : ${qchem_setup}
 Method       : ${method}
 Basis        : ${basis}
 Dispersion   : ${disp_mode}
 Grid         : ${grid}
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
# PROCESS ONE MOLECULE
# ============================================================================
process_molecule() {
    local xyz_file=$1 tag=$2

    local workdir="${tag}/01_gas_opt"
    mkdir -p "$workdir"

    local inp_file="${workdir}/${tag}.inp"
    local slurm_file="${workdir}/${tag}.slurm"

    read_xyz_header "$xyz_file"
    write_qchem_input "$inp_file" "$xyz_file"

    if $dry_run; then
        log "[${tag}] dry run — input written to ${inp_file}"
        echo "--- Q-Chem input: ${inp_file} ---"
        cat "$inp_file"
        write_slurm "$slurm_file" "$tag" "$workdir"
        echo ""
        echo "--- SLURM script: ${slurm_file} ---"
        cat "$slurm_file"
        echo ""
        return
    fi

    write_slurm "$slurm_file" "$tag" "$workdir"
    log "[${tag}] submitting to SLURM (partition=${partition})"
    sbatch "$slurm_file"
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
            local tag xyz_path
            if [[ $entry == *.xyz ]]; then
                xyz_path=$entry
                tag=$(basename "$entry" .xyz)
            else
                xyz_path="${XYZ_DIR}/${entry}.xyz"
                tag=$entry
            fi
            if [[ ! -f $xyz_path ]]; then
                warn "'${xyz_path}' not found — skipping"
                continue
            fi
            process_molecule "$xyz_path" "$tag"
        done < "$list_file"
    else
        local tag
        tag=$(basename "$xyz_arg" .xyz)
        process_molecule "$xyz_arg" "$tag"
    fi

    log "Stage 1 complete."
}

main "$@"
