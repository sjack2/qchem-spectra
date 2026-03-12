#!/usr/bin/env bash
# ============================================================================
# 1-qchem-init-opt.sh -- Gas-phase geometry optimization (Stage 1)
# ============================================================================
#
# OVERVIEW
#   For each supplied molecule this script:
#     1. Reads charge and multiplicity from the XYZ comment line,
#     2. Writes a Q-Chem geometry-optimization input file,
#     3. Either runs Q-Chem directly (--local) or submits a SLURM job.
#
#   Execution mode is chosen automatically: if sbatch is available and
#   --local is not set, SLURM mode is used. Otherwise, Q-Chem is run
#   directly in the current shell.
#
# Usage:
#   1-qchem-init-opt.sh TAG
#   1-qchem-init-opt.sh path/to/foo.xyz --method PBE0
#   1-qchem-init-opt.sh --list molecules.txt --dry-run
#   1-qchem-init-opt.sh --local --cpus 4 aspirin
#
# Flags:
#   -m | --method METHOD           DFT functional               [B3LYP]
#   -b | --basis  BASIS            Basis set                    [def2-SVP]
#        --disp {auto,none,D3BJ}   Dispersion correction        [auto]
#        --max-iter N              SCF iteration limit           [150]
#   -c | --cpus N                  CPU cores (threads + SLURM)  [12]
#   -g | --grid {SG1,SG2,SG3}     Integration grid              [SG3]
#        --mem-per-cpu MB          Memory per core in MB         [2048]
#        --qchem-setup PATH        Path to Q-Chem setup script   [auto]
#        --list FILE               File of TAGs or XYZ paths
#        --local                   Run Q-Chem directly (no SLURM)
#        --dry-run                 Write inputs but do not run
#   -h | --help                    Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --partition NAME          SLURM partition               [general]
#        --time HH:MM:SS           Wall-clock limit              [02:00:00]
#
# Cluster configuration (cluster.cfg):
#   Create a file named cluster.cfg in the same directory as this script
#   to set site-specific defaults without using flags every time.
#   See cluster.cfg.example for the full list of supported variables.
#   Variables set by cluster.cfg are overridden by explicit CLI flags.
#
# XYZ charge/multiplicity convention:
#   Line 2 of each .xyz file may specify charge and multiplicity in
#   either of two formats:
#     charge=0 mult=1       (key=value tokens, any order)
#     0 1                   (bare integers: charge multiplicity)
#   If neither is found, the defaults charge=0 mult=1 are used.
#
# Directory layout produced:
#   <TAG>/
#   -- 01_gas_opt/
#       |-- <TAG>.inp          Q-Chem input
#       |-- <TAG>.out          Q-Chem output  (after execution)
#       -- <TAG>.slurm        SLURM script (HPC mode only)
#
# Examples:
#   # Local workstation, 4 cores, default method
#   1-qchem-init-opt.sh --local --cpus 4 aspirin
#
#   # HPC cluster, PBE0 functional, submit to SLURM
#   1-qchem-init-opt.sh --method PBE0 --partition gpu aspirin
#
#   # Dry run: inspect inputs without running anything
#   1-qchem-init-opt.sh --dry-run --local --list molecules.txt
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
DEFAULT_MAX_ITER=150
DEFAULT_CPUS=12
DEFAULT_GRID="SG3"
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="general"
DEFAULT_WALL="02:00:00"

XYZ_DIR="pre_xyz"           # where starting geometries live

# ============================================================================
# CLUSTER CONFIG
# ============================================================================
# Loads cluster.cfg from the script's own directory if present.
# Supported variables (all optional):
#   QCHEM_SETUP              -- path to Q-Chem environment setup script
#   CLUSTER_PARTITION         -- default SLURM partition
#   CLUSTER_WALL              -- default SLURM wall-clock limit (HH:MM:SS)
#   CLUSTER_LD_LIBRARY_PATH   -- extra LD_LIBRARY_PATH for compute nodes
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
# Q-CHEM SETUP RESOLUTION
# ============================================================================
# Priority: --qchem-setup flag  >  QCHEM_SETUP env var (incl. cluster.cfg)
#           >  qchem already on PATH (no setup needed)
# The variable qchem_setup is set during CLI parsing and finalized in main().
resolve_qchem_setup() {
    # $1 = value from --qchem-setup flag (may be empty)
    local flag_val=$1

    if [[ -n $flag_val ]]; then
        printf '%s' "$flag_val"
        return
    fi

    if [[ -n ${QCHEM_SETUP:-} ]]; then
        printf '%s' "$QCHEM_SETUP"
        return
    fi

    # if qchem is already on PATH, no setup script needed
    if command -v qchem >/dev/null 2>&1; then
        printf ''
        return
    fi

    # nothing found -- caller will decide whether this is fatal
    printf ''
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
die()  { printf 'Error: %s\n' "$*" >&2; exit 1; }
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
warn() { printf '[%s] Warning: %s\n' "$(date '+%F %T')" "$*" >&2; }

require_file() { [[ -f $1 ]] || die "File '$1' not found"; }

show_help() {
    sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# \{0,1\}//'
    exit 0
}

# ============================================================================
# EXECUTION-MODE DETECTION
# ============================================================================
# Returns "local" or "slurm". --local flag always wins; otherwise we
# check whether sbatch is available.
detect_mode() {
    local force_local=$1
    if $force_local; then
        printf 'local'
    elif command -v sbatch >/dev/null 2>&1; then
        printf 'slurm'
    else
        printf 'local'
    fi
}

# ============================================================================
# CHARGE / MULTIPLICITY PARSER
# ============================================================================
# Reads line 2 of an XYZ file and sets the global variables `charge`
# and `mult`. Supports two formats:
#   charge=0 mult=1     (key=value, any surrounding text)
#   0 1                 (bare integers)
read_xyz_header() {
    local xyz_file=$1
    local header_line
    header_line=$(sed -n '2p' "$xyz_file")

    # try key=value format first
    if [[ $header_line =~ charge[[:space:]]*=[[:space:]]*([+-]?[0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
    else
        charge=0
    fi

    if [[ $header_line =~ mult[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
        mult=${BASH_REMATCH[1]}
    # fall back to bare "int int" format
    elif [[ $header_line =~ ^[[:space:]]*([+-]?[0-9]+)[[:space:]]+([0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
        mult=${BASH_REMATCH[2]}
    else
        mult=1
    fi
}

# ============================================================================
# GRID LOGIC
# ============================================================================
# Maps user-facing grid names (SG1/SG2/SG3) to Q-Chem XC_GRID integer.
# Also accepts bare integers 1-3 for convenience.
grid_number() {
    case ${grid^^} in
        SG1|1) printf '1' ;;
        SG2|2) printf '2' ;;
        SG3|3) printf '3' ;;
        *)     die "--grid must be SG1, SG2, or SG3 (got '${grid}')" ;;
    esac
}

# ============================================================================
# DISPERSION LOGIC
# ============================================================================
# Returns the Q-Chem $rem line for DFT-D3(BJ) dispersion, or an empty
# string if none is needed. In "auto" mode the script avoids doubling
# up on functionals that already include dispersion.
disp_line() {
    local method_upper=${method^^}
    case $disp_mode in
        none|NONE)
            printf '' ;;
        D3BJ|d3bj)
            printf '  DFT_D              D3_BJ' ;;
        auto|AUTO)
            # skip if the method string already contains a dispersion suffix
            if [[ $method =~ (-D[0-9]?|-D3BJ|-D3ZERO|-D4)($|[[:space:]]) ]]; then
                printf ''; return
            fi
            # skip for functionals with built-in dispersion
            case $method_upper in
                WB97X-D|WB97X-D3|WB97X-D4|WB97X-V|WB97XD|WB97M-V|B97-D|B97-D3)
                    printf ''; return ;;
            esac
            # default: add D3BJ
            printf '  DFT_D              D3_BJ'
            ;;
        *)
            die "--disp must be auto, none, or D3BJ (got '$disp_mode')" ;;
    esac
}

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    method=$DEFAULT_METHOD
    basis=$DEFAULT_BASIS
    disp_mode=$DEFAULT_DISP
    max_iter=$DEFAULT_MAX_ITER
    cpus=$DEFAULT_CPUS
    grid=$DEFAULT_GRID
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=${CLUSTER_PARTITION:-$DEFAULT_PARTITION}
    wall=${CLUSTER_WALL:-$DEFAULT_WALL}
    dry_run=false
    force_local=false
    list_file=""
    xyz_arg=""
    qchem_setup_flag=""        # raw value from --qchem-setup

    local opts
    opts=$(getopt -o hb:m:c:g: \
        --long help,basis:,method:,disp:,max-iter:,cpus:,grid:,\
mem-per-cpu:,partition:,time:,list:,dry-run,local,qchem-setup: -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)      method=$2;             shift 2 ;;
            -b|--basis)       basis=$2;              shift 2 ;;
            --disp)           disp_mode=$2;          shift 2 ;;
            --max-iter)       max_iter=$2;           shift 2 ;;
            -c|--cpus)        cpus=$2;               shift 2 ;;
            -g|--grid)        grid=$2;               shift 2 ;;
            --mem-per-cpu)    mem_mb=$2;             shift 2 ;;
            --partition)      partition=$2;           shift 2 ;;
            --time)           wall=$2;               shift 2 ;;
            --list)           list_file=$2;          shift 2 ;;
            --qchem-setup)    qchem_setup_flag=$2;   shift 2 ;;
            --local)          force_local=true;      shift ;;
            --dry-run)        dry_run=true;          shift ;;
            -h|--help)        show_help ;;
            --)               shift; break ;;
            *)                die "Unknown option '$1'" ;;
        esac
    done

    # positional: either --list was given, or exactly one TAG/XYZ
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
    local disp grid_num
    disp=$(disp_line)
    grid_num=$(grid_number)

    cat >"$inp_file" <<EOF
\$rem
  JOBTYPE              OPT
  METHOD               ${method}
  BASIS                ${basis}
  SCF_ALGORITHM        DIIS_GDM
  MAX_SCF_CYCLES       ${max_iter}
  XC_GRID              ${grid_num}
  XC_SMART_GRID        TRUE
  GEN_SCFMAN           TRUE
  MEM_TOTAL            $(( cpus * mem_mb ))
EOF

    # append dispersion line only if the function returned one
    if [[ -n $disp ]]; then
        echo "${disp}" >>"$inp_file"
    fi

    cat >>"$inp_file" <<EOF
\$end

\$molecule
  ${charge} ${mult}
EOF
    tail -n +3 "$xyz_file" >>"$inp_file"
    echo '$end' >>"$inp_file"
}

# ============================================================================
# SLURM SCRIPT WRITER  (HPC mode only)
# ============================================================================
write_slurm() {
    local slurm_file=$1 tag=$2 workdir=$3

    # resolve absolute working directory for the SLURM environment block
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

# ---- Q-Chem environment ----
source ${qchem_setup}
${CLUSTER_LD_LIBRARY_PATH:+export LD_LIBRARY_PATH=${CLUSTER_LD_LIBRARY_PATH}:\$LD_LIBRARY_PATH}
export QCSCRATCH=/tmp/\$SLURM_JOB_ID

qchem -nt ${cpus} ${abs_workdir}/${tag}.inp ${abs_workdir}/${tag}.out
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# PROCESS ONE MOLECULE
# ============================================================================
process_molecule() {
    local xyz_file=$1 tag=$2

    local workdir="${tag}/01_gas_opt"
    mkdir -p "$workdir"

    local inp_file="${workdir}/${tag}.inp"
    local out_file="${workdir}/${tag}.out"

    # parse charge / multiplicity from the XYZ comment line
    read_xyz_header "$xyz_file"

    # write the Q-Chem input
    write_qchem_input "$inp_file" "$xyz_file"

    # ---- dry run: show what would happen and return --------------------
    if $dry_run; then
        log "[${tag}] dry run (${exec_mode} mode)"
        echo "--- Q-Chem input: ${inp_file} ---"
        cat "$inp_file"
        if [[ $exec_mode == slurm ]]; then
            local slurm_file="${workdir}/${tag}.slurm"
            write_slurm "$slurm_file" "$tag" "$workdir"
            echo ""
            echo "--- SLURM script: ${slurm_file} ---"
            cat "$slurm_file"
        else
            echo ""
            echo "--- Would run: qchem -nt ${cpus} ${inp_file} ${out_file} ---"
        fi
        echo ""
        return
    fi

    # ---- live execution ------------------------------------------------
    if [[ $exec_mode == slurm ]]; then
        local slurm_file="${workdir}/${tag}.slurm"
        write_slurm "$slurm_file" "$tag" "$workdir"
        log "[${tag}] submitting to SLURM (partition=${partition})"
        (cd "$workdir" && sbatch "$(basename "$slurm_file")")
    else
        log "[${tag}] running Q-Chem locally (${cpus} cores)"
        if [[ -n $qchem_setup ]]; then
            source "$qchem_setup"
        fi
        qchem -nt "$cpus" "$inp_file" "$out_file" 2>&1
        if grep -q "OPTIMIZATION CONVERGED\|Have a nice day" "$out_file" 2>/dev/null; then
            log "[${tag}] optimization converged"
        else
            warn "[${tag}] convergence string not found -- check ${out_file}"
        fi
    fi
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    local disp
    disp=$(disp_line)
    [[ -n $disp ]] && disp="D3_BJ" || disp="none"
    cat >&2 <<EOF
=============================================================
 Stage 1: Gas-Phase Geometry Optimization
-------------------------------------------------------------
 Mode        : ${exec_mode}
 Q-Chem      : ${qchem_setup:-"(qchem on PATH)"}
 Method      : ${method}
 Basis       : ${basis}
 Dispersion  : ${disp}
 Grid        : ${grid}
 SCF MaxIter : ${max_iter}
 Cores       : ${cpus}
 Memory      : ${mem_mb} MB/core
 Dry run     : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"

    # resolve execution mode and Q-Chem setup path
    exec_mode=$(detect_mode "$force_local")
    qchem_setup=$(resolve_qchem_setup "$qchem_setup_flag")

    # in non-dry-run SLURM mode, a setup script is required
    if ! $dry_run && [[ $exec_mode == slurm && -z $qchem_setup ]]; then
        die "Q-Chem setup script not found. Set QCHEM_SETUP in cluster.cfg, use --qchem-setup, or add qchem to PATH."
    fi
    # verify the resolved setup script is a readable file (if one was given)
    if [[ -n $qchem_setup && ! -f $qchem_setup ]]; then
        if $dry_run; then
            warn "Q-Chem setup script '${qchem_setup}' not found -- dry-run will proceed"
        else
            die "Q-Chem setup script '${qchem_setup}' does not exist or is not readable."
        fi
    fi
    # in local mode without a setup script, qchem must be on PATH
    if ! $dry_run && [[ $exec_mode == local && -z $qchem_setup ]]; then
        command -v qchem >/dev/null 2>&1 \
            || die "qchem not found on PATH. Provide --qchem-setup or set QCHEM_SETUP in cluster.cfg."
    fi

    print_banner

    # process molecules from --list or single positional argument
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
                warn "'${xyz_path}' not found -- skipping"
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
