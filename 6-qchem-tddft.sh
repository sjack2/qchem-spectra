#!/usr/bin/env bash
# ============================================================================
# 6-qchem-tddft.sh -- TD-DFT excited-state calculations (Stage 6)
# ============================================================================
#
# OVERVIEW
#   For every conformer that passed the Boltzmann filter (Stage 5), this
#   script builds a Q-Chem TD-DFT input with implicit solvent (SMD) and
#   either runs it locally or submits a SLURM job.
#
#   The resulting outputs contain UV-Vis oscillator strengths and ECD
#   rotatory strengths for all requested excited states. Q-Chem prints
#   both automatically for every TD-DFT calculation.
#
# Usage:
#   6-qchem-tddft.sh TAG
#   6-qchem-tddft.sh --local --cpus 4 --method CAM-B3LYP pna
#   6-qchem-tddft.sh --list mols.txt --dry-run
#
# Flags:
#   -m | --method NAME         TD-DFT functional               [CAM-B3LYP]
#   -b | --basis NAME          Basis set                        [def2-SVP]
#        --roots N             Number of excited states          [30]
#        --solvent NAME        SMD solvent keyword               [water]
#        --max-scf N           Max SCF cycles                    [150]
#   -c | --cpus N              CPU cores                         [12]
#        --mem-per-cpu MB      Memory per core in MB             [2048]
#        --qchem-setup PATH    Path to Q-Chem setup script       [auto]
#        --list FILE           File of molecule TAGs
#        --local               Run Q-Chem directly (no SLURM)
#        --dry-run             Write inputs but do not run
#   -h | --help                Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --partition NAME      SLURM partition                   [general]
#        --time HH:MM:SS       Wall-clock limit                  [04:00:00]
#
# Cluster configuration (cluster.cfg):
#   Create a file named cluster.cfg in the same directory as this script
#   to set site-specific defaults without using flags every time.
#   See cluster.cfg.example for the full list of supported variables.
#   Variables set by cluster.cfg are overridden by explicit CLI flags.
#
# Directory layout:
#   <TAG>/
#   |-- 03_solvent_opt/<CID>/<CID>.out      <- Stage 4 geometry source
#   |-- 04_boltzmann/<TAG>_bw_labels.dat    <- conformer list (from Stage 5)
#   -- 05_tddft/
#       |-- <CID>/
#       |   |-- <CID>.inp                Q-Chem TD-DFT input
#       |   |-- <CID>.out                Q-Chem output
#       |   -- <CID>.slurm             (HPC mode only)
#       -- ...
#
# Examples:
#   # Local: compute ECD for all populated conformers of pna
#   6-qchem-tddft.sh --local --cpus 4 --roots 30 pna
#
#   # HPC: submit with a range-separated hybrid
#   6-qchem-tddft.sh --method CAM-B3LYP --partition general ephedrine
#
#   # Dry run: inspect inputs
#   6-qchem-tddft.sh --dry-run --local --list molecules.txt
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
# HELPERS
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
# Q-CHEM SETUP RESOLUTION
# ============================================================================
resolve_qchem_setup() {
    local flag_val=$1
    if [[ -n $flag_val ]]; then printf '%s' "$flag_val"; return; fi
    if [[ -n ${QCHEM_SETUP:-} ]]; then printf '%s' "$QCHEM_SETUP"; return; fi
    if command -v qchem >/dev/null 2>&1; then printf ''; return; fi
    printf ''
}

detect_mode() {
    local force_local=$1
    if $force_local; then printf 'local'
    elif command -v sbatch >/dev/null 2>&1; then printf 'slurm'
    else printf 'local'; fi
}

# ============================================================================
# Q-CHEM OUTPUT -> XYZ EXTRACTOR (last geometry in file)
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
    dry_run=false
    force_local=false
    list_file=""
    single_tag=""
    qchem_setup_flag=""

    local opts
    opts=$(getopt -o hb:m:c: \
        --long help,method:,basis:,roots:,solvent:,max-scf:,cpus:,\
mem-per-cpu:,partition:,time:,list:,qchem-setup:,local,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)      method=$2;             shift 2 ;;
            -b|--basis)       basis=$2;              shift 2 ;;
            --roots)          roots=$2;              shift 2 ;;
            --solvent)        solvent=$2;            shift 2 ;;
            --max-scf)        max_scf=$2;            shift 2 ;;
            -c|--cpus)        cpus=$2;               shift 2 ;;
            --mem-per-cpu)    mem_mb=$2;             shift 2 ;;
            --partition)      partition=$2;          shift 2 ;;
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

    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional args not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide exactly one molecule TAG"
        single_tag=$1
    fi
}

# ============================================================================
# Q-CHEM INPUT WRITER
# ============================================================================
write_qchem_input() {
    local cid=$1 xyz_file=$2 inp_file=$3

    cat >"$inp_file" <<EOF
\$comment
TD-DFT UV-Vis + ECD (SMD ${solvent}) -- generated by 6-qchem-tddft.sh
\$end

\$molecule
${charge} ${mult}
$(tail -n +3 "$xyz_file")
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
# SLURM WRITER
# ============================================================================
write_slurm() {
    local cid=$1 inp_file=$2 slurm_file=$3 workdir=$4
    local abs_workdir
    abs_workdir=$(cd "$workdir" && pwd)

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=tddft_${cid}
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

qchem -nt ${cpus} ${abs_workdir}/${cid}.inp ${abs_workdir}/${cid}.out
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# PROCESS ONE CONFORMER
# ============================================================================
process_conformer() {
    local tag=$1 cid=$2 outfile=$3

    local dir="${tag}/${OUT_SUBDIR}/${cid}"
    mkdir -p "$dir"

    # extract geometry from Stage 4 .out file into a temp XYZ
    local xyz_tmp
    xyz_tmp=$(mktemp --suffix=.xyz)
    qcout2xyz "$outfile" >"$xyz_tmp" || {
        warn "  [${cid}] cannot extract geometry from ${outfile}"
        rm -f "$xyz_tmp"; return
    }

    local inp_file="${dir}/${cid}.inp"
    local log_file="${dir}/${cid}.out"

    write_qchem_input "$cid" "$xyz_tmp" "$inp_file"
    rm -f "$xyz_tmp"

    if $dry_run; then
        log "  [${cid}] dry run -- input written to ${inp_file}"
        return
    fi

    if [[ $exec_mode == slurm ]]; then
        local slurm_file="${dir}/${cid}.slurm"
        write_slurm "$cid" "$inp_file" "$slurm_file" "$dir"
        log "  [${cid}] submitting (partition=${partition})"
        (cd "$dir" && sbatch "$(basename "$slurm_file")")
    else
        log "  [${cid}] running TD-DFT (${roots} roots, solvent=${solvent})"
        if [[ -n $qchem_setup ]]; then
            # shellcheck source=/dev/null
            source "$qchem_setup"
        fi
        qchem -nt "$cpus" "$inp_file" "$log_file" 2>&1
        if grep -q "TDDFT/TDA Excitation Energies\|Excited state" "$log_file" 2>/dev/null; then
            log "  [${cid}] TD-DFT completed successfully"
        else
            warn "  [${cid}] excited-state output not found -- check ${log_file}"
        fi
    fi
}

# ============================================================================
# PROCESS ONE MOLECULE (all its Boltzmann-filtered conformers)
# ============================================================================
process_tag() {
    local tag=$1

    # find the label file from Stage 5
    local lab_file="${tag}/${BW_SUBDIR}/${tag}_bw_labels.dat"
    if [[ ! -f $lab_file ]]; then
        warn "[${tag}] label file not found: ${lab_file}"
        return
    fi

    # read conformer IDs
    local cids=()
    mapfile -t cids < <(grep -v '^[[:space:]]*$' "$lab_file")
    if (( ${#cids[@]} == 0 )); then
        warn "[${tag}] no conformers in label file"
        return
    fi

    # read charge/multiplicity from the original XYZ
    charge=0; mult=1
    local pre_xyz="${XYZ_DIR}/${tag}.xyz"
    if [[ -f $pre_xyz ]]; then
        read_xyz_header "$pre_xyz"
    else
        warn "[${tag}] ${pre_xyz} not found -- using charge=0 mult=1"
    fi

    log "[${tag}] ${#cids[@]} conformers to process (charge=${charge}, mult=${mult})"

    for cid in "${cids[@]}"; do
        cid=$(echo "$cid" | xargs)  # trim whitespace
        [[ -z $cid ]] && continue

        # look for the Stage 4 output file
        local outpath="${tag}/${SOLV_OUT_DIR}/${cid}/${cid}.out"
        if [[ ! -f $outpath ]]; then
            warn "  [${cid}] solvent optimization output not found -- skipping"
            continue
        fi
        process_conformer "$tag" "$cid" "$outpath"
    done
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 6: TD-DFT Excited-State Calculation (UV-Vis + ECD)
-------------------------------------------------------------
 Mode         : ${exec_mode}
 Q-Chem setup : ${qchem_setup:-<qchem in PATH>}
 Method       : ${method}
 Basis        : ${basis}
 Roots        : ${roots}
 Solvent      : ${solvent} (SMD)
 Max SCF      : ${max_scf}
 Cores        : ${cpus}
 Mem/core     : ${mem_mb} MB
 Dry run      : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"

    exec_mode=$(detect_mode "$force_local")
    qchem_setup=$(resolve_qchem_setup "$qchem_setup_flag")

    if ! $dry_run && [[ -z $qchem_setup ]]; then
        if ! command -v qchem >/dev/null 2>&1; then
            die "Q-Chem not found. Set QCHEM_SETUP in cluster.cfg, use --qchem-setup, or add qchem to PATH."
        fi
    fi
    if ! $dry_run && [[ -n $qchem_setup && ! -f $qchem_setup ]]; then
        die "Q-Chem setup script '${qchem_setup}' not found. Check QCHEM_SETUP in cluster.cfg or --qchem-setup."
    fi
    if [[ -z $qchem_setup ]] && $dry_run; then
        warn "Q-Chem setup not found -- dry-run will proceed without it"
    fi

    print_banner

    if [[ -n $list_file ]]; then
        while IFS= read -r tag || [[ -n $tag ]]; do
            [[ -z $tag || $tag == \#* || $tag == \;* ]] && continue
            process_tag "$tag"
        done < "$list_file"
    else
        process_tag "$single_tag"
    fi

    log "Stage 6-tddft complete."
}

main "$@"
