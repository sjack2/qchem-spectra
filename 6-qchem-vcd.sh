#!/usr/bin/env bash
# ============================================================================
# 6-qchem-vcd.sh -- IR + VCD frequency calculations (Stage 6-vcd)
# ============================================================================
#
# OVERVIEW
#   For every conformer that passed the Boltzmann filter (Stage 5), this
#   script builds a Q-Chem frequency input requesting VCD and Raman with
#   implicit solvent (SMD) and either runs it locally or submits a SLURM
#   job array.
#
#   The resulting outputs contain IR absorption intensities and VCD
#   rotatory strengths for all normal modes. Post-process with
#   qc_vcd_ir_tools.py to produce Boltzmann-weighted broadened spectra.
#
#   This script is a parallel alternative to 6-qchem-tddft.sh (TD-DFT).
#   Both read from the same Stage 5 output and Stage 4 geometries;
#   they write to separate directories (05_vcd/ vs 05_tddft/).
#
# Usage:
#   6-qchem-vcd.sh TAG
#   6-qchem-vcd.sh --local --cpus 4 --method B3LYP methyloxirane
#   6-qchem-vcd.sh --list mols.txt --dry-run
#
# Flags:
#   -m | --method NAME         DFT functional                   [B3LYP]
#   -b | --basis NAME          Basis set                        [6-31+G(d)]
#        --solvent NAME        SMD solvent keyword               [water]
#        --max-scf N           Max SCF cycles                    [150]
#   -c | --cpus N              CPU cores                         [12]
#        --mem-per-cpu MB      Memory per core in MB             [2048]
#        --max-running N       Max simultaneous SLURM array tasks [10]
#        --qchem-setup PATH    Path to Q-Chem setup script       [auto]
#        --list FILE           File of molecule TAGs
#        --local               Run Q-Chem directly (no SLURM)
#        --dry-run             Write inputs but do not run
#   -h | --help                Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --partition NAME      SLURM partition                   [general]
#        --time HH:MM:SS       Wall-clock limit                  [06:00:00]
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
#   -- 05_vcd/
#       |-- <TAG>_conf_list.txt          (SLURM mode: list of working dirs)
#       |-- <TAG>_array.slurm            (SLURM mode: array job script)
#       |-- <CID>/
#       |   |-- <CID>.inp                Q-Chem FREQ+VCD input
#       |   -- <CID>.out                Q-Chem output
#       -- ...
#
# Post-processing:
#   qc_vcd_ir_tools.py --bw <TAG>/04_boltzmann/<TAG>_energies.dat \
#       --prefix <TAG>_vib --stick <TAG>/05_vcd/**/*.out
#
# Examples:
#   # Local: compute IR/VCD for all populated conformers of methyloxirane
#   6-qchem-vcd.sh --local --cpus 4 --solvent water methyloxirane
#
#   # HPC: submit with B3LYP and a larger basis
#   6-qchem-vcd.sh --method B3LYP --basis 6-311+G(d,p) aspirin
#
#   # Dry run: inspect inputs
#   6-qchem-vcd.sh --dry-run --local --list molecules.txt
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_METHOD="B3LYP"
DEFAULT_BASIS="6-31+G(d)"
DEFAULT_SOLVENT="water"
DEFAULT_MAX_SCF=150
DEFAULT_CPUS=4
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="general"
DEFAULT_WALL="06:00:00"
DEFAULT_MAX_RUNNING=10

XYZ_DIR="pre_xyz"
SOLV_OUT_DIR="03_solvent_opt"
BW_SUBDIR="04_boltzmann"
OUT_SUBDIR="05_vcd"

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
    solvent=$DEFAULT_SOLVENT
    max_scf=$DEFAULT_MAX_SCF
    cpus=$DEFAULT_CPUS
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=${CLUSTER_PARTITION:-$DEFAULT_PARTITION}
    wall=${CLUSTER_WALL:-$DEFAULT_WALL}
    max_running=${CLUSTER_MAX_RUNNING:-$DEFAULT_MAX_RUNNING}
    dry_run=false
    force_local=false
    list_file=""
    single_tag=""
    qchem_setup_flag=""

    local opts
    opts=$(getopt -o hb:m:c: \
        --long help,method:,basis:,solvent:,max-scf:,cpus:,\
mem-per-cpu:,max-running:,partition:,time:,list:,qchem-setup:,local,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)      method=$2;             shift 2 ;;
            -b|--basis)       basis=$2;              shift 2 ;;
            --solvent)        solvent=$2;            shift 2 ;;
            --max-scf)        max_scf=$2;            shift 2 ;;
            -c|--cpus)        cpus=$2;               shift 2 ;;
            --mem-per-cpu)    mem_mb=$2;             shift 2 ;;
            --max-running)    max_running=$2;        shift 2 ;;
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
    local mem_total=$(( cpus * mem_mb ))

    cat >"$inp_file" <<EOF
\$comment
IR + VCD (SMD ${solvent}) -- generated by 6-qchem-vcd.sh
\$end

\$molecule
${charge} ${mult}
$(tail -n +3 "$xyz_file")
\$end

\$rem
  JOB_TYPE            FREQ
  METHOD              ${method}
  BASIS               ${basis}
  VCD                 TRUE
  RAMAN               TRUE
  SOLVENT_METHOD      SMD
  MAX_SCF_CYCLES      ${max_scf}
  SCF_CONVERGENCE     8
  SYM_IGNORE          TRUE
  XC_GRID             3
  MEM_TOTAL           ${mem_total}
  MEM_STATIC          500
\$end
EOF

    # Q-Chem uses water as the default SMD solvent, so only add $smx
    # for non-water solvents
    if [[ ${solvent,,} != water ]]; then
        cat >>"$inp_file" <<EOF

\$smx
  solvent ${solvent}
\$end
EOF
    fi
}

# ============================================================================
# SLURM ARRAY WRITER
# ============================================================================
write_array_slurm() {
    local tag=$1 conf_list=$2 slurm_file=$3 n=$4
    local abs_list out_dir
    abs_list=$(cd "$(dirname "$conf_list")" && pwd)/$(basename "$conf_list")
    out_dir=$(dirname "$abs_list")

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=vcd_${tag}
#SBATCH --partition=${partition}
#SBATCH --array=1-${n}%${max_running}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem-per-cpu=${mem_mb}
#SBATCH --time=${wall}
#SBATCH --output=${out_dir}/slurm-%A_%a.out
#SBATCH --error=${out_dir}/slurm-%A_%a.err

# ---- Q-Chem environment ----
${qchem_setup:+source ${qchem_setup}}
${CLUSTER_LD_LIBRARY_PATH:+export LD_LIBRARY_PATH=${CLUSTER_LD_LIBRARY_PATH}:\$LD_LIBRARY_PATH}
export QCSCRATCH=/tmp/\$SLURM_JOB_ID

WORKDIR=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${abs_list}")
CID=\$(basename "\$WORKDIR")
cd "\$WORKDIR"
qchem -nt ${cpus} "\${CID}.inp" "\${CID}.out"
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# PROCESS ONE CONFORMER (local mode only)
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

    log "  [${cid}] running FREQ (VCD, solvent=${solvent})"
    if [[ -n $qchem_setup ]]; then
        # shellcheck source=/dev/null
        source "$qchem_setup"
    fi
    qchem -nt "$cpus" "$inp_file" "$log_file" 2>&1
    if grep -q "VIBRATIONAL ANALYSIS\|Mode:" "$log_file" 2>/dev/null; then
        log "  [${cid}] frequency calculation completed successfully"
    else
        warn "  [${cid}] expected output blocks not found -- check ${log_file}"
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

    # ---- local mode: run each conformer directly -------------------------
    if [[ $exec_mode != slurm ]]; then
        for cid in "${cids[@]}"; do
            cid=$(echo "$cid" | xargs)
            [[ -z $cid ]] && continue
            local outpath="${tag}/${SOLV_OUT_DIR}/${cid}/${cid}.out"
            if [[ ! -f $outpath ]]; then
                warn "  [${cid}] solvent optimization output not found -- skipping"
                continue
            fi
            process_conformer "$tag" "$cid" "$outpath"
        done
        return
    fi

    # ---- SLURM mode: write all inputs, then submit one array job --------
    local out_dir="${tag}/${OUT_SUBDIR}"
    mkdir -p "$out_dir"

    local conf_dirs=()
    local cid outpath dir xyz_tmp inp_file abs_dir
    for cid in "${cids[@]}"; do
        cid=$(echo "$cid" | xargs)
        [[ -z $cid ]] && continue

        outpath="${tag}/${SOLV_OUT_DIR}/${cid}/${cid}.out"
        if [[ ! -f $outpath ]]; then
            warn "  [${cid}] solvent optimization output not found -- skipping"
            continue
        fi

        dir="${tag}/${OUT_SUBDIR}/${cid}"
        mkdir -p "$dir"

        xyz_tmp=$(mktemp --suffix=.xyz)
        if ! qcout2xyz "$outpath" >"$xyz_tmp" 2>/dev/null; then
            warn "  [${cid}] cannot extract geometry from ${outpath} -- skipping"
            rm -f "$xyz_tmp"
            continue
        fi

        inp_file="${dir}/${cid}.inp"
        write_qchem_input "$cid" "$xyz_tmp" "$inp_file"
        rm -f "$xyz_tmp"

        if $dry_run; then
            log "  [${cid}] dry run -- input written to ${inp_file}"
            continue
        fi

        abs_dir=$(cd "$dir" && pwd)
        conf_dirs+=("$abs_dir")
    done

    $dry_run && return

    local n=${#conf_dirs[@]}
    if (( n == 0 )); then
        warn "[${tag}] no valid conformer inputs written -- nothing to submit"
        return
    fi

    local conf_list="${out_dir}/${tag}_conf_list.txt"
    printf '%s\n' "${conf_dirs[@]}" >"$conf_list"

    local array_slurm="${out_dir}/${tag}_array.slurm"
    write_array_slurm "$tag" "$conf_list" "$array_slurm" "$n"

    local jid
    jid=$(sbatch --parsable "$array_slurm")
    log "[${tag}] submitted array job ${jid} (${n} conformers, max ${max_running} concurrent)"
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 6-vcd: Frequency Calculation (IR + VCD)
-------------------------------------------------------------
 Mode         : ${exec_mode}
 Q-Chem setup : ${qchem_setup:-<qchem in PATH>}
 Method       : ${method}
 Basis        : ${basis}
 Solvent      : ${solvent} (SMD)
 Max SCF      : ${max_scf}
 Cores        : ${cpus}
 Mem/core     : ${mem_mb} MB
 Max running  : ${max_running} (SLURM array throttle)
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

    log "Stage 6-vcd complete."
}

main "$@"
