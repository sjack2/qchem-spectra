#!/usr/bin/env bash
# ============================================================================
# 2b-qchem-crest-conf-search.sh -- Conformer enumeration with CREST (Stage 2b)
# ============================================================================
#
# OVERVIEW
#   Alternative to Stage 2 (Confab). Uses CREST's GFN2-xTB metadynamics
#   to explore the conformational landscape more thoroughly than force-field
#   methods. Better for flexible molecules; slower but more reliable.
#
#   The optimized geometry is extracted automatically from the Q-Chem Stage 1
#   output file (<TAG>/01_gas_opt/<TAG>.out). If the geometry was already
#   extracted by a previous run it is reused directly.
#
#   Execution mode:
#     --local   -> runs CREST directly in the current shell (blocking)
#     (default) -> submits a SLURM batch job if sbatch is available
#
# Usage:
#   2b-qchem-crest-conf-search.sh TAG
#   2b-qchem-crest-conf-search.sh --local --cpus 4 ephedrine
#   2b-qchem-crest-conf-search.sh --list mols.txt --dry-run
#
# Flags:
#   -c | --cpus N         CPU threads for CREST            [4]
#        --list FILE      File of TAGs or XYZ paths
#        --local          Run CREST directly (no SLURM)
#        --dry-run        Echo actions without running
#        --pre-xyz        Also look in pre_xyz/ for geometries
#   -h | --help           Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --mem MB         Memory for SLURM job (MB)         [4096]
#        --partition NAME SLURM partition                    [general]
#        --time HH:MM:SS  SLURM wall time                    [06:00:00]
#
# Directory layout:
#   <TAG>/
#   |-- 01_gas_opt/<TAG>.out           <- Q-Chem Stage 1 output (geometry source)
#   -- 02_conf_search/
#       |-- <TAG>.xyz               extracted geometry
#       |-- crest_conformers.xyz    CREST output (multi-structure XYZ)
#       |-- crest.log               CREST log file
#       -- crest.energies          parsed energy table
#
#   After running this, use 3b-qchem-crest-conf-split.sh to split conformers
#   into individual XYZ files.
#
# Dependencies: crest, python3
#
# Examples:
#   2b-qchem-crest-conf-search.sh --local --cpus 8 ephedrine
#   2b-qchem-crest-conf-search.sh --list molecules.txt --partition general
#   2b-qchem-crest-conf-search.sh --dry-run --local aspirin
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_CPUS=4
DEFAULT_MEM_MB=4096
DEFAULT_PARTITION="general"
DEFAULT_WALL="06:00:00"
XYZ_DIR="pre_xyz"

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
# Q-CHEM OUTPUT -> XYZ EXTRACTOR
# ============================================================================
# Creates a temporary Python script that extracts the last optimized
# geometry from a Q-Chem .out file. Written once and reused across all
# molecules in a batch run.
make_qxyz_py() {
    QXYZ_PY=$(mktemp)
    cat >"$QXYZ_PY" <<'PY'
#!/usr/bin/env python3
"""Extract the last optimized geometry from a Q-Chem output file."""
import sys, re, pathlib

if len(sys.argv) != 2:
    sys.exit("Usage: qxyz.py file.out")

lines = pathlib.Path(sys.argv[1]).read_text(errors="ignore").splitlines()
hdr = re.compile(
    r'(Standard Nuclear Orientation|Input orientation|Coordinates \(Angstroms\))',
    re.I
)

start = None
for idx in range(len(lines) - 1, -1, -1):
    if hdr.search(lines[idx]):
        start = idx + 3   # skip header line + column titles + dashes
        break

if start is None:
    sys.exit(f"qxyz: no coordinate block found in {sys.argv[1]}")

geom = []
for line in lines[start:]:
    if not line.strip() or line.lstrip().startswith('-'):
        break
    parts = line.split()
    if len(parts) >= 5:       # index  symbol  x  y  z
        symbol, x, y, z = parts[1], parts[2], parts[3], parts[4]
    elif len(parts) == 4:     # symbol  x  y  z
        symbol, x, y, z = parts
    else:
        continue
    geom.append(f"{symbol} {x} {y} {z}")

if not geom:
    sys.exit(f"qxyz: coordinate parsing failed in {sys.argv[1]}")

print(len(geom))
print()
print("\n".join(geom))
PY
    chmod +x "$QXYZ_PY"
}

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    cpus=$DEFAULT_CPUS
    mem_mb=$DEFAULT_MEM_MB
    partition=$DEFAULT_PARTITION
    wall=$DEFAULT_WALL
    dry_run=false
    force_local=false
    pre_xyz_fallback=false
    list_file=""
    worker_tag=""
    positional=()

    local opts
    opts=$(getopt -o hc: \
        --long help,cpus:,mem:,partition:,time:,list:,local,dry-run,pre-xyz,worker: -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -c|--cpus)    cpus=$2;               shift 2 ;;
            --mem)        mem_mb=$2;             shift 2 ;;
            --partition)  partition=$2;           shift 2 ;;
            --time)       wall=$2;               shift 2 ;;
            --list)       list_file=$2;          shift 2 ;;
            --local)      force_local=true;      shift ;;
            --dry-run)    dry_run=true;          shift ;;
            --pre-xyz)    pre_xyz_fallback=true; shift ;;
            --worker)     worker_tag=$2;         shift 2 ;;
            -h|--help)    show_help ;;
            --)           shift; break ;;
            *)            die "Unknown option '$1'" ;;
        esac
    done
    positional+=("$@")

    # --worker is an internal flag used by SLURM jobs to re-invoke this script
    if [[ -n $worker_tag ]]; then
        run_worker "$worker_tag"
        exit 0
    fi

    if [[ -n $list_file ]]; then
        [[ ${#positional[@]} -eq 0 ]] || die "Positional args invalid with --list"
        require_file "$list_file"
    else
        [[ ${#positional[@]} -gt 0 ]] || die "Provide at least one TAG"
    fi

    command -v crest >/dev/null 2>&1 || \
        { $force_local && die "crest not found in PATH"; } || \
        warn "crest not found in PATH -- SLURM jobs will need it on compute nodes"
    command -v python3 >/dev/null 2>&1 || die "python3 not found in PATH"
}

# ============================================================================
# FIND OR EXTRACT OPTIMIZED GEOMETRY
# ============================================================================
# Returns the path to a usable XYZ file. If only a Q-Chem .out is available,
# extracts the geometry on the fly and writes it to 02_conf_search/<TAG>.xyz.
find_optimized_xyz() {
    local tag=$1

    # already extracted from a prior run
    local probe="${tag}/02_conf_search/${tag}.xyz"
    [[ -f $probe ]] && { printf '%s' "$probe"; return; }

    # Stage 1 Q-Chem output -- extract geometry on the fly
    probe="${tag}/01_gas_opt/${tag}.out"
    if [[ -f $probe ]]; then
        local out_dir="${tag}/02_conf_search"
        mkdir -p "$out_dir"
        local xyz_out="${out_dir}/${tag}.xyz"
        if python3 "$QXYZ_PY" "$probe" > "$xyz_out"; then
            printf '%s' "$xyz_out"
            return
        else
            warn "[${tag}] geometry extraction failed from ${probe}"
            rm -f "$xyz_out"
            printf ''
            return
        fi
    fi

    # pre_xyz fallback (if --pre-xyz flag supplied)
    if $pre_xyz_fallback; then
        probe="${XYZ_DIR}/${tag}.xyz"
        [[ -f $probe ]] && { printf '%s' "$probe"; return; }
    fi

    printf ''
}

# ============================================================================
# WORKER: runs inside the job directory (called directly or via SLURM)
# ============================================================================
run_worker() {
    local tag=$1
    local xyz="${tag}.xyz"

    [[ -f $xyz ]] || die "[worker] ${xyz} not found in $(pwd)"
    command -v crest >/dev/null 2>&1 || die "[worker] crest not in PATH"

    log "[${tag}] launching CREST with ${cpus} threads"
    crest "$xyz" --gfn2 --T "$cpus" > crest.log 2>&1

    [[ -s crest_conformers.xyz ]] || die "[${tag}] crest_conformers.xyz not produced"

    # parse a quick energy summary
    if grep -q 'TOTAL' crest.log; then
        awk '/TOTAL/{printf "%17.10f %8.3f %s\n",$4,$6,$7}' crest.log > crest.energies
    fi

    local n_conf
    n_conf=$(grep -c '^[[:space:]]*[0-9]' crest_conformers.xyz 2>/dev/null || echo '?')
    log "[${tag}] CREST finished -- ${n_conf} structures"
}

# ============================================================================
# PROCESS ONE TAG
# ============================================================================
process_tag() {
    local tag=$1 xyz_path=$2

    local job_dir="${tag}/02_conf_search"
    mkdir -p "$job_dir"
    # xyz_path may already be inside job_dir (extracted on the fly); skip if same file
    [[ "$xyz_path" -ef "${job_dir}/${tag}.xyz" ]] || cp -f "$xyz_path" "${job_dir}/${tag}.xyz"

    # ---- dry run -------------------------------------------------------
    if $dry_run; then
        if [[ $exec_mode == slurm ]]; then
            log "[${tag}] (dry run) would submit CREST SLURM job in ${job_dir}"
        else
            log "[${tag}] (dry run) would run: crest ${tag}.xyz --gfn2 --T ${cpus}"
        fi
        return
    fi

    # ---- SLURM mode: submit batch job ----------------------------------
    if [[ $exec_mode == slurm ]]; then
        local self_path
        self_path=$(realpath "$0")
        local jobid
        jobid=$(sbatch \
            --chdir="$(realpath "$job_dir")" \
            --job-name="${tag}_crest" \
            --partition="$partition" \
            --cpus-per-task="$cpus" \
            --mem="${mem_mb}" \
            --time="$wall" \
            --output="slurm-%j.out" \
            --error="slurm-%j.err" \
            --wrap="bash ${self_path} --worker ${tag} --cpus ${cpus}" \
            | awk '{print $4}')
        log "[${tag}] submitted as SLURM job ${jobid}"
        return
    fi

    # ---- local mode: run directly --------------------------------------
    log "[${tag}] running CREST locally (${cpus} threads)"
    (cd "$job_dir" && run_worker "$tag")
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 2b: Conformer Search (CREST / GFN2-xTB)
-------------------------------------------------------------
 Mode        : ${exec_mode}
 CPU threads : ${cpus}
 Dry run     : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    parse_cli "$@"
    make_qxyz_py

    exec_mode=$(detect_mode "$force_local")
    print_banner

    local entries=()
    if [[ -n $list_file ]]; then
        mapfile -t entries < <(grep -vE '^[[:space:]]*(#|;|$)' "$list_file")
    else
        entries=("${positional[@]}")
    fi

    for entry in "${entries[@]}"; do
        local tag xyz_path
        if [[ $entry == *.xyz ]]; then
            tag=$(basename "$entry" .xyz)
            xyz_path=$entry
            require_file "$xyz_path"
        else
            tag=$entry
            xyz_path=$(find_optimized_xyz "$tag")
            if [[ -z $xyz_path ]]; then
                warn "[${tag}] optimized geometry not found (checked 01_gas_opt/${tag}.out and 02_conf_search/${tag}.xyz) -- skipping"
                continue
            fi
        fi
        process_tag "$tag" "$xyz_path"
    done

    rm -f "$QXYZ_PY"
    log "Stage 2b complete."
}

main "$@"
