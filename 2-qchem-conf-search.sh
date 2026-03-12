#!/usr/bin/env bash
# ============================================================================
# 2-qchem-conf-search.sh -- Conformer enumeration with Open Babel Confab
# ============================================================================
#
# OVERVIEW
#   Stage 2 of the Q-Chem workflow. For each molecule it:
#     1. Retrieves the optimized geometry from the Stage 1 Q-Chem output,
#     2. Extracts coordinates to XYZ format using an embedded Python script,
#     3. Converts that geometry to SDF format,
#     4. Runs Open Babel Confab to generate conformers within an energy window.
#
#   Q-Chem does not write a standalone .xyz file, so the geometry must be
#   extracted from the .out file. The embedded Python extractor searches
#   for the last "Standard Nuclear Orientation" block and parses it into
#   standard XYZ format.
#
#   Confab operates by systematically rotating torsion angles and evaluating
#   each resulting geometry with a force field. Structures above the energy
#   cutoff are discarded. For rigid molecules Confab may return only a single
#   conformer; for flexible molecules it can return hundreds.
#
# Usage:
#   2-qchem-conf-search.sh TAG
#   2-qchem-conf-search.sh path/to/foo.xyz --conf 500
#   2-qchem-conf-search.sh --list molecules.txt --ecut 4 --dry-run
#
# Flags:
#   --conf N       Maximum number of conformers              [1000]
#   --ecut E       Energy cutoff in kcal/mol                 [5]
#   --list FILE    Text file of TAGs or XYZ paths
#   --dry-run      Echo commands without executing them
#   -h | --help    Show this help and exit
#
# Directory layout (reads from Stage 1, writes to):
#   <TAG>/
#   |-- 01_gas_opt/
#   |   -- <TAG>.out              <- input (from Stage 1)
#   -- 02_conf_search/
#       |-- <TAG>.xyz              extracted optimized geometry
#       |-- <TAG>.sdf              SDF conversion
#       -- <TAG>_combined.sdf     all conformers in one file (run Stage 3 to split)
#
# Dependencies: obabel (Open Babel >= 3.1), python3
#
# Examples:
#   2-qchem-conf-search.sh --ecut 10 --conf 500 ephedrine
#   2-qchem-conf-search.sh --list molecules.txt --dry-run
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_CONF_COUNT=1000
DEFAULT_ECUT_KCAL=5
XYZ_DIR="pre_xyz"

# ============================================================================
# CLUSTER CONFIG (sourced for consistency; no SLURM used in this stage)
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
die()  { printf 'Error: %s\n' "$*" >&2; exit 1; }
log()  { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
warn() { printf '[%s] Warning: %s\n' "$(date '+%F %T')" "$*" >&2; }

require_file() { [[ -f $1 ]] || die "File '$1' not found"; }

show_help() {
    sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# \{0,1\}//'
    exit 0
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
    conf_count=$DEFAULT_CONF_COUNT
    ecut_kcal=$DEFAULT_ECUT_KCAL
    dry_run=false
    list_file=""
    positional=()

    local opts
    opts=$(getopt -o h --long help,list:,conf:,ecut:,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            --list)    list_file=$2; shift 2 ;;
            --conf)    conf_count=$2; shift 2 ;;
            --ecut)    ecut_kcal=$2; shift 2 ;;
            --dry-run) dry_run=true; shift ;;
            -h|--help) show_help ;;
            --)        shift; break ;;
            *)         die "Unknown option '$1'" ;;
        esac
    done
    positional+=("$@")

    if [[ -n $list_file ]]; then
        [[ ${#positional[@]} -eq 0 ]] || die "Positional args invalid with --list"
        require_file "$list_file"
    else
        [[ ${#positional[@]} -gt 0 ]] || die "Provide at least one TAG or XYZ path"
    fi

    command -v obabel >/dev/null 2>&1 || die "obabel not found in PATH"
    command -v python3 >/dev/null 2>&1 || die "python3 not found in PATH"
}

# ============================================================================
# PATH RESOLUTION
# ============================================================================
# Look for the Q-Chem .out from Stage 1 in the expected location.
find_optimized_out() {
    local tag=$1
    local probe="${tag}/01_gas_opt/${tag}.out"
    if [[ -f $probe ]]; then
        printf '%s' "$probe"
    else
        printf ''
    fi
}

# ============================================================================
# CORE: CONFAB CONFORMER GENERATION
# ============================================================================
process_tag() {
    local tag=$1

    # locate the Stage 1 Q-Chem output
    local qcout
    qcout=$(find_optimized_out "$tag")
    if [[ -z $qcout ]]; then
        warn "[${tag}] Stage 1 output not found (${tag}/01_gas_opt/${tag}.out) -- skipping"
        return
    fi

    local out_dir="${tag}/02_conf_search"
    mkdir -p "$out_dir"

    # extract optimized geometry from Q-Chem output
    local xyz_copy="${out_dir}/${tag}.xyz"
    if ! python3 "$QXYZ_PY" "$qcout" >"$xyz_copy"; then
        warn "[${tag}] geometry extraction failed from ${qcout}"
        return
    fi

    local sdf="${out_dir}/${tag}.sdf"
    local combined="${out_dir}/${tag}_combined.sdf"

    if $dry_run; then
        log "[${tag}] (dry run) python3 qxyz.py ${qcout} > ${xyz_copy}"
        log "[${tag}] (dry run) obabel -ixyz ${xyz_copy} -osdf -O ${sdf}"
        log "[${tag}] (dry run) obabel ${sdf} -O ${combined} --confab --original --conf ${conf_count} --ecutoff ${ecut_kcal}"
        log "[${tag}] (dry run) run Stage 3 to split conformers"
        return
    fi

    # convert XYZ -> SDF (Confab requires SDF input)
    log "[${tag}] converting to SDF"
    obabel -ixyz "$xyz_copy" -osdf -O "$sdf" 2>/dev/null
    if [[ ! -s $sdf ]]; then
        warn "[${tag}] SDF conversion produced empty file -- skipping"
        return
    fi

    # run Confab
    log "[${tag}] running Confab (ecutoff=${ecut_kcal} kcal/mol, max=${conf_count})"
    obabel "$sdf" -O "$combined" \
        --confab --original --conf "$conf_count" --ecutoff "$ecut_kcal" 2>/dev/null

    if [[ ! -s $combined ]]; then
        warn "[${tag}] Confab produced no conformers"
        return
    fi

    local count
    count=$(grep -c '^\$\$\$\$' "$combined" 2>/dev/null || echo 0)
    log "[${tag}] Confab generated ${count} conformer(s) -> ${combined} (run Stage 3 to split)"
}

# ============================================================================
# PROCESS A USER-SUPPLIED XYZ (bypasses Stage 1 extraction)
# ============================================================================
process_xyz() {
    local xyz_path=$1
    local tag
    tag=$(basename "$xyz_path" .xyz)

    local out_dir="${tag}/02_conf_search"
    mkdir -p "$out_dir"

    # copy the user-supplied XYZ into the conformer directory
    local xyz_copy="${out_dir}/${tag}.xyz"
    cp -f "$xyz_path" "$xyz_copy"

    local sdf="${out_dir}/${tag}.sdf"
    local combined="${out_dir}/${tag}_combined.sdf"

    if $dry_run; then
        log "[${tag}] (dry run) obabel -ixyz ${xyz_copy} -osdf -O ${sdf}"
        log "[${tag}] (dry run) obabel ${sdf} -O ${combined} --confab --original --conf ${conf_count} --ecutoff ${ecut_kcal}"
        log "[${tag}] (dry run) run Stage 3 to split conformers"
        return
    fi

    # convert XYZ -> SDF (Confab requires SDF input)
    log "[${tag}] converting to SDF"
    obabel -ixyz "$xyz_copy" -osdf -O "$sdf" 2>/dev/null
    if [[ ! -s $sdf ]]; then
        warn "[${tag}] SDF conversion produced empty file -- skipping"
        return
    fi

    # run Confab
    log "[${tag}] running Confab (ecutoff=${ecut_kcal} kcal/mol, max=${conf_count})"
    obabel "$sdf" -O "$combined" \
        --confab --original --conf "$conf_count" --ecutoff "$ecut_kcal" 2>/dev/null

    if [[ ! -s $combined ]]; then
        warn "[${tag}] Confab produced no conformers"
        return
    fi

    local count
    count=$(grep -c '^\$\$\$\$' "$combined" 2>/dev/null || echo 0)
    log "[${tag}] Confab generated ${count} conformer(s) -> ${combined} (run Stage 3 to split)"
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 2: Conformer Generation (Confab)
-------------------------------------------------------------
 Energy cutoff : ${ecut_kcal} kcal/mol
 Max conformers: ${conf_count}
 Dry run       : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"
    make_qxyz_py
    print_banner

    local entries=()
    if [[ -n $list_file ]]; then
        mapfile -t entries < <(grep -vE '^[[:space:]]*(#|;|$)' "$list_file")
    else
        entries=("${positional[@]}")
    fi

    for entry in "${entries[@]}"; do
        if [[ $entry == *.xyz ]]; then
            require_file "$entry"
            process_xyz "$entry"
        else
            process_tag "$entry"
        fi
    done

    rm -f "$QXYZ_PY"
    log "Stage 2 complete."
}

main "$@"
