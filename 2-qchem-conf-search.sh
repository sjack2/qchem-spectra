#!/usr/bin/env bash
# ============================================================================
# 2-qchem-conf-search.sh — Conformer ensemble generation (Stage 2)
# ============================================================================
#
# OVERVIEW
#   Extracts the optimised geometry from the Stage 1 Q-Chem output, converts
#   it to SDF, and runs Open Babel's Confab to generate a conformer ensemble.
#
# Usage:
#   2-qchem-conf-search.sh TAG
#   2-qchem-conf-search.sh --list molecules.txt --conf 500 --ecut 3
#
# Flags:
#   --conf N       Max conformers to generate             [1000]
#   --ecut KCAL    Energy cutoff above minimum (kcal/mol)  [5]
#   --list FILE    File of TAGs or XYZ paths
#   --dry-run      Echo commands without running
#   -h | --help    Show this help and exit
#
# Directory layout:
#   <TAG>/
#   ├── 01_gas_opt/<TAG>.out         ← Stage 1 output (input here)
#   └── 02_conf_search/
#       ├── <TAG>.xyz                extracted optimised geometry
#       ├── <TAG>.sdf                single-conformer SDF
#       └── <TAG>_combined.sdf       all conformers  → Stage 3 input
#
# Dependencies: obabel (Open Babel >= 3.1), python3
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
DEFAULT_CONF_COUNT=1000
DEFAULT_ECUT_KCAL=5

GAS_OPT_SUBDIR="01_gas_opt"
CONF_SUBDIR="02_conf_search"

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
die()          { printf 'Error: %s\n' "$*" >&2; exit 1; }
log()          { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
warn()         { printf '[%s] Warning: %s\n' "$(date '+%F %T')" "$*" >&2; }
require_file() { [[ -f $1 ]] || die "File '$1' not found"; }

show_help() {
    sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# \{0,1\}//'
    exit 0
}

# ============================================================================
# Q-CHEM OUTPUT → XYZ EXTRACTOR
# Written to a temp file once and reused across all molecules in a batch run.
# ============================================================================
make_qxyz_py() {
    QXYZ_PY=$(mktemp)
    cat >"$QXYZ_PY" <<'PY'
#!/usr/bin/env python3
"""Extract the last optimised geometry from a Q-Chem output file."""
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
    ecut=$DEFAULT_ECUT_KCAL
    dry_run=false
    list_file=""
    positional=()

    local opts
    opts=$(getopt -o h --long help,conf:,ecut:,list:,dry-run \
        -- "$@") || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case "$1" in
            --conf)     conf_count=$2; shift 2 ;;
            --ecut)     ecut=$2;       shift 2 ;;
            --list)     list_file=$2;  shift 2 ;;
            --dry-run)  dry_run=true;  shift ;;
            -h|--help)  show_help ;;
            --)         shift; break ;;
            *)          die "Unknown option '$1'" ;;
        esac
    done
    positional+=("$@")

    if [[ -n $list_file ]]; then
        [[ ${#positional[@]} -eq 0 ]] || die "Positional args not allowed with --list"
        require_file "$list_file"
    else
        [[ ${#positional[@]} -ge 1 ]] || die "Provide at least one TAG or XYZ path"
    fi

    command -v obabel  >/dev/null 2>&1 || die "obabel not in PATH"
    command -v python3 >/dev/null 2>&1 || die "python3 not in PATH"
}

# ============================================================================
# HELPERS
# ============================================================================
find_qchem_out() {
    local tag=$1
    local f="${tag}/${GAS_OPT_SUBDIR}/${tag}.out"
    [[ -f $f ]] && { printf '%s' "$f"; return; }
    printf ''
}

confab_run() {
    local tag=$1 xyz=$2 sdf=$3 combined=$4

    if $dry_run; then
        log "[${tag}] obabel -ixyz ${xyz} -osdf -O ${sdf}"
    else
        obabel -ixyz "$xyz" -osdf -O "$sdf"
    fi
    [[ -s $sdf ]] || { warn "[${tag}] SDF empty after conversion — skipping Confab"; return; }

    local args=(--confab --original -xC --conf "$conf_count" --ecutoff "$ecut" --errorlevel 1)
    if $dry_run; then
        log "[${tag}] obabel ${sdf} -O ${combined} ${args[*]}"
    else
        obabel "$sdf" -O "$combined" "${args[@]}"
    fi
    [[ -s $combined ]] || warn "[${tag}] Confab produced no output"
}

# ============================================================================
# PER-MOLECULE PROCESSING
# ============================================================================
process_tag() {
    local tag=$1
    local qcout
    qcout=$(find_qchem_out "$tag")
    if [[ -z $qcout ]]; then
        warn "[${tag}] Stage 1 output not found (${tag}/${GAS_OPT_SUBDIR}/${tag}.out) — skipping"
        return
    fi

    local out_dir="${tag}/${CONF_SUBDIR}"
    mkdir -p "$out_dir"

    local xyz="${out_dir}/${tag}.xyz"
    if ! "$QXYZ_PY" "$qcout" >"$xyz"; then
        warn "[${tag}] geometry extraction failed from ${qcout}"
        return
    fi

    confab_run "$tag" "$xyz" "${out_dir}/${tag}.sdf" "${out_dir}/${tag}_combined.sdf"
    log "[${tag}] conformer search complete → ${out_dir}/${tag}_combined.sdf"
}

process_xyz() {
    local xyz_path=$1
    local stem
    stem=$(basename "${xyz_path%.xyz}")
    local out_dir="${stem}/${CONF_SUBDIR}"
    mkdir -p "$out_dir"
    cp "$xyz_path" "${out_dir}/${stem}.xyz"
    confab_run "$stem" "${out_dir}/${stem}.xyz" \
        "${out_dir}/${stem}.sdf" "${out_dir}/${stem}_combined.sdf"
    log "[${stem}] conformer search complete → ${out_dir}/${stem}_combined.sdf"
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"
    make_qxyz_py

    local entries=()
    if [[ -n $list_file ]]; then
        mapfile -t entries < <(grep -vE '^[[:space:]]*(#|;|$)' "$list_file")
    else
        entries=("${positional[@]}")
    fi

    for entry in "${entries[@]}"; do
        if [[ $entry == *.xyz ]]; then
            process_xyz "$entry"
        else
            process_tag "$entry"
        fi
    done

    rm -f "$QXYZ_PY"
    $dry_run && log "Dry run complete — no files created" || log "Stage 2 complete."
}

main "$@"
