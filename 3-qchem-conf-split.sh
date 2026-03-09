#!/usr/bin/env bash
# ============================================================================
# 3-qchem-conf-split.sh — Split conformer SDF into per-conformer files (Stage 3)
# ============================================================================
#
# OVERVIEW
#   Splits the multi-conformer SDF produced by Stage 2 into individual SDF
#   and XYZ files, one per conformer.
#
# Usage:
#   3-qchem-conf-split.sh TAG
#   3-qchem-conf-split.sh --list molecules.txt --dry-run
#   3-qchem-conf-split.sh path/to/molecule_combined.sdf
#
# Flags:
#   --list FILE    File of TAGs or combined-SDF paths
#   --dry-run      Echo commands without running
#   -h | --help    Show this help and exit
#
# Directory layout:
#   <TAG>/
#   ├── 02_conf_search/<TAG>_combined.sdf    ← Stage 2 output (input here)
#   └── 02_conf_search/
#       ├── split_sdf/<TAG>_N.sdf            one SDF per conformer
#       └── split_xyz/<TAG>_N.xyz            one XYZ per conformer → Stage 4
#
# Dependencies: obabel (Open Babel >= 3.1)
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS
# ============================================================================
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
# CLI PARSER
# ============================================================================
parse_cli() {
    dry_run=false
    list_file=""
    positional=()

    local opts
    opts=$(getopt -o h --long help,list:,dry-run \
        -- "$@") || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case "$1" in
            --list)     list_file=$2; shift 2 ;;
            --dry-run)  dry_run=true; shift ;;
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
        [[ ${#positional[@]} -ge 1 ]] || die "Provide at least one TAG or SDF path"
    fi

    command -v obabel >/dev/null 2>&1 || die "obabel not in PATH"
}

# ============================================================================
# PATH HELPERS
# ============================================================================
find_combined_sdf() {
    local tag=$1
    local f="${tag}/${CONF_SUBDIR}/${tag}_combined.sdf"
    [[ -f $f ]] && printf '%s' "$f"
}

tag_from_sdf() { basename "$1" | sed -E 's/_combined\.sdf$//'; }

# ============================================================================
# CORE PROCESSING
# ============================================================================
split_and_convert() {
    local tag=$1 sdf=$2

    local base_dir="${tag}/${CONF_SUBDIR}"
    local split_sdf="${base_dir}/split_sdf"
    local split_xyz="${base_dir}/split_xyz"
    mkdir -p "$split_sdf" "$split_xyz"

    log "[${tag}] splitting ${sdf} → individual SDF files"
    if $dry_run; then
        log "[${tag}] obabel ${sdf} -O ${split_sdf}/${tag}_.sdf -m"
    else
        obabel "$sdf" -O "${split_sdf}/${tag}_".sdf -m --errorlevel 1
    fi

    local count
    count=$(printf '%s\n' "${split_sdf}/${tag}_"*.sdf 2>/dev/null | wc -l)
    [[ $count -eq 0 ]] && { warn "[${tag}] no split SDF files found"; return; }

    log "[${tag}] converting ${count} SDF → XYZ"
    for s in "${split_sdf}/${tag}_"*.sdf; do
        local stem
        stem=$(basename "${s%.sdf}")
        if $dry_run; then
            log "[${tag}] obabel ${s} -O ${split_xyz}/${stem}.xyz"
        else
            obabel "$s" -O "${split_xyz}/${stem}.xyz" --errorlevel 1
        fi
    done

    log "[${tag}] split complete → ${split_xyz}/"
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"

    local entries=()
    if [[ -n $list_file ]]; then
        mapfile -t entries < <(grep -vE '^[[:space:]]*(#|;|$)' "$list_file")
    else
        entries=("${positional[@]}")
    fi

    for entry in "${entries[@]}"; do
        if [[ $entry == *.sdf ]]; then
            require_file "$entry"
            local tag
            tag=$(tag_from_sdf "$entry")
            split_and_convert "$tag" "$entry"
        else
            local tag=$entry
            local combined
            combined=$(find_combined_sdf "$tag" || true)
            if [[ -z $combined ]]; then
                warn "[${tag}] combined SDF not found — skipping"
                continue
            fi
            split_and_convert "$tag" "$combined"
        fi
    done

    $dry_run && log "Dry run complete — no files created" || log "Stage 3 complete."
}

main "$@"
