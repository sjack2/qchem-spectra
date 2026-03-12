#!/usr/bin/env bash
# ============================================================================
# 3-qchem-conf-split.sh -- Split Confab output into per-conformer XYZ (Stage 3)
# ============================================================================
#
# OVERVIEW
#   Reads the combined SDF produced by Stage 2 (Confab), splits it into
#   individual SDF files, and converts each to XYZ format. The resulting
#   XYZ files in split_xyz/ are the input for Stage 4 (solvent optimization).
#
# Usage:
#   3-qchem-conf-split.sh TAG
#   3-qchem-conf-split.sh --list molecules.txt --dry-run
#
# Flags:
#   --list FILE    Text file of TAGs or SDF paths
#   --dry-run      Echo actions without executing
#   -h | --help    Show this help and exit
#
# Directory layout:
#   <TAG>/02_conf_search/
#   |-- <TAG>_combined.sdf     <- input (from Stage 2)
#   |-- split_sdf/             individual SDF files
#   -- split_xyz/             individual XYZ files -> Stage 4 input
#
# Examples:
#   3-qchem-conf-split.sh ephedrine
#   3-qchem-conf-split.sh --list molecules.txt
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

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
# CLI PARSER
# ============================================================================
parse_cli() {
    dry_run=false
    list_file=""
    positional=()

    local opts
    opts=$(getopt -o h --long help,list:,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            --list)    list_file=$2; shift 2 ;;
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
        [[ ${#positional[@]} -gt 0 ]] || die "Provide at least one TAG or SDF path"
    fi

    command -v obabel >/dev/null 2>&1 || die "obabel not found in PATH"
}

# ============================================================================
# CORE PROCESSING
# ============================================================================
find_combined_sdf() {
    local tag=$1
    local probe="${tag}/02_conf_search/${tag}_combined.sdf"
    [[ -f $probe ]] && printf '%s' "$probe" || printf ''
}

split_and_convert() {
    local tag=$1 sdf_file=$2
    local base_dir="${tag}/02_conf_search"
    local split_sdf="${base_dir}/split_sdf"
    local split_xyz="${base_dir}/split_xyz"
    mkdir -p "$split_sdf" "$split_xyz" pre_xyz

    # copy reference geometry to pre_xyz/ if absent (needed for charge/mult)
    if [[ ! -f "pre_xyz/${tag}.xyz" && -f "${base_dir}/${tag}.xyz" ]]; then
        cp -f "${base_dir}/${tag}.xyz" "pre_xyz/"
        log "[${tag}] copied reference geometry to pre_xyz/"
    fi

    if $dry_run; then
        log "[${tag}] (dry run) would split ${sdf_file} into ${split_xyz}/"
        return
    fi

    log "[${tag}] splitting ${sdf_file}"
    obabel "$sdf_file" -O "${split_sdf}/${tag}_".sdf -m 2>/dev/null

    shopt -s nullglob
    local sdf_files=("${split_sdf}/${tag}_"*.sdf)
    shopt -u nullglob

    local count=${#sdf_files[@]}
    if (( count == 0 )); then
        warn "[${tag}] no split SDF files produced"
        return
    fi

    log "[${tag}] converting ${count} conformers to XYZ"
    for s in "${sdf_files[@]}"; do
        local base
        base=$(basename "${s%.sdf}")
        obabel "$s" -O "${split_xyz}/${base}.xyz" 2>/dev/null
    done

    log "[${tag}] ${count} conformers ready in ${split_xyz}/"
}

# ============================================================================
# MAIN
# ============================================================================
main() {
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
            tag=$(basename "$entry" | sed -E 's/_combined\.sdf$//')
            split_and_convert "$tag" "$entry"
        else
            local tag=$entry
            local sdf_path
            sdf_path=$(find_combined_sdf "$tag")
            if [[ -z $sdf_path ]]; then
                warn "[${tag}] combined SDF not found -- skipping"
                continue
            fi
            split_and_convert "$tag" "$sdf_path"
        fi
    done

    log "Stage 3 complete."
}

main "$@"
