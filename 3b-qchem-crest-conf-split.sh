#!/usr/bin/env bash
# ============================================================================
# 3b-qchem-crest-conf-split.sh -- Split CREST output into per-conformer XYZ (Stage 3b)
# ============================================================================
#
# OVERVIEW
#   Reads crest_conformers.xyz produced by Stage 2b (CREST) and splits it
#   into individual numbered XYZ files. Also copies the reference geometry
#   to pre_xyz/ if not already present (needed by later stages for
#   charge/multiplicity parsing).
#
# Usage:
#   3b-qchem-crest-conf-split.sh TAG
#   3b-qchem-crest-conf-split.sh --list mols.txt --dry-run
#
# Flags:
#   --list FILE    Text file of TAGs (one per line)
#   --dry-run      Echo actions without creating files
#   -h | --help    Show this help and exit
#
# Directory layout:
#   <TAG>/02_conf_search/
#   |-- crest_conformers.xyz       <- input (from Stage 2b)
#   -- split_xyz/
#       |-- <TAG>_001.xyz          individual conformers (zero-padded)
#       |-- <TAG>_002.xyz          -> Stage 4 input
#       -- ...
#
# Examples:
#   3b-qchem-crest-conf-split.sh ephedrine
#   3b-qchem-crest-conf-split.sh --list molecules.txt
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
    tags=()

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
    tags+=("$@")

    if [[ -n $list_file ]]; then
        mapfile -t tags < <(grep -vE '^[[:space:]]*(#|;|$)' "$list_file")
    fi
    (( ${#tags[@]} > 0 )) || die "No TAGs supplied"
}

# ============================================================================
# SPLITTER
# ============================================================================
split_crest_xyz() {
    local tag=$1
    local src="${tag}/02_conf_search/crest_conformers.xyz"
    local outdir="${tag}/02_conf_search/split_xyz"

    if [[ ! -f $src ]]; then
        warn "[${tag}] crest_conformers.xyz not found -- skipping"
        return
    fi

    mkdir -p "$outdir" pre_xyz

    # copy reference geometry to pre_xyz/ if absent (needed for charge/mult)
    if [[ ! -f "pre_xyz/${tag}.xyz" && -f "${tag}/02_conf_search/${tag}.xyz" ]]; then
        cp -f "${tag}/02_conf_search/${tag}.xyz" "pre_xyz/"
    fi

    if $dry_run; then
        log "[${tag}] (dry run) would split ${src} into ${outdir}/${tag}_NNN.xyz"
        return
    fi

    # awk-based splitter for multi-structure XYZ format:
    #   <atom_count>
    #   <comment line>
    #   <coord line 1>
    #   ...
    #   <atom_count>       <- next structure starts
    #   ...
    awk -v tag="$tag" -v outdir="$outdir" '
    BEGIN { idx = 0 }
    {
        if (newblock) {
            comment = $0
            idx++
            cid = sprintf("%s_%03d", tag, idx)
            fname = outdir "/" cid ".xyz"
            print atoms > fname
            print comment >> fname
            for (i = 1; i <= atoms; i++) {
                getline line
                split(line, a)
                printf "%s %s %s %s\n", a[1], a[2], a[3], a[4] >> fname
            }
            newblock = 0
        } else {
            atoms = $1
            newblock = 1
        }
    }' "$src"

    shopt -s nullglob
    local count
    count=$(find "$outdir" -maxdepth 1 -name "${tag}_*.xyz" -type f | wc -l)
    shopt -u nullglob

    log "[${tag}] wrote ${count} conformers to ${outdir}/"
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    parse_cli "$@"

    for tag in "${tags[@]}"; do
        split_crest_xyz "$tag"
    done

    log "Stage 3b complete."
}

main "$@"
