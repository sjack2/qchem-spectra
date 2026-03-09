#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 3-qchem-conf-split.sh – Split multi-conformer SDF files (from step 2) into
# individual SDF and XYZ files.
#
# INPUT  (default path for tag “caffeine”)
#   caffeine/qchem_opt_conf/caffeine_combined.sdf
#
# OUTPUT
#   caffeine/qchem_opt_conf/split_sdf/caffeine_N.sdf
#   caffeine/qchem_opt_conf/split_xyz/caffeine_N.xyz
#
# COMMON FLAGS (accepted for uniform CLI but ignored here)
#   -m/--method  -b/--basis  --disp  --grid  --cpus  --mem-per-cpu
#
# STEP-SPECIFIC FLAGS
#   --list FILE   List of tags or SDF paths
#   --dry-run     Echo commands, do not run Open Babel
#   -h, --help    Show this help and exit
#
# DEPENDENCIES
#   Bash ≥ 4, Open Babel ≥ 3.1  (obabel must be in PATH)
# -----------------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

# --------------------------- utilities ---------------------------------------
help()  { grep -E '^(#|$)' "$0" | head -n 120 | sed 's/^# \{0,1\}//'; }
die()   { printf 'Error: %s\n' "$*" >&2; exit 1; }
note()  { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# --------------------------- CLI parsing -------------------------------------
parse_cli() {
  DRY_RUN=false
  LIST_FILE=""
  POSITIONAL=()

  local opts
  opts=$(getopt -o hb:m: --long help,basis:,method:,disp:,grid:,cpus:,mem-per-cpu:,list:,dry-run -- "$@")
  eval set -- "$opts"
  while true; do
    case "$1" in
      --list)        LIST_FILE=$2; shift 2 ;;
      --dry-run)     DRY_RUN=true; shift ;;
      -h|--help)     help; exit 0 ;;
      # shared (ignored) flags ------------------------------------------------
      -m|--method|-b|--basis|--disp|--grid|--cpus|--mem-per-cpu)
                     shift 2 ;;
      --) shift; break ;;
      *) die "Unknown option $1" ;;
    esac
  done
  POSITIONAL+=("$@")

  if [[ -n $LIST_FILE ]]; then
    [[ ${#POSITIONAL[@]} -eq 0 ]] || die 'Positional args not allowed with --list'
    [[ -f $LIST_FILE ]] || die "List file '$LIST_FILE' not found"
  else
    [[ ${#POSITIONAL[@]} -ge 1 ]] || die 'Provide at least one tag or SDF path'
  fi

  command -v obabel >/dev/null 2>&1 || die 'obabel not in PATH'
}

# ------------------------- path helpers --------------------------------------
find_combined_sdf() {          # arg1 = tag
  local tag=$1
  local file="${tag}/qchem_opt_conf/${tag}_combined.sdf"
  [[ -f $file ]] && printf '%s' "$file"
}

tag_from_sdf() { basename "$1" | sed -E 's/_combined\.sdf$//' ; }

# ---------------------- core processing --------------------------------------
split_and_convert() {
  local tag=$1 sdf=$2

  local base_dir="${tag}/qchem_opt_conf"
  local split_sdf="${base_dir}/split_sdf"
  local split_xyz="${base_dir}/split_xyz"
  mkdir -p "$split_sdf" "$split_xyz"

  note "[$tag] Splitting SDF → individual files"
  if $DRY_RUN; then
    note "[$tag] obabel $sdf -O ${split_sdf}/${tag}_ .sdf -m"
  else
    obabel "$sdf" -O "${split_sdf}/${tag}_".sdf -m --errorlevel 1
  fi

  local count
  count=$(printf '%s\n' "${split_sdf}/${tag}_"*.sdf 2>/dev/null | wc -l)
  [[ $count -eq 0 ]] && { note "[$tag] No split SDF files found"; return; }

  note "[$tag] Converting $count SDF → XYZ"
  for s in "${split_sdf}/${tag}_"*.sdf; do
    local stem; stem=$(basename "${s%.sdf}")
    if $DRY_RUN; then
      note "[$tag] obabel $s -O ${split_xyz}/${stem}.xyz"
    else
      obabel "$s" -O "${split_xyz}/${stem}.xyz" --errorlevel 1
    fi
  done
}

# -------------------------------- main ---------------------------------------
main() {
  parse_cli "$@"

  local entries=()
  if [[ -n $LIST_FILE ]]; then
    mapfile -t entries < <(grep -vE '^[[:space:]]*(#|;|$)' "$LIST_FILE")
  else
    entries=("${POSITIONAL[@]}")
  fi

  for entry in "${entries[@]}"; do
    if [[ $entry == *.sdf ]]; then                     # explicit SDF path
      [[ -f $entry ]] || { note "SDF '$entry' not found – skipping"; continue; }
      local tag; tag=$(tag_from_sdf "$entry")
      split_and_convert "$tag" "$entry"
    else                                               # treat as molecule tag
      local tag=$entry
      local combined; combined=$(find_combined_sdf "$tag" || true)
      if [[ -z $combined ]]; then
        note "Combined SDF not found for tag '$tag' – skipping"
        continue
      fi
      split_and_convert "$tag" "$combined"
    fi
  done

  $DRY_RUN && note "Dry run complete – no files created"
  note "All molecules processed"
}

main "$@"
