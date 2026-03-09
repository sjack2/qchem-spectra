#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 2-qchem-conf-search.sh  -  Generate conformer ensembles from step-1 Q-Chem
# gas-phase optimisations or from user-supplied XYZ files.
# -----------------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

DEFAULT_CONF_COUNT=1000
DEFAULT_ECUT_KCAL=5

# --------------------------- utilities ---------------------------------------
show_help() { grep -E '^(#|$)' "$0" | head -n 120 | sed 's/^# \{0,1\}//'; }
die()       { printf 'Error: %s\n' "$*" >&2; exit 1; }
note()      { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# --------------------------- CLI parsing -------------------------------------
parse_cli() {
  CONF_COUNT=$DEFAULT_CONF_COUNT
  ECUT=$DEFAULT_ECUT_KCAL
  DRY_RUN=false
  LIST_FILE=""
  POSITIONAL=()

  local opts
  opts=$(getopt -o hb:m: --long help,basis:,method:,disp:,grid:,cpus:,mem-per-cpu:,list:,conf:,ecut:,dry-run -- "$@")
  eval set -- "$opts"
  while true; do
    case "$1" in
      --list)        LIST_FILE=$2; shift 2 ;;
      --conf)        CONF_COUNT=$2; shift 2 ;;
      --ecut)        ECUT=$2; shift 2 ;;
      --dry-run)     DRY_RUN=true; shift ;;
      -h|--help)     show_help; exit 0 ;;
      -m|--method|-b|--basis|--disp|--grid|--cpus|--mem-per-cpu)
                     shift 2 ;;  # parsed but ignored
      --) shift; break ;;
      *) die "Unknown option $1" ;;
    esac
  done
  POSITIONAL+=("$@")

  if [[ -n $LIST_FILE ]]; then
    [[ ${#POSITIONAL[@]} -eq 0 ]] || die "Positional args not allowed with --list"
    [[ -f $LIST_FILE ]] || die "List file '$LIST_FILE' not found"
  else
    [[ ${#POSITIONAL[@]} -ge 1 ]] || die "Provide at least one tag or XYZ when not using --list"
  fi

  command -v obabel  >/dev/null 2>&1 || die "obabel not in PATH"
  command -v python3 >/dev/null 2>&1 || die "python3 not in PATH"
}

# -------------- robust Q-Chem → XYZ extractor (backwards scanning) -----------
make_qxyz_py() {
  QXYZ_PY=$(mktemp)
  cat >"$QXYZ_PY" <<'PY'
#!/usr/bin/env python3
import sys, re, pathlib

if len(sys.argv) != 2:
    sys.exit("Usage: qxyz.py file.out")

lines = pathlib.Path(sys.argv[1]).read_text(errors="ignore").splitlines()
hdr = re.compile(r'(Standard Nuclear Orientation|Input orientation|Coordinates \(Angstroms\))', re.I)

start = None
for idx in range(len(lines) - 1, -1, -1):
    if hdr.search(lines[idx]):
        start = idx + 3           # skip header + dashed line + column titles
        break
if start is None:
    sys.exit("No coordinate block found")

geom = []
for line in lines[start:]:
    if not line.strip() or line.lstrip().startswith('-'):
        break
    parts = line.split()
    if len(parts) >= 5:                 # '1  O  x  y  z'
        symbol, x, y, z = parts[1], parts[2], parts[3], parts[4]
    elif len(parts) == 4:               # 'O  x  y  z'
        symbol, x, y, z = parts
    else:
        continue
    geom.append(f"{symbol} {x} {y} {z}")

if not geom:
    sys.exit("Coordinate parsing failed")

print(len(geom))
print()
print("\n".join(geom))
PY
  chmod +x "$QXYZ_PY"
}

# ------------------------- helpers ------------------------------------------
find_qchem_out() {
  local tag=$1
  local new="${tag}/${tag}_qchem_opt/${tag}.out"
  local old="${tag}/init_opt/output/${tag}_qchem_init_opt.out"
  [[ -f $new ]] && { printf '%s' "$new"; return; }
  [[ -f $old ]] && { printf '%s' "$old"; return; }
  printf ''
}

confab_run() {
  local tag=$1 xyz=$2 sdf=$3 combined=$4
  if $DRY_RUN; then
    note "[$tag] obabel -ixyz $xyz -osdf -O $sdf"
  else
    obabel -ixyz "$xyz" -osdf -O "$sdf"
  fi
  [[ -s $sdf ]] || { note "[$tag] SDF empty, skipping"; return; }

  local args=( --confab --original -xC --conf "$CONF_COUNT" --ecutoff "$ECUT" --errorlevel 1 )
  if $DRY_RUN; then
    note "[$tag] obabel $sdf -O $combined ${args[*]}"
  else
    obabel "$sdf" -O "$combined" "${args[@]}"
  fi
  [[ -s $combined ]] || note "[$tag] Confab produced no output"
}

process_tag() {
  local tag=$1
  local qcout; qcout=$(find_qchem_out "$tag")
  if [[ -z $qcout ]]; then
    note "[$tag] Q-Chem output not found, skipping"
    return
  fi
  local out_dir="${tag}/qchem_opt_conf"; mkdir -p "$out_dir"
  local xyz="$out_dir/${tag}.xyz"
  if ! "$QXYZ_PY" "$qcout" >"$xyz"; then
    note "[$tag] Coordinate extraction failed"
    return
  fi
  confab_run "$tag" "$xyz" "$out_dir/${tag}.sdf" "$out_dir/${tag}_combined.sdf"
}

process_xyz() {
  local xyz_path=$1
  local stem; stem=$(basename "${xyz_path%.xyz}")
  local out_dir="${stem}/qchem_opt_conf"; mkdir -p "$out_dir"
  cp "$xyz_path" "$out_dir/${stem}.xyz"
  confab_run "$stem" "$out_dir/${stem}.xyz" "$out_dir/${stem}.sdf" "$out_dir/${stem}_combined.sdf"
}

# -------------------------------- main ---------------------------------------
main() {
  parse_cli "$@"; make_qxyz_py

  local entries=()
  if [[ -n $LIST_FILE ]]; then
    mapfile -t entries < <(grep -vE '^[[:space:]]*(#|;|$)' "$LIST_FILE")
  else
    entries=("${POSITIONAL[@]}")
  fi

  for entry in "${entries[@]}"; do
    if [[ $entry == *.xyz ]]; then
      process_xyz "$entry"
    else
      process_tag "$entry"
    fi
  done

  $DRY_RUN && note "Dry run complete - no files created"
  rm -f "$QXYZ_PY"
  note "All molecules processed"
}

main "$@"
