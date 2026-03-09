#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 5-qchem-boltzmann-weight.sh - Read final SCF energies from solvent-phase
# Q-Chem outputs (step 4), compute Boltzmann populations, and write:
#
#   <tag>/bw_results/<tag>_energies.dat   # sorted table (CID, E, dE, p)
#   <tag>/bw_results/<tag>_bw_labels.dat  # CIDs with p >= cutoff
#
# CLI superset:  -m/--method  -b/--basis  --disp  --grid  --cpus  --mem-per-cpu
# Step flags:    --temp K   --p-cut P   --list FILE   --dry-run
# -----------------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

# --------------------------- constants ---------------------------------------
H2KCAL=627.509474          # Hartree -> kcal mol-1
R_J=8.314462618            # J mol-1 K-1
DEFAULT_T=298.15           # K
DEFAULT_P_CUT=0.01         # probability cutoff

# --------------------------- helpers -----------------------------------------
die()  { printf 'Error: %s\n' "$*" >&2; exit 1; }
note() { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

show_help() {
cat <<EOF
Usage:
  $(basename "$0") MoleculeTag [OPTIONS]
  $(basename "$0") --list tags.txt [OPTIONS]

Common flags (accepted for parity, ignored here):
  -m/--method, -b/--basis, --disp, --grid, --cpus, --mem-per-cpu

Step-specific:
  --temp K        Temperature in Kelvin     [$DEFAULT_T]
  --p-cut P       Probability cutoff        [$DEFAULT_P_CUT]
  --list FILE     Molecule tags (blank/#/; lines ignored)
  --dry-run       Echo actions only
  -h, --help      Show this help
EOF
}

# ----------------------------- CLI -------------------------------------------
parse_cli() {
  T=$DEFAULT_T; P_CUT=$DEFAULT_P_CUT; DRY=false; LIST=""; SINGLE=""

  local opts
  opts=$(getopt -o hb:m: --long help,basis:,method:,disp:,grid:,cpus:,mem-per-cpu:,temp:,p-cut:,list:,dry-run -- "$@")
  eval set -- "$opts"
  while true; do
    case "$1" in
      --temp)          T=$2; shift 2;;
      --p-cut)         P_CUT=$2; shift 2;;
      --list)          LIST=$2; shift 2;;
      --dry-run)       DRY=true; shift;;
      -h|--help)       show_help; exit 0;;
      # shared (ignored) flags
      -m|--method|-b|--basis|--disp|--grid|--cpus|--mem-per-cpu)
                       shift 2;;
      --) shift; break;;
      *) die "Unknown option $1";;
    esac
  done

  if [[ -n $LIST ]]; then
    [[ $# -eq 0 ]] || die "Positional args not allowed with --list"
    [[ -f $LIST ]] || die "List file '$LIST' not found"
  else
    [[ $# -eq 1 ]] || die "Provide one molecule tag"
    SINGLE=$1
  fi

  command -v python3 >/dev/null 2>&1 || die "python3 not in PATH"
}

# --------------------------- Python parser -----------------------------------
make_parser() {
  PARSER=$(mktemp)
  cat >"$PARSER" <<'PY'
#!/usr/bin/env python3
import sys, re, pathlib
pat = re.compile(r'Final energy is\s+(-?\d+\.\d+)')
for f in sys.argv[1:]:
    txt = pathlib.Path(f).read_text(errors='ignore')
    m = pat.search(txt)
    if m:
        print(f"{f}\t{m.group(1)}")
PY
  chmod +x "$PARSER"
}

# ------------------------- per-molecule routine ------------------------------
process_tag() {
  local TAG=$1
  local OPT_DIR="${TAG}/solvent_opt"
  [[ -d $OPT_DIR ]] || { note "[$TAG] solvent_opt dir missing"; return; }

  mapfile -t OUTS < <(find "$OPT_DIR" -name '*.out' | sort)
  [[ ${#OUTS[@]} -gt 0 ]] || { note "[$TAG] no .out files"; return; }

  local RAW; RAW=$(python3 "$PARSER" "${OUTS[@]}") || true
  [[ -n $RAW ]] || { note "[$TAG] energies not found"; return; }

  python3 - <<PY
import math, pathlib
H2KCAL = $H2KCAL
R_J    = $R_J
T      = float("$T")
p_cut  = float("$P_CUT")
tag    = "$TAG"

rows = []
for line in """$RAW""".splitlines():
    fn, e = line.split('\t')
    cid   = pathlib.Path(fn).stem
    rows.append((cid, float(e)))

e_min = min(e for _, e in rows)
rows = [(cid,
         e,
         (e - e_min) * H2KCAL,
         math.exp(-(e - e_min)*H2KCAL*4184/(R_J*T)))
        for cid, e in rows]

Z = sum(r[3] for r in rows)
rows = [(cid, e, dE, b/Z) for cid, e, dE, b in rows]
rows.sort(key=lambda r: r[3], reverse=True)

out_dir = pathlib.Path(tag) / "bw_results"
out_dir.mkdir(parents=True, exist_ok=True)
dat = out_dir / f"{tag}_energies.dat"
lab = out_dir / f"{tag}_bw_labels.dat"

if "$DRY" == "True":
    print(f"[{tag}] would create {dat.name} and label list")
else:
    with dat.open('w') as f:
        f.write("#CID  E(Ha)        dE(kcal)     p\n")
        for cid, e, dE, p in rows:
            f.write(f"{cid}  {e:12.6f}  {dE:10.3f}  {p:10.5f}\n")
    with lab.open('w') as f:
        for cid, _, _, p in rows:
            if p >= p_cut:
                f.write(cid + "\n")
    kept = sum(1 for _ in lab.open())
    print(f"[{tag}] wrote {dat.name} ({kept} kept)")
PY
}

# -------------------------------- main ---------------------------------------
main() {
  parse_cli "$@"; make_parser

  if [[ -n $LIST ]]; then
    while IFS= read -r line || [[ -n $line ]]; do
      [[ $line =~ ^[[:space:]]*$ ]] && continue
      [[ $line =~ '^[[:space:]]*[#;]' ]] && continue
      process_tag "$line"
    done <"$LIST"
  else
    process_tag "$SINGLE"
  fi

  rm -f "$PARSER"
}

main "$@"
