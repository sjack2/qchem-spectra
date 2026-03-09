#!/usr/bin/env bash
# ============================================================================
# 5-qchem-boltzmann-weight.sh — Boltzmann weighting & filtering (Stage 5)
# ============================================================================
#
# OVERVIEW
#   Collects Q-Chem final SCF energies from the solvent-phase optimisation
#   outputs (Stage 4), computes relative energies and Boltzmann populations
#   at a given temperature, and writes two files:
#
#     <TAG>/04_boltzmann/<TAG>_energies.dat    full table sorted by probability
#     <TAG>/04_boltzmann/<TAG>_bw_labels.dat   conformer IDs above the cutoff
#
#   The label file is read by Stage 6 to determine which conformers receive
#   excited-state or frequency calculations.
#
# Usage:
#   5-qchem-boltzmann-weight.sh TAG
#   5-qchem-boltzmann-weight.sh --temp 310 --p-cut 0.02 aspirin
#   5-qchem-boltzmann-weight.sh --list mols.txt --dry-run
#
# Flags:
#   --temp K       Temperature in Kelvin               [298.15]
#   --p-cut VAL    Probability cutoff (e.g. 0.01 = 1%) [0.01]
#   --list FILE    Text file of molecule TAGs
#   --dry-run      Show what would be computed without writing
#   -h | --help    Show this help and exit
#
# Directory layout:
#   <TAG>/
#   ├── 03_solvent_opt/
#   │   ├── <TAG>_1/<TAG>_1.out      ← Stage 4 outputs (input here)
#   │   └── <TAG>_2/<TAG>_2.out
#   └── 04_boltzmann/
#       ├── <TAG>_energies.dat       full table (CID  E  dE  p)
#       └── <TAG>_bw_labels.dat      filtered labels → Stage 6
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS & CONSTANTS
# ============================================================================
H2KCAL=627.509474       # Hartree → kcal/mol
R_J=8.314462618         # J mol⁻¹ K⁻¹
DEFAULT_TEMP=298.15
DEFAULT_P_CUT=0.01

SOLV_OPT_SUBDIR="03_solvent_opt"
BW_SUBDIR="04_boltzmann"

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
    temp=$DEFAULT_TEMP
    p_cut=$DEFAULT_P_CUT
    dry_run=false
    list_file=""
    single=""

    local opts
    opts=$(getopt -o h --long help,temp:,p-cut:,list:,dry-run \
        -- "$@") || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case "$1" in
            --temp)     temp=$2;       shift 2 ;;
            --p-cut)    p_cut=$2;      shift 2 ;;
            --list)     list_file=$2;  shift 2 ;;
            --dry-run)  dry_run=true;  shift ;;
            -h|--help)  show_help ;;
            --)         shift; break ;;
            *)          die "Unknown option '$1'" ;;
        esac
    done

    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional arguments not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide one molecule TAG"
        single=$1
    fi

    command -v python3 >/dev/null 2>&1 || die "python3 not in PATH"
}

# ============================================================================
# ENERGY EXTRACTOR — finds "Final energy is" lines in Q-Chem output files
# ============================================================================
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

# ============================================================================
# PER-MOLECULE BOLTZMANN CALCULATION
# ============================================================================
process_tag() {
    local tag=$1
    local opt_dir="${tag}/${SOLV_OPT_SUBDIR}"
    if [[ ! -d $opt_dir ]]; then
        warn "[${tag}] ${SOLV_OPT_SUBDIR} directory not found — skipping"
        return
    fi

    mapfile -t outs < <(find "$opt_dir" -name '*.out' | sort)
    if [[ ${#outs[@]} -eq 0 ]]; then
        warn "[${tag}] no .out files found in ${opt_dir}"
        return
    fi

    local raw
    raw=$(python3 "$PARSER" "${outs[@]}") || true
    if [[ -z $raw ]]; then
        warn "[${tag}] no energies found (check that Stage 4 jobs completed)"
        return
    fi

    python3 - <<PY
import math, pathlib
H2KCAL = ${H2KCAL}
R_J    = ${R_J}
T      = float("${temp}")
p_cut  = float("${p_cut}")
tag    = "${tag}"
bw_dir = "${BW_SUBDIR}"

rows = []
for line in """${raw}""".splitlines():
    fn, e = line.split('\t')
    cid   = pathlib.Path(fn).stem
    rows.append((cid, float(e)))

e_min = min(e for _, e in rows)
rows = [
    (cid, e,
     (e - e_min) * H2KCAL,
     math.exp(-(e - e_min) * H2KCAL * 4184 / (R_J * T)))
    for cid, e in rows
]

Z = sum(r[3] for r in rows)
rows = [(cid, e, dE, b / Z) for cid, e, dE, b in rows]
rows.sort(key=lambda r: r[3], reverse=True)

out_dir = pathlib.Path(tag) / bw_dir
out_dir.mkdir(parents=True, exist_ok=True)
dat = out_dir / f"{tag}_energies.dat"
lab = out_dir / f"{tag}_bw_labels.dat"

if "${dry_run}" == "true":
    print(f"[{tag}] would write {dat.name} and {lab.name} ({len(rows)} conformers)")
else:
    with dat.open('w') as f:
        f.write("#CID  E(Ha)        dE(kcal)     p\n")
        for cid, e, dE, p in rows:
            f.write(f"{cid}  {e:12.6f}  {dE:10.3f}  {p:10.5f}\n")
    kept = [cid for cid, _, _, p in rows if p >= p_cut]
    with lab.open('w') as f:
        for cid in kept:
            f.write(cid + "\n")
    print(f"[{tag}] wrote {dat.name} — {len(kept)}/{len(rows)} conformers kept (p >= {p_cut})")
PY
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"
    make_parser

    if [[ -n $list_file ]]; then
        while IFS= read -r entry || [[ -n $entry ]]; do
            [[ -z $entry || $entry == \#* || $entry == \;* ]] && continue
            process_tag "$entry"
        done < "$list_file"
    else
        process_tag "$single"
    fi

    rm -f "$PARSER"
    $dry_run && log "Dry run complete." || log "Stage 5 complete."
}

main "$@"
