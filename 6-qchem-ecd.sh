#!/usr/bin/env bash
# 6-qchem-ecd.sh – EOM-EE-CCSD ECD/UV-Vis for Boltzmann-kept conformers
# -----------------------------------------------------------------------------
set -euo pipefail
IFS=$'\n\t'

# ---------------- defaults ---------------------------------------------------
DEFAULT_METHOD="ccsd"
DEFAULT_BASIS="aug-cc-pVDZ"
DEFAULT_ROOTS=10
DEFAULT_SOLVENT="water"
DEFAULT_CPUS=12
DEFAULT_MEM_MB=3800
DEFAULT_PART="circe"
DEFAULT_WALL="08:00:00"
DEFAULT_MAX_SCF=200

SOLV_OUT_DIR="solvent_opt"          # <tag>/solvent_opt/<CID>/<CID>.out
PRE_DIR="pre_xyz"                   # charge / multiplicity source

# ---------------- helper funcs ----------------------------------------------
die(){ printf 'Error: %s\n' "$*" >&2; exit 1; }
note(){ printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
pre_xyz(){ printf '%s/%s.xyz' "$PRE_DIR" "$1"; }
charge_mult(){ awk 'NR==2{print $1,$2}' "$1"; }

# ---- on-the-fly extractor (same logic as script 2, but returns to stdout) ---
qcout2xyz() {                         # $1 = path to .out
  python3 - "$1" <<'PY'
import pathlib, re, sys, textwrap; f=sys.argv[1]
lines=pathlib.Path(f).read_text(errors="ignore").splitlines()
hdr=re.compile(r'(Standard Nuclear Orientation|Input orientation|Coordinates \(Angstroms\))',re.I)
start=None
for i in range(len(lines)-1,-1,-1):
    if hdr.search(lines[i]): start=i+3; break
if start is None: sys.exit(1)
geom=[]
for l in lines[start:]:
    if not l.strip() or l.lstrip().startswith('-'): break
    p=l.split()
    if len(p)>=5: sym,x,y,z=p[1:5]
    elif len(p)==4: sym,x,y,z=p
    else: continue
    geom.append(f"{sym} {x} {y} {z}")
print(len(geom));print();print("\n".join(geom))
PY
}

# ---------------- CLI parsing ------------------------------------------------
parse_cli(){
  METHOD=$DEFAULT_METHOD; BASIS=$DEFAULT_BASIS; ROOTS=$DEFAULT_ROOTS
  SOLVENT=$DEFAULT_SOLVENT; CPUS=$DEFAULT_CPUS; MEM_MB=$DEFAULT_MEM_MB
  PART=$DEFAULT_PART; MAX_SCF=$DEFAULT_MAX_SCF
  LIST=""; SINGLE=""; DRY=false

  local opts
  opts=$(getopt -o h --long help,method:,basis:,roots:,solvent:,cpus:,mem-per-cpu:,partition:,max-scf-cycles:,list:,dry-run -- "$@")
  eval set -- "$opts"
  while true; do
    case "$1" in
      --method) METHOD=$2; shift 2;;
      --basis)  BASIS=$2; shift 2;;
      --roots)  ROOTS=$2; shift 2;;
      --solvent) SOLVENT=$2; shift 2;;
      --cpus)   CPUS=$2; shift 2;;
      --mem-per-cpu) MEM_MB=$2; shift 2;;
      --partition) PART=$2; shift 2;;
      --max-scf-cycles) MAX_SCF=$2; shift 2;;
      --list) LIST=$2; shift 2;;
      --dry-run) DRY=true; shift;;
      -h|--help)
        sed -n '2,40p' "$0"; exit 0;;
      --) shift; break;;
      *) die "Unknown option $1";;
    esac
  done
  if [[ -n $LIST ]]; then
    [[ $# -eq 0 ]] || die "Positional args not allowed with --list"
    [[ -f $LIST ]] || die "List file '$LIST' not found"
  else
    [[ $# -eq 1 ]] || die "Provide a molecule tag"
    SINGLE=$1
  fi
  command -v sbatch >/dev/null 2>&1 || die "sbatch not in PATH"
}

# ---------------- writers ----------------------------------------------------
write_input(){
  local TAG=$1 CID=$2 XYZ=$3 OUT=$4
  local CH=0 ML=1
  local P=$(pre_xyz "$TAG"); [[ -f $P ]] && read -r CH ML < <(charge_mult "$P")
  local MEM_TOTAL=$((CPUS * MEM_MB))

cat >"$OUT" <<EOF
\$comment
EOM-EE-CCSD ECD / UV-Vis (SMD $SOLVENT)
\$end

\$molecule
$CH $ML
$(tail -n +3 "$XYZ")
\$end

\$rem
  JOB_TYPE                 SP
  METHOD                   $METHOD
  BASIS                    $BASIS
  EE_SINGLETS              [$ROOTS]
  EE_TRIPLETS              [0]
  CC_TRANS_PROP            1
  CC_EOM_ECD               TRUE
  EOM_DAVIDSON_CONVERGENCE 5
  GEN_SCFMAN               TRUE
  SCF_ALGORITHM            DIIS_GDM
  SOLVENT_METHOD           SMD
  MAX_SCF_CYCLES           $MAX_SCF
  CC_CONVERGENCE           6
  SCF_CONVERGENCE          8
  THRESH                   12
  SYM_IGNORE               TRUE
  MEM_TOTAL                $MEM_TOTAL
  MEM_STATIC               500

\$end
EOF

  if [[ ${SOLVENT,,} != water ]]; then
cat >>"$OUT" <<EOF

\$smx
  solvent $SOLVENT
  StateSpecific
\$end

\$gauge_origin
  center_of_mass
\$end
EOF
  fi
}

write_slurm() {
  local CID=$1 INP=$2 OUT=$3 SLM=$4

  cat >"$SLM" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=$CID
#SBATCH --partition=$PART
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem-per-cpu=$MEM_MB
#SBATCH --time=$DEFAULT_WALL
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

module purge
source /shares/chem_hlw/qchem/7.0.0-qmmm-nt-intel-pre/setqc
export QCSCRATCH=/tmp/\$SLURM_JOB_ID

qchem -nt $CPUS $(basename "$INP") $(basename "$OUT")
EOF
  chmod +x "$SLM"
}


# ---------------- main workhorses -------------------------------------------
submit_one(){
  local TAG=$1 CID=$2 OUTFILE=$3
  local DIR="$TAG/ecd_cc/$CID"; mkdir -p "$DIR"

  local XYZ_TMP; XYZ_TMP=$(mktemp --suffix=.xyz)
  qcout2xyz "$OUTFILE" > "$XYZ_TMP" || { note "[$TAG] cannot parse $CID"; rm -f "$XYZ_TMP"; return; }

  local INP="$DIR/$CID.inp" OUT="$DIR/$CID.out" SLM="$DIR/$CID.slurm"
  write_input "$TAG" "$CID" "$XYZ_TMP" "$INP"
  write_slurm "$CID" "$INP" "$OUT" "$SLM"
  rm -f "$XYZ_TMP"

  if $DRY; then
    note "[$TAG] sbatch --chdir=$DIR $(basename "$SLM")"
  else
    sbatch --chdir="$DIR" "$(basename "$SLM")"
  fi
}

process_tag(){
  local TAG=$1
  local LAB="$TAG/bw_results/${TAG}_bw_labels.dat"
  [[ -f $LAB ]] || { note "[$TAG] label list missing"; return; }
  mapfile -t CIDS < <(grep -v '^[[:space:]]*$' "$LAB")
  [[ ${#CIDS[@]} -eq 0 ]] && { note "[$TAG] label list empty"; return; }

  for CID in "${CIDS[@]}"; do
    local OUTPATH="$TAG/$SOLV_OUT_DIR/$CID/$CID.out"
    [[ -f $OUTPATH ]] || { note "[$TAG] missing $CID.out"; continue; }
    submit_one "$TAG" "$CID" "$OUTPATH"
  done
}

# -------------------------------- main ---------------------------------------
main(){
  parse_cli "$@"

  if [[ -n $LIST ]]; then
    grep -vE '^[[:space:]]*(#|;|$)' "$LIST" | while read -r tag; do
      process_tag "$tag"
    done
  else
    process_tag "$SINGLE"
  fi

  $DRY && note "Dry run complete - no jobs submitted"
}

main "$@"
