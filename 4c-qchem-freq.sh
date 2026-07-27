#!/usr/bin/env bash
# ============================================================================
# 4c-qchem-freq.sh -- Harmonic frequencies for Gibbs corrections (Stage 4c)
# ============================================================================
#
# OVERVIEW
#   Optional stage between Stage 4b (dedup) and Stage 5 (Boltzmann weighting).
#   Runs a Q-Chem FREQ job on each Stage-4 optimized geometry so that
#   4c-qchem-thermo.py can turn the electronic energies into free energies and
#   Stage 5 --gibbs can weight conformers on G instead of E.
#
#   Why bother: weighting on E alone systematically over-populates
#   intramolecularly hydrogen-bonded conformers, which are enthalpically
#   favoured but entropically penalized. For ECD, where folded and extended
#   conformers can contribute opposite-sign rotatory strengths, that bias can
#   change band signs rather than merely rescale them.
#
#   Second payoff: an imaginary frequency means Stage 4 converged to a saddle
#   point rather than a minimum. Without this stage nothing in the pipeline
#   notices, and the bad structure is weighted and passed to Stage 6 as though
#   it were real.
#
# LEVEL OF THEORY -- READ THIS
#   A Hessian is only meaningful at the geometry's own level of theory. Run
#   this with the SAME --method/--basis/--grid/--solvent/--disp you gave
#   Stage 4. If they differ, the geometry is not a stationary point for this
#   Hamiltonian and you will get spurious imaginary modes. The defaults here
#   deliberately mirror 4-qchem-solvent-opt.sh so that plain invocations agree.
#
# Usage:
#   4c-qchem-freq.sh TAG
#   4c-qchem-freq.sh --solvent methanol --grid fine ephedrine
#   4c-qchem-freq.sh --all --local --cpus 4 ephedrine
#
# Conformer selection (the two-pass workflow):
#   By default this reads the Stage-5 label file, so the intended sequence is
#
#     5-qchem-boltzmann-weight.sh --p-cut 0 --e-window 5 TAG   # loose prefilter
#     4c-qchem-freq.sh TAG                                     # freq on survivors
#     4c-qchem-thermo.py TAG                                   # G corrections
#     5-qchem-boltzmann-weight.sh --gibbs qrrho TAG            # final populations
#
#   That keeps the frequency cost off conformers that are far too high in
#   energy to matter. Widen --e-window if you would rather be safe than fast:
#   a conformer dropped in the first pass cannot be rescued by the second.
#
# Flags:
#   -m | --method NAME         DFT functional                  [B3LYP]
#   -b | --basis NAME          Basis set                       [def2-TZVP]
#        --disp KW             Dispersion: auto|none|D3BJ       [auto]
#        --ri KW               Density fitting: none|j|jk       [none]
#   -g | --grid LEVEL          Integration grid: coarse|default|fine|ultrafine [default]
#        --scf-conv LEVEL      SCF convergence: default|verytight|extreme [default]
#        --solvent NAME        SMD solvent keyword; 'none' = gas phase       [water]
#        --solvent-model NAME  Implicit model: smd|cpcm|iefpcm|cosmo  [smd]
#        --max-scf N           Max SCF cycles                  [150]
#        --all                 Use the Stage-4b keep-list (or every conformer)
#                              instead of the Stage-5 label file
#        --labels FILE         Explicit conformer-ID list
#   -c | --cpus N              CPU cores                       [4]
#        --mem-per-cpu MB      Memory per core in MB            [2048]
#        --max-running N       Max simultaneous SLURM array tasks [10]
#        --qchem-setup PATH    Path to Q-Chem setup script      [auto]
#        --variant LABEL       Suffix the output dir -> 03b_freq_LABEL [none]
#        --list FILE           File of molecule TAGs
#        --local               Run Q-Chem directly (no SLURM)
#        --dry-run             Write inputs but do not run
#   -h | --help                Show this help and exit
#
#   SLURM-only flags (ignored in --local mode):
#        --partition NAME      SLURM partition                  [general]
#        --time HH:MM:SS       Wall-clock limit                 [06:00:00]
#
# Directory layout:
#   <TAG>/
#   |-- 03_solvent_opt/<CID>/<CID>.out      <- Stage 4 geometry source
#   |-- 04_boltzmann/<TAG>_bw_labels.dat    <- conformer list (default source)
#   -- 03b_freq/
#       |-- <TAG>_conf_list.txt          (SLURM mode: list of working dirs)
#       |-- <TAG>_array.slurm            (SLURM mode: array job script)
#       |-- <CID>/<CID>.inp / .out
#       -- <TAG>_thermo.dat              (written later by 4c-qchem-thermo.py)
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS  (mirror 4-qchem-solvent-opt.sh -- see LEVEL OF THEORY above)
# ============================================================================
DEFAULT_METHOD="B3LYP"
DEFAULT_BASIS="def2-TZVP"
DEFAULT_DISP="auto"
DEFAULT_RI="none"
DEFAULT_GRID="default"
DEFAULT_SCF_CONV="default"
DEFAULT_SOLVENT="water"
DEFAULT_SOLVENT_MODEL="smd"
DEFAULT_MAX_SCF=150
DEFAULT_CPUS=4
DEFAULT_MEM_PER_CPU=2048
DEFAULT_PARTITION="general"
DEFAULT_WALL="06:00:00"
DEFAULT_MAX_RUNNING=10

XYZ_DIR="pre_xyz"
SOLV_OUT_DIR="03_solvent_opt"
BW_SUBDIR="04_boltzmann"
OUT_SUBDIR="03b_freq"

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
# CLUSTER CONFIG
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

resolve_qchem_setup() {
    local flag_val=$1
    if [[ -n $flag_val ]]; then printf '%s' "$flag_val"; return; fi
    if [[ -n ${QCHEM_SETUP:-} ]]; then printf '%s' "$QCHEM_SETUP"; return; fi
    printf ''
}

detect_mode() {
    local force_local=$1
    if $force_local; then printf 'local'
    elif command -v sbatch >/dev/null 2>&1; then printf 'slurm'
    else printf 'local'; fi
}

# ============================================================================
# Q-CHEM OUTPUT -> XYZ EXTRACTOR (last geometry in file)
# ============================================================================
qcout2xyz() {
    python3 - "$1" <<'PY'
import pathlib, re, sys
lines = pathlib.Path(sys.argv[1]).read_text(errors="ignore").splitlines()
hdr = re.compile(
    r'(Standard Nuclear Orientation|Input orientation|Coordinates \(Angstroms\))', re.I
)
start = None
for i in range(len(lines) - 1, -1, -1):
    if hdr.search(lines[i]):
        start = i + 3; break
if start is None:
    sys.exit(f"qcout2xyz: no coordinate block found in {sys.argv[1]}")
geom = []
for l in lines[start:]:
    if not l.strip() or l.lstrip().startswith('-'): break
    p = l.split()
    if len(p) >= 5:   sym, x, y, z = p[1], p[2], p[3], p[4]
    elif len(p) == 4: sym, x, y, z = p
    else: continue
    geom.append(f"{sym} {x} {y} {z}")
if not geom:
    sys.exit(f"qcout2xyz: coordinate parsing failed in {sys.argv[1]}")
print(len(geom)); print(); print("\n".join(geom))
PY
}

# ============================================================================
# CHARGE / MULTIPLICITY PARSER
# ============================================================================
read_xyz_header() {
    local xyz_file=$1
    local header_line
    header_line=$(sed -n '2p' "$xyz_file")

    if [[ $header_line =~ charge[[:space:]]*=[[:space:]]*([+-]?[0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
    else
        charge=0
    fi

    if [[ $header_line =~ mult[[:space:]]*=[[:space:]]*([0-9]+) ]]; then
        mult=${BASH_REMATCH[1]}
    elif [[ $header_line =~ ^[[:space:]]*([+-]?[0-9]+)[[:space:]]+([0-9]+) ]]; then
        charge=${BASH_REMATCH[1]}
        mult=${BASH_REMATCH[2]}
    else
        mult=1
    fi
}

# ============================================================================
# RI / DISPERSION / GRID LOGIC  (identical to Stages 4 and 6)
# ============================================================================
ri_lines() {
    case $ri_mode in
        none|NONE) printf '' ;;
        j|J)
            [[ ${basis,,} == def2-* ]] || \
                die "--ri j needs a def2-* basis (got '${basis}')"
            printf '\n  AUX_BASIS_J         RIJ-%s' "$basis" ;;
        jk|JK)
            [[ ${basis,,} == def2-* ]] || \
                die "--ri jk needs a def2-* basis (got '${basis}')"
            printf '\n  AUX_BASIS_J         RIJK-%s\n  AUX_BASIS_K         RIJK-%s' "$basis" "$basis" ;;
        *) die "--ri must be none, j, or jk" ;;
    esac
}

disp_line() {
    local method_upper=${method^^}
    case $disp_mode in
        none|NONE) printf '' ;;
        D3BJ|d3bj) printf '\n  DFT_D               D3_BJ' ;;
        auto|AUTO)
            if [[ $method =~ (-D[0-9]?|-D3BJ|-D3ZERO|-D4)($|[[:space:]]) ]]; then
                printf ''; return; fi
            case $method_upper in
                WB97X-D|WB97X-D3|WB97X-D4|WB97X-V|WB97XD|WB97M-V|B97-D|B97-D3)
                    printf ''; return ;; esac
            printf '\n  DFT_D               D3_BJ' ;;
        *) die "--disp must be auto, none, or D3BJ" ;;
    esac
}

xc_grid_value() {
    case ${grid,,} in
        coarse)    printf '1' ;;
        default)   printf '3' ;;
        fine)      printf '000099000590' ;;
        ultrafine) printf '000099000974' ;;
        sg1|1)     printf '1' ;;
        sg2|2)     printf '2' ;;
        sg3|3)     printf '3' ;;
        *)         printf '%s' "$grid" ;;
    esac
}

scf_conv_value() {
    case ${scf_conv,,} in
        default)   printf '8' ;;
        verytight) printf '10' ;;
        extreme)   printf '11' ;;
        *)         printf '8' ;;
    esac
}

# ============================================================================
# CLI PARSER
# ============================================================================
parse_cli() {
    method=$DEFAULT_METHOD
    basis=$DEFAULT_BASIS
    disp_mode=$DEFAULT_DISP
    ri_mode=$DEFAULT_RI
    grid=$DEFAULT_GRID
    scf_conv=$DEFAULT_SCF_CONV
    solvent=$DEFAULT_SOLVENT
    solvent_model=$DEFAULT_SOLVENT_MODEL
    max_scf=$DEFAULT_MAX_SCF
    cpus=$DEFAULT_CPUS
    mem_mb=$DEFAULT_MEM_PER_CPU
    partition=${CLUSTER_PARTITION:-$DEFAULT_PARTITION}
    wall=${CLUSTER_WALL:-$DEFAULT_WALL}
    max_running=${CLUSTER_MAX_RUNNING:-$DEFAULT_MAX_RUNNING}
    dry_run=false
    force_local=false
    use_all=false
    labels_file=""
    list_file=""
    single_tag=""
    qchem_setup_flag=""
    variant=""

    local opts
    opts=$(getopt -o hb:m:c:g: \
        --long help,method:,basis:,disp:,ri:,grid:,scf-conv:,variant:,solvent:,solvent-model:,\
max-scf:,cpus:,mem-per-cpu:,max-running:,partition:,time:,list:,labels:,all,qchem-setup:,local,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            -m|--method)      method=$2;            shift 2 ;;
            -b|--basis)       basis=$2;             shift 2 ;;
            --disp)           disp_mode=$2;         shift 2 ;;
            --ri)             ri_mode=$2;           shift 2 ;;
            -g|--grid)        grid=$2;              shift 2 ;;
            --scf-conv)       scf_conv=$2;          shift 2 ;;
            --variant)        variant=$2;           shift 2 ;;
            --solvent)        solvent=$2;           shift 2 ;;
            --solvent-model)  solvent_model=$2;     shift 2 ;;
            --max-scf)        max_scf=$2;           shift 2 ;;
            --all)            use_all=true;         shift ;;
            --labels)         labels_file=$2;       shift 2 ;;
            -c|--cpus)        cpus=$2;              shift 2 ;;
            --mem-per-cpu)    mem_mb=$2;            shift 2 ;;
            --max-running)    max_running=$2;       shift 2 ;;
            --partition)      partition=$2;         shift 2 ;;
            --time)           wall=$2;              shift 2 ;;
            --list)           list_file=$2;         shift 2 ;;
            --qchem-setup)    qchem_setup_flag=$2;  shift 2 ;;
            --local)          force_local=true;     shift ;;
            --dry-run)        dry_run=true;         shift ;;
            -h|--help)        show_help ;;
            --)               shift; break ;;
            *)                die "Unknown option '$1'" ;;
        esac
    done

    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional args not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide exactly one molecule TAG"
        single_tag=$1
    fi

    $use_all && [[ -n $labels_file ]] && die "--all and --labels are mutually exclusive"
    [[ -n $labels_file ]] && require_file "$labels_file"

    if [[ -n $variant ]]; then
        [[ $variant =~ ^[A-Za-z0-9._-]+$ ]] || die "--variant must contain only letters, digits, . _ - (got '${variant}')"
    fi

    case ${grid,,} in
        coarse|default|fine|ultrafine|sg1|sg2|sg3|1|2|3) ;;
        *) [[ $grid =~ ^[0-9]{12}$ ]] || die "--grid must be coarse|default|fine|ultrafine, SG1/SG2/SG3, 1/2/3, or a 12-digit XC_GRID code (got '${grid}')" ;;
    esac
    case ${scf_conv,,} in
        default|verytight|extreme) ;;
        *) die "--scf-conv must be default, verytight, or extreme (got '${scf_conv}')" ;;
    esac
    case $solvent_model in
        smd|cpcm|iefpcm|cosmo) ;;
        *) die "--solvent-model must be smd, cpcm, iefpcm, or cosmo (got '${solvent_model}')" ;;
    esac
}

# ============================================================================
# Q-CHEM INPUT WRITER
# ============================================================================
write_qchem_input() {
    local cid=$1 xyz_file=$2 inp_file=$3
    local mem_total=$(( cpus * mem_mb ))
    local ri disp grid_val scf_val thresh_line
    ri=$(ri_lines)
    disp=$(disp_line)
    grid_val=$(xc_grid_value)
    scf_val=$(scf_conv_value)
    [[ ${scf_conv,,} == default ]] && thresh_line="" || thresh_line=$'\n  THRESH              14'

    local solvent_rem solvent_blocks solvent_desc pcm_theory
    if [[ ${solvent,,} == none ]]; then
        solvent_rem=""; solvent_blocks=""; solvent_desc="gas phase"
    elif [[ $solvent_model == smd ]]; then
        solvent_rem=$'\n  SOLVENT_METHOD      SMD'
        printf -v solvent_blocks '\n\n$smx\n  solvent %s\n$end' "$solvent"
        solvent_desc="SMD ${solvent}"
    elif [[ $solvent_model == cosmo ]]; then
        solvent_rem=$'\n  SOLVENT_METHOD      COSMO'
        printf -v solvent_blocks '\n\n$solvent\n  SolventName %s\n$end' "$solvent"
        solvent_desc="COSMO ${solvent}"
    else
        [[ $solvent_model == iefpcm ]] && pcm_theory=IEFPCM || pcm_theory=CPCM
        solvent_rem=$'\n  SOLVENT_METHOD      PCM'
        printf -v solvent_blocks '\n\n$pcm\n  Theory %s\n$end\n\n$solvent\n  SolventName %s\n$end' "$pcm_theory" "$solvent"
        solvent_desc="${pcm_theory} ${solvent}"
    fi

    cat >"$inp_file" <<EOF
\$comment
Harmonic frequencies for Gibbs corrections (${solvent_desc}) -- generated by 4c-qchem-freq.sh
Level of theory must match Stage 4 or the Hessian is not at a stationary point.
\$end

\$molecule
${charge} ${mult}
$(tail -n +3 "$xyz_file")
\$end

\$rem
  JOB_TYPE            FREQ
  METHOD              ${method}
  BASIS               ${basis}${solvent_rem}
  SCF_CONVERGENCE     ${scf_val}${thresh_line}
  MAX_SCF_CYCLES      ${max_scf}
  SYM_IGNORE          TRUE
  XC_GRID             ${grid_val}${disp}${ri}
  MEM_TOTAL           ${mem_total}
  MEM_STATIC          500
\$end${solvent_blocks}
EOF
}

# ============================================================================
# SLURM ARRAY WRITER
# ============================================================================
write_array_slurm() {
    local tag=$1 conf_list=$2 slurm_file=$3 n=$4
    local abs_list out_dir
    abs_list=$(cd "$(dirname "$conf_list")" && pwd)/$(basename "$conf_list")
    out_dir=$(dirname "$abs_list")

    cat >"$slurm_file" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=freq_${tag}
#SBATCH --partition=${partition}
#SBATCH --array=1-${n}%${max_running}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem-per-cpu=${mem_mb}
#SBATCH --time=${wall}
#SBATCH --output=${out_dir}/slurm-%A_%a.out
#SBATCH --error=${out_dir}/slurm-%A_%a.err

# ---- Q-Chem environment ----
${qchem_setup:+source ${qchem_setup}}
${CLUSTER_LD_LIBRARY_PATH:+export LD_LIBRARY_PATH=${CLUSTER_LD_LIBRARY_PATH}:\$LD_LIBRARY_PATH}
export QCSCRATCH=/tmp/\$SLURM_JOB_ID

WORKDIR=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${abs_list}")
CID=\$(basename "\$WORKDIR")
cd "\$WORKDIR"
qchem -nt ${cpus} "\${CID}.inp" "\${CID}.out"
EOF
    chmod +x "$slurm_file"
}

# ============================================================================
# CONFORMER LIST RESOLUTION
# ============================================================================
# Default source is the Stage-5 label file, which makes the two-pass workflow
# the path of least resistance. --all falls back to the Stage-4b keep-list, and
# then to every conformer directory present.
resolve_cids() {
    local tag=$1
    local src=""

    if [[ -n $labels_file ]]; then
        src=$labels_file
    elif $use_all; then
        local keep="${tag}/${SOLV_OUT_DIR}/${tag}_unique.dat"
        if [[ -f $keep ]]; then
            src=$keep
        else
            # no dedup keep-list: use every Stage-4 conformer directory
            find "${tag}/${SOLV_OUT_DIR}" -mindepth 1 -maxdepth 1 -type d \
                -printf '%f\n' 2>/dev/null | sort
            return
        fi
    else
        src="${tag}/${BW_SUBDIR}/${tag}_bw_labels.dat"
        if [[ ! -f $src ]]; then
            warn "[${tag}] ${src} not found."
            warn "[${tag}] Run Stage 5 first for a prefiltered list, or pass --all."
            return
        fi
    fi

    tr -d '\r' < "$src" | grep -v '^[[:space:]]*$' || true
}

# ============================================================================
# PROCESS ONE MOLECULE
# ============================================================================
process_tag() {
    local tag=$1

    local cids=()
    mapfile -t cids < <(resolve_cids "$tag")
    if (( ${#cids[@]} == 0 )); then
        warn "[${tag}] no conformers selected -- skipping"
        return
    fi

    charge=0; mult=1
    local pre_xyz="${XYZ_DIR}/${tag}.xyz"
    if [[ -f $pre_xyz ]]; then
        read_xyz_header "$pre_xyz"
    else
        warn "[${tag}] ${pre_xyz} not found -- using charge=0 mult=1"
    fi

    log "[${tag}] ${#cids[@]} conformers to process (charge=${charge}, mult=${mult})"

    local out_dir="${tag}/${OUT_SUBDIR}"
    mkdir -p "$out_dir"

    local conf_dirs=()
    local cid outpath dir xyz_tmp inp_file abs_dir
    for cid in "${cids[@]}"; do
        cid=$(echo "$cid" | xargs)
        [[ -z $cid ]] && continue

        outpath="${tag}/${SOLV_OUT_DIR}/${cid}/${cid}.out"
        if [[ ! -f $outpath ]]; then
            warn "  [${cid}] solvent optimization output not found -- skipping"
            continue
        fi

        dir="${out_dir}/${cid}"
        mkdir -p "$dir"

        xyz_tmp=$(mktemp --suffix=.xyz)
        if ! qcout2xyz "$outpath" >"$xyz_tmp" 2>/dev/null; then
            warn "  [${cid}] cannot extract geometry from ${outpath} -- skipping"
            rm -f "$xyz_tmp"
            continue
        fi

        inp_file="${dir}/${cid}.inp"
        write_qchem_input "$cid" "$xyz_tmp" "$inp_file"
        rm -f "$xyz_tmp"

        if $dry_run; then
            log "  [${cid}] dry run -- input written to ${inp_file}"
            continue
        fi

        abs_dir=$(cd "$dir" && pwd)
        conf_dirs+=("$abs_dir")
    done

    $dry_run && return

    local n=${#conf_dirs[@]}
    if (( n == 0 )); then
        warn "[${tag}] no valid conformer inputs written -- nothing to run"
        return
    fi

    # ---- local mode --------------------------------------------------------
    if [[ $exec_mode != slurm ]]; then
        if [[ -n $qchem_setup ]]; then
            # shellcheck source=/dev/null
            source "$qchem_setup"
        fi
        for abs_dir in "${conf_dirs[@]}"; do
            cid=$(basename "$abs_dir")
            log "  [${cid}] running FREQ"
            ( cd "$abs_dir" && qchem -nt "$cpus" "${cid}.inp" "${cid}.out" ) 2>&1
            if grep -q 'Frequency:' "${abs_dir}/${cid}.out" 2>/dev/null; then
                log "  [${cid}] frequencies complete"
            else
                warn "  [${cid}] no frequencies found -- check ${abs_dir}/${cid}.out"
            fi
        done
        log "[${tag}] next: 4c-qchem-thermo.py ${tag}"
        return
    fi

    # ---- SLURM mode --------------------------------------------------------
    local conf_list="${out_dir}/${tag}_conf_list.txt"
    printf '%s\n' "${conf_dirs[@]}" >"$conf_list"

    local array_slurm="${out_dir}/${tag}_array.slurm"
    write_array_slurm "$tag" "$conf_list" "$array_slurm" "$n"

    local jid
    jid=$(sbatch --parsable "$array_slurm")
    log "[${tag}] submitted array job ${jid} (${n} conformers, max ${max_running} concurrent)"
    log "[${tag}] when it finishes: 4c-qchem-thermo.py ${tag}"
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    local src_desc
    if [[ -n $labels_file ]]; then src_desc="--labels ${labels_file}"
    elif $use_all; then src_desc="all (Stage-4b keep-list)"
    else src_desc="Stage-5 label file (prefiltered)"; fi

    cat >&2 <<EOF
=============================================================
 Stage 4c: Harmonic Frequencies (for Gibbs corrections)
-------------------------------------------------------------
 Mode         : ${exec_mode}
 Q-Chem setup : ${qchem_setup:-<qchem in PATH>}
 Method       : ${method}
 Basis        : ${basis}
 Dispersion   : ${disp_mode}
 Grid         : ${grid} (XC_GRID $(xc_grid_value))
 SCF conv     : ${scf_conv} (SCF_CONVERGENCE $(scf_conv_value))
 Solvent      : ${solvent}$([[ ${solvent,,} == none ]] && printf ' (gas phase)' || printf ' (%s)' "$solvent_model")
 Conformers   : ${src_desc}
 Cores        : ${cpus}
 Mem/core     : ${mem_mb} MB
 Output dir   : <TAG>/${OUT_SUBDIR}
 Max running  : ${max_running} (SLURM array throttle)
 Dry run      : ${dry_run}
-------------------------------------------------------------
 Level of theory MUST match Stage 4, or the Hessian is not at
 a stationary point and you will see spurious imaginary modes.
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    source_cluster_cfg
    parse_cli "$@"

    [[ -n $variant ]] && OUT_SUBDIR="${OUT_SUBDIR}_${variant}"

    exec_mode=$(detect_mode "$force_local")
    qchem_setup=$(resolve_qchem_setup "$qchem_setup_flag")

    if ! $dry_run && [[ -z $qchem_setup ]]; then
        if ! command -v qchem >/dev/null 2>&1; then
            die "Q-Chem not found. Set QCHEM_SETUP in cluster.cfg, use --qchem-setup, or add qchem to PATH."
        fi
    fi
    if ! $dry_run && [[ -n $qchem_setup && ! -f $qchem_setup ]]; then
        die "Q-Chem setup script '${qchem_setup}' not found."
    fi
    if [[ -z $qchem_setup ]] && $dry_run; then
        warn "Q-Chem setup not found -- dry-run will proceed without it"
    fi

    print_banner

    if [[ -n $list_file ]]; then
        while IFS= read -r tag || [[ -n $tag ]]; do
            [[ -z $tag || $tag == \#* || $tag == \;* ]] && continue
            process_tag "$tag"
        done < "$list_file"
    else
        process_tag "$single_tag"
    fi

    log "Stage 4c complete."
}

main "$@"
