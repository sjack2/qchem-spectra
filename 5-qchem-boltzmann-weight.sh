#!/usr/bin/env bash
# ============================================================================
# 5-qchem-boltzmann-weight.sh -- Boltzmann weighting & filtering (Stage 5)
# ============================================================================
#
# OVERVIEW
#   Collects Q-Chem final SCF energies from the solvent-phase optimization
#   outputs (Stage 4), computes relative energies and Boltzmann populations
#   at a given temperature, and writes two files:
#
#     <TAG>/04_boltzmann/<TAG>_energies.dat    full table sorted by probability
#     <TAG>/04_boltzmann/<TAG>_bw_labels.dat   conformer IDs above the cutoff
#
#   The label file is read by Stage 6 to determine which conformers get
#   TD-DFT or frequency calculations.
#
# Usage:
#   5-qchem-boltzmann-weight.sh TAG
#   5-qchem-boltzmann-weight.sh --temp 310 --p-cut 0.02 aspirin
#   5-qchem-boltzmann-weight.sh --list mols.txt --dry-run
#
# Flags:
#   --temp K         Temperature in Kelvin                    [298.15]
#   --p-cut VAL      Probability cutoff (e.g. 0.01 = 1%)     [0.01]
#   --list FILE      Text file of molecule TAGs
#   --dry-run        Show what would be computed without writing
#   -h | --help      Show this help and exit
#
# Directory layout:
#   <TAG>/
#   |-- 03_solvent_opt/
#   |   |-- <TAG>_001/<TAG>_001.out      <- input (Stage 4 outputs)
#   |   -- <TAG>_002/<TAG>_002.out
#   -- 04_boltzmann/
#       |-- <TAG>_energies.dat           full table (CID E dE p)
#       -- <TAG>_bw_labels.dat          filtered labels -> Stage 6
#
# Examples:
#   5-qchem-boltzmann-weight.sh ephedrine
#   5-qchem-boltzmann-weight.sh --temp 310 --p-cut 0.05 --list molecules.txt
#
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ============================================================================
# DEFAULTS & CONSTANTS
# ============================================================================
H2KCAL=627.509474       # Hartree -> kcal/mol
R_J=8.314462618          # J mol-1 K-1
DEFAULT_TEMP=298.15
DEFAULT_P_CUT=0.01

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
    temp_k=$DEFAULT_TEMP
    p_cut=$DEFAULT_P_CUT
    dry_run=false
    list_file=""
    single_tag=""

    local opts
    opts=$(getopt -o h --long help,temp:,p-cut:,list:,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            --temp)    temp_k=$2;   shift 2 ;;
            --p-cut)   p_cut=$2;    shift 2 ;;
            --list)    list_file=$2; shift 2 ;;
            --dry-run) dry_run=true; shift ;;
            -h|--help) show_help ;;
            --)        shift; break ;;
            *)         die "Unknown option '$1'" ;;
        esac
    done

    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional args not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide exactly one molecule TAG"
        single_tag=$1
    fi
}

# ============================================================================
# ENERGY EXTRACTION
# ============================================================================
# Scans all .out files under <TAG>/03_solvent_opt/ for the Q-Chem
# energy string. Outputs "conformer_id  energy_hartree" lines.
extract_energies() {
    local tag=$1
    local opt_dir="${tag}/03_solvent_opt"

    if [[ ! -d $opt_dir ]]; then
        warn "[${tag}] 03_solvent_opt/ directory not found"
        return
    fi

    # find all Q-Chem output files and extract the final energy from each
    while IFS= read -r outfile; do
        local energy_line cid energy
        energy_line=$(grep 'Final energy is' "$outfile" | tail -n1) || true
        [[ -n $energy_line ]] || continue
        energy=$(awk '{print $NF}' <<< "$energy_line")
        cid=$(basename "${outfile%.*}")
        printf '%s\t%s\n' "$cid" "$energy"
    done < <(find "$opt_dir" -type f -name '*.out' | sort)
}

# ============================================================================
# PROCESS ONE MOLECULE
# ============================================================================
process_tag() {
    local tag=$1

    local raw
    raw=$(extract_energies "$tag") || true
    if [[ -z $raw ]]; then
        warn "[${tag}] no energies extracted -- skipping"
        return
    fi

    local n_confs
    n_confs=$(wc -l <<< "$raw")
    log "[${tag}] extracted energies for ${n_confs} conformers"

    if $dry_run; then
        log "[${tag}] (dry run) would compute Boltzmann weights at ${temp_k} K, cutoff ${p_cut}"
        echo "$raw" | column -t >&2
        return
    fi

    local out_dir="${tag}/04_boltzmann"
    local dat_file="${out_dir}/${tag}_energies.dat"
    local lab_file="${out_dir}/${tag}_bw_labels.dat"
    mkdir -p "$out_dir"

    # remove any previous label file (awk appends)
    rm -f "$lab_file"

    # compute Boltzmann weights using awk (no Python dependency)
    awk -v T="$temp_k" -v R="$R_J" -v H2K="$H2KCAL" -v pc="$p_cut" \
        -v dat="$dat_file" -v lab="$lab_file" '
    BEGIN { FS = OFS = "\t" }
    {
        cid[NR] = $1
        E[NR]   = $2
        n = NR
        if (NR == 1 || $2 < Emin) Emin = $2
    }
    END {
        # compute partition function
        Z = 0
        for (i = 1; i <= n; i++) {
            dE_kcal = (E[i] - Emin) * H2K
            boltz   = exp(-dE_kcal * 4184 / (R * T))
            B[i]    = boltz
            Z      += boltz
            dE[i]   = dE_kcal
        }

        # write full table (will be sorted afterwards)
        printf "#%-24s  %-16s  %-10s  %-8s\n", "CID", "E(Ha)", "dE(kcal)", "p" > dat
        for (i = 1; i <= n; i++) {
            p = B[i] / Z
            printf "%-25s  %16.8f  %10.3f  %8.5f\n", cid[i], E[i], dE[i], p >> dat
            if (p >= pc) print cid[i] >> lab
        }
    }' <<< "$raw"

    # sort the energies file by probability (descending)
    if [[ -f $dat_file ]]; then
        local header body
        header=$(head -n1 "$dat_file")
        body=$(tail -n +2 "$dat_file" | sort -k4 -nr)
        { echo "$header"; echo "$body"; } > "$dat_file"
    fi

    local n_keep=0
    [[ -f $lab_file ]] && n_keep=$(wc -l < "$lab_file")
    log "[${tag}] ${n_keep}/${n_confs} conformers above p >= ${p_cut} threshold"
    log "[${tag}] wrote ${dat_file} and ${lab_file}"
}

# ============================================================================
# SUMMARY BANNER
# ============================================================================
print_banner() {
    cat >&2 <<EOF
=============================================================
 Stage 5: Boltzmann Weighting & Filtering
-------------------------------------------------------------
 Temperature : ${temp_k} K
 P cutoff    : ${p_cut}
 Dry run     : ${dry_run}
=============================================================
EOF
}

# ============================================================================
# MAIN
# ============================================================================
main() {
    parse_cli "$@"
    print_banner

    if [[ -n $list_file ]]; then
        while IFS= read -r tag || [[ -n $tag ]]; do
            [[ -z $tag || $tag == \#* || $tag == \;* ]] && continue
            process_tag "$tag"
        done < "$list_file"
    else
        process_tag "$single_tag"
    fi

    log "Stage 5 complete."
}

main "$@"
