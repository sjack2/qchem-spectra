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
#   --e-window KCAL  Also require dE (or dG) <= KCAL          [none]
#   --gibbs MODE     Weight on free energy: none|rrho|qrrho   [none]
#   --keep-saddles   Keep conformers with an imaginary mode   [drop them]
#   --freq-dir D     Stage-4c output subdir                   [03b_freq]
#   --list FILE      Text file of molecule TAGs
#   --dry-run        Show what would be computed without writing
#   -h | --help      Show this help and exit
#
# FREE-ENERGY WEIGHTING (--gibbs)
#   By default this stage weights on the Stage-4 SCF electronic energy. That
#   systematically over-populates intramolecularly hydrogen-bonded conformers,
#   which are enthalpically favoured but entropically penalized. --gibbs adds
#   the thermal correction from 4c-qchem-thermo.py so the weighting is done on
#   G instead:
#
#       G = E_elec (Stage 4)  +  G_corr (<TAG>/03b_freq/<TAG>_thermo.dat)
#
#   MODE picks which correction column to use:
#     rrho    textbook harmonic. Vibrational entropy diverges as the frequency
#             goes to zero, so for a flexible molecule the low torsional modes
#             dominate the entropy and are the least trustworthy numbers here.
#     qrrho   Grimme's quasi-harmonic damping of those modes. Preferred for
#             flexible molecules; compare against rrho to see how much the low
#             modes are actually moving your populations.
#
#   Conformers whose frequency job found a significant imaginary mode are
#   saddle points, not minima, and are dropped unless --keep-saddles is given.
#
# TWO-PASS WORKFLOW
#   Frequency jobs are expensive, so run this stage twice and put the cost only
#   on conformers that could plausibly matter:
#
#     5-qchem-boltzmann-weight.sh --p-cut 0 --e-window 5 TAG   # loose prefilter
#     4c-qchem-freq.sh TAG                                     # freq on survivors
#     4c-qchem-thermo.py TAG                                   # G corrections
#     5-qchem-boltzmann-weight.sh --gibbs qrrho TAG            # final populations
#
#   A conformer cut in the first pass cannot come back in the second, so keep
#   the first window generous.
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
DEFAULT_GIBBS="none"
DEFAULT_FREQ_DIR="03b_freq"

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
    gibbs_mode=$DEFAULT_GIBBS
    freq_dir=$DEFAULT_FREQ_DIR
    e_window=""
    keep_saddles=false
    dry_run=false
    list_file=""
    single_tag=""

    local opts
    opts=$(getopt -o h --long help,temp:,p-cut:,e-window:,gibbs:,freq-dir:,keep-saddles,list:,dry-run -- "$@") \
        || die "Failed to parse options (try --help)"
    eval set -- "$opts"

    while true; do
        case $1 in
            --temp)         temp_k=$2;         shift 2 ;;
            --p-cut)        p_cut=$2;          shift 2 ;;
            --e-window)     e_window=$2;       shift 2 ;;
            --gibbs)        gibbs_mode=$2;     shift 2 ;;
            --freq-dir)     freq_dir=$2;       shift 2 ;;
            --keep-saddles) keep_saddles=true; shift ;;
            --list)         list_file=$2;      shift 2 ;;
            --dry-run)      dry_run=true;      shift ;;
            -h|--help)      show_help ;;
            --)             shift; break ;;
            *)              die "Unknown option '$1'" ;;
        esac
    done

    if [[ -n $list_file ]]; then
        [[ $# -eq 0 ]] || die "Positional args not allowed with --list"
        require_file "$list_file"
    else
        [[ $# -eq 1 ]] || die "Provide exactly one molecule TAG"
        single_tag=$1
    fi

    case ${gibbs_mode,,} in
        none|rrho|qrrho) gibbs_mode=${gibbs_mode,,} ;;
        *) die "--gibbs must be none, rrho, or qrrho (got '${gibbs_mode}')" ;;
    esac
    if [[ -n $e_window ]]; then
        [[ $e_window =~ ^[0-9]*\.?[0-9]+$ ]] || \
            die "--e-window must be a non-negative number in kcal/mol (got '${e_window}')"
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

    # honor the Stage-4b dedup keep-list (unique conformers) when present
    local keep_file="${opt_dir}/${tag}_unique.dat"
    [[ -f $keep_file ]] && \
        log "[${tag}] honoring dedup keep-list ($(wc -l < "$keep_file") unique conformers)"

    # find all Q-Chem output files and extract the final energy from each
    while IFS= read -r outfile; do
        local energy_line cid energy
        cid=$(basename "${outfile%.*}")
        if [[ -f $keep_file ]] && ! tr -d '\r' < "$keep_file" | grep -qxF "$cid"; then
            continue
        fi
        energy_line=$(grep 'Final energy is' "$outfile" | tail -n1) || true
        [[ -n $energy_line ]] || continue
        energy=$(awk '{print $NF}' <<< "$energy_line")
        printf '%s\t%s\n' "$cid" "$energy"
    done < <(find "$opt_dir" -type f -name '*.out' | sort)
}

# ============================================================================
# GIBBS CORRECTION JOIN
# ============================================================================
# Adds G_corr from the Stage-4c thermo table to each electronic energy, so the
# downstream Boltzmann math is unchanged -- it just receives G instead of E.
# Conformers with no thermo row, or with an unusable one, are dropped rather
# than silently falling back to E: a table that mixes G for some conformers and
# E for others would produce populations that mean nothing.
apply_gibbs() {
    local tag=$1 raw=$2
    local thermo="${tag}/${freq_dir}/${tag}_thermo.dat"

    if [[ ! -f $thermo ]]; then
        die "[${tag}] --gibbs ${gibbs_mode} needs ${thermo}, which does not exist.
       Run 4c-qchem-freq.sh then 4c-qchem-thermo.py first, or drop --gibbs."
    fi

    # column 2 = RRHO correction, column 3 = qRRHO correction
    local gcol
    [[ $gibbs_mode == rrho ]] && gcol=2 || gcol=3

    awk -v gcol="$gcol" -v keep="$keep_saddles" -v want_t="$temp_k" -v tag="$tag" '
    BEGIN { FS = "[ \t]+"; OFS = "\t"; ttab = "" }
    # ---- first file: the thermo table ----
    NR == FNR {
        if ($1 ~ /^#/) {
            if ($2 == "temperature_K") ttab = $3
            next
        }
        if (NF < 9) next
        g[$1]  = $(gcol)
        st[$1] = $9
        next
    }
    # ---- second file: cid <tab> E_elec ----
    FNR == 1 && ttab != "" {
        # a G computed at one temperature and Boltzmann-weighted at another is
        # simply wrong, so this is worth shouting about
        dt = ttab - want_t
        if (dt < 0) dt = -dt
        if (dt > 0.01)
            printf "[%s] Warning: thermo table is at %s K but weighting at %s K\n",
                   tag, ttab, want_t > "/dev/stderr"
    }
    {
        cid = $1
        if (!(cid in g)) { miss++; next }
        if (g[cid] == "nan") { bad++; next }
        if (st[cid] == "saddle" && keep != "true") { sad++; next }
        printf "%s\t%.10f\n", cid, $2 + g[cid]
    }
    END {
        if (miss) printf "[%s] Warning: %d conformer(s) had no thermo row -- dropped\n", tag, miss > "/dev/stderr"
        if (bad)  printf "[%s] Warning: %d conformer(s) had unusable thermochemistry -- dropped\n", tag, bad > "/dev/stderr"
        if (sad)  printf "[%s] %d conformer(s) dropped as saddle points (--keep-saddles overrides)\n", tag, sad > "/dev/stderr"
    }' "$thermo" <(printf '%s\n' "$raw")
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

    if [[ $gibbs_mode != none ]]; then
        local n_before
        n_before=$(wc -l <<< "$raw")
        raw=$(apply_gibbs "$tag" "$raw")
        if [[ -z $raw ]]; then
            warn "[${tag}] no conformers survived the Gibbs join -- skipping"
            return
        fi
        log "[${tag}] free-energy weighting (${gibbs_mode}): $(wc -l <<< "$raw")/${n_before} conformers"
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

    # column labels track whether we are weighting on E or on G
    local e_label d_label
    if [[ $gibbs_mode == none ]]; then
        e_label="E(Ha)";  d_label="dE(kcal)"
    else
        e_label="G(Ha)";  d_label="dG(kcal)"
    fi

    # compute Boltzmann weights using awk (no Python dependency)
    awk -v T="$temp_k" -v R="$R_J" -v H2K="$H2KCAL" -v pc="$p_cut" \
        -v ewin="$e_window" -v elab="$e_label" -v dlab="$d_label" \
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
        printf "#%-24s  %-16s  %-10s  %-8s\n", "CID", elab, dlab, "p" > dat
        for (i = 1; i <= n; i++) {
            p = B[i] / Z
            printf "%-25s  %16.8f  %10.3f  %8.5f\n", cid[i], E[i], dE[i], p >> dat
            # --e-window is an AND with --p-cut, so `--p-cut 0 --e-window 5`
            # selects purely on the energy window (the prefilter pass) while
            # the default `--p-cut 0.01` with no window behaves as before
            keep = (p >= pc)
            if (ewin != "" && dE[i] > ewin + 0) keep = 0
            if (keep) print cid[i] >> lab
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
    local crit="p >= ${p_cut}"
    [[ -n $e_window ]] && crit="${crit} and ${d_label} <= ${e_window}"
    log "[${tag}] ${n_keep}/${n_confs} conformers pass ${crit}"
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
 Energy win. : ${e_window:-<none>}$([[ -n $e_window ]] && printf ' kcal/mol')
 Weighting   : $([[ $gibbs_mode == none ]] && printf 'electronic energy (Stage 4 SCF)' || printf 'free energy G, %s (%s/<TAG>_thermo.dat)' "$gibbs_mode" "$freq_dir")
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
