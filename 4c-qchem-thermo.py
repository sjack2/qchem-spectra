#!/usr/bin/env python3
"""
4c-qchem-thermo.py -- Gibbs free-energy corrections from Stage-4c frequencies.

Stage 5 weights conformers on the Stage-4 SCF electronic energy alone. That
systematically over-populates intramolecularly hydrogen-bonded conformers,
which are enthalpically favoured but entropically penalized (the hydrogen bond
freezes a torsion). This tool computes the thermal correction that turns E into
G, so Stage 5 --gibbs can weight on free energy instead.

    G = E_elec + G_corr        <- E_elec from Stage 4, G_corr from here

Two entropy treatments are reported side by side for every conformer:

  RRHO   textbook rigid-rotor/harmonic-oscillator.
  qRRHO  Grimme's quasi-harmonic treatment (Chem. Eur. J. 2012, 18, 9955).
         Harmonic vibrational entropy diverges as the frequency goes to zero,
         so for a flexible molecule the 20-60 cm-1 torsional modes -- exactly
         the ones that distinguish folded from extended conformers -- dominate
         the entropy and are the least reliable numbers in the calculation.
         qRRHO interpolates each mode between the harmonic and free-rotor
         limits with a damping function centred on --qrrho-cutoff, which keeps
         those modes from running away. Applied to the ENTROPY only, which is
         Grimme's original formulation.

Comparing the two columns tells you how much the low-frequency modes are
actually moving the populations. If they disagree badly, the harmonic numbers
are the ones to distrust.

Reads:   <TAG>/03b_freq/<CID>/<CID>.out       (Stage 4c frequency jobs)
Writes:  <TAG>/03b_freq/<TAG>_thermo.dat      (read by Stage 5 --gibbs)

Usage:
    4c-qchem-thermo.py TAG
    4c-qchem-thermo.py --temp 310 ephedrine
    4c-qchem-thermo.py --list molecules.txt
    4c-qchem-thermo.py --dry-run ephedrine       # report only, write nothing

Flags:
    --temp K            Temperature in Kelvin                      [298.15]
    --pressure ATM      Pressure in atmospheres                    [1.0]
    --qrrho-cutoff CM   qRRHO damping centre, cm-1                 [100.0]
    --imag-tol CM       |imaginary| below this is treated as noise [50.0]
    --freq-dir D        Stage-4c output subdir                     [03b_freq]
    --list FILE         Text file of TAGs (one per line)
    --dry-run           Print the report but do not write the table
    -h, --help          Show this help

Symmetry number:
    Fixed at 1. Every conformer of a chiral flexible molecule is C1, and the
    pipeline runs with SYM_IGNORE TRUE throughout. If you ever point this at a
    symmetric species the rotational entropy will be too high by R*ln(sigma).

Engine:
    Q-Chem only. The ORCA mirror prints its frequencies and geometry in a
    different format and needs its own parser.
"""
import argparse
import glob
import math
import os
import re
import sys

# ---------------------------------------------------------------------------
# CODATA 2018 constants (SI unless noted)
# ---------------------------------------------------------------------------
H_PLANCK = 6.62607015e-34        # J s
K_B = 1.380649e-23               # J / K
C_LIGHT = 2.99792458e10          # cm / s   (note: cm, to pair with cm-1)
N_A = 6.02214076e23              # 1 / mol
R_GAS = 8.314462618              # J / (mol K)
AMU = 1.66053906660e-27          # kg
HARTREE_J = 4.3597447222071e-18  # J
ATM_PA = 101325.0                # Pa
H2KCAL = 627.509474              # Hartree -> kcal/mol

# Grimme's fixed average molecular moment of inertia for the free-rotor limit
B_AV = 1.0e-44                   # kg m^2
QRRHO_ALPHA = 4                  # damping-function exponent

# Most-abundant-isotope masses, matching Q-Chem's default thermochemistry.
MASSES = {
    "H": 1.00782503, "D": 2.01410178, "He": 4.00260325,
    "Li": 7.0160034, "Be": 9.0121831, "B": 11.0093054,
    "C": 12.0, "N": 14.0030740, "O": 15.9949146, "F": 18.9984032,
    "Ne": 19.9924402, "Na": 22.9897693, "Mg": 23.9850417,
    "Al": 26.9815385, "Si": 27.9769265, "P": 30.9737620,
    "S": 31.9720712, "Cl": 34.9688527, "Ar": 39.9623831,
    "K": 38.9637065, "Ca": 39.9625909, "Br": 78.9183376,
    "I": 126.9044719,
}

_FREQ_RE = re.compile(r"^\s*Frequency:\s+(.+)$")


def parse_floats(s):
    """Whitespace-separated floats, tolerating Fortran 'D' exponents."""
    vals = []
    for tok in s.split():
        try:
            vals.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    return vals


def read_frequencies(path):
    """All harmonic frequencies (cm-1) printed by a Q-Chem FREQ job.

    Q-Chem prints 3N-6 vibrational modes in blocks of three, having already
    projected out translations and rotations, so every value here is a genuine
    vibrational mode. Imaginary modes appear as negative numbers.
    """
    freqs = []
    with open(path, errors="ignore") as fh:
        for line in fh:
            m = _FREQ_RE.match(line)
            if m:
                freqs.extend(parse_floats(m.group(1)))
    return freqs


def read_geometry(path):
    """(symbols, coords_angstrom) from the last standard orientation block.

    Same block and same walk as 4b-qchem-dedup.py, kept deliberately identical
    so the two tools cannot disagree about which geometry a conformer has.
    """
    txt = open(path, errors="ignore").read()
    if "Standard Nuclear Orientation (Angstroms)" not in txt:
        return [], []
    blk = txt.split("Standard Nuclear Orientation (Angstroms)")[-1].splitlines()
    els, xyz, dash, on = [], [], 0, False
    for ln in blk:
        s = ln.strip()
        if s and set(s) == {"-"}:
            dash += 1
            on = (dash == 1)
            if dash == 2:
                break
            continue
        if on:
            p = ln.split()
            if len(p) >= 5 and p[0].isdigit():
                els.append(p[1])
                xyz.append([float(p[2]), float(p[3]), float(p[4])])
    return els, xyz


def job_finished(path):
    """True if Q-Chem reported a clean exit."""
    txt = open(path, errors="ignore").read()
    return "Thank you very much for using Q-Chem" in txt


def principal_moments(els, xyz):
    """Principal moments of inertia in kg m^2, ascending."""
    masses = []
    for e in els:
        sym = e.strip().capitalize()
        if sym not in MASSES:
            raise KeyError("no mass for element '%s'" % e)
        masses.append(MASSES[sym] * AMU)
    coords = [[v * 1.0e-10 for v in row] for row in xyz]  # Angstrom -> m

    mtot = sum(masses)
    com = [sum(masses[i] * coords[i][k] for i in range(len(masses))) / mtot
           for k in range(3)]
    rel = [[coords[i][k] - com[k] for k in range(3)] for i in range(len(masses))]

    # inertia tensor
    ixx = iyy = izz = ixy = ixz = iyz = 0.0
    for i, m in enumerate(masses):
        x, y, z = rel[i]
        ixx += m * (y * y + z * z)
        iyy += m * (x * x + z * z)
        izz += m * (x * x + y * y)
        ixy -= m * x * y
        ixz -= m * x * z
        iyz -= m * y * z
    tensor = [[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]]

    return sorted(jacobi_eigenvalues(tensor)), mtot


def jacobi_eigenvalues(a):
    """Eigenvalues of a real symmetric 3x3 matrix (cyclic Jacobi).

    Hand-rolled rather than numpy so this tool has no hard third-party
    dependency -- it may need to run on a bare cluster login node.

    The matrix is normalized to unit scale before diagonalizing. Inertia
    tensors in SI are ~1e-45 kg m^2, so every convergence test here has to be
    relative; an absolute off-diagonal threshold silently skips all rotations
    at that magnitude and hands back the undiagonalized matrix.
    """
    scale = max(abs(a[i][j]) for i in range(3) for j in range(3))
    if scale <= 0.0:
        return [0.0, 0.0, 0.0]
    m = [[v / scale for v in row] for row in a]

    for _ in range(100):
        off = abs(m[0][1]) + abs(m[0][2]) + abs(m[1][2])
        if off < 1.0e-15:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if abs(m[p][q]) < 1.0e-18:
                continue
            theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q])
            t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
            c = 1.0 / math.sqrt(t * t + 1.0)
            s = t * c
            for k in range(3):
                mkp, mkq = m[k][p], m[k][q]
                m[k][p] = c * mkp - s * mkq
                m[k][q] = s * mkp + c * mkq
            for k in range(3):
                mpk, mqk = m[p][k], m[q][k]
                m[p][k] = c * mpk - s * mqk
                m[q][k] = s * mpk + c * mqk
    return [m[0][0] * scale, m[1][1] * scale, m[2][2] * scale]


def s_translational(mass_kg, temp, pressure_pa):
    """Sackur-Tetrode, J/(mol K). Identical for every conformer, so it cancels
    in the relative G -- computed anyway so the absolute G_corr can be checked
    against Q-Chem's own printed thermochemistry."""
    q = ((2.0 * math.pi * mass_kg * K_B * temp) / (H_PLANCK ** 2)) ** 1.5 \
        * (K_B * temp / pressure_pa)
    return R_GAS * (math.log(q) + 2.5)


def s_rotational(moments, temp, sigma=1):
    """Classical rigid-rotor entropy, J/(mol K). Does NOT cancel between
    conformers -- folded and extended structures have genuinely different
    moments of inertia."""
    ia, ib, ic = moments
    if ia <= 0.0:
        return 0.0
    q = (math.sqrt(math.pi) / sigma) \
        * ((8.0 * math.pi ** 2 * K_B * temp / (H_PLANCK ** 2)) ** 1.5) \
        * math.sqrt(ia * ib * ic)
    return R_GAS * (math.log(q) + 1.5)


def vib_terms(freqs_cm, temp):
    """(ZPE, thermal vibrational energy, harmonic S_vib per mode) in J/mol."""
    zpe = 0.0
    e_th = 0.0
    s_modes = []
    for nu in freqs_cm:
        # J/mol for half a quantum
        zpe += 0.5 * H_PLANCK * C_LIGHT * nu * N_A
        x = H_PLANCK * C_LIGHT * nu / (K_B * temp)
        if x > 700.0:
            # exp() would overflow; the mode is completely frozen out and
            # contributes nothing beyond its ZPE. Append in place rather than
            # skipping, so s_modes stays index-aligned with freqs_cm for qRRHO.
            s_modes.append(0.0)
            continue
        e_th += R_GAS * temp * x / (math.exp(x) - 1.0)
        s_modes.append(R_GAS * (x / (math.exp(x) - 1.0) - math.log(1.0 - math.exp(-x))))
    return zpe, e_th, s_modes


def s_free_rotor(nu_cm, temp):
    """Grimme's free-rotor entropy for one mode, J/(mol K)."""
    nu_si = nu_cm * C_LIGHT                      # s^-1
    mu = H_PLANCK / (8.0 * math.pi ** 2 * nu_si)  # kg m^2
    mu_eff = mu * B_AV / (mu + B_AV)
    return R_GAS * (0.5 + math.log(
        math.sqrt(8.0 * math.pi ** 3 * mu_eff * K_B * temp / (H_PLANCK ** 2))))


def s_vib_qrrho(freqs_cm, s_harmonic, temp, cutoff):
    """Damped interpolation between harmonic and free-rotor entropy."""
    total = 0.0
    for nu, s_ho in zip(freqs_cm, s_harmonic):
        w = 1.0 / (1.0 + (cutoff / nu) ** QRRHO_ALPHA)
        total += w * s_ho + (1.0 - w) * s_free_rotor(nu, temp)
    return total


def thermo_for_conformer(out_path, temp, pressure_pa, cutoff, imag_tol):
    """Returns a result dict, or a dict with 'status' set to a failure reason."""
    cid = os.path.basename(os.path.dirname(os.path.abspath(out_path)))
    res = {"cid": cid, "n_imag": 0, "nu_min": float("nan"),
           "g_rrho": None, "g_qrrho": None, "h_corr": None,
           "ts_rrho": None, "ts_qrrho": None, "status": "ok"}

    if not job_finished(out_path):
        res["status"] = "incomplete"
        return res

    freqs = read_frequencies(out_path)
    if not freqs:
        res["status"] = "no_freqs"
        return res

    els, xyz = read_geometry(out_path)
    if not els:
        res["status"] = "no_geom"
        return res

    expected = 3 * len(els) - 6
    if len(freqs) != expected:
        res["status"] = "nmode_mismatch(%d!=%d)" % (len(freqs), expected)
        return res

    imag = [f for f in freqs if f < 0.0]
    real = [f for f in freqs if f > 0.0]
    res["n_imag"] = len(imag)
    res["nu_min"] = min(real) if real else float("nan")

    if imag:
        worst = min(imag)
        # A large imaginary mode means Stage 4 landed on a saddle point, not a
        # minimum. Small ones are grid/convergence noise and are tolerated.
        res["status"] = "saddle" if abs(worst) > imag_tol else "small_imag"

    try:
        moments, mass_kg = principal_moments(els, xyz)
    except KeyError as exc:
        res["status"] = "mass_error:%s" % exc
        return res

    zpe, e_vib, s_modes = vib_terms(real, temp)
    s_tr = s_translational(mass_kg, temp, pressure_pa)
    s_rot = s_rotational(moments, temp)
    s_vib_h = sum(s_modes)
    s_vib_q = s_vib_qrrho(real, s_modes, temp, cutoff)

    # H = U + RT; U = E_elec + ZPE + E_vib + E_trans(3/2 RT) + E_rot(3/2 RT)
    h_corr = zpe + e_vib + 4.0 * R_GAS * temp          # J/mol
    ts_rrho = temp * (s_tr + s_rot + s_vib_h)          # J/mol
    ts_qrrho = temp * (s_tr + s_rot + s_vib_q)         # J/mol

    j2ha = 1.0 / (HARTREE_J * N_A)
    res["h_corr"] = h_corr * j2ha
    res["ts_rrho"] = ts_rrho * j2ha
    res["ts_qrrho"] = ts_qrrho * j2ha
    res["g_rrho"] = (h_corr - ts_rrho) * j2ha
    res["g_qrrho"] = (h_corr - ts_qrrho) * j2ha
    return res


def process_tag(tag, args):
    freq_dir = os.path.join(tag, args.freq_dir)
    if not os.path.isdir(freq_dir):
        sys.stderr.write("Warning: [%s] %s/ not found -- skipping\n" % (tag, freq_dir))
        return

    outs = sorted(glob.glob(os.path.join(freq_dir, "*", "*.out")))
    if not outs:
        sys.stderr.write("Warning: [%s] no frequency outputs under %s/\n" % (tag, freq_dir))
        return

    results = []
    for path in outs:
        results.append(thermo_for_conformer(
            path, args.temp, args.pressure * ATM_PA, args.qrrho_cutoff, args.imag_tol))

    good = [r for r in results if r["g_rrho"] is not None]
    if not good:
        sys.stderr.write("Warning: [%s] no usable thermochemistry\n" % tag)
        for r in results:
            sys.stderr.write("    %-28s %s\n" % (r["cid"], r["status"]))
        return

    # relative numbers, for the on-screen report only
    ref_r = min(r["g_rrho"] for r in good)
    ref_q = min(r["g_qrrho"] for r in good)

    sys.stderr.write("\n[%s] thermochemistry at %.2f K (%d conformers)\n"
                     % (tag, args.temp, len(results)))
    sys.stderr.write("  %-28s %10s %10s %6s %9s  %s\n"
                     % ("CID", "dG_RRHO", "dG_qRRHO", "n_imag", "nu_min", "status"))
    for r in sorted(results, key=lambda x: (x["g_qrrho"] is None,
                                            x["g_qrrho"] if x["g_qrrho"] is not None else 0.0)):
        if r["g_rrho"] is None:
            sys.stderr.write("  %-28s %10s %10s %6d %9s  %s\n"
                             % (r["cid"], "-", "-", r["n_imag"], "-", r["status"]))
            continue
        sys.stderr.write("  %-28s %10.3f %10.3f %6d %9.1f  %s\n"
                         % (r["cid"],
                            (r["g_rrho"] - ref_r) * H2KCAL,
                            (r["g_qrrho"] - ref_q) * H2KCAL,
                            r["n_imag"], r["nu_min"], r["status"]))

    n_saddle = sum(1 for r in results if r["status"] == "saddle")
    if n_saddle:
        sys.stderr.write("  !! %d conformer(s) have an imaginary mode above %.0f cm-1 --\n"
                         "     these are saddle points, not minima. Stage 5 --gibbs skips\n"
                         "     them unless you pass --keep-saddles.\n" % (n_saddle, args.imag_tol))

    # how much the low modes actually matter, on the ranking that will be used
    spread = max(abs((r["g_rrho"] - ref_r) - (r["g_qrrho"] - ref_q)) * H2KCAL
                 for r in good)
    sys.stderr.write("  RRHO vs qRRHO: max disagreement in relative G = %.2f kcal/mol\n"
                     % spread)

    if args.dry_run:
        sys.stderr.write("  (dry run) no table written\n")
        return

    table = os.path.join(freq_dir, "%s_thermo.dat" % tag)
    with open(table, "w") as fh:
        fh.write("# Gibbs corrections from 4c-qchem-thermo.py\n")
        fh.write("# temperature_K %.4f\n" % args.temp)
        fh.write("# pressure_atm %.4f\n" % args.pressure)
        fh.write("# qrrho_cutoff_cm %.2f\n" % args.qrrho_cutoff)
        fh.write("# imag_tol_cm %.2f\n" % args.imag_tol)
        fh.write("# columns: CID G_corr_rrho(Ha) G_corr_qrrho(Ha) "
                 "H_corr(Ha) TS_rrho(Ha) TS_qrrho(Ha) n_imag nu_min(cm-1) status\n")
        for r in sorted(results, key=lambda x: x["cid"]):
            if r["g_rrho"] is None:
                fh.write("%-28s %16s %16s %16s %16s %16s %4d %10s %s\n"
                         % (r["cid"], "nan", "nan", "nan", "nan", "nan",
                            r["n_imag"], "nan", r["status"]))
                continue
            fh.write("%-28s %16.9f %16.9f %16.9f %16.9f %16.9f %4d %10.2f %s\n"
                     % (r["cid"], r["g_rrho"], r["g_qrrho"], r["h_corr"],
                        r["ts_rrho"], r["ts_qrrho"], r["n_imag"], r["nu_min"],
                        r["status"]))
    sys.stderr.write("  wrote %s\n" % table)


def main():
    ap = argparse.ArgumentParser(
        description="Gibbs corrections from Q-Chem frequency jobs (Stage 4c).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tags", nargs="*", help="molecule TAG(s)")
    ap.add_argument("--temp", type=float, default=298.15, help="temperature, K")
    ap.add_argument("--pressure", type=float, default=1.0, help="pressure, atm")
    ap.add_argument("--qrrho-cutoff", type=float, default=100.0,
                    help="qRRHO damping centre, cm-1")
    ap.add_argument("--imag-tol", type=float, default=50.0,
                    help="|imaginary| below this is treated as numerical noise")
    ap.add_argument("--freq-dir", default="03b_freq", help="Stage-4c output subdir")
    ap.add_argument("--list", dest="list_file", help="text file of TAGs")
    ap.add_argument("--dry-run", action="store_true", help="report only")
    args = ap.parse_args()

    tags = list(args.tags)
    if args.list_file:
        if tags:
            ap.error("positional TAGs not allowed with --list")
        with open(args.list_file) as fh:
            for ln in fh:
                ln = ln.strip()
                if ln and not ln.startswith("#") and not ln.startswith(";"):
                    tags.append(ln)
    if not tags:
        ap.error("provide at least one TAG (or --list FILE)")

    if args.temp <= 0.0:
        ap.error("--temp must be positive")
    if args.qrrho_cutoff <= 0.0:
        ap.error("--qrrho-cutoff must be positive")

    for tag in tags:
        process_tag(tag, args)


if __name__ == "__main__":
    main()
