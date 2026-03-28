#!/usr/bin/env python3
# ============================================================================
# qc_vcd_ir_tools.py - Parse Q-Chem FREQ .out files and build
# Boltzmann-weighted IR absorption and VCD spectra.
#
# OVERVIEW
#   - Accepts one or more *.out files, glob patterns or directories.
#   - Optionally reads a Boltzmann-weight table produced by the stage-5
#     script and scales intensities accordingly.
#   - Writes two CSV tables (*_ir.csv, *_vcd.csv) and generates Gaussian-
#     broadened spectra as PNG + PDF.
#
# Q-Chem output parsing
#   The script reads the VIBRATIONAL ANALYSIS section produced by
#   JOB_TYPE FREQ with VCD TRUE. Modes are reported in blocks of N columns:
#
#     Mode:            1      2      3  ...
#     Frequency:    XXX.X  XXX.X  XXX.X
#     Force Cnst:   X.XXX  X.XXX  X.XXX
#     Red. Mass:    X.XXX  X.XXX  X.XXX
#     IR Activ:     X.XXX  X.XXX  X.XXX      (km/mol)
#     Raman Activ:  X.XXX  X.XXX  X.XXX
#     Rot. Str.(L): X.XXX  X.XXX  X.XXX      (10^-44 esu^2 cm^2)
#     Rot. Str.(V): X.XXX  X.XXX  X.XXX
#
#   Only modes with positive frequencies (real vibrations) are kept.
#   VCD rotational strengths: length gauge preferred; velocity gauge used
#   as fallback if only one gauge is present.
#
# Quick examples
#   # Plot a single conformer with default settings
#   qc_vcd_ir_tools.py conf01.out --stick
#
#   # Recursively gather outputs, apply Boltzmann weights and invert IR axis
#   qc_vcd_ir_tools.py --bw l-glu/04_boltzmann/l-glu_energies.dat \
#       --prefix l-glu --invert_ir l-glu/06_vcd/**/*.out
#
#   # Custom FWHM and axis limits
#   qc_vcd_ir_tools.py *.out --ir_fwhm 12 --vcd_fwhm 8 \
#       --xlim 1000 3500 --ir_ylim -0.01 0.08 --vcd_ylim -4 4
#
# Flags (see --help for full list)
#   --bw PATH          Boltzmann weight file
#   --prefix STR       Prefix for all outputs               [vib]
#   --stick            Overlay stick spectra
#   --invert_ir        Plot IR absorption peaks downward
#   --ir_fwhm CM-1     Gaussian FWHM for IR   (cm-1)        [10]
#   --vcd_fwhm CM-1    Gaussian FWHM for VCD  (cm-1)        [6]
#
# ============================================================================

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------- #
#                               CONSTANTS                                #
# ---------------------------------------------------------------------- #

SQRT2PI = np.sqrt(2.0 * np.pi)

_VIB_SECTION_RE = re.compile(r"VIBRATIONAL ANALYSIS", re.IGNORECASE)
_MODE_RE = re.compile(r"^\s*Mode:\s+(\d+(?:\s+\d+)*)\s*$")
_FREQ_RE = re.compile(r"^\s*Frequency:\s+(.+)$")
_IR_RE = re.compile(r"^\s*IR Activ(?:ities)?(?:\s*\([^)]*\))?:\s+(.+)$", re.IGNORECASE)
# Length gauge: "Rot. Str.(L)", "Rot. Str. (L)", "Rot. Str. (Len)"
_ROT_L_RE = re.compile(
    r"^\s*Rot(?:atory)?\.?\s*Str(?:ength)?\.?\s*\(L(?:en(?:gth)?)?\):\s+(.+)$",
    re.IGNORECASE,
)
# Velocity gauge: "Rot. Str.(V)", "Rot. Str. (Vel)"
_ROT_V_RE = re.compile(
    r"^\s*Rot(?:atory)?\.?\s*Str(?:ength)?\.?\s*\(V(?:el(?:ocity)?)?\):\s+(.+)$",
    re.IGNORECASE,
)
# Unlabelled / single-gauge: "Rot. Str.:" without a gauge qualifier
_ROT_ANY_RE = re.compile(
    r"^\s*Rot(?:atory)?\.?\s*Str(?:ength)?\.?:\s+(.+)$",
    re.IGNORECASE,
)


def _parse_floats(s: str) -> List[float]:
    """Parse a whitespace-separated string of floats (handle 'D' exponents)."""
    vals = []
    for tok in s.split():
        try:
            vals.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    return vals


# ---------------------------------------------------------------------- #
#                           PARSING UTILITIES                            #
# ---------------------------------------------------------------------- #


def parse_qchem_vib(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a Q-Chem FREQ (+ VCD TRUE) output file.

    Returns
    -------
    ir_df  : DataFrame with columns nu_cm, intensity  (IR km/mol)
    vcd_df : DataFrame with columns nu_cm, intensity  (rotatory strength,
             length gauge preferred, velocity as fallback)
    """
    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    n = len(lines)

    # Find the VIBRATIONAL ANALYSIS section
    vib_start: Optional[int] = None
    for i, line in enumerate(lines):
        if _VIB_SECTION_RE.search(line):
            vib_start = i
            break

    if vib_start is None:
        return pd.DataFrame(columns=["nu_cm", "intensity"]), pd.DataFrame(
            columns=["nu_cm", "intensity"]
        )

    # ------------------------------------------------------------------ #
    # Split the VIBRATIONAL ANALYSIS section into per-block groups.       #
    # Each group starts with a "Mode:" line.                              #
    # ------------------------------------------------------------------ #
    # Collect (frequencies, ir_activs, rot_len, rot_vel) per block
    ir_rows: List[Dict[str, float]] = []
    vcd_rows: List[Dict[str, float]] = []

    # State within a block
    freqs: List[float] = []
    ir_activ: List[float] = []
    rot_len: List[float] = []
    rot_vel: List[float] = []
    rot_any: List[float] = []  # for single-gauge output without L/V label
    in_block = False

    def _flush_block():
        nonlocal freqs, ir_activ, rot_len, rot_vel, rot_any
        if not freqs:
            return
        # Prefer length gauge; fall back to velocity; fall back to unlabeled
        rot = rot_len if rot_len else rot_vel if rot_vel else rot_any
        for k, freq in enumerate(freqs):
            if freq <= 0.0:  # skip imaginary and translation/rotation (~0)
                continue
            if k < len(ir_activ):
                ir_rows.append({"nu_cm": freq, "intensity": ir_activ[k]})
            if k < len(rot):
                vcd_rows.append({"nu_cm": freq, "intensity": rot[k]})
        freqs = []
        ir_activ = []
        rot_len = []
        rot_vel = []
        rot_any = []

    for i in range(vib_start, n):
        line = lines[i]

        # A new Mode: line starts a new block
        m_mode = _MODE_RE.match(line)
        if m_mode:
            _flush_block()
            in_block = True
            continue

        if not in_block:
            continue

        # Frequency row
        m = _FREQ_RE.match(line)
        if m:
            freqs = _parse_floats(m.group(1))
            continue

        # IR activity row
        m = _IR_RE.match(line)
        if m:
            ir_activ = _parse_floats(m.group(1))
            continue

        # VCD rotatory strength -- length gauge
        m = _ROT_L_RE.match(line)
        if m:
            rot_len = _parse_floats(m.group(1))
            continue

        # VCD rotatory strength -- velocity gauge
        m = _ROT_V_RE.match(line)
        if m:
            rot_vel = _parse_floats(m.group(1))
            continue

        # VCD rotatory strength -- unlabeled (single gauge)
        # Only match if neither length nor velocity gauge seen yet
        if not rot_len and not rot_vel:
            m = _ROT_ANY_RE.match(line)
            if m:
                rot_any = _parse_floats(m.group(1))
                continue

        # A section header that is clearly outside VIBRATIONAL ANALYSIS
        if re.match(r"^\s*\*{10,}", line) and i > vib_start + 3:
            # Check if this is a NEW section header (e.g. Thermochemistry)
            # rather than the decorative banner that opens VIBRATIONAL ANALYSIS
            peek = "".join(lines[i : i + 3])
            if "THERMOCHEMICAL" in peek.upper() or "ZERO POINT" in peek.upper():
                break

    _flush_block()  # don't forget the last block

    return (
        pd.DataFrame(ir_rows, columns=["nu_cm", "intensity"]) if ir_rows
        else pd.DataFrame(columns=["nu_cm", "intensity"]),
        pd.DataFrame(vcd_rows, columns=["nu_cm", "intensity"]) if vcd_rows
        else pd.DataFrame(columns=["nu_cm", "intensity"]),
    )


# ---------------------------------------------------------------------- #
#                        BROADENING & PLOTTING                           #
# ---------------------------------------------------------------------- #


def _sigma_from_fwhm(fwhm: float) -> float:
    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def broaden(df: pd.DataFrame, *, sigma: float, nu_grid: np.ndarray) -> np.ndarray:
    """
    Convolve stick spectrum with a Gaussian of width *sigma* (cm-1).

    The Gaussians are normalized by 1/(sigma*sqrt(2*pi)) so peak height equals stick
    intensity.
    """
    if sigma <= 0:
        raise ValueError("sigma must be positive.")
    curve = np.zeros_like(nu_grid, dtype=float)
    pref = 1.0 / (sigma * SQRT2PI)
    for nu0, inten in zip(df["nu_cm"], df["intensity"]):
        curve += inten * pref * np.exp(-0.5 * ((nu_grid - nu0) / sigma) ** 2)
    return curve


def _stick_spectrum(
    ax,
    df: pd.DataFrame,
    *,
    signed: bool,
    invert: bool = False,
    nu_min: float = 0.0,
    nu_max: float = float("inf"),
):
    lo, hi = min(nu_min, nu_max), max(nu_min, nu_max)
    for _, row in df.iterrows():
        if not (lo <= row["nu_cm"] <= hi):
            continue
        height = -row["intensity"] if invert and not signed else row["intensity"]
        color = (
            "royalblue"
            if signed and height >= 0
            else "crimson"
            if signed
            else "grey"
        )
        ax.vlines(row["nu_cm"], 0, height, color=color, alpha=0.4, lw=0.8)


def _plot(
    df: pd.DataFrame,
    *,
    sigma: float,
    nu_min: float,
    nu_max: float,
    png: str,
    pdf: str,
    signed: bool,
    sticks: bool,
    ylim: Optional[Tuple[float, float]],
    title: Optional[str],
    invert: bool,
):
    if df.empty:
        print(f"[SKIP] {png}: no transitions")
        return

    nu_grid = np.arange(nu_min, nu_max, 1.0)
    curve = broaden(df, sigma=sigma, nu_grid=nu_grid)
    if invert and not signed:
        curve = -curve

    fig, ax = plt.subplots(figsize=(4.75, 3.25))
    ax.plot(nu_grid, curve, color="black", lw=1.1)
    if sticks:
        _stick_spectrum(ax, df, signed=signed, invert=invert, nu_min=nu_min, nu_max=nu_max)
    if signed:
        ax.axhline(0, color="grey", lw=0.8)

    ax.set_xlim(nu_max, nu_min)  # IR axis decreases left-to-right
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(r"Wavenumber / cm$^{-1}$")
    ax.set_ylabel("Rotatory strength (arb.)" if signed else "Intensity (arb.)")
    if title:
        ax.set_title(title, fontsize=10)

    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)


# ---------------------------------------------------------------------- #
#                        BOLTZMANN WEIGHTS                               #
# ---------------------------------------------------------------------- #


def load_weights(path: str) -> Dict[str, float]:
    df = pd.read_csv(
        path,
        comment="#",
        sep=r"\s+",
        engine="python",
        header=None,
        names=["cid", "E", "dE", "p"],
    )
    return dict(zip(df["cid"].astype(str), df["p"].astype(float)))


# ---------------------------------------------------------------------- #
#                                 CLI                                    #
# ---------------------------------------------------------------------- #


def _cli() -> None:
    parser = argparse.ArgumentParser(
        prog="qc_vcd_ir_tools.py",
        description="Parse Q-Chem FREQ .out files and plot Boltzmann-weighted IR & VCD spectra.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("logs", nargs="+", help="Q-Chem .out files, glob patterns, or directories")
    parser.add_argument("--bw", dest="bw_file", help="Boltzmann weight file")
    parser.add_argument(
        "--outdir",
        metavar="TAG",
        help="Molecule directory (e.g. methyloxirane). Outputs go to <TAG>/06_spectra/ "
             "with prefix <TAG>/06_spectra/<TAG>. Overridden by --prefix.",
    )
    parser.add_argument("--prefix", default=None, help="Explicit output prefix [vib]")
    parser.add_argument("--stick", action="store_true", help="Overlay stick spectra")
    parser.add_argument("--invert_ir", action="store_true", help="Plot IR peaks downward")
    parser.add_argument("--ir_fwhm", type=float, default=10.0, help="IR FWHM / cm-1")
    parser.add_argument("--vcd_fwhm", type=float, default=6.0, help="VCD FWHM / cm-1")
    parser.add_argument(
        "--xlim",
        nargs=2,
        metavar=("MIN", "MAX"),
        type=float,
        help="nu range / cm-1",
    )
    parser.add_argument("--ir_ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("--vcd_ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("--no_title", action="store_true")
    args = parser.parse_args()

    # -------------------------------------------------------------- #
    #   Resolve output prefix                                         #
    # -------------------------------------------------------------- #
    if args.prefix is None:
        if args.outdir:
            tag = Path(args.outdir).name
            out_dir = Path(args.outdir) / "06_spectra"
            out_dir.mkdir(parents=True, exist_ok=True)
            prefix = str(out_dir / tag)
        else:
            prefix = "vib"
    else:
        prefix = args.prefix
        Path(prefix).parent.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------- #
    #   Collect output paths                                          #
    # -------------------------------------------------------------- #
    log_files: List[str] = []
    for token in args.logs:
        if os.path.isdir(token):
            log_files.extend(
                glob.glob(os.path.join(token, "**", "*.out"), recursive=True)
            )
        else:
            log_files.extend(glob.glob(token))
    # Drop SLURM output files (slurm-XXXXXXXX.out) -- not Q-Chem outputs
    log_files = [f for f in log_files if not Path(f).stem.startswith("slurm-")]
    if not log_files:
        raise SystemExit("No .out files found.")

    # -------------------------------------------------------------- #
    #   Load Boltzmann weights                                        #
    # -------------------------------------------------------------- #
    weights: Dict[str, float] = {Path(f).stem: 1.0 for f in log_files}
    if args.bw_file:
        weights.update(load_weights(args.bw_file))

    # -------------------------------------------------------------- #
    #   Aggregate transitions                                         #
    # -------------------------------------------------------------- #
    ir_all, vcd_all = [], []
    for f in log_files:
        ir_df, vcd_df = parse_qchem_vib(f)
        cid = Path(f).stem
        w = weights.get(cid, 1.0)
        if not ir_df.empty:
            ir_df["intensity"] *= w
        if not vcd_df.empty:
            vcd_df["intensity"] *= w
        ir_all.append(ir_df)
        vcd_all.append(vcd_df)

    ir_df = pd.concat(ir_all, ignore_index=True)
    vcd_df = pd.concat(vcd_all, ignore_index=True)

    ir_df.to_csv(f"{prefix}_ir.csv", index=False)
    vcd_df.to_csv(f"{prefix}_vcd.csv", index=False)

    nu_min, nu_max = (0.0, 4000.0) if not args.xlim else tuple(map(float, args.xlim))
    sigma_ir = _sigma_from_fwhm(args.ir_fwhm)
    sigma_vcd = _sigma_from_fwhm(args.vcd_fwhm)

    base = Path(prefix).name.replace("_", " ")
    title_ir = None if args.no_title else f"{base} IR Spectrum"
    title_vcd = None if args.no_title else f"{base} VCD Spectrum"

    _plot(
        ir_df,
        sigma=sigma_ir,
        nu_min=nu_min,
        nu_max=nu_max,
        png=f"{prefix}_ir.png",
        pdf=f"{prefix}_ir.pdf",
        signed=False,
        sticks=args.stick,
        ylim=tuple(args.ir_ylim) if args.ir_ylim else None,
        title=title_ir,
        invert=args.invert_ir,
    )

    _plot(
        vcd_df,
        sigma=sigma_vcd,
        nu_min=nu_min,
        nu_max=nu_max,
        png=f"{prefix}_vcd.png",
        pdf=f"{prefix}_vcd.pdf",
        signed=True,
        sticks=args.stick,
        ylim=tuple(args.vcd_ylim) if args.vcd_ylim else None,
        title=title_vcd,
        invert=False,
    )

    print(
        f"{len(log_files)} file(s) -> {len(ir_df)} IR and {len(vcd_df)} VCD transitions "
        f"(sigma_ir={sigma_ir:.1f} cm-1, sigma_vcd={sigma_vcd:.1f} cm-1)"
    )


if __name__ == "__main__":
    _cli()
