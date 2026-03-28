#!/usr/bin/env python3
# ============================================================================
# qc_ecd_uvvis_tools.py - Parse Q-Chem TD-DFT .out files and build
# Boltzmann-weighted, Gaussian-broadened UV/Vis and ECD spectra.
#
# OVERVIEW
#   1. Collects *.out files (positional patterns, wild-cards, or directories),
#   2. Extracts stick data from the full "TDDFT Excitation Energies" block
#      (not the TDA pre-step) and the "Electronic Circular Dichroism (ECD)
#      Spectrum" table,
#   3. Applies optional Boltzmann weights (stage-5 output),
#   4. Writes CSV tables and plots broadened spectra (PNG + PDF).
#
# Quick examples
#   # simple: average of outputs in the current folder
#   qc_ecd_uvvis_tools.py *.out
#
#   # recurse through folders, supply Boltzmann weights, 0.35 eV Gaussian FWHM
#   qc_ecd_uvvis_tools.py --bw aspirin/04_boltzmann/aspirin_energies.dat \
#       --uv_fwhm 0.35 --ecd_fwhm 0.25 ./aspirin/05_tddft/**
#
#   # publish-quality plot with sticks and custom limits
#   qc_ecd_uvvis_tools.py *.out --stick --xlim 190 350 --uv_ylim 0 0.05 \
#       --ecd_ylim -2 2 --prefix aspirin_spectra
#
# Flags (see --help for full list)
#   --bw PATH          Boltzmann weight file (stage-5 output, energies.dat)
#   --prefix STR       Prefix for all outputs              [spectra]
#   --uv_fwhm EV       Gaussian FWHM for UV/Vis (eV)       [0.35]
#   --ecd_fwhm EV      Gaussian FWHM for ECD   (eV)        [0.25]
#   --stick            Overlay stick spectrum
#   --flip_x           Plot lambda increasing left -> right
#   --scale FACTOR     Multiply intensities after weighting [1.0]
#   --ecd_gauge STR    ECD rotatory strength: length|velocity [length]
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

HC_OVER_EV_NM = 1239.84193  # Planck*c in eV*nm
SQRT2PI = np.sqrt(2.0 * np.pi)

# Matches the full-TDDFT header line (NOT "TDDFT/TDA Excitation Energies")
_TDDFT_HEADER_RE = re.compile(r"^\s*TDDFT Excitation Energies\s*$")

# Per-state line patterns inside the TDDFT block
_STATE_RE = re.compile(
    r"Excited state\s+\d+:\s+excitation energy \(eV\)\s*=\s*([+-]?\d+\.\d+)"
)
_MULT_RE = re.compile(r"Multiplicity:\s*(\w+)")
_STRENGTH_RE = re.compile(
    r"Strength\s*:\s*([+-]?\d+\.\d+(?:[eEdD][+-]?\d+)?)"
)

_ECD_HEADER = "Electronic Circular Dichroism (ECD) Spectrum"

# ---------------------------------------------------------------------- #
#                           PARSING UTILITIES                            #
# ---------------------------------------------------------------------- #


def parse_qchem_out(
    path: str, *, ecd_gauge: str = "length"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract stick UV/Vis and ECD data from a Q-Chem TD-DFT .out file.

    UV-Vis intensities are oscillator strengths from the *full* TDDFT block
    (the last "TDDFT Excitation Energies" section; Q-Chem prints a TDA
    pre-step first which is intentionally skipped).

    ECD rotatory strengths are taken from the tabular
    "Electronic Circular Dichroism (ECD) Spectrum" section.

    Parameters
    ----------
    path : str
        Path to the Q-Chem .out file.
    ecd_gauge : str
        'length'   --> use R(length)   column  (default, more reliable)
        'velocity' --> use R(velocity) column

    Returns
    -------
    uv_df : DataFrame with columns energy_eV, wavelength_nm, intensity
    ecd_df : DataFrame with columns energy_eV, wavelength_nm, intensity
    """
    gauge_col = 8 if ecd_gauge == "length" else 9

    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()
    n = len(lines)

    # ------------------------------------------------------------------ #
    #  Step 1: locate the last "TDDFT Excitation Energies" block start    #
    #          and the ECD table start (used as an upper bound).          #
    # ------------------------------------------------------------------ #
    last_tddft_start: Optional[int] = None
    ecd_line: Optional[int] = None

    for i, line in enumerate(lines):
        if _TDDFT_HEADER_RE.match(line):
            last_tddft_start = i
        if ecd_line is None and _ECD_HEADER in line:
            ecd_line = i

    # ------------------------------------------------------------------ #
    #  Step 2: parse UV-Vis (oscillator strengths) from the TDDFT block   #
    # ------------------------------------------------------------------ #
    uv_rows: List[Dict[str, float]] = []

    if last_tddft_start is not None:
        # Stop at the ECD table (or end of file) so we don't pick up stray
        # "Strength" lines from other sections.
        parse_end = ecd_line if ecd_line is not None else n
        cur_energy: Optional[float] = None
        cur_singlet = False

        for i in range(last_tddft_start + 1, parse_end):
            line = lines[i]

            m = _STATE_RE.search(line)
            if m:
                cur_energy = float(m.group(1))
                cur_singlet = False
                continue

            mm = _MULT_RE.search(line)
            if mm:
                cur_singlet = mm.group(1).lower() == "singlet"
                continue

            ms = _STRENGTH_RE.search(line)
            if ms and cur_energy is not None and cur_singlet:
                fosc = float(ms.group(1).replace("D", "E").replace("d", "e"))
                wl = HC_OVER_EV_NM / cur_energy
                uv_rows.append(
                    {"energy_eV": cur_energy, "wavelength_nm": wl, "intensity": fosc}
                )
                cur_energy = None
                cur_singlet = False

    # ------------------------------------------------------------------ #
    #  Step 3: parse ECD rotatory strengths from the tabular section      #
    #                                                                      #
    #  Table layout (5 header lines before data):                         #
    #    ECD Spectrum header line                                          #
    #    ----  (dashes)                                                    #
    #    State  Energy  PX  PY  PZ  MX  MY  MZ  R(len)  R(vel)           #
    #    (eV)   ...                                                        #
    #    ----  (dashes)                                                    #
    #    data rows ...                                                     #
    #    ----  (closing dashes)                                            #
    # ------------------------------------------------------------------ #
    ecd_rows: List[Dict[str, float]] = []

    if ecd_line is not None:
        data_start = ecd_line + 5  # skip 5 header lines
        for i in range(data_start, n):
            raw = lines[i].strip()
            if not raw or raw.startswith("-"):
                break  # closing separator or blank -- end of table
            toks = raw.split()
            if len(toks) <= gauge_col:
                continue
            try:
                e = float(toks[1])
                r = float(toks[gauge_col])
            except (ValueError, IndexError):
                continue
            wl = HC_OVER_EV_NM / e
            ecd_rows.append(
                {"energy_eV": e, "wavelength_nm": wl, "intensity": r}
            )

    _cols = ["energy_eV", "wavelength_nm", "intensity"]
    uv = pd.DataFrame(uv_rows, columns=_cols) if uv_rows else pd.DataFrame(columns=_cols)
    ecd = pd.DataFrame(ecd_rows, columns=_cols) if ecd_rows else pd.DataFrame(columns=_cols)
    return uv, ecd


# ---------------------------------------------------------------------- #
#                      BROADENING & PLOTTING TOOLS                       #
# ---------------------------------------------------------------------- #


def _energy_from_lambda_nm(lam_nm: np.ndarray) -> np.ndarray:
    """Convert wavelength (nm) to energy (eV)."""
    return HC_OVER_EV_NM / lam_nm


def broaden(
    df: pd.DataFrame,
    *,
    sigma_eV: float,
    lam_nm: np.ndarray,
    energy_grid_eV: np.ndarray,
    jacobian: bool,
) -> np.ndarray:
    """
    Return broadened curve on the supplied lambda grid.

    Each Gaussian is normalized with 1/(sigma*sqrt(2*pi)) so its peak equals
    the stick intensity. If *jacobian* is True (ECD), multiply by
    |dE/dlambda| = E/lambda so area is preserved after the change of variable.
    """
    if sigma_eV <= 0:
        raise ValueError("sigma must be positive.")
    curve = np.zeros_like(lam_nm, dtype=float)
    prefactor = 1.0 / (sigma_eV * SQRT2PI)
    for e, inten in zip(df["energy_eV"], df["intensity"]):
        gaussian = np.exp(-0.5 * ((energy_grid_eV - e) / sigma_eV) ** 2)
        if jacobian:
            gaussian *= e / lam_nm
        curve += inten * prefactor * gaussian
    return curve


def _stick_spectrum(
    ax, df: pd.DataFrame, *, is_ecd: bool, scale: float, lam_min: float, lam_max: float
):
    """Draw vertical lines representing individual transitions (clipped to plot range)."""
    lo, hi = min(lam_min, lam_max), max(lam_min, lam_max)
    for _, row in df.iterrows():
        if not (lo <= row["wavelength_nm"] <= hi):
            continue
        color = (
            "royalblue"
            if is_ecd and row["intensity"] >= 0
            else "crimson"
            if is_ecd
            else "grey"
        )
        ax.vlines(
            row["wavelength_nm"],
            0,
            row["intensity"] * scale,
            color=color,
            alpha=0.4,
            lw=0.8,
        )


def _make_plot(
    df: pd.DataFrame,
    *,
    sigma_eV: float,
    lam_min: float,
    lam_max: float,
    outfile_png: str,
    outfile_pdf: str,
    y_limits: Optional[Tuple[float, float]],
    sticks: bool,
    is_ecd: bool,
    title: Optional[str],
    flip_x: bool,
    scale: float,
):
    """Create and save one spectrum plot."""
    if df.empty:
        print(f"[SKIP] {outfile_png}: no transitions")
        return

    lam_nm = np.arange(lam_min, lam_max, 0.1)
    energy_grid_eV = _energy_from_lambda_nm(lam_nm)

    curve = broaden(
        df,
        sigma_eV=sigma_eV,
        lam_nm=lam_nm,
        energy_grid_eV=energy_grid_eV,
        jacobian=is_ecd,
    )
    curve *= scale

    fig, ax = plt.subplots(figsize=(4.75, 3.25))
    ax.plot(lam_nm, curve, color="black", lw=1.1)

    if sticks:
        _stick_spectrum(ax, df, is_ecd=is_ecd, scale=scale, lam_min=lam_min, lam_max=lam_max)

    if is_ecd:
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_ylabel("Rotatory strength (arb.)")
    else:
        ax.set_ylabel("Relative intensity (arb.)")

    ax.set_xlabel("Wavelength / nm")
    if flip_x:
        ax.set_xlim(lam_min, lam_max)
    else:
        ax.set_xlim(lam_max, lam_min)

    if y_limits:
        ax.set_ylim(*y_limits)
    if title:
        ax.set_title(title, fontsize=10)

    fig.tight_layout()
    fig.savefig(outfile_png, dpi=300)
    fig.savefig(outfile_pdf)
    plt.close(fig)


# ---------------------------------------------------------------------- #
#                      BOLTZMANN WEIGHT HANDLING                         #
# ---------------------------------------------------------------------- #


def _load_weights(path: str) -> Dict[str, float]:
    """Read stage-5 energies.dat; return {CID: Boltzmann_probability}."""
    df = pd.read_csv(
        path,
        comment="#",
        sep=r"\s+",
        engine="python",
        header=None,
        names=["cid", "E_Ha", "dE_kcal", "p"],
    )
    return dict(zip(df["cid"].astype(str), df["p"].astype(float)))


# ---------------------------------------------------------------------- #
#                              CLI ENTRY                                 #
# ---------------------------------------------------------------------- #


def _cli() -> None:
    parser = argparse.ArgumentParser(
        prog="qc_ecd_uvvis_tools.py",
        description="Parse Q-Chem TD-DFT .out files and generate UV/Vis & ECD spectra.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "logs",
        nargs="+",
        help="Q-Chem .out files, glob patterns, or directories (recursive).",
    )
    parser.add_argument("--bw", dest="bw_file", metavar="PATH", help="Boltzmann weight file")
    parser.add_argument(
        "--outdir",
        metavar="TAG",
        help="Molecule directory (e.g. pna). Outputs go to <TAG>/06_spectra/ "
             "with prefix <TAG>/06_spectra/<TAG>. Overridden by --prefix.",
    )
    parser.add_argument("--prefix", default=None, help="Explicit output prefix [spectra]")
    parser.add_argument("--no_title", action="store_true", help="Suppress figure titles")
    parser.add_argument("--stick", action="store_true", help="Overlay stick spectra")
    parser.add_argument(
        "--flip_x",
        action="store_true",
        help="Plot short -> long wavelength (e.g. 190 -> 400 nm)",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        metavar="FAC",
        help="Multiply intensities after Boltzmann weighting",
    )
    parser.add_argument(
        "--ecd_gauge",
        choices=["length", "velocity"],
        default="length",
        help="ECD rotatory strength gauge to use",
    )
    parser.add_argument("--uv_fwhm", type=float, default=0.35, help="UV-Vis Gaussian FWHM / eV")
    parser.add_argument("--ecd_fwhm", type=float, default=0.25, help="ECD Gaussian FWHM / eV")
    parser.add_argument(
        "--xlim",
        nargs=2,
        metavar=("MIN", "MAX"),
        type=float,
        help="Wavelength range / nm",
    )
    parser.add_argument("--uv_ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("--ecd_ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    args = parser.parse_args()

    # ------------------------------------------------------------------ #
    #   Resolve output prefix                                             #
    # ------------------------------------------------------------------ #
    if args.prefix is None:
        if args.outdir:
            tag = Path(args.outdir).name
            out_dir = Path(args.outdir) / "06_spectra"
            out_dir.mkdir(parents=True, exist_ok=True)
            prefix = str(out_dir / tag)
        else:
            prefix = "spectra"
    else:
        prefix = args.prefix
        Path(prefix).parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    #   Locate .out files (skip SLURM log files)                         #
    # ------------------------------------------------------------------ #
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

    # ------------------------------------------------------------------ #
    #   Boltzmann weights                                                 #
    # ------------------------------------------------------------------ #
    weights: Dict[str, float] = {Path(f).stem: 1.0 for f in log_files}
    if args.bw_file:
        bw_map = _load_weights(args.bw_file)
        for cid in weights:
            if cid in bw_map:
                weights[cid] = bw_map[cid]

    # ------------------------------------------------------------------ #
    #   Aggregate all transitions                                         #
    # ------------------------------------------------------------------ #
    uv_all, ecd_all = [], []
    for f in log_files:
        uv_df, ecd_df = parse_qchem_out(f, ecd_gauge=args.ecd_gauge)
        cid = Path(f).stem
        weight = weights.get(cid, 1.0) * args.scale
        print(
            f"  {cid}: {len(uv_df)} UV transitions, {len(ecd_df)} ECD transitions"
        )
        if "intensity" not in uv_df.columns:
            print(f"    [WARNING] UV DataFrame missing 'intensity'; columns: {list(uv_df.columns)}")
        else:
            uv_df["intensity"] *= weight
        if "intensity" not in ecd_df.columns:
            print(f"    [WARNING] ECD DataFrame missing 'intensity'; columns: {list(ecd_df.columns)}")
        else:
            ecd_df["intensity"] *= weight
        uv_all.append(uv_df)
        ecd_all.append(ecd_df)

    _out_cols = ["energy_eV", "wavelength_nm", "intensity"]
    _uv_valid = [df for df in uv_all if not df.empty]
    _ecd_valid = [df for df in ecd_all if not df.empty]
    uv_df = (pd.concat(_uv_valid, ignore_index=True)
             if _uv_valid else pd.DataFrame(columns=_out_cols))
    ecd_df = (pd.concat(_ecd_valid, ignore_index=True)
              if _ecd_valid else pd.DataFrame(columns=_out_cols))

    uv_df.to_csv(f"{prefix}_uvvis.csv", index=False)
    ecd_df.to_csv(f"{prefix}_ecd.csv", index=False)

    lam_min, lam_max = (160.0, 400.0) if not args.xlim else sorted(args.xlim)
    fwhm_to_sigma = lambda f: f / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sigma_uv = fwhm_to_sigma(args.uv_fwhm)
    sigma_ecd = fwhm_to_sigma(args.ecd_fwhm)
    title = None if args.no_title else Path(prefix).name.replace("_", " ")

    _make_plot(
        uv_df,
        sigma_eV=sigma_uv,
        lam_min=lam_min,
        lam_max=lam_max,
        outfile_png=f"{prefix}_uvvis.png",
        outfile_pdf=f"{prefix}_uvvis.pdf",
        y_limits=tuple(args.uv_ylim) if args.uv_ylim else None,
        sticks=args.stick,
        is_ecd=False,
        title=title,
        flip_x=args.flip_x,
        scale=1.0,
    )

    _make_plot(
        ecd_df,
        sigma_eV=sigma_ecd,
        lam_min=lam_min,
        lam_max=lam_max,
        outfile_png=f"{prefix}_ecd.png",
        outfile_pdf=f"{prefix}_ecd.pdf",
        y_limits=tuple(args.ecd_ylim) if args.ecd_ylim else None,
        sticks=args.stick,
        is_ecd=True,
        title=title,
        flip_x=args.flip_x,
        scale=1.0,
    )

    print(
        f"{len(log_files)} file(s) -> {len(uv_df)} UV and {len(ecd_df)} ECD transitions "
        f"(sigma_uv={sigma_uv:.3f} eV, sigma_ecd={sigma_ecd:.3f} eV, "
        f"gauge={args.ecd_gauge}, scale={args.scale})"
    )


if __name__ == "__main__":
    _cli()
