#!/usr/bin/env python3
# ============================================================================
# or_vcd_ir_tools.py - Parse ORCA *AnFreq* logs and build Boltzmann‑weighted
# IR absorption and VCD spectra (stage 6‑alt post‑processing).
#
# OVERVIEW
#   - Accepts one or more *.log* files, glob patterns or directories.  
#   - Optionally reads a Boltzmann‑weight table produced by
#     *5‑orca‑boltzmann‑weight.sh* and scales intensities accordingly.  
#   - Writes two CSV tables (`*_ir.csv`, `*_vcd.csv`) and generates Gaussian‑
#     broadened spectra as PNG + PDF.  
#
# Quick examples
#   # Plot a single conformer with default settings
#   or_vcd_ir_tools.py conf01.log --stick
#
#   # Recursively gather logs, apply Boltzmann weights and invert IR axis
#   or_vcd_ir_tools.py --bw l-glu/bw_results/l-glu_energies.dat \
#       --prefix l-glu --invert_ir **/solvent_opt/**/*.log
#
#   # Custom FWHM and axis limits
#   or_vcd_ir_tools.py *.log --ir_fwhm 12 --vcd_fwhm 8 \
#       --xlim 1000 3500 --ir_ylim -0.01 0.08 --vcd_ylim -4 4
#
# Flags (see --help for full list)
#   --bw PATH          Boltzmann weight file
#   --prefix STR       Prefix for all outputs               [vib]
#   --stick            Overlay stick spectra
#   --invert_ir        Plot IR absorption peaks downward
#   --ir_fwhm CM-1     Gaussian FWHM for IR   (cm‑1)        [10]
#   --vcd_fwhm CM-1    Gaussian FWHM for VCD  (cm‑1)        [6]
#
# ============================================================================

from __future__ import annotations

import argparse
import glob
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------- #
#                               CONSTANTS                                #
# ---------------------------------------------------------------------- #

_IR_RE = re.compile(r"^IR SPECTRUM")
_VCD_RE = re.compile(r"^VCD SPECTRUM")

SQRT2PI = np.sqrt(2.0 * np.pi)

# ---------------------------------------------------------------------- #
#                           PARSING UTILITIES                            #
# ---------------------------------------------------------------------- #


def _iter_blocks(lines: Sequence[str], header_re: re.Pattern[str]) -> List[List[str]]:
    """Return list of blocks that start with *header_re* until next header."""
    blocks, block = [], []
    inside = False
    for line in lines:
        if header_re.match(line.lstrip()):
            if block:
                blocks.append(block)
                block = []
            inside = True
            continue
        if inside:
            if _IR_RE.match(line.lstrip()) or _VCD_RE.match(line.lstrip()):
                blocks.append(block)
                block = []
                inside = False
            else:
                block.append(line)
    if block:
        blocks.append(block)
    return blocks


def _parse_block(block: List[str], signed: bool) -> List[Dict[str, float]]:
    """Parse one IR or VCD block into rows with ν (cm‑1) and intensity."""
    # skip until numeric lines begin
    idx = 0
    while idx < len(block) and not re.match(r"\s*\d+", block[idx]):
        idx += 1
    rows: List[Dict[str, float]] = []
    for line in block[idx:]:
        if not re.match(r"\s*\d+", line):
            break
        toks = line.split()
        if len(toks) < 3:
            continue
        nu_cm = float(toks[1])
        inten = float(toks[2] if signed else toks[3] if len(toks) >= 4 else toks[2])
        rows.append({"nu_cm": nu_cm, "intensity": inten})
    return rows


def parse_orca_vib(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse ORCA -AnFreq- log.

    Returns
    -------
    ir_df   : DataFrame with ν (cm‑1) and intensity
    vcd_df  : DataFrame with ν (cm‑1) and rotatory strength
    """
    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    ir_rows: List[Dict[str, float]] = []
    for block in _iter_blocks(lines, _IR_RE):
        ir_rows.extend(_parse_block(block, signed=False))

    vcd_rows: List[Dict[str, float]] = []
    for block in _iter_blocks(lines, _VCD_RE):
        vcd_rows.extend(_parse_block(block, signed=True))

    return pd.DataFrame(ir_rows), pd.DataFrame(vcd_rows)


# ---------------------------------------------------------------------- #
#                        BROADENING & PLOTTING                           #
# ---------------------------------------------------------------------- #


def _sigma_from_fwhm(fwhm: float) -> float:
    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def broaden(df: pd.DataFrame, *, sigma: float, nu_grid: np.ndarray) -> np.ndarray:
    """
    Convolve stick spectrum with a Gaussian of width *sigma* (cm‑1).

    The Gaussians are *not* normalised by 1/(σ√2π) so peak height equals stick.
    """
    if sigma <= 0:
        raise ValueError("σ must be positive.")
    curve = np.zeros_like(nu_grid)
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
):
    for _, row in df.iterrows():
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
        _stick_spectrum(ax, df, signed=signed, invert=invert)
    if signed:
        ax.axhline(0, color="grey", lw=0.8)

    ax.set_xlim(nu_max, nu_min)  # IR axis decreases
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
        prog="or_vcd_ir_tools.py",
        description="Parse ORCA AnFreq logs and plot Boltzmann‑weighted IR & VCD spectra.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("logs", nargs="+", help="Log files, glob patterns, or directories")
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
    parser.add_argument("--ir_fwhm", type=float, default=10.0, help="IR FWHM / cm‑1")
    parser.add_argument("--vcd_fwhm", type=float, default=6.0, help="VCD FWHM / cm‑1")
    parser.add_argument(
        "--xlim",
        nargs=2,
        metavar=("MIN", "MAX"),
        type=float,
        help="ν range / cm‑1",
    )
    parser.add_argument("--ir_ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("--vcd_ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("--no_title", action="store_true")
    args = parser.parse_args()

    # -------------------------------------------------------------- #
    #   Resolve output prefix (--outdir takes effect if --prefix absent) #
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

    # -------------------------------------------------------------- #
    #   Collect log paths                                            #
    # -------------------------------------------------------------- #
    log_files: List[str] = []
    for token in args.logs:
        if os.path.isdir(token):
            log_files.extend(glob.glob(os.path.join(token, "**", "*.log"), recursive=True))
        else:
            log_files.extend(glob.glob(token))
    if not log_files:
        raise SystemExit("No log files found.")

    # -------------------------------------------------------------- #
    #   Load Boltzmann weights                                       #
    # -------------------------------------------------------------- #
    weights: Dict[str, float] = {Path(f).stem: 1.0 for f in log_files}
    if args.bw_file:
        weights.update(load_weights(args.bw_file))

    # -------------------------------------------------------------- #
    #   Aggregate transitions                                        #
    # -------------------------------------------------------------- #
    ir_all, vcd_all = [], []
    for f in log_files:
        ir_df, vcd_df = parse_orca_vib(f)
        cid = Path(f).stem
        w = weights.get(cid, 1.0)
        ir_df["intensity"] *= w
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
        f"{len(log_files)} log(s) → {len(ir_df)} IR and {len(vcd_df)} VCD transitions "
        f"(σ_ir={sigma_ir:.1f} cm⁻¹, σ_vcd={sigma_vcd:.1f} cm⁻¹)"
    )


if __name__ == "__main__":
    _cli()
