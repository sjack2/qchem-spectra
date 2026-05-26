#!/usr/bin/env python3
# ============================================================================
# qc_raman_tools.py - Parse Q-Chem FREQ .out files and build Boltzmann-weighted
# Raman scattering spectra.
#
# OVERVIEW
#   - Accepts one or more *.out files, glob patterns or directories.
#   - Optionally reads a Boltzmann-weight table produced by the stage-5 script
#     and scales activities accordingly. When --bw is given it also acts as the
#     post-dedup keep-list, restricting the file set to its conformers.
#   - Writes a CSV table (*_raman.csv) and a Gaussian-broadened spectrum (PNG+PDF).
#
# Q-Chem output parsing
#   Reads the per-mode "Raman Intens:" lines in the VIBRATIONAL ANALYSIS section
#   produced by JOB_TYPE FREQ with RAMAN TRUE (the default in 6-qchem-vcd.sh).
#   Only modes with positive frequencies (real vibrations) are kept.
#
# Convention
#   Raman shift is plotted INCREASING left-to-right (low -> high), the opposite of
#   the IR/VCD convention, matching how experimental Raman spectra are displayed.
#
# Quick examples
#   # Plot a single conformer with stick overlay
#   qc_raman_tools.py conf01.out --stick
#
#   # Recursively gather outputs and apply Boltzmann weights
#   qc_raman_tools.py --bw l-glu/04_boltzmann/l-glu_energies.dat \
#       --outdir l-glu l-glu/05_vcd
#
# Flags (see --help for full list)
#   --bw PATH          Boltzmann weight / keep-list file
#   --outdir TAG       write <TAG>/06_spectra/<TAG>_raman.*
#   --prefix STR       explicit output prefix (overrides --outdir)      [vib]
#   --stick            overlay stick spectrum
#   --raman_fwhm CM-1  Gaussian FWHM for Raman (cm-1)                   [10]
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

SQRT2PI = np.sqrt(2.0 * np.pi)

_VIB_SECTION_RE = re.compile(r"VIBRATIONAL ANALYSIS", re.IGNORECASE)
_MODE_RE = re.compile(r"^\s*Mode:\s+\d")
_FREQ_RE = re.compile(r"^\s*Frequency:\s+(.+)$")
# Q-Chem prints per-mode "Raman Intens:" (the boolean "Raman Active:" line is skipped)
_RAMAN_RE = re.compile(r"^\s*Raman\s+Intens(?:ity|ities)?:\s+(.+)$", re.IGNORECASE)


def _parse_floats(s: str) -> List[float]:
    vals = []
    for tok in s.split():
        try:
            vals.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            pass
    return vals


def parse_qchem_raman(path: str) -> pd.DataFrame:
    """Parse Q-Chem FREQ (RAMAN TRUE) output -> DataFrame[nu_cm, intensity]."""
    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    start: Optional[int] = None
    for i, line in enumerate(lines):
        if _VIB_SECTION_RE.search(line):
            start = i
            break
    if start is None:
        return pd.DataFrame(columns=["nu_cm", "intensity"])

    rows: List[Dict[str, float]] = []

    def flush(freqs, raman):
        if not freqs or not raman:
            return
        for k, f in enumerate(freqs):
            if f > 0.0 and k < len(raman):
                rows.append({"nu_cm": f, "intensity": raman[k]})

    freqs, raman = None, None
    for line in lines[start:]:
        if _MODE_RE.match(line):
            flush(freqs, raman)
            freqs, raman = [], []
            continue
        if freqs is None:
            continue
        m = _FREQ_RE.match(line)
        if m:
            freqs = _parse_floats(m.group(1))
            continue
        m = _RAMAN_RE.match(line)
        if m:
            raman = _parse_floats(m.group(1))
            continue
        if "THERMOCHEMICAL" in line.upper():
            break
    flush(freqs, raman)

    return (pd.DataFrame(rows, columns=["nu_cm", "intensity"]) if rows
            else pd.DataFrame(columns=["nu_cm", "intensity"]))


def _sigma_from_fwhm(fwhm: float) -> float:
    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def broaden(df: pd.DataFrame, *, sigma: float, nu_grid: np.ndarray) -> np.ndarray:
    if sigma <= 0:
        raise ValueError("sigma must be positive.")
    curve = np.zeros_like(nu_grid, dtype=float)
    pref = 1.0 / (sigma * SQRT2PI)
    for nu0, inten in zip(df["nu_cm"], df["intensity"]):
        curve += inten * pref * np.exp(-0.5 * ((nu_grid - nu0) / sigma) ** 2)
    return curve


def _stick_spectrum(ax, df, *, nu_min, nu_max):
    lo, hi = min(nu_min, nu_max), max(nu_min, nu_max)
    for _, row in df.iterrows():
        if lo <= row["nu_cm"] <= hi:
            ax.vlines(row["nu_cm"], 0, row["intensity"], color="grey", alpha=0.4, lw=0.8)


def _plot(df, *, sigma, nu_min, nu_max, png, pdf, sticks, ylim, title):
    if df.empty:
        print(f"[SKIP] {png}: no transitions")
        return
    nu_grid = np.arange(nu_min, nu_max, 1.0)
    curve = broaden(df, sigma=sigma, nu_grid=nu_grid)

    fig, ax = plt.subplots(figsize=(4.75, 3.25))
    ax.plot(nu_grid, curve, color="black", lw=1.1)
    if sticks:
        _stick_spectrum(ax, df, nu_min=nu_min, nu_max=nu_max)

    ax.set_xlim(nu_min, nu_max)  # Raman convention: low -> high (opposite IR/VCD)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(r"Raman shift / cm$^{-1}$")
    ax.set_ylabel("Raman activity (arb.)")
    if title:
        ax.set_title(title, fontsize=10, wrap=True)

    fig.tight_layout()
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)


def load_weights(path: str) -> Dict[str, float]:
    df = pd.read_csv(path, comment="#", sep=r"\s+", engine="python",
                     header=None, names=["cid", "E", "dE", "p"])
    return dict(zip(df["cid"].astype(str), df["p"].astype(float)))


def _cli() -> None:
    parser = argparse.ArgumentParser(
        prog="qc_raman_tools.py",
        description="Parse Q-Chem FREQ .out files and plot Boltzmann-weighted Raman spectra.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("logs", nargs="+", help="Q-Chem .out files, glob patterns, or directories")
    parser.add_argument("--bw", dest="bw_file", help="Boltzmann weight file")
    parser.add_argument("--outdir", metavar="TAG",
                        help="Molecule directory; outputs go to <TAG>/06_spectra/<TAG>. "
                             "Overridden by --prefix.")
    parser.add_argument("--prefix", default=None, help="Explicit output prefix [vib]")
    parser.add_argument("--stick", action="store_true", help="Overlay stick spectrum")
    parser.add_argument("--raman_fwhm", type=float, default=10.0, help="Raman FWHM / cm-1")
    parser.add_argument("--xlim", nargs=2, metavar=("MIN", "MAX"), type=float, help="nu range / cm-1")
    parser.add_argument("--ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("--no_title", action="store_true")
    parser.add_argument("--title", default=None, metavar="STR",
                        help="Custom figure title base (overrides prefix-derived title)")
    args = parser.parse_args()

    # ---- resolve output prefix ----
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

    # ---- collect output files ----
    log_files: List[str] = []
    for token in args.logs:
        if os.path.isdir(token):
            log_files.extend(glob.glob(os.path.join(token, "**", "*.out"), recursive=True))
        else:
            log_files.extend(glob.glob(token))
    log_files = [f for f in log_files if not Path(f).stem.startswith("slurm-")]
    if not log_files:
        raise SystemExit("No .out files found.")

    # ---- Boltzmann weights (when given, also the post-dedup keep-list) ----
    if args.bw_file:
        weights = load_weights(args.bw_file)
        log_files = [f for f in log_files if Path(f).stem in weights]
        if not log_files:
            raise SystemExit("No .out files match the Boltzmann keep-list.")
    else:
        weights = {Path(f).stem: 1.0 for f in log_files}

    # ---- aggregate ----
    raman_all = []
    for f in log_files:
        df = parse_qchem_raman(f)
        w = weights.get(Path(f).stem, 1.0)
        if not df.empty:
            df["intensity"] *= w
        raman_all.append(df)
    raman_df = pd.concat(raman_all, ignore_index=True)
    raman_df.to_csv(f"{prefix}_raman.csv", index=False)

    nu_min, nu_max = (0.0, 4000.0) if not args.xlim else tuple(map(float, args.xlim))
    base = args.title if args.title else Path(prefix).name.replace("_", " ")
    title = None if args.no_title else f"{base} Raman Spectrum"

    _plot(raman_df, sigma=_sigma_from_fwhm(args.raman_fwhm), nu_min=nu_min, nu_max=nu_max,
          png=f"{prefix}_raman.png", pdf=f"{prefix}_raman.pdf", sticks=args.stick,
          ylim=tuple(args.ylim) if args.ylim else None, title=title)

    print(f"{len(log_files)} file(s) -> {len(raman_df)} Raman transitions "
          f"(sigma_raman={_sigma_from_fwhm(args.raman_fwhm):.1f} cm-1)")


if __name__ == "__main__":
    _cli()
