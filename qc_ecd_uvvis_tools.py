#!/usr/bin/env python3
# ============================================================================
# or_ecd_uvvis_tools.py - Parse ORCA TD‑DFT log files and build
# Boltzmann‑weighted, Gaussian‑broadened UV/Vis and ECD spectra.
#
# OVERVIEW
#   1. Collects *.log files (positional patterns, wild‑cards, or directories),
#   2. Extracts stick data from “ABSORPTION” and “CD” spectrum blocks,
#   3. Applies optional Boltzmann weights (5‑orca‑boltzmann‑weight output),
#   4. Writes CSV tables and plots broadened spectra (PNG + PDF).
#
# Quick examples
#   # simple: average of logs in the current folder
#   or_ecd_uvvis_tools.py *.log
#
#   # recurse through folders, supply Boltzmann weights, 0.35 eV Gaussian FWHM
#   or_ecd_uvvis_tools.py --bw aspirin/bw_results/aspirin_energies.dat \
#       --uv_fwhm 0.35 --ecd_fwhm 0.25 ./aspirin/**
#
#   # publish‑quality plot with sticks and custom limits
#   or_ecd_uvvis_tools.py *.log --stick --xlim 190 350 --uv_ylim 0 0.05 \
#       --ecd_ylim -2 2 --prefix aspirin_spectra
#
# Flags (see --help for full list)
#   --bw PATH          Boltzmann weight file (stage‑5 output, energies.dat)
#   --prefix STR       Prefix for all outputs              [spectra]
#   --uv_fwhm EV       Gaussian FWHM for UV/Vis (eV)       [0.35]
#   --ecd_fwhm EV      Gaussian FWHM for ECD   (eV)        [0.25]
#   --stick            Overlay stick spectrum
#   --flip_x           Plot lambda increasing left -> right
#   --scale FACTOR     Multiply intensities after weighting [1.0]
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

HC_OVER_EV_NM = 1239.84193  # Planck·c in eV·nm
SQRT2PI = np.sqrt(2.0 * np.pi)

_BLOCK_HEADER_RE = re.compile(
    r"^(ABSORPTION|CD) SPECTRUM VIA TRANSITION (ELECTRIC|VELOCITY) DIPOLE MOMENTS"
)

# ---------------------------------------------------------------------- #
#                           PARSING UTILITIES                            #
# ---------------------------------------------------------------------- #


def _canonical_tokens(tokens: Sequence[str]) -> List[str]:
    """Remove leading 'X -> Y' so numeric columns line up."""
    return tokens[3:] if len(tokens) > 2 and tokens[1] == "->" else list(tokens)


def _safe_float(s: str) -> Optional[float]:
    """Return float value; accept FORTRAN 'D' exponents."""
    try:
        return float(s.replace("D", "E"))
    except ValueError:
        return None


def parse_orca_log(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract stick data from a single ORCA TD‑DFT log.

    Returns
    -------
    uv_df : DataFrame
        Columns: energy_eV, wavelength_nm, intensity
    ecd_df : DataFrame
        Same as uv_df but intensity is rotatory strength (R)
    """
    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    uv_rows: List[Dict[str, float]] = []
    ecd_rows: List[Dict[str, float]] = []

    i = 0
    n_lines = len(lines)
    while i < n_lines:
        if not _BLOCK_HEADER_RE.match(lines[i].lstrip()):
            i += 1
            continue

        is_uv = lines[i].lstrip().startswith("ABSORPTION")
        block: List[str] = []
        j = i + 1
        while j < n_lines and not _BLOCK_HEADER_RE.match(lines[j].lstrip()):
            block.append(lines[j])
            j += 1

        try:
            unit_idx = next(idx for idx, l in enumerate(block) if "(nm)" in l)
        except StopIteration:
            i = j
            continue

        hdr_tokens = block[unit_idx - 1].split()
        unit_tokens = block[unit_idx].split()

        # ORCA ≤5 used "(eV)" energy and "X -> Y" state notation;
        # _canonical_tokens strips those 3 tokens so column indices need
        # no offset.  ORCA 6 uses "(cm-1)" and a bare state index, which
        # is NOT stripped, shifting every data column right by +1.
        if "(eV)" in unit_tokens:
            col_offset = 0
        elif "(cm-1)" in unit_tokens:
            col_offset = 1
        else:
            i = j
            continue

        wl_col = unit_tokens.index("(nm)") + col_offset

        if is_uv:
            inten_token = next(
                (t for t in ("fosc(D2)", "fosc(P2)", "fosc") if t in hdr_tokens),
                None,
            )
        else:
            inten_token = "R" if "R" in hdr_tokens else None

        if inten_token is None:
            i = j
            continue

        inten_col = hdr_tokens.index(inten_token) - 1 + col_offset

        for raw in block[unit_idx + 1 :]:
            if not raw.strip() or raw.lstrip().startswith("-"):
                continue
            toks = _canonical_tokens(raw.split())
            if len(toks) <= max(wl_col, inten_col):
                continue
            w = _safe_float(toks[wl_col])
            inten = _safe_float(toks[inten_col])
            if None in (w, inten):
                continue
            e = HC_OVER_EV_NM / w  # derive energy from wavelength; avoids cm-1/eV ambiguity
            (uv_rows if is_uv else ecd_rows).append(
                {"energy_eV": e, "wavelength_nm": w, "intensity": inten}
            )

        i = j

    return pd.DataFrame(uv_rows), pd.DataFrame(ecd_rows)


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
    Return broadened curve on the supplied λ grid.

    Each Gaussian is normalised with 1/(σ√(2π)) so its peak equals the stick
    intensity.  If *jacobian* is True (ECD), multiply by |dE/dλ| = E/λ so area
    is preserved after change of variable.
    """
    if sigma_eV <= 0:
        raise ValueError("σ must be positive.")

    curve = np.zeros_like(lam_nm, dtype=float)
    prefactor = 1.0 / (sigma_eV * SQRT2PI)
    for e, inten in zip(df["energy_eV"], df["intensity"]):
        gaussian = np.exp(-0.5 * ((energy_grid_eV - e) / sigma_eV) ** 2)
        if jacobian:
            gaussian *= e / lam_nm
        curve += inten * prefactor * gaussian
    return curve


def _stick_spectrum(ax, df: pd.DataFrame, *, is_ecd: bool, scale: float):
    """Draw vertical lines representing individual transitions."""
    for _, row in df.iterrows():
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
        _stick_spectrum(ax, df, is_ecd=is_ecd, scale=scale)

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
    """Read 5-orca-boltzmann-weight energies.dat; return {CID: p}."""
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
        prog="or_ecd_uvvis_tools.py",
        description="Parse ORCA TD-DFT logs and generate UV/Vis & ECD spectra.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "logs",
        nargs="+",
        help="ORCA .log files, glob patterns, or directories (recursive).",
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
        help="Plot long → short wavelength (e.g. 350 → 190 nm)",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        metavar="FAC",
        help="Multiply intensities after Boltzmann weighting",
    )
    parser.add_argument("--uv_fwhm", type=float, default=0.35, help="UV‑Vis Gaussian FWHM / eV")
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
    #   Resolve output prefix (--outdir takes effect if --prefix absent) #
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

    # ------------------------------------------------------------------ #
    #   Locate .log files                                                #
    # ------------------------------------------------------------------ #
    log_files: List[str] = []
    for token in args.logs:
        if os.path.isdir(token):
            log_files.extend(glob.glob(os.path.join(token, "**", "*.log"), recursive=True))
        else:
            log_files.extend(glob.glob(token))
    if not log_files:
        raise SystemExit("No .log files found.")

    # ------------------------------------------------------------------ #
    #   Boltzmann weights                                                #
    # ------------------------------------------------------------------ #
    weights: Dict[str, float] = {
        Path(f).stem: 1.0 for f in log_files
    }
    if args.bw_file:
        bw_map = _load_weights(args.bw_file)
        for cid in weights:
            if cid in bw_map:
                weights[cid] = bw_map[cid]

    # ------------------------------------------------------------------ #
    #   Aggregate all transitions                                        #
    # ------------------------------------------------------------------ #
    uv_all, ecd_all = [], []
    for f in log_files:
        uv_df, ecd_df = parse_orca_log(f)
        cid = Path(f).stem
        weight = weights.get(cid, 1.0) * args.scale
        uv_df["intensity"] *= weight
        ecd_df["intensity"] *= weight
        uv_all.append(uv_df)
        ecd_all.append(ecd_df)

    uv_df = pd.concat(uv_all, ignore_index=True)
    ecd_df = pd.concat(ecd_all, ignore_index=True)

    uv_df.to_csv(f"{prefix}_uvvis.csv", index=False)
    ecd_df.to_csv(f"{prefix}_ecd.csv", index=False)

    lam_min, lam_max = (160.0, 400.0) if not args.xlim else sorted(args.xlim)
    fwhm_to_sigma = lambda f: f / (2 * np.sqrt(2 * np.log(2)))
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
        f"{len(log_files)} log(s) → {len(uv_df)} UV and {len(ecd_df)} ECD transitions "
        f"(σ_uv={sigma_uv:.3f} eV, σ_ecd={sigma_ecd:.3f} eV, scale={args.scale})"
    )


if __name__ == "__main__":
    _cli()
