#!/usr/bin/env python3
"""
conformer_dedup.py -- deduplicate Stage-4 conformers before Boltzmann weighting.

Multiple Confab seeds frequently relax to the SAME minimum during the Stage-4
DFT optimization. Left in place, each redundant copy is counted again in the
Stage-5 Boltzmann sum, distorting the populations (and the Stage-6 spectra).

This tool clusters the optimized conformers by heavy-atom RMSD (with an energy
guard) and writes the lowest-energy representative of each cluster to:

    <TAG>/03_solvent_opt/<TAG>_unique.dat   (one conformer ID per line)

Stage 5 reads that file (when present) and weights only the unique conformers.

Engine is auto-detected per conformer directory:
    *.out present  -> Q-Chem  (geometry from last "Standard Nuclear Orientation")
    *.log present  -> ORCA    (geometry from <cid>.xyz; energy from the .log)

Usage:
    conformer_dedup.py TAG
    conformer_dedup.py --list molecules.txt
    conformer_dedup.py --rmsd 0.10 --ecut 0.05 ephedrine
    conformer_dedup.py --dry-run ephedrine        # report only, write nothing

Flags:
    --rmsd A      heavy-atom RMSD cluster threshold, Angstrom   [0.10]
    --ecut KCAL   energy guard for merging, kcal/mol            [0.05]
    --opt-dir D   stage-4 output subdir                         [03_solvent_opt]
    --list FILE   text file of TAGs (one per line)
    --dry-run     print the report but do not write the keep-list
    -h, --help    show this help
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

H2KCAL = 627.5094740631


def read_orca_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        els.append(p[0]); xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, np.array(xyz)


def read_qchem_geom(path):
    txt = open(path, errors="ignore").read()
    blk = txt.split("Standard Nuclear Orientation (Angstroms)")[-1].splitlines()
    els, xyz, dash, on = [], [], 0, False
    for ln in blk:
        s = ln.strip()
        if s and set(s) == {"-"}:
            dash += 1; on = (dash == 1)
            if dash == 2:
                break
            continue
        if on:
            p = ln.split()
            if len(p) >= 5 and p[0].isdigit():
                els.append(p[1]); xyz.append([float(p[2]), float(p[3]), float(p[4])])
    return els, np.array(xyz)


def qchem_energy(path):
    m = re.findall(r"Final energy is\s+(-?\d+\.\d+)", open(path, errors="ignore").read())
    return float(m[-1]) if m else None


def orca_energy(path):
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", open(path, errors="ignore").read())
    return float(m[-1]) if m else None


def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))]) @ Wt
    return np.sqrt(((Pc @ R - Qc) ** 2).sum() / len(P))


def heavy(els, xyz):
    return xyz[[i for i, e in enumerate(els) if e != "H"]]


def load_conformer(cdir):
    """Return dict(cid, els, xyz, E) for one conformer dir, or None."""
    cid = os.path.basename(os.path.normpath(cdir))
    out = glob.glob(os.path.join(cdir, "*.out"))
    log = glob.glob(os.path.join(cdir, "*.log"))
    if out:                                     # Q-Chem
        els, xyz = read_qchem_geom(out[0]); E = qchem_energy(out[0])
    elif log:                                   # ORCA
        xs = [x for x in glob.glob(os.path.join(cdir, "*.xyz")) if not x.endswith("_trj.xyz")]
        if not xs:
            return None
        els, xyz = read_orca_xyz(xs[0]); E = orca_energy(log[0])
    else:
        return None
    if E is None or len(xyz) == 0:
        return None
    return {"cid": cid, "els": els, "xyz": xyz, "E": E}


def cluster(confs, rmsd_thresh, ecut):
    """Lowest-energy-first greedy clustering. Returns list of clusters (each a
    list of conformer dicts, representative first)."""
    clusters = []
    for c in sorted(confs, key=lambda x: x["E"]):
        for cl in clusters:
            rep = cl[0]
            if rep["els"] != c["els"]:
                continue
            r = kabsch_rmsd(heavy(rep["els"], rep["xyz"]), heavy(c["els"], c["xyz"]))
            if r < rmsd_thresh and abs(c["E"] - rep["E"]) * H2KCAL < ecut:
                cl.append(c)
                break
        else:
            clusters.append([c])
    return clusters


def process(tag, opt_subdir, rmsd_thresh, ecut, dry_run):
    opt_dir = os.path.join(tag, opt_subdir)
    if not os.path.isdir(opt_dir):
        print(f"[{tag}] {opt_dir} not found -- skipping", file=sys.stderr)
        return
    confs = [c for c in (load_conformer(d) for d in sorted(glob.glob(opt_dir + "/*/"))) if c]
    if not confs:
        print(f"[{tag}] no parseable conformers in {opt_dir}", file=sys.stderr)
        return
    clusters = cluster(confs, rmsd_thresh, ecut)
    reps = [cl[0]["cid"] for cl in clusters]
    print(f"[{tag}] {len(confs)} conformers -> {len(reps)} unique")
    for cl in clusters:
        if len(cl) > 1:
            members = ", ".join(c["cid"] for c in cl[1:])
            print(f"    keep {cl[0]['cid']}  (merged: {members})")
    keep_file = os.path.join(opt_dir, f"{tag}_unique.dat")
    if dry_run:
        print(f"    (dry run) would write {keep_file}", file=sys.stderr)
    else:
        with open(keep_file, "w", newline="\n") as fh:   # force LF (Unix) endings
            fh.write("\n".join(reps) + "\n")
        print(f"    wrote {keep_file}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("tag", nargs="?")
    ap.add_argument("--rmsd", type=float, default=0.10)
    ap.add_argument("--ecut", type=float, default=0.05)
    ap.add_argument("--opt-dir", default="03_solvent_opt")
    ap.add_argument("--list", dest="list_file")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-h", "--help", action="store_true")
    a = ap.parse_args()
    if a.help or (not a.tag and not a.list_file):
        print(__doc__)
        sys.exit(0)
    tags = []
    if a.list_file:
        tags = [t.strip() for t in open(a.list_file)
                if t.strip() and not t.lstrip().startswith(("#", ";"))]
    else:
        tags = [a.tag]
    for t in tags:
        process(t, a.opt_dir, a.rmsd, a.ecut, a.dry_run)


if __name__ == "__main__":
    main()
