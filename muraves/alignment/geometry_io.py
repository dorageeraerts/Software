#!/usr/bin/env python3
"""
Helpers of the alignment workflow (align_batch.py):

* run selection   : run range + list of good runs -> files
* geometry export : aligned telescope<DET>.cfg and aligned reco_parameters.json, rebuilt from the
                    configuration stored in the METADATA of the ROOT files
* .gitignore      : make sure the image directory is not tracked

Mapping between the telescope cfg (module list) and the alignment planes
------------------------------------------------------------------------
cfg column n_station  1..3  -> plane number k      (station 4 is not aligned yet and is left untouched)
cfg type  x1 / x2           -> Z plane of that station (measures z; two modules = two halves of one plane)
cfg type  y1 / y2           -> Y plane of that station
cfg z                       -> position along the detector axis = -x_pos of the ROOT/track frame
                               (NERO: cfg z = 0.26, 0, -0.262  <->  x_pos = -0.26, 0, 0.262)
cfg disallinemento_misurato -> offset along the measured coordinate = z_add (x-type) / y_add (y-type)

    new disallinemento = offset used in the reconstruction of the batch + du
    new z              = -(x_pos used in the reconstruction + dw)

The *offset used in the reconstruction* is taken from reco_parameters.json in the METADATA (z_add,
y_add, x_pos): these are the numbers that reproduce the stored residuals / chi2 exactly.  They differ
slightly from the cfg numbers (e.g. NERO z_add 0.293 vs cfg 0.289), so ALL modules of stations 1-3 are
rewritten (fixed planes with du = dw = 0), which keeps the written geometry self-consistent.
The tilts (dalpha, dbeta, dgamma) have no column in the cfg and are NOT exported.
"""

from __future__ import annotations

import copy
import re
import subprocess
from pathlib import Path

import numpy as np

import muraves_alignment as ma


# ----------------------------------------------------------------------------
# run selection
# ----------------------------------------------------------------------------
def parse_run_list(path):
    """Run numbers from a text file: integers separated by whitespace / commas / semicolons / newlines,
    ranges 'a-b' allowed, '#' starts a comment.  Other tokens (headers, words) are ignored and counted."""
    runs, ignored = set(), []
    for line in Path(path).expanduser().read_text().splitlines():
        for tok in re.split(r"[,\s;]+", line.split("#")[0].strip()):
            if not tok:
                continue
            m = re.fullmatch(r"(\d+)-(\d+)", tok)
            if m:
                a, b = int(m.group(1)), int(m.group(2))
                runs.update(range(min(a, b), max(a, b) + 1))
            elif tok.isdigit():
                runs.add(int(tok))
            else:
                ignored.append(tok)
    return runs, ignored


def run_of(path):
    m = re.findall(r"run_?(\d+)", Path(path).stem, flags=re.I)
    return int(m[-1]) if m else None


def select_runs(data_dir, first, last, good=None, pattern="*.root"):
    """files of the runs first..last (inclusive) that are in `good` (if given).
    Returns (list of (run, path) sorted by run, report dict)."""
    d = Path(data_dir).expanduser()
    if not d.is_dir():
        raise SystemExit(f"data directory not found: {d}")
    by_run, unparsable = {}, []
    for p in sorted(d.glob(pattern)):
        r = run_of(p)
        if r is None:
            unparsable.append(p.name)
        elif r in by_run:
            raise SystemExit(f"two files for run {r}: {by_run[r].name} and {p.name}")
        else:
            by_run[r] = p
    in_range = {r: p for r, p in by_run.items() if first <= r <= last}
    chosen = {r: p for r, p in in_range.items() if good is None or r in good}
    report = {"n_files_dir": len(by_run), "n_in_range": len(in_range), "n_selected": len(chosen),
              "not_good": sorted(set(in_range) - set(chosen)),
              "good_without_file": sorted(r for r in (good or set()) if first <= r <= last and r not in by_run),
              "unparsable": unparsable}
    return sorted(chosen.items()), report


# ----------------------------------------------------------------------------
# configuration stored in the METADATA
# ----------------------------------------------------------------------------
def read_configs(root_file, detector):
    """(telescope cfg rows, reco_parameters dict) from the METADATA of a ROOT file"""
    meta = ma.read_metadata(root_file)
    if meta is None:
        raise SystemExit(f"no METADATA in {root_file}: cannot rebuild the configuration files")
    cfg = reco = None
    for c in meta.get("detector_config_files", []):
        name = c["path"].split("/")[-1]
        # also accept the aligned copies (telescope<DET>_aligned*.cfg, reco_parameters_aligned*.json)
        if re.fullmatch(rf"telescope{detector}(_\w*)?\.cfg", name, flags=re.I):
            cfg = c["parsed_content"]
        elif name.startswith("reco_parameters") and name.endswith(".json"):
            reco = c["parsed_content"]
    if cfg is None or reco is None:
        raise SystemExit(f"telescope{detector}.cfg and/or reco_parameters.json not found in the METADATA of {root_file}")
    return cfg, reco


def plane_of_module(row):
    """('Z'|'Y', station index 0..3) of a cfg module row"""
    letter = "Z" if row["type"].lower().startswith("x") else "Y"
    return letter, int(row["n_station"]) - 1


# ----------------------------------------------------------------------------
# aligned cfg / reco_parameters
# ----------------------------------------------------------------------------
def build_aligned_cfg(cfg_rows, geom, params):
    """new cfg rows (list of dicts, same keys/order) + list of (module, column, old, new) changes."""
    out, changes, warns = [], [], []
    for row in cfg_rows:
        new = dict(row)
        L, s = plane_of_module(row)
        if s <= 2:
            du, dw = params[f"{L}{s + 1}"][:2]
            base_off = (geom.z_add if L == "Z" else geom.y_add)[s]
            z_new = -(geom.x_reco[L][s] + dw)
            d_new = base_off + du
            if abs(float(row["z"]) + geom.x_reco[L][s]) > 1e-6:
                warns.append(f"module {row['n_mod']} ({row['type']}, station {s + 1}): cfg z = {row['z']} but "
                             f"-x_pos = {-geom.x_reco[L][s]:.4f}; written z is based on x_pos")
            new["z"] = f"{z_new:.6f}"
            new["disallinemento_misurato"] = f"{d_new:.6f}"
            for col in ("z", "disallinemento_misurato"):
                if abs(float(new[col]) - float(row[col])) > 5e-7:
                    changes.append((row["n_mod"], row["type"], s + 1, col, row[col], new[col]))
        out.append(new)
    return out, changes, warns


def write_cfg(rows, path):
    keys = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in keys) + "\n")


def read_cfg(path):
    lines = [l for l in Path(path).read_text().splitlines() if l.strip() and not l.lstrip().startswith("#")]
    keys = lines[0].split()
    return [dict(zip(keys, l.split())) for l in lines[1:]]


def implied_alignment(cfg_rows, geom):
    """(du, dw) per plane implied by a cfg relative to the geometry used for the reconstruction
    (round-trip check of build_aligned_cfg)."""
    out = {}
    for row in cfg_rows:
        L, s = plane_of_module(row)
        if s > 2:
            continue
        base_off = (geom.z_add if L == "Z" else geom.y_add)[s]
        du = float(row["disallinemento_misurato"]) - base_off
        dw = -(float(row["z"]) + geom.x_reco[L][s])
        prev = out.setdefault(f"{L}{s + 1}", (du, dw))
        if abs(prev[0] - du) > 1e-9 or abs(prev[1] - dw) > 1e-9:
            raise ValueError(f"the two modules of plane {L}{s + 1} do not agree")
    return out


def build_aligned_reco_parameters(reco, geom, params):
    """reco_parameters.json with z_add / y_add shifted by du and the plane positions shifted by dw:
    x_pos (Z view) and x_pos_y (Y view, read by the reco instead of x_pos for the Y tracks)."""
    out = copy.deepcopy(reco)
    g = out["detector_geometry"][geom.detector]
    g["x_pos_y"] = list(g.get("x_pos_y", g["x_pos"]))
    for k in range(3):
        g["z_add"][k] = round(float(geom.z_add[k] + params[f"Z{k + 1}"][0]), 6)
        out["reconstruction"]["y_add"][k] = round(float(geom.y_add[k] + params[f"Y{k + 1}"][0]), 6)
        g["x_pos"][k] = round(float(geom.x_reco["Z"][k] + params[f"Z{k + 1}"][1]), 6)
        g["x_pos_y"][k] = round(float(geom.x_reco["Y"][k] + params[f"Y{k + 1}"][1]), 6)
    return out


# ----------------------------------------------------------------------------
# .gitignore
# ----------------------------------------------------------------------------
def ensure_gitignore(entry, gitignore=Path(".gitignore")):
    """append `entry` (e.g. 'images/') to ./.gitignore if no equivalent line exists.
    Returns (added, already_tracked_files)"""
    norm = lambda s: s.strip().lstrip("/").rstrip("/")  # noqa: E731
    lines = gitignore.read_text().splitlines() if gitignore.exists() else []
    added = False
    if norm(entry) not in {norm(l) for l in lines if l.strip() and not l.lstrip().startswith("#")}:
        with open(gitignore, "a") as f:
            if lines and not gitignore.read_text().endswith("\n"):
                f.write("\n")
            f.write(f"# alignment plots (generated)\n{norm(entry)}/\n")
        added = True
    tracked = []
    try:
        res = subprocess.run(["git", "ls-files", norm(entry)], capture_output=True, text=True, timeout=10)
        tracked = [l for l in res.stdout.splitlines() if l]
    except (OSError, subprocess.SubprocessError):
        pass
    return added, tracked