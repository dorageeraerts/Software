#!/usr/bin/env python3
"""Append-only log of alignment results (one row per fit / per evaluation of fixed parameters).

The CSV is meant to be used to follow the alignment over a run period:
  * mode = fit    : parameters fitted on that batch (value + stat. error), residuals before/after
  * mode = apply  : fixed parameters (e.g. from a reference fit) evaluated on that batch;
                    mean/rms residuals before/after show whether they still correct that run.
Units in the log: shifts in mm, angles in mrad, residual means in um, rms in mm.
A JSON-lines twin (<log>.jsonl) keeps the full information (incl. the reco-equivalent offsets)."""

from __future__ import annotations

import csv
import datetime as dt
import json
from pathlib import Path

import numpy as np

import muraves_alignment as ma

UNIT = {"du": "mm", "dw": "mm", "dalpha": "mrad", "dbeta": "mrad", "dgamma": "mrad"}
META_COLS = ["timestamp", "mode", "file", "run_min", "run_max", "n_runs", "first_event", "max_events",
             "n_tracks", "detector", "solver", "free_planes", "free_params", "converged", "n_iter",
             "params_in", "tag", "chi2_xz_before", "chi2_xz_after", "chi2_xy_before", "chi2_xy_after"]
PARAM_COLS = [f"{p}_{n}_{UNIT[n]}{s}" for p in ma.PLANES for n in ma.PARAMS for s in ("", "_err")]
RES_COLS = [f"{p}_{q}" for p in ma.PLANES for q in ("mean_before_um", "mean_after_um", "rms_before_mm", "rms_after_mm")]
COLUMNS = META_COLS + PARAM_COLS + RES_COLS


def make_row(geom, files, sample, result, fits0, fits1, mode, solver, first_event, max_events,
             params_in=None, tag=""):
    runs = np.unique(sample.run)
    row = {"timestamp": dt.datetime.now().isoformat(timespec="seconds"), "mode": mode,
           "file": (";".join(Path(f).name for f in files) if len(files) <= 3 else
                    f"{Path(files[0]).name} .. {Path(files[-1]).name} ({len(files)} files)"), "run_min": int(runs.min()), "run_max": int(runs.max()),
           "n_runs": len(runs), "first_event": first_event, "max_events": "" if max_events is None else max_events,
           "n_tracks": len(sample), "detector": geom.detector, "solver": solver,
           "free_planes": " ".join(result["free_planes"]), "free_params": " ".join(result["free_params"]),
           "converged": "" if mode == "apply" else int(bool(result["converged"])),
           "n_iter": "" if mode == "apply" else len(result["history"]) - 1,
           "params_in": "" if params_in is None else str(params_in), "tag": tag,
           "chi2_xz_before": f"{fits0['Z']['chi2'].mean():.5g}", "chi2_xz_after": f"{fits1['Z']['chi2'].mean():.5g}",
           "chi2_xy_before": f"{fits0['Y']['chi2'].mean():.5g}", "chi2_xy_after": f"{fits1['Y']['chi2'].mean():.5g}"}
    for p in ma.PLANES:
        v = result["params"][p]
        e = result["errors"].get(p)
        for i, n in enumerate(ma.PARAMS):
            row[f"{p}_{n}_{UNIT[n]}"] = f"{v[i] * 1e3:.6g}"
            row[f"{p}_{n}_{UNIT[n]}_err"] = "nan" if e is None or not np.isfinite(e[i]) else f"{e[i] * 1e3:.4g}"
    for L in ma.LETTERS:
        for k in range(3):
            r0, r1 = fits0[L]["r"][:, k], fits1[L]["r"][:, k]
            p = f"{L}{k + 1}"
            row[f"{p}_mean_before_um"] = f"{r0.mean() * 1e6:.4g}"
            row[f"{p}_mean_after_um"] = f"{r1.mean() * 1e6:.4g}"
            row[f"{p}_rms_before_mm"] = f"{r0.std() * 1e3:.5g}"
            row[f"{p}_rms_after_mm"] = f"{r1.std() * 1e3:.5g}"
    return row


def append_row(csv_path, row, extra=None):
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    new = not csv_path.exists()
    if not new:
        with open(csv_path, newline="") as f:
            header = next(csv.reader(f), None)
        if header != COLUMNS:
            raise SystemExit(f"{csv_path} has a different column layout than this version of the scripts; "
                             f"use another --log file (or move the old one).")
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        if new:
            w.writeheader()
        w.writerow(row)
    with open(csv_path.with_suffix(".jsonl"), "a") as f:   # full record incl. reco offsets
        rec = dict(row)
        if extra:
            rec.update(extra)
        f.write(json.dumps(rec) + "\n")


def parse_row(row):
    """dict of strings -> dict with floats where possible"""
    out = {}
    for k, v in row.items():
        try:
            out[k] = float(v)
        except (TypeError, ValueError):
            out[k] = v
    return out


def read_log(csv_path):
    """list of dict rows with numbers converted to float where possible"""
    with open(csv_path, newline="") as f:
        return [parse_row(r) for r in csv.DictReader(f)]


def period_summary(rows):
    """How well do FIXED parameters work run by run?  rows = parsed mode=apply rows (one per run).
    Per plane, over the runs:
      after_wmean   inverse-variance weighted mean of the per-run mean residual [um] (+- its error)
      after_spread  rms over runs of the per-run mean residual [um]
      typ_stat      typical statistical error of one run's mean (median rms/sqrt(N)) [um]
      intrinsic     sqrt(max(spread^2 - <stat^2>, 0)): run-to-run variation beyond statistics [um]
      n_z3          number of runs with |mean| > 3 stat. errors after alignment
      chi2_N        sum z^2 / N_runs for 'mean residual = 0' (about 1 if the parameters are good in every run)
    """
    out = []
    n = np.array([r["n_tracks"] for r in rows], float)
    for p in ma.PLANES:
        mb = np.array([r[f"{p}_mean_before_um"] for r in rows], float)
        ma_ = np.array([r[f"{p}_mean_after_um"] for r in rows], float)
        sem = np.array([r[f"{p}_rms_after_mm"] for r in rows], float) * 1e3 / np.sqrt(n)
        w = 1.0 / sem ** 2
        z = ma_ / sem
        spread = ma_.std(ddof=1) if len(rows) > 1 else float("nan")
        out.append({"plane": p, "n_runs": len(rows), "before_mean_um": mb.mean(),
                    "before_spread_um": mb.std(ddof=1) if len(rows) > 1 else float("nan"),
                    "after_wmean_um": np.sum(w * ma_) / w.sum(), "after_wmean_err_um": 1 / np.sqrt(w.sum()),
                    "after_spread_um": spread, "typ_stat_um": float(np.median(sem)),
                    "intrinsic_um": float(np.sqrt(max(spread ** 2 - np.mean(sem ** 2), 0))) if len(rows) > 1 else float("nan"),
                    "n_z3": int(np.sum(np.abs(z) > 3)), "chi2_N": float(np.mean(z ** 2))})
    return out


def print_period_summary(summ, csv_out=None):
    print("\nFixed parameters applied run by run (all numbers in um, mean residual per run):")
    print(f"{'plane':5s}{'runs':>5s}{'before':>9s}{'spread':>8s}{'after':>9s}{'+-':>7s}{'spread':>8s}"
          f"{'stat':>7s}{'intrinsic':>10s}{'|z|>3':>7s}{'chi2/N':>8s}")
    for s in summ:
        print(f"{s['plane']:5s}{s['n_runs']:5d}{s['before_mean_um']:9.1f}{s['before_spread_um']:8.1f}"
              f"{s['after_wmean_um']:9.1f}{s['after_wmean_err_um']:7.1f}{s['after_spread_um']:8.1f}"
              f"{s['typ_stat_um']:7.1f}{s['intrinsic_um']:10.1f}{s['n_z3']:7d}{s['chi2_N']:8.2f}")
    print("  after = weighted mean over runs; spread = rms over runs; stat = typical error of one run's mean;\n"
          "  intrinsic = run-to-run variation beyond statistics (the parameters are 'enough' if ~0 and 'after' ~ 0)")
    if csv_out:
        with open(csv_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summ[0].keys()), lineterminator="\n")
            w.writeheader()
            for s in summ:
                w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in s.items()})