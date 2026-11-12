#!/usr/bin/env python3
"""
Plot the alignment log written by run_alignment.py (alignment_log.csv) over a run period.

    python plot_alignment_log.py period.csv                 # everything in the log
    python plot_alignment_log.py period.csv --mode fit      # only independent per-run fits
    python plot_alignment_log.py period.csv --mode apply --tag "ref 5736"

Figures:
  period_parameters.png      fitted parameters vs run (mode=fit rows), +- stat. error, with the weighted
                             mean and the chi2/ndf of "constant in time" -> steps / drifts of the geometry
  period_residual_means.png  mean residual per plane before / after alignment vs run (the 'are the residual
                             plots fixed over the whole period?' check; use rows with mode=apply to test one
                             fixed set of parameters, mode=fit for run-by-run fits)
  period_chi2.png            mean 3p track chi2 before / after vs run
  period_check_summary.csv   (+ printed table) per plane: weighted mean / spread / intrinsic run-to-run
                             variation of the mean residual after applying FIXED parameters
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import alignment_log as alog
import muraves_alignment as ma


def get(rows, key):
    return np.array([r[key] for r in rows], float)


def xaxis(rows, kind):
    if kind == "index":
        return np.arange(len(rows), dtype=float), "row in log"
    return get(rows, "run_min"), "run"


def plot_parameters(rows, kind, out):
    fits = [r for r in rows if r["mode"] == "fit"]
    if not fits:
        print("  no mode=fit rows: skipping parameter plot")
        return
    x, xl = xaxis(fits, kind)
    planes = [p for p in ma.PLANES
              if any(np.isfinite(r[f"{p}_du_mm_err"]) or np.isfinite(r[f"{p}_dw_mm_err"]) for r in fits)]
    fig, axs = plt.subplots(len(ma.PARAMS), len(planes), figsize=(5 * len(planes), 2.2 * len(ma.PARAMS) + 0.8),
                            sharex=True, squeeze=False)
    for c, p in enumerate(planes):
        for r_, n in enumerate(ma.PARAMS):
            col = f"{p}_{n}_{alog.UNIT[n]}"
            y, e = get(fits, col), get(fits, col + "_err")
            ok = np.isfinite(e) & (e > 0)
            ax = axs[r_, c]
            if ok.any():
                ax.errorbar(x[ok], y[ok], yerr=e[ok], fmt="o", ms=3.5, lw=0.8, color="C0")
                w = 1 / e[ok] ** 2
                mean = np.sum(w * y[ok]) / w.sum()
                chi2 = np.sum(w * (y[ok] - mean) ** 2)
                ax.axhline(mean, color="0.4", ls="--", lw=0.8)
                ax.set_title(f"mean {mean:+.3f}, $\\chi^2$/ndf = {chi2:.1f}/{ok.sum() - 1}", loc="right", fontsize=8)
            ax.grid(True, alpha=0.3, ls="--")
            ax.set_ylabel(ma.PARAM_LABEL[n] if hasattr(ma, "PARAM_LABEL") else n, fontsize=9)
            if r_ == 0:
                ax.set_title(f"plane {p}", loc="left", fontsize=10, fontweight="bold")
            if r_ == len(ma.PARAMS) - 1:
                ax.set_xlabel(xl)
    fig.suptitle("Fitted alignment parameters per batch (error bars: stat.)")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out / "period_parameters.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_residual_means(rows, kind, out):
    """mean residual per plane vs run.  A pooled fit over several runs (mode=fit, n_runs>1) is shown as a
    grey band = the training runs; all other rows as points."""
    train = [r for r in rows if r["mode"] == "fit" and r["n_runs"] > 1]
    rows = [r for r in rows if not (r["mode"] == "fit" and r["n_runs"] > 1)]
    x, xl = xaxis(rows, kind)
    n = get(rows, "n_tracks")
    modes = sorted({r["mode"] for r in rows})
    mk = {"fit": "o", "apply": "s"}
    fig, axs = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True, squeeze=False)
    for i, L in enumerate(ma.LETTERS):
        for k in range(3):
            p = f"{L}{k + 1}"
            ax = axs[i, k]
            if kind == "run":
                for j, t in enumerate(train):
                    ax.axvspan(t["run_min"], t["run_max"], color="0.8", alpha=0.7, lw=0,
                               label="training runs (pooled fit)" if j == 0 else None)
            rms = get(rows, f"{p}_rms_before_mm") * 1e3          # um
            ax.errorbar(x, get(rows, f"{p}_mean_before_um"), yerr=rms / np.sqrt(n), fmt="o", ms=3.5, lw=0.8,
                        color="tab:blue", label="before")
            for m in modes:
                sel = np.array([r["mode"] == m for r in rows])
                rms_a = get(rows, f"{p}_rms_after_mm") * 1e3
                ax.errorbar(x[sel], get(rows, f"{p}_mean_after_um")[sel], yerr=(rms_a / np.sqrt(n))[sel],
                            fmt=mk.get(m, "o"), ms=3.5, lw=0.8, color="tab:orange", label=f"after ({m})")
            ax.axhline(0, color="black", lw=0.8)
            ax.set_title(p, fontsize=10)
            ax.grid(True, alpha=0.3, ls="--")
            if k == 0:
                ax.set_ylabel("mean residual [$\\mu$m]")
            if i == 1:
                ax.set_xlabel(xl)
            if i == 0 and k == 0:
                ax.legend(fontsize=8)
    fig.suptitle("Mean 3p residual (measured - fit) per run, before / after alignment (error bars: rms/$\\sqrt{N}$)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out / "period_residual_means.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_chi2(rows, kind, out):
    x, xl = xaxis(rows, kind)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for ax, proj in zip(axs, ("xz", "xy")):
        ax.plot(x, get(rows, f"chi2_{proj}_before"), "o", ms=4, color="tab:blue", label="before")
        ax.plot(x, get(rows, f"chi2_{proj}_after"), "o", ms=4, color="tab:orange", label="after")
        ax.set_xlabel(xl)
        ax.set_ylabel(f"<chi2> 3p {proj}")
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend()
    fig.tight_layout()
    fig.savefig(out / "period_chi2.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", type=Path)
    ap.add_argument("--mode", choices=["fit", "apply", "all"], default="all")
    ap.add_argument("--tag", default=None, help="only rows whose tag equals this")
    ap.add_argument("--x", choices=["run", "index"], default="run", help="x axis (index = row order in the log)")
    ap.add_argument("--output-dir", type=Path, default=None)
    args = ap.parse_args()
    rows = alog.read_log(args.log)
    if args.mode != "all":
        rows = [r for r in rows if r["mode"] == args.mode]
    if args.tag is not None:
        rows = [r for r in rows if r["tag"] == args.tag]
    if not rows:
        raise SystemExit("no rows selected")
    if args.x == "run":
        rows.sort(key=lambda r: (r["run_min"], r["timestamp"]))
    out = args.output_dir or args.log.parent
    out.mkdir(parents=True, exist_ok=True)
    print(f"{len(rows)} rows ({', '.join(sorted({r['mode'] for r in rows}))})")
    plot_parameters(rows, args.x, out)
    plot_residual_means(rows, args.x, out)
    plot_chi2([r for r in rows if not (r["mode"] == "fit" and r["n_runs"] > 1)], args.x, out)
    applied = [r for r in rows if r["mode"] == "apply"]
    if applied:
        alog.print_period_summary(alog.period_summary(applied), out / "period_check_summary.csv")
    print(f"plots written to {out.resolve()}")


if __name__ == "__main__":
    main()