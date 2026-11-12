#!/usr/bin/env python3
"""
How many tracks does the alignment need?  Statistical scan + stability over time.

Pools the 3D 3p tracks of a run range (same selection as align_batch.py), then

 1. reference fit on ALL tracks.
 2. subsample scan: for each N, split the tracks into K DISJOINT chunks of N tracks and fit each chunk
      * random  chunks (tracks shuffled): pure statistics -> spread = A / sqrt(N).  A is the
                  statistical constant of the estimator used, so the tracks needed for a target precision is
                  N_needed = (A / target)^2.
      * ordered chunks (consecutive tracks in run order): statistics + real changes over time.  When this curve
                  stops falling with N and stays above the random one, time variation (or other systematics), not
                  statistics, sets the precision.
    The spread of the chunks is measured directly, so it is also right for --estimator median / gauss, where the
    error printed by the fit (least-squares formula, rms of the residuals) does not apply.
 3. run-group scan: one fit per group of --runs-per-group consecutive runs.  Compares the scatter between the
    groups with the statistical error A / sqrt(N_group): chi2/ndf of "constant parameter" and the extra
    (systematic) spread  sigma_sys = sqrt(max(0, weighted variance - <stat err^2>)).

Outputs (in --output-dir):  stat_scan.csv, run_groups.csv, summary.txt,
                             stat_scan_<plane>.png (spread vs N), run_groups_<plane>.png (parameter vs run)

Example
-------
python track_stat_scan.py --data-dir ~/MURAVES/test_data/alignment/RECONSTRUCTED/ROSSO/v1 \\
        --runs 4000 4100 --estimator gauss --target 0.2 --runs-per-group 5
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2 as chi2_dist

import geometry_io as gio
import muraves_alignment as ma

C_RANDOM, C_ORDERED = "#2A9D5C", "#D1495B"
LENGTH = ("du", "dw")


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", default="~/MURAVES/test_data/alignment/RECONSTRUCTED/ROSSO/v1")
    p.add_argument("--pattern", default="MURAVES_AnalyzedData_run*.root")
    p.add_argument("--runs", nargs=2, type=int, metavar=("FIRST", "LAST"), default=[4000, 4100])
    p.add_argument("--good-runs", type=Path, help="good-run list (as in align_batch.py)")
    p.add_argument("--detector", choices=["NERO", "ROSSO", "BLU"])
    p.add_argument("--free", nargs="+", default=["Z2", "Y2"], metavar="PLANE")
    p.add_argument("--params", nargs="+", default=["du", "dw"], choices=ma.PARAMS)
    p.add_argument("--leverage", action="store_true", help="faster-converging update (same result)")
    p.add_argument("--unbiased", action="store_true")
    p.add_argument("--estimator", choices=["mean", "median", "gauss"], default="mean")
    p.add_argument("--gauss-window", choices=["quantile", "iterative"], default="quantile")
    p.add_argument("--max-events", type=int, default=None, help="events per file (default: all)")
    p.add_argument("--chi2-max", type=float, default=None)
    p.add_argument("--corr", action="store_true", help="use the shoulder-corrected *_corr branches")
    p.add_argument("--n-values", nargs="+", type=int,
                   help="chunk sizes to scan (default: 8 log-spaced values from 1000 to N_total/4)")
    p.add_argument("--max-chunks", type=int, default=40,
                   help="at most this many chunks per N (more = better estimate of the spread, slower)")
    p.add_argument("--runs-per-group", type=int, default=1, help="consecutive runs per group in the run scan")
    p.add_argument("--min-group-tracks", type=int, default=2000, help="skip run groups with fewer tracks")
    p.add_argument("--target", type=float, default=0.2, help="wanted precision on du, dw [mm]")
    p.add_argument("--target-angle", type=float, default=0.5, help="wanted precision on the tilts [mrad]")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--output-dir", type=Path, default=Path("./stat_scan"))
    return p.parse_args()


def fit(geom, sample, args, free):
    """one alignment fit -> (params dict, converged)"""
    res = ma.run_alignment(geom, sample, free_planes=free, free_params=args.params, leverage=args.leverage,
                           unbiased=args.unbiased, estimator=args.estimator, gauss_window=args.gauss_window,
                           verbose=False)
    return res["params"], res["converged"], res["errors"]


def keys(free, params):
    return [(p, n) for p in free for n in params]


def target_of(name, args):
    return args.target if name in LENGTH else args.target_angle


def scan_chunks(geom, sample, args, free, n_values, rng):
    """rows: one per (kind, N, plane, param) with the spread of the chunk fits"""
    n_tot = len(sample)
    order = np.lexsort((sample.event, sample.run))          # run order (ordered chunks)
    perm = rng.permutation(n_tot)                           # random chunks
    rows = []
    for N in n_values:
        k_max = min(args.max_chunks, n_tot // N)
        if k_max < 3:
            print(f"  N={N:,}: fewer than 3 chunks possible, skipped")
            continue
        for kind, idx in (("random", perm), ("ordered", order)):
            if kind == "random":
                starts = np.arange(k_max) * N
            else:   # spread the ordered chunks over the whole period
                starts = np.linspace(0, n_tot - N, k_max).astype(int)
                starts = starts[np.r_[True, np.diff(starts) >= N]]   # keep them disjoint
            vals, n_bad = [], 0
            for s in starts:
                try:
                    prm, conv, _ = fit(geom, sample.subset(idx[s:s + N]), args, free)
                except (RuntimeError, ValueError):    # Gaussian core fit failed (low statistics)
                    conv = False
                if not conv:
                    n_bad += 1
                    continue
                vals.append([prm[p][ma.PARAMS.index(n)] for p, n in keys(free, args.params)])
            vals = np.array(vals) * 1e3   # mm / mrad
            K = len(vals)
            print(f"  N={N:>9,} {kind:7s}: {K} chunks" + (f" ({n_bad} not converged, dropped)" if n_bad else ""))
            if K < 3:
                continue
            for j, (p, n) in enumerate(keys(free, args.params)):
                sd = vals[:, j].std(ddof=1)
                rows.append({"kind": kind, "N": N, "plane": p, "param": n, "n_chunks": K,
                             "mean": vals[:, j].mean(), "spread": sd, "spread_err": sd / np.sqrt(2 * (K - 1)),
                             "n_not_converged": n_bad})
    return rows


def stat_constant(rows, p, n):
    """A from the random chunks: spread = A / sqrt(N), weighted fit (weights from spread_err)"""
    r = [x for x in rows if x["kind"] == "random" and x["plane"] == p and x["param"] == n]
    if not r:
        return np.nan, np.nan
    a = np.array([x["spread"] * np.sqrt(x["N"]) for x in r])
    ae = np.array([x["spread_err"] * np.sqrt(x["N"]) for x in r])
    w = 1 / ae ** 2
    return float((w * a).sum() / w.sum()), float(1 / np.sqrt(w.sum()))


def scan_run_groups(geom, sample, args, free):
    runs = np.unique(sample.run)
    groups = [runs[i:i + args.runs_per_group] for i in range(0, len(runs), args.runs_per_group)]
    rows = []
    for g in groups:
        sub = sample.subset(np.isin(sample.run, g))
        if len(sub) < args.min_group_tracks:
            print(f"  runs {g[0]}-{g[-1]}: {len(sub)} tracks < {args.min_group_tracks}, skipped")
            continue
        prm, conv, err = fit(geom, sub, args, free)
        row = {"run_min": int(g[0]), "run_max": int(g[-1]), "n_tracks": len(sub), "converged": conv}
        for p, n in keys(free, args.params):
            i = ma.PARAMS.index(n)
            row[f"{p}_{n}"] = prm[p][i] * 1e3
            row[f"{p}_{n}_fit_err"] = err[p][i] * 1e3
        rows.append(row)
    print(f"  {len(rows)} run groups fitted")
    return rows


def group_stats(rows, p, n, A):
    """constant-parameter test of the run groups with stat errors A/sqrt(N_group)"""
    r = [x for x in rows if x["converged"]]
    if len(r) < 2 or not np.isfinite(A):
        return None
    v = np.array([x[f"{p}_{n}"] for x in r])
    e = A / np.sqrt([x["n_tracks"] for x in r])
    w = 1 / e ** 2
    mean = (w * v).sum() / w.sum()
    chi2 = float((((v - mean) / e) ** 2).sum())
    var = (w * (v - mean) ** 2).sum() / w.sum() * len(v) / (len(v) - 1)
    sys = float(np.sqrt(max(0.0, var - np.mean(e ** 2))))
    return {"mean": float(mean), "chi2": chi2, "ndf": len(v) - 1, "p": float(chi2_dist.sf(chi2, len(v) - 1)),
            "sigma_sys": sys, "n_groups": len(v)}


def plot_scan(rows, p, args, ref, out, title):
    names = args.params
    fig, axs = plt.subplots(1, len(names), figsize=(5.2 * len(names), 4.3), squeeze=False, layout="constrained")
    for ax, n in zip(axs[0], names):
        for kind, col in (("random", C_RANDOM), ("ordered", C_ORDERED)):
            r = sorted((x for x in rows if x["kind"] == kind and x["plane"] == p and x["param"] == n),
                       key=lambda x: x["N"])
            if r:
                ax.errorbar([x["N"] for x in r], [x["spread"] for x in r], [x["spread_err"] for x in r],
                            fmt="o-", color=col, ms=4, capsize=2,
                            label=f"{kind} chunks" + (" (pure stat.)" if kind == "random" else " (run order)"))
        ax.plot(ref["N"], ref["err"][(p, n)], "k*", ms=10, label=f"fit error, all {ref['N']:,} tracks")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("tracks per fit")
        ax.set_ylabel(f"spread of {ma.PARAM_LABEL[n]}")
        ax.set_title(f"{p} {n}")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7)
    fig.suptitle(title, fontsize=10)
    fig.savefig(out / f"stat_scan_{p}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_groups(grows, p, args, A, ref, gst, out, title):
    names = args.params
    r = [x for x in grows if x["converged"]]
    if not r:
        return
    fig, axs = plt.subplots(len(names), 1, figsize=(10, 2.8 * len(names)), squeeze=False, layout="constrained",
                            sharex=True)
    x = np.array([0.5 * (g["run_min"] + g["run_max"]) for g in r])
    for ax, n in zip(axs[:, 0], names):
        if args.runs_per_group > 1:   # alternating bands mark the runs of each group
            for i, g in enumerate(r):
                if i % 2 == 0:
                    ax.axvspan(g["run_min"] - 0.5, g["run_max"] + 0.5, color="0.92", lw=0, zorder=0)
        v = np.array([g[f"{p}_{n}"] for g in r])
        e = A[(p, n)][0] / np.sqrt([g["n_tracks"] for g in r])
        ax.errorbar(x, v, e, fmt="o", ms=3, color=C_ORDERED, capsize=0, label="run group")
        ax.axhline(ref["val"][(p, n)], color="k", lw=1, label="fit on all tracks")
        g = gst.get((p, n))
        if g:
            ax.text(0.01, 0.04, f"$\\chi^2$/ndf = {g['chi2']:.1f}/{g['ndf']}",
                    transform=ax.transAxes, va="bottom", fontsize=8)
        ax.set_ylabel(ma.PARAM_LABEL[n])
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="upper right")
    axs[-1, 0].set_xlabel("run")
    fig.suptitle(f"{p}: {title}", fontsize=10)
    fig.savefig(out / f"run_groups_{p}.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def main():
    args = parse()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    first, last = sorted(args.runs)
    free = [p.upper() for p in args.free]
    try:
        ma.check_free_planes(free, args.leverage, args.unbiased)
    except ValueError as exc:
        raise SystemExit(f"Error: {exc}")

    good = gio.parse_run_list(args.good_runs)[0] if args.good_runs else None
    sel, _ = gio.select_runs(args.data_dir, first, last, good, args.pattern)
    files = [f for _, f in sel]
    if not files:
        raise SystemExit("no run selected")
    geom = ma.load_geometry(files[0], args.detector)
    for f in files[1:]:
        if not ma.same_geometry(geom, ma.load_geometry(f, args.detector)):
            raise SystemExit(f"{f}: different geometry in METADATA than {files[0]}; cannot pool")
    print(geom.describe())
    sample = ma.read_samples(files, geom, max_events=args.max_events, chi2_max=args.chi2_max, corr=args.corr)
    n_tot = len(sample)
    print(f"Total: {n_tot:,} 3D 3p tracks;  estimator={args.estimator}, free={free}, params={args.params}")

    # 1. reference fit on everything
    prm, conv, err = fit(geom, sample, args, free)
    ref = {"N": n_tot, "val": {(p, n): prm[p][ma.PARAMS.index(n)] * 1e3 for p, n in keys(free, args.params)},
           "err": {(p, n): err[p][ma.PARAMS.index(n)] * 1e3 for p, n in keys(free, args.params)}}
    print("reference fit " + ("converged" if conv else "NOT CONVERGED"))

    # 2. chunk scan
    n_values = args.n_values or sorted(set(int(round(x, -2)) for x in np.geomspace(1000, n_tot / 4, 8)))
    print(f"\nSubsample scan, N = {n_values}")
    rows = scan_chunks(geom, sample, args, free, n_values, rng)
    A = {(p, n): stat_constant(rows, p, n) for p, n in keys(free, args.params)}

    # 3. run groups
    print(f"\nRun-group scan ({args.runs_per_group} run(s) per group)")
    grows = scan_run_groups(geom, sample, args, free)
    gst = {(p, n): group_stats(grows, p, n, A[(p, n)][0]) for p, n in keys(free, args.params)}

    # summary
    lines = [f"{geom.detector} runs {first}-{last}: {n_tot:,} tracks, estimator={args.estimator}"
             + (f" ({args.gauss_window})" if args.estimator == "gauss" else "")
             + f", leverage={args.leverage}, unbiased={args.unbiased}, chi2_max={args.chi2_max}, corr={args.corr}",
             "",
             f"{'':10s}{'value':>10s}{'fit err':>10s}{'A':>12s}{'A/sqrtN':>10s}{'target':>8s}{'N needed':>12s}"
             f"{'chi2/ndf':>12s}{'p':>8s}{'sigma_sys':>11s}  (mm / mrad; A, A/sqrtN for all tracks from the random chunks)"]
    for p, n in keys(free, args.params):
        a, ae = A[(p, n)]
        t = target_of(n, args)
        g = gst[(p, n)]
        lines.append(f"{p + ' ' + n:10s}{ref['val'][(p, n)]:10.4f}{ref['err'][(p, n)]:10.4f}"
                     f"{a:8.3g}±{ae:<3.1g}{a / np.sqrt(n_tot):10.4f}{t:8.3g}{(a / t) ** 2:12,.0f}"
                     + (f"{g['chi2']:7.1f}/{g['ndf']:<4d}{g['p']:8.2g}{g['sigma_sys']:11.4f}" if g
                        else f"{'-':>12s}{'-':>8s}{'-':>11s}"))
        if g and g["p"] < 0.05 and g["sigma_sys"] > t:
            lines.append(f"{'':10s}-> sigma_sys > target: more tracks will not reach {t:g}; run-to-run variation "
                         f"dominates (consider time-dependent constants)")
    lines += ["",
              "A        : stat. constant, spread of the chunk fits = A / sqrt(N).  N needed = (A / target)^2.",
              "fit err  : error printed by the alignment (least-squares formula, uses the rms of the residuals);",
              "           compare with A/sqrtN - for median/gauss estimators only the empirical one is right.",
              "chi2/ndf : run groups vs a constant, stat. errors A/sqrt(N_group).  ~1 = stable within statistics;",
              "p        : its p-value.  p < 0.05: significant variation between run groups.",
              "sigma_sys: spread between run groups beyond statistics; the precision floor of a single constant set.",
              "           Only meaningful when p is small (otherwise it is a statistical fluctuation of the estimate)."]
    txt = "\n".join(lines)
    print("\n" + txt)
    (out / "summary.txt").write_text(txt + "\n")
    write_csv(out / "stat_scan.csv", rows)
    write_csv(out / "run_groups.csv", grows)

    title = (f"{geom.detector} runs {first}-{last}, {n_tot:,} tracks, estimator {args.estimator}"
             + (", leverage" if args.leverage else "") + (", unbiased" if args.unbiased else ""))
    for p in free:
        plot_scan(rows, p, args, ref, out, title)
        plot_groups(grows, p, args, A, ref, gst, out, title)
    print(f"\nOutputs in {out.resolve()}")


if __name__ == "__main__":
    main()
