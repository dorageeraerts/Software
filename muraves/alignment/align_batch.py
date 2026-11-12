#!/usr/bin/env python3
"""
Alignment of one batch of runs  ->  ONE parameter set  ->  residuals before / after  ->  geometry files.

Workflow
--------
 1. reconstruct the batch (one ROOT file per run, MURAVES_AnalyzedData_run<N>.root)
 2. python align_batch.py --data-dir ~/MURAVES/test_data/alignment/RECONSTRUCTED/ROSSO/v1 \\
                          --runs 4000 4100 --good-runs good_runs.txt
       * selects the files with FIRST <= run <= LAST that are in the good-run list
       * pools all their 3D 3p tracks and fits ONE set of parameters (planes Z2, Y2; Z1, Z3, Y1, Y3 are the
         references; see muraves_alignment.py for why)
       * recomputes the 3p residuals with the aligned hit positions and a NEW straight-line fit through them
         ("before" = nominal geometry, which reproduces the stored Residue_Track3p_* exactly)
       * writes
            ./images/    alignment_<DET>_runs<A>-<B>_*.png            (./images is added to .gitignore)
            ./geometry/  alignment_<DET>_runs<A>-<B>_parameters.json/.csv, _residual_summary.csv, _runs_used.txt
                         telescope<DET>_aligned_runs<A>-<B>.cfg         <- module list, same columns as the input cfg
                         reco_parameters_aligned_runs<A>-<B>.json       <- same file with z_add / y_add / x_pos updated + x_pos_y
                         alignment_log.csv / .jsonl                     <- one row per fit (see plot_alignment_log.py)
 3. use the cfg / json in the track reconstruction and re-reconstruct (at least a few runs) to confirm
 4. check the parameters on the whole period (run by run):
       python run_alignment.py <all runs>/*.root --per-run --params-in geometry/alignment_<...>_parameters.json \\
              --log geometry/alignment_log.csv
       python plot_alignment_log.py geometry/alignment_log.csv --mode apply --output-dir images

Train / test split (--split ordered|random):
    the selected runs are split BY RUN into a training set (--train-fraction, default 0.5) and a test set.
    ordered: the first runs (in run order) train, the later runs test  -> also probes time stability
    random : runs drawn at random (--seed)                              -> pure statistical cross-check
    The alignment is fitted on the training set only (geometry files as above); the fitted parameters are then
    applied unchanged to the pooled test set (no fitting): residuals before / after with the prefix
    alignment_<tag>_test_*, a log row with mode 'apply', and a train-vs-test table <prefix>train_vs_test.csv.
    Default tag: <DET>_runs<A>-<B>_<split>, so the two splits do not overwrite each other.

--params:  default 'du dw' = exactly what the cfg can express (shift along the measured coordinate and the
position along the detector axis).  The tilts dalpha/dbeta/dgamma cannot be written to the cfg; if you fit them
(--params du dw dalpha dbeta dgamma) they are only stored in the parameters json and the plots then show
residuals that the exported geometry does NOT fully reproduce.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import os
from pathlib import Path

import numpy as np

import geometry_io as gio
import muraves_alignment as ma
import run_alignment as ra


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", default="~/MURAVES/test_data/alignment/RECONSTRUCTED/ROSSO/v1",
                   help="directory with one ROOT file per run")
    p.add_argument("--pattern", default="MURAVES_AnalyzedData_run*.root",
                   help="file pattern inside --data-dir (default: analyzed-data ROOT files)")
    p.add_argument("--runs", nargs=2, type=int, metavar=("FIRST", "LAST"), default=[4000, 4100],
                   help="run range, both ends included (default 4000 4100)")
    p.add_argument("--good-runs", type=Path,
                   help="text file with the good runs (integers separated by whitespace/commas/newlines, "
                        "'a-b' ranges and '#' comments allowed). Without it ALL runs in the range are used.")
    p.add_argument("--detector", choices=["NERO", "ROSSO", "BLU"], help="default: from the METADATA")
    p.add_argument("--free", nargs="+", default=["Z2", "Y2"], metavar="PLANE", help="planes to align (default Z2 Y2)")
    p.add_argument("--params", nargs="+", default=["du", "dw"], choices=ma.PARAMS,
                   help="parameters to fit (default: du dw, the ones the cfg can express)")
    p.add_argument("--leverage", action="store_true", help="faster-converging update (same result)")
    p.add_argument("--unbiased", action="store_true",
                   help="unbiased residuals (aligned plane's hit not in the track fit) in the update and the plots")
    p.add_argument("--estimator", choices=["mean", "median", "gauss"], default="mean",
                   help="align to zero mean residual (least squares, default), zero median residual "
                        "(least absolute deviations) or zero mean of a Gaussian fitted to the residual core")
    p.add_argument("--gauss-window", choices=["quantile", "iterative"], default="quantile",
                   help="fit window of --estimator gauss: central 68%% of the residuals (default) or iterative "
                        "mu +- sigma of the fit")
    p.add_argument("--max-events", type=int, default=None, help="events per file (default: all)")
    p.add_argument("--chi2-max", type=float, default=None, help="optional cut on the stored 3p track chi2")
    p.add_argument("--corr", action="store_true",
                   help="use the shoulder-corrected *_corr branches (hit positions, residuals, best track, chi2)")
    p.add_argument("--forward", action="store_true",
                   help="only forward tracks: 90 < Phi_3p < 270 deg (fit, residuals, plots, test set)")
    p.add_argument("--split", choices=["ordered", "random"],
                   help="split the selected runs into a training set (fit) and a test set (apply only)")
    p.add_argument("--train-fraction", type=float, default=0.5, help="fraction of the runs in the training set")
    p.add_argument("--seed", type=int, default=12345, help="random seed of --split random")
    p.add_argument("--images-dir", type=Path, default=Path("./images"))
    p.add_argument("--geometry-dir", type=Path, default=Path("./geometry"))
    p.add_argument("--tag", default="", help="name of the batch (default <DET>_runs<FIRST>-<LAST>)")
    p.add_argument("--no-gitignore", action="store_true", help="do not touch ./.gitignore")
    return p.parse_args()


def fmt_runs(runs, n=12):
    runs = list(runs)
    if not runs:
        return "-"
    return " ".join(map(str, runs[:n])) + (f" ... ({len(runs)} total)" if len(runs) > n else "")


def main():
    args = parse()
    first, last = sorted(args.runs)
    print("=" * 78)
    print(f"Alignment batch: runs {first}-{last}   data: {args.data_dir}")

    # ---- 1. which runs
    good = None
    if args.good_runs:
        good, ignored = gio.parse_run_list(args.good_runs)
        n_in = sum(first <= r <= last for r in good)
        print(f"good-run list {args.good_runs}: {len(good)} runs, {n_in} inside {first}-{last}")
        if ignored:
            print(f"  [info] ignored non-numeric tokens: {ignored[:6]}{' ...' if len(ignored) > 6 else ''}")
        if n_in == 0:
            raise SystemExit("none of the listed good runs is inside the requested range - wrong file/format?")
    else:
        print("[warn] no --good-runs list: ALL runs in the range are used")
    sel, rep = gio.select_runs(args.data_dir, first, last, good, args.pattern)
    print(f"files in directory: {rep['n_files_dir']},  in range: {rep['n_in_range']},  selected: {rep['n_selected']}")
    if rep["not_good"]:
        print(f"  in range but NOT in the good-run list (skipped): {fmt_runs(rep['not_good'])}")
    if rep["good_without_file"]:
        print(f"  [warn] good runs in range WITHOUT a file: {fmt_runs(rep['good_without_file'])}")
    if rep["unparsable"]:
        print(f"  [warn] files whose run number could not be read from the name: {fmt_runs(rep['unparsable'], 4)}")
    if not sel:
        raise SystemExit("no run selected")
    test = []
    if args.split:
        sel, test = split_runs(sel, args.split, args.train_fraction, args.seed)
        print(f"{args.split} split (train fraction {args.train_fraction}"
              + (f", seed {args.seed}" if args.split == "random" else "") + f"): "
              f"{len(sel)} training runs, {len(test)} test runs")
        print(f"  train: {fmt_runs(r for r, _ in sel)}\n  test : {fmt_runs(r for r, _ in test)}")
    files = [p for _, p in sel]

    # ---- 2. .gitignore for the images
    if not args.no_gitignore:
        try:
            rel = os.path.normpath(str(args.images_dir.resolve().relative_to(Path.cwd().resolve())))
        except ValueError:
            rel = None
            print(f"[info] {args.images_dir} is outside the current directory: .gitignore not touched")
        if rel and rel != ".":
            added, tracked = gio.ensure_gitignore(rel)
            if added:
                print(f"added '{rel}/' to .gitignore")
            if tracked:
                print(f"[warn] {len(tracked)} file(s) under {rel}/ are already tracked by git; .gitignore does not "
                      f"affect them. To stop tracking them:  git rm -r --cached {rel}")

    # ---- 3. alignment on the pooled batch
    geom = ma.load_geometry(files[0], args.detector)
    tag = (args.tag or f"{geom.detector}_runs{first}-{last}" + (f"_{args.split}" if args.split else "")
           + ("_forward" if args.forward else ""))
    prefix = f"alignment_{tag}_"
    ns = ra.parse_args([str(files[0])])
    ns.input, ns.detector, ns.max_events, ns.chi2_max = files, args.detector, args.max_events, args.chi2_max
    ns.corr, ns.estimator, ns.gauss_window = args.corr, args.estimator, args.gauss_window
    ns.forward = args.forward
    ns.free, ns.params, ns.leverage, ns.unbiased = args.free, args.params, args.leverage, args.unbiased
    ns.tag = f"batch {first}-{last}" + (f" {args.split} split train" if args.split else "")
    free = [p.upper() for p in args.free]
    bad = [p for p in free if p not in ma.PLANES]
    if bad:
        raise SystemExit(f"Unknown plane(s) {bad}; allowed: {ma.PLANES}")
    try:
        ma.check_free_planes(free, ns.leverage, ns.unbiased)
    except ValueError as exc:
        raise SystemExit(f"Error: {exc}")
    ns.free = free
    args.geometry_dir.mkdir(parents=True, exist_ok=True)
    args.images_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.geometry_dir / "alignment_log.csv"
    ctx = {}
    row = ra.process_batch(files, ns, args.geometry_dir, None, True, log_path,
                           img_dir=args.images_dir, prefix=prefix, ctx=ctx)
    if row is None:
        raise SystemExit("alignment batch skipped (too few tracks)")
    result, params = ctx["result"], ctx["result"]["params"]

    # ---- 4. geometry files for the reconstruction
    cfg_rows, reco = gio.read_configs(files[0], geom.detector)
    new_rows, changes, warns = gio.build_aligned_cfg(cfg_rows, geom, params)
    span = f"runs{first}-{last}" + (f"_{args.split}_train" if args.split else "")
    cfg_path = args.geometry_dir / f"telescope{geom.detector}_aligned_{span}.cfg"
    reco_path = args.geometry_dir / f"reco_parameters_aligned_{span}.json"
    gio.write_cfg(new_rows, cfg_path)
    reco_path.write_text(json.dumps(gio.build_aligned_reco_parameters(reco, geom, params), indent=2) + "\n")
    (args.geometry_dir / f"{prefix}runs_used.txt").write_text(
        "# runs used for the alignment\n" + "\n".join(str(r) for r, _ in sel) + "\n"
        + ("# in range but not in the good-run list\n" + "\n".join(f"# {r}" for r in rep["not_good"]) + "\n"
           if rep["not_good"] else ""))

    # round trip: the cfg written to disk must imply exactly the fitted (du, dw)
    implied = gio.implied_alignment(gio.read_cfg(cfg_path), geom)
    worst = max(max(abs(implied[p][0] - params[p][0]), abs(implied[p][1] - params[p][1])) for p in implied)
    print("\n" + "-" * 78)
    print(f"Geometry written to {args.geometry_dir.resolve()}")
    print(f"  {cfg_path.name}\n  {reco_path.name}")
    print(f"  cfg round-trip check: implied (du, dw) vs fitted: max difference {worst * 1e6:.2f} um (rounding only)")
    for w in warns:
        print("  [warn]", w)
    print("\n  changes with respect to the telescope cfg stored in the METADATA (stations 1-3):")
    print(f"  {'module':>6s} {'type':>4s} {'stat':>4s}  {'column':24s}{'old':>10s}{'new':>12s}")
    for mod, typ, st, col, old, new in changes:
        print(f"  {mod:>6s} {typ:>4s} {st:>4d}  {col:24s}{old:>10s}{new:>12s}")
    print("  NOTE: for the planes that are not aligned the numbers also change slightly, because the batch was "
          "reconstructed with the z_add / x_pos of reco_parameters.json (the numbers that reproduce the stored "
          "residuals), which differ from the cfg values; the written geometry is the one the fit assumed.")
    tilts = [n for n in args.params if n in ("dalpha", "dbeta", "dgamma")]
    if tilts:
        print(f"  [warn] fitted {tilts} cannot be expressed in the cfg / reco_parameters: they are only in "
              f"{prefix}parameters.json; re-reconstructing with the exported files reproduces the offset part only.")
    print("  The reco reads the geometry from reco_parameters.json only (the cfg gives just the module mapping): "
          "du -> z_add / y_add, dw -> x_pos (Z view) / x_pos_y (Y view).")
    print(f"\nPlots in {args.images_dir.resolve()} (prefix {prefix})")
    print(f"Log row appended to {log_path.resolve()}")

    if test:
        run_test(args, ns, test, prefix, log_path, ctx)


def split_runs(sel, mode, frac, seed):
    """(train, test) lists of (run, path); `sel` is sorted by run"""
    if not 0 < frac < 1:
        raise SystemExit("--train-fraction must be between 0 and 1")
    n_train = int(round(frac * len(sel)))
    if n_train == 0 or n_train == len(sel):
        raise SystemExit(f"cannot split {len(sel)} run(s) with --train-fraction {frac}")
    if mode == "ordered":
        return sel[:n_train], sel[n_train:]
    idx = np.random.default_rng(seed).permutation(len(sel))
    return [sel[i] for i in sorted(idx[:n_train])], [sel[i] for i in sorted(idx[n_train:])]


def centre(r, estimator, window):
    """(centre, error) of a residual distribution in the sense of the alignment estimator"""
    if estimator == "gauss":
        g = ma.gauss_core_fit(r, window)
        return g["mu"], g["mu_err"]
    if estimator == "median":
        return np.median(r), 1.2533 * r.std() / np.sqrt(len(r))
    return r.mean(), r.std() / np.sqrt(len(r))


def run_test(args, ns, test, prefix, log_path, train_ctx):
    """apply the parameters fitted on the training runs, unchanged, to the pooled test runs"""
    params_path = args.geometry_dir / f"{prefix}parameters.json"
    (args.geometry_dir / f"{prefix}runs_test.txt").write_text(
        "# test runs (alignment parameters applied, not fitted)\n" + "\n".join(str(r) for r, _ in test) + "\n")
    tns = copy.copy(ns)
    tns.params_in = params_path
    tns.tag = ns.tag.replace("train", "test")
    print("\n" + "=" * 78)
    print(f"TEST: applying {params_path.name} to {len(test)} test runs")
    ctx = {}
    row = ra.process_batch([p for _, p in test], tns, args.geometry_dir, json.loads(params_path.read_text()),
                           True, log_path, img_dir=args.images_dir, prefix=f"{prefix}test_", ctx=ctx)
    if row is None:
        raise SystemExit("test set skipped (too few tracks)")

    kind = {"mean": "mean", "median": "median", "gauss": f"gauss core mu ({args.gauss_window})"}[args.estimator]
    print("\n" + "-" * 78)
    print(f"Train vs test: {kind} of the {'unbiased' if args.unbiased else 'biased'} residuals [um] "
          f"(train: {len(train_ctx['sample']):,} tracks, test: {len(ctx['sample']):,} tracks)")
    print(f"{'plane':6s}{'train before':>16s}{'train after':>20s}{'test before':>16s}{'test after':>20s}")
    rows = []
    for L in ma.LETTERS:
        for k in range(3):
            vals = []
            for c in (train_ctx, ctx):
                for f in ("fits0", "fits1"):
                    vals += centre(c[f][L]["r"][:, k], args.estimator, args.gauss_window)
            v = np.array(vals) * 1e6
            print(f"{L}{k + 1:<5d}{v[0]:16.1f}{v[2]:11.1f} +- {v[3]:5.1f}{v[4]:16.1f}{v[6]:11.1f} +- {v[7]:5.1f}")
            rows.append([f"{L}{k + 1}"] + [f"{x:.3f}" for x in v])
    with open(args.geometry_dir / f"{prefix}train_vs_test.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["plane"] + [f"{s}_{ba}_{q}_um" for s in ("train", "test") for ba in ("before", "after")
                                for q in ("centre", "err")])
        w.writerows(rows)
    print(f"<chi2>  test  xz: {ctx['fits0']['Z']['chi2'].mean():.4f} -> {ctx['fits1']['Z']['chi2'].mean():.4f}   "
          f"xy: {ctx['fits0']['Y']['chi2'].mean():.4f} -> {ctx['fits1']['Y']['chi2'].mean():.4f}")
    print(f"Test outputs: prefix {prefix}test_ in {args.geometry_dir.resolve()} and {args.images_dir.resolve()}")


if __name__ == "__main__":
    main()