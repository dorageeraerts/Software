#!/usr/bin/env python3
"""
MURAVES track-based alignment of planes 1-3 (Karimaki et al., physics/0306034).

Examples
--------
# small test batch: first 20000 events of one run, align the middle planes Z2/Y2
python run_alignment.py MURAVES_AnalyzedData_run5736.root --max-events 20000

# several runs, faster-converging variant of the update, only fit du and dw
python run_alignment.py /pnfs/.../MURAVES_AnalyzedData_run57{36,37,38}.root \
        --leverage --params du dw

# apply previously found parameters to another batch (no fitting): residuals before / after
python run_alignment.py other_run.root --params-in alignment_output/alignment_parameters.json

# train on a small batch (pooled fit over these files, ONE parameter set), then check on the whole period
python run_alignment.py train_run{5736..5745}.root --leverage \
        --check /pnfs/.../MURAVES_AnalyzedData_run*.root --log period.csv

# follow the whole period: ONE independent fit per file, one row per run in the log
python run_alignment.py /pnfs/.../MURAVES_AnalyzedData_run*.root --per-run --leverage --log period.csv
# ... and check that fixed reference parameters still work in every run (no fitting)
python run_alignment.py /pnfs/.../MURAVES_AnalyzedData_run*.root --per-run \
        --params-in alignment_output/alignment_parameters.json --log period.csv
python plot_alignment_log.py period.csv

Which planes are aligned?  See the module docstring of muraves_alignment.py: with 3 planes
per projection only one combination of the positions is observable, so reference planes
have to be fixed.  Default: Z1, Z3, Y1, Y3 fixed; Z2, Y2 aligned (--free to change).
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path

import numpy as np

import alignment_plots as ap
import alignment_log as alog
import muraves_alignment as ma
import plot_alignment_log as pal


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="+", help="MURAVES_AnalyzedData_run*.root file(s)")
    p.add_argument("--detector", choices=["NERO", "ROSSO", "BLU"],
                   help="default: taken from the METADATA of the first file")
    p.add_argument("--max-events", type=int, default=None,
                   help="use only the first N events of every file (small test batch)")
    p.add_argument("--first-event", type=int, default=0,
                   help="skip the first N events of every file (e.g. to validate on events not used for the fit)")
    p.add_argument("--chi2-max", type=float, default=None,
                   help="optional cut on the stored 3p track chi2 (both projections) for the fit sample")
    p.add_argument("--corr", action="store_true",
                   help="use the shoulder-corrected *_corr branches (hit positions, residuals, best track, chi2)")
    p.add_argument("--forward", action="store_true",
                   help="only forward tracks: 90 < Phi_3p < 270 deg (fit sample, residuals and plots)")
    p.add_argument("--free", nargs="+", default=["Z2", "Y2"], metavar="PLANE",
                   help="planes to align (Z1 Z2 Z3 Y1 Y2 Y3 or 'all'); the others are references. "
                        "Default: Z2 Y2")
    p.add_argument("--params", nargs="+", default=ma.PARAMS, choices=ma.PARAMS,
                   help="subset of parameters to fit (others stay 0). Default: all five")
    p.add_argument("--leverage", action="store_true",
                   help="weight the normal matrix with (1-h_kk): same minimum, converges in 2-3 "
                        "iterations instead of ~6-10")
    p.add_argument("--unbiased", action="store_true",
                   help="use unbiased residuals e = r/(1-h_kk) (aligned plane's hit not in the track fit) in "
                        "the update AND in all residual tables/plots/log; same result as --leverage")
    p.add_argument("--estimator", choices=["mean", "median", "gauss"], default="mean",
                   help="mean: least-squares update (zero mean residual, default); median: least-absolute-"
                        "deviations update (zero median residual, less sensitive to tails/outliers); gauss: zero "
                        "mean of a Gaussian fitted to the residual core (see --gauss-window)")
    p.add_argument("--gauss-window", choices=["quantile", "iterative"], default="quantile",
                   help="fit window of --estimator gauss: quantile = central 68%% of the residuals (default); "
                        "iterative = refit in mu +- sigma of the previous fit until stable")
    p.add_argument("--rigid-station", nargs="+", type=int, default=[], choices=[1, 2, 3], metavar="K",
                   help="treat station K as one rigid body: its Z and Y layer share the shifts (dx, dy, dz) "
                        "(Z: du=dz, dw=dx; Y: du=dy, dw=dx); tilts, if fitted, stay per layer. "
                        "Z<K> and Y<K> must be free. Only with --estimator mean")
    p.add_argument("--yz-spacing", type=float, default=ma.YZ_SPACING, metavar="M",
                   help="x(Y layer) - x(Z layer) within a station [m], added to the reco x_pos_y in the alignment "
                        "model only (outputs for the reco keep the reco convention). Default: %(default)s. "
                        "Does not change the residuals, only v (enters the tilt terms)")
    p.add_argument("--max-iter", type=int, default=30)
    p.add_argument("--tol-length", type=float, default=2e-6, help="convergence: max step of du,dw [m]")
    p.add_argument("--tol-angle", type=float, default=2e-6, help="convergence: max step of angles [rad]")
    p.add_argument("--params-in", type=Path, help="JSON from a previous run: only evaluate, do not fit")
    p.add_argument("--per-run", action="store_true",
                   help="treat every input file as its own batch: one fit (or one evaluation with "
                        "--params-in) per file, each with the geometry of ITS OWN METADATA; "
                        "one log row per file; plots only with --plots-per-run")
    p.add_argument("--check", nargs="+", type=Path, metavar="FILE",
                   help="after the (pooled) fit on the input files, apply the fitted parameters unchanged to each of "
                        "these files separately (e.g. the whole period) and log/plot/summarise the residuals run by run")
    p.add_argument("--check-max-events", type=int, default=None,
                   help="events per file for the --check runs (default: all; --max-events only limits the training batch)")
    p.add_argument("--plots-per-run", action="store_true", help="with --per-run: also make the full set of plots per file")
    p.add_argument("--log", type=Path, default=None,
                   help="CSV log to append to (default: <output-dir>/alignment_log.csv). One row per fit/evaluation; "
                        "plot with plot_alignment_log.py")
    p.add_argument("--tag", default="", help="free text stored in the log row (e.g. 'ref fit Z2Y2 leverage')")
    p.add_argument("--output-dir", type=Path, default=Path("./alignment_output"))
    return p.parse_args(argv)


def solver_name(args):
    med = {"mean": "", "median": " (median)", "gauss": f" (gauss core, {args.gauss_window})"}[args.estimator]
    if args.rigid_station:
        med += f" (rigid station {', '.join(map(str, args.rigid_station))})"
    if args.unbiased:
        return "unbiased residuals" + med
    return ("paper + leverage weights" if args.leverage else "paper") + med


def mad_sigma(x):
    """robust width: 1.4826 * median absolute deviation"""
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def write_outputs(out_dir, geom, args, sample, files, result, fits0, fits1, prefix=""):
    params, errors = result["params"], result["errors"]
    n = len(sample)
    meta = {
        "detector": geom.detector, "geometry": geom.describe(), "files": [str(f) for f in files],
        "n_tracks": n, "runs": [int(x) for x in np.unique(sample.run)], "max_events_per_file": args.max_events, "chi2_max": args.chi2_max,
        "free_planes": result["free_planes"], "free_params": result["free_params"],
        "rigid_stations": result.get("rigid_stations", []), "yz_spacing_m": geom.yz_spacing,
        "solver": solver_name(args),
        "residuals": "unbiased (r / (1 - h_kk))" if args.unbiased else "biased (hit included in the track fit)",
        "corr": bool(args.corr), "branches": ma._needed_branches(args.corr),
        "forward_only": bool(args.forward),
        "converged": result["converged"], "iterations": len(result["history"]) - 1,
        "units": "du, dw in metres; dalpha, dbeta, dgamma in radians",
        "convention": "corrected hit: x' = x_plane + dw + dbeta*u + dalpha*v ; u' = u_meas + du - dgamma*v ; "
                      "u,v relative to the plane centre; v = coordinate along the strips taken from the other projection",
    }
    meta["reco_equivalent_offsets_m"] = ma.reco_equivalent_offsets(geom, params)
    f_json = f"{prefix}parameters.json" if prefix else "alignment_parameters.json"
    f_csv = f"{prefix}parameters.csv" if prefix else "alignment_parameters.csv"
    with open(out_dir / f_json, "w") as f:
        json.dump({"meta": meta, "parameters": ma.params_to_json(params, errors)}, f, indent=2)

    with open(out_dir / f_csv, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["plane", "role"] + [f"{nme}{u}" for nme, u in zip(
            ma.PARAMS, ["[mm]", "[mm]", "[mrad]", "[mrad]", "[mrad]"])] + [f"{nme}_err" for nme in ma.PARAMS])
        for plane in ma.PLANES:
            v = params[plane]
            e = errors.get(plane, np.full(len(ma.PARAMS), np.nan))
            w.writerow([plane, "aligned" if plane in result["free_planes"] else "reference"]
                       + [f"{x * 1e3:.5g}" for x in v] + [f"{x * 1e3:.3g}" for x in e])

    with open(out_dir / f"{prefix}residual_summary.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["plane", "n_tracks", "mean_before_um", "mean_after_um", "mean_err_um",
                    "rms_before_mm", "rms_after_mm", "sigmaMAD_before_mm", "sigmaMAD_after_mm"])
        for L in ma.LETTERS:
            for k in range(3):
                r0, r1 = fits0[L]["r"][:, k], fits1[L]["r"][:, k]
                w.writerow([f"{L}{k + 1}", len(r0), f"{r0.mean() * 1e6:.3f}", f"{r1.mean() * 1e6:.3f}",
                            f"{r1.std() / np.sqrt(len(r1)) * 1e6:.3f}", f"{r0.std() * 1e3:.4f}",
                            f"{r1.std() * 1e3:.4f}", f"{mad_sigma(r0) * 1e3:.4f}", f"{mad_sigma(r1) * 1e3:.4f}"])


def print_table(result):
    print("\nAlignment parameters (value +- stat. error), units: mm / mrad")
    print(f"{'plane':6s}" + "".join(f"{n:>20s}" for n in ma.PARAMS))
    for plane in ma.PLANES:
        v = result["params"][plane]
        e = result["errors"].get(plane)
        row = f"{plane:6s}"
        for i in range(len(ma.PARAMS)):
            if e is not None and np.isfinite(e[i]):
                row += f"{v[i] * 1e3:>9.4f} +- {e[i] * 1e3:<7.4f}"
            else:
                row += f"{v[i] * 1e3:>9.4f} (fixed)  "
        print(row + ("   <- aligned" if plane in result["free_planes"] else "   (reference)"))


def process_batch(files, args, out, ref, make_plots, log_path, quiet=False, img_dir=None, prefix="", ctx=None):
    """One batch of tracks (one file in --per-run mode, otherwise all files together):
    fit (or evaluate the reference parameters), write outputs, append one row to the log.
    Returns the logged row (dict) or None if the batch was skipped."""
    out.mkdir(parents=True, exist_ok=True)
    img_dir = out if img_dir is None else Path(img_dir)
    img_dir.mkdir(parents=True, exist_ok=True)
    say = (lambda *a, **k: None) if quiet else print
    geom = ma.load_geometry(files[0], args.detector, args.yz_spacing)
    say(geom.describe())
    for f in files[1:]:
        if not ma.same_geometry(geom, ma.load_geometry(f, args.detector, args.yz_spacing)):
            raise SystemExit(f"{f}: the geometry in its METADATA differs from {files[0]}. Files with different "
                             f"geometries must not be pooled - use --per-run.")
    sample = ma.read_samples(files, geom, max_events=args.max_events, chi2_max=args.chi2_max,
                             first_event=args.first_event, corr=args.corr, forward=args.forward)
    say(f"Total: {len(sample):,} 3D 3p tracks" + (" (forward only, 90 < Phi_3p < 270)" if args.forward else ""))
    if len(sample) < 500:
        print(f"[skip] {', '.join(Path(f).name for f in files)}: only {len(sample)} tracks (need >= 500)")
        return None

    # --- reference point: nominal geometry; must reproduce the numbers stored in the file
    fits0 = ma.evaluate(geom, sample, ma.zero_params())
    dres, dchi = ma.check_against_stored(geom, sample, fits0)
    msg = f"Check of the track model vs the ROOT file: max|dres| = {dres:.1e} m, max|dchi2| = {dchi:.1e}"
    if dres >= 1e-9:
        msg += "   <-- MISMATCH: geometry/hit definition differs from the reconstruction!"
    if dres >= 1e-9 or not quiet:
        print(msg)

    if ref is not None:
        params = ma.params_from_json(ref["parameters"])
        free = ref["meta"]["free_planes"]
        fits1 = ma.evaluate(geom, sample, params)
        result = {"params": params,
                  "errors": {p: np.array([ref["parameters"][p].get(f"{n}_err", np.nan) for n in ma.PARAMS])
                             for p in free if p in ref["parameters"]},
                  "history": [ma.snapshot(0, params, fits1, geom)], "free_planes": free,
                  "free_params": ref["meta"]["free_params"], "converged": True}
        if ref["meta"].get("geometry") and ref["meta"]["geometry"] != geom.describe():
            print("[warn] the reconstruction geometry of this batch differs from the one the parameters were "
                  "derived with - the parameters are corrections ON TOP of the reco geometry, so they may not transfer.")
        say(f"Applying parameters from {args.params_in}")
        label, mode, solver = f"applied {args.params_in.name}", "apply", "-"
    else:
        free = ma.PLANES if args.free == ["all"] else [p.upper() for p in args.free]
        say(f"\nAligning {free}  (references: {[p for p in ma.PLANES if p not in free]}), "
            f"parameters: {args.params}")
        result = ma.run_alignment(geom, sample, free_planes=free, free_params=args.params,
                                  max_iter=args.max_iter, tol_len=args.tol_length, tol_ang=args.tol_angle,
                                  leverage=args.leverage, unbiased=args.unbiased, estimator=args.estimator,
                                  gauss_window=args.gauss_window, rigid_stations=args.rigid_station,
                                  verbose=not quiet)
        fits1 = result["fits"]
        solver = solver_name(args)
        label, mode = ("unbiased residuals" if args.unbiased else
                       "paper update" + (" + leverage weights" if args.leverage else "")), "fit"
        label += f", {args.estimator} residual"
        if not result["converged"]:
            print(f"[warn] {Path(files[0]).name}: alignment did not converge")

    raw0, raw1 = fits0, fits1          # biased residuals, for the slope plot
    if args.unbiased:
        fits0, fits1 = ma.unbiased_residuals(fits0), ma.unbiased_residuals(fits1)

    if not quiet:
        print_table(result)
        print(f"\nResiduals (measured - fit, {'unbiased' if args.unbiased else 'biased'}); "
              f"mean +- error of the mean, width = rms and robust sigma (1.4826*MAD):")
        print(f"{'plane':6s}{'mean before':>13s}{'mean after':>20s}{'rms before':>12s}{'rms after':>11s}"
              f"{'MAD before':>12s}{'MAD after':>11s}")
        for L in ma.LETTERS:
            for k in range(3):
                r0, r1 = fits0[L]["r"][:, k], fits1[L]["r"][:, k]
                print(f"{L}{k + 1:<5d}{r0.mean() * 1e6:10.1f} um{r1.mean() * 1e6:11.1f} +- {r1.std() / np.sqrt(len(r1)) * 1e6:4.1f} um"
                      f"{r0.std() * 1e3:9.4f} mm{r1.std() * 1e3:8.4f} mm"
                      f"{mad_sigma(r0) * 1e3:9.4f} mm{mad_sigma(r1) * 1e3:8.4f} mm")
        if args.estimator == "gauss":
            print(f"Gaussian core mean ({args.gauss_window} window): " + "  ".join(
                f"{L}{k + 1}: {ma.gauss_core_fit(fits0[L]['r'][:, k], args.gauss_window)['mu'] * 1e6:.1f} -> "
                f"{ma.gauss_core_fit(fits1[L]['r'][:, k], args.gauss_window)['mu'] * 1e6:.1f} um"
                for L in ma.LETTERS for k in range(3)))
        print(f"<chi2>  xz: {fits0['Z']['chi2'].mean():.4f} -> {fits1['Z']['chi2'].mean():.4f}   "
              f"xy: {fits0['Y']['chi2'].mean():.4f} -> {fits1['Y']['chi2'].mean():.4f}")

    write_outputs(out, geom, args, sample, files, result, fits0, fits1, prefix=prefix)
    row = alog.make_row(geom, files, sample, result, fits0, fits1, mode, solver, args.first_event,
                        args.max_events, params_in=args.params_in, tag=args.tag)
    alog.append_row(log_path, row, extra={"geometry": geom.describe(),
                                          "reco_equivalent_offsets_m": ma.reco_equivalent_offsets(geom, result["params"])})
    if quiet:
        worst = max(abs(fits1[L]["r"][:, k].mean()) for L in ma.LETTERS for k in range(3)) * 1e6
        print(f"run {row['run_min']}: {len(sample):>7,} tracks  {mode:5s} "
              + ("converged" if result["converged"] else "NOT CONVERGED")
              + f"  max |mean residual| before -> after: "
              f"{max(abs(fits0[L]['r'][:, k].mean()) for L in ma.LETTERS for k in range(3)) * 1e6:6.1f} -> {worst:6.1f} um")

    if make_plots:
        runs = np.unique(sample.run)
        run_txt = f"run {runs[0]}" if len(runs) == 1 else f"{len(runs)} runs ({runs.min()}-{runs.max()})"
        title = (f"{geom.detector}, {run_txt}, {len(sample):,} 3D 3p {'forward ' if args.forward else ''}tracks, aligned: "
                 f"{', '.join(result['free_planes'])} ({label})")
        ap.plot_residuals_before_after(fits0, fits1, img_dir / f"{prefix}residuals_3p_before_after.png", title)
        ap.plot_profiles(geom, fits0, fits1, result["free_planes"], img_dir, title, prefix=prefix)
        ap.plot_residual_vs_slope(raw0, raw1, img_dir / f"{prefix}residual_vs_slope_Z2_Y2.png", title,
                                  gauss_window=args.gauss_window)
        ap.plot_track_chi2(fits0, fits1, img_dir / f"{prefix}track_chi2_before_after.png", title)
        if args.estimator == "gauss":
            for f, tag in ((fits0, "before"), (fits1, "after")):
                ap.plot_gauss_fits(f, img_dir / f"{prefix}gauss_core_fits_{tag}.png", f"{title}, {tag} alignment",
                                   args.gauss_window)
        if ref is None:
            ap.plot_parameter_convergence(result["history"], result["free_planes"], result["free_params"],
                                          result["errors"], img_dir / f"{prefix}convergence_parameters.png", title)
            ap.plot_convergence_quality(result["history"], img_dir / f"{prefix}convergence_quality.png", title)
    if ctx is not None:
        ctx.update(geom=geom, result=result, sample=sample, fits0=fits0, fits1=fits1, row=row)
    return row


def main():
    args = parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    log_path = args.log or (out / "alignment_log.csv")
    print("=" * 78)
    print("MURAVES track-based alignment, planes 1-3 (physics/0306034)")

    ref = None
    if args.params_in:
        ref = json.loads(args.params_in.read_text())
    else:
        free = ma.PLANES if args.free == ["all"] else [p.upper() for p in args.free]
        bad = [p for p in free if p not in ma.PLANES]
        if bad:
            raise SystemExit(f"Unknown plane(s) {bad}; allowed: {ma.PLANES}")
        try:
            ma.check_free_planes(free, args.leverage, args.unbiased)
            ma.check_rigid_stations(args.rigid_station, free, args.estimator)
        except ValueError as exc:
            raise SystemExit(f"Error: {exc}")

    if args.check and args.per_run:
        raise SystemExit("--check trains on the pooled input files; it cannot be combined with --per-run")
    if args.check and args.params_in:
        raise SystemExit("--check fits first; to only apply existing parameters use --params-in with --per-run")

    if args.per_run:
        n_ok = 0
        for f in args.input:
            sub = out / "per_run" / Path(f).stem
            n_ok += process_batch([f], args, sub, ref, args.plots_per_run, log_path, quiet=not args.plots_per_run) is not None
        print(f"\n{n_ok}/{len(args.input)} files logged to {log_path.resolve()}")
    else:
        train_row = process_batch(args.input, args, out, ref, True, log_path)
        print(f"\nOutput written to {out.resolve()}")
        print(f"Log row appended to {log_path.resolve()}")
        if args.check:
            if train_row is None:
                raise SystemExit("Training batch was skipped; nothing to check.")
            run_check(args, out, log_path, train_row)


def run_check(args, out, log_path, train_row):
    """apply the parameters just fitted on the pooled training batch to every --check file separately"""
    ref_path = out / "alignment_parameters.json"
    ref = json.loads(ref_path.read_text())
    cargs = copy.copy(args)
    cargs.params_in = ref_path
    cargs.max_events = args.check_max_events
    cargs.first_event = 0
    cargs.tag = (args.tag + " " if args.tag else "") + "check"
    print("\n" + "=" * 78)
    print(f"Checking the fixed parameters on {len(args.check)} files, run by run "
          f"(training: {train_row['n_runs']:.0f} run(s), run {train_row['run_min']:.0f}-{train_row['run_max']:.0f})")
    rows = [train_row]
    for f in args.check:
        row = process_batch([f], cargs, out / "check_runs" / Path(f).stem, ref, args.plots_per_run,
                            log_path, quiet=not args.plots_per_run)
        if row is not None:
            rows.append(row)
    parsed = [alog.parse_row(r) for r in rows]
    applied = [r for r in parsed if r["mode"] == "apply"]
    if not applied:
        raise SystemExit("No check run could be evaluated.")
    alog.print_period_summary(alog.period_summary(applied), out / "period_check_summary.csv")
    ap_out = out / "period_plots"
    ap_out.mkdir(exist_ok=True)
    pal.plot_residual_means(parsed, "run", ap_out)
    pal.plot_chi2([r for r in parsed if r["mode"] == "apply"], "run", ap_out)
    print(f"\n{len(applied)} runs checked. Plots and period_check_summary.csv in {ap_out.resolve()} / {out.resolve()}")


if __name__ == "__main__":
    main()