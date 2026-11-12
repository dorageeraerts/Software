#!/usr/bin/env python3
"""
Effect of a tilt (or shift) of ONE plane on the reconstructed 3p tracks, theta/phi and the image,
when the reconstruction ignores it (the reco cfg has no tilt columns, see geometry_io.py).

Each plane in --plane is tilted on its own (the others stay flat) by --angle about --param:
    dalpha / dbeta / dgamma = rotation about u / v / w of that plane (muraves_alignment.py conventions)
    du / dw                 = shift along the measured coordinate / along the telescope axis

    Delta = (reco with flat planes) - (correct reconstruction)

data mode (--runs):  stored 3p tracks of the ROOT files.  "Correct" = the same hits refitted with the
                     tilt applied through the linearised model of muraves_alignment.corrected_hits
                     (small angles only).
toy mode (--toy N):  straight tracks (flat in slope, inside the acceptance), hits computed with the
                     EXACT rotated plane (any angle).  "Correct" = the same tracks with the flat plane.
                     The "linear" rows refit the tilted toy hits with the linearised model; they must
                     agree with the exact rows for small angles.

Lever arm: with planes at x ~ (-0.265, 0, 0.262) m and equal weights, a shift d of plane k changes the
slope by (x_k - <x>) d / Sxx.  For the middle plane x_2 - <x> ~ 0: a Z2/Y2 tilt only shows up in the
residuals / chi2, not in the angles.  Tilts of planes 1 and 3 move theta/phi.

Delta phi is taken modulo 180 deg (= Delta atan(slope_xy)), so a track near theta = 0 whose up/down
side flips does not show up as a 180 deg jump.

    python tilt_effect.py --plane all --param dbeta --angle 7.1275 --unit mrad --runs 4000 4100
    python tilt_effect.py --plane Z1 Z2 Z3 --param dbeta --angle 1 --unit deg --toy 200000
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np

import geometry_io as gio
import muraves_alignment as ma
from alignment_plots import _profile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tracks_reconstruction" / "4_MainRecSof"))
from Evaluate_Angular_Coordinates import TrackAngularCoordinates  # noqa: E402

UNITS = {"rad": 1.0, "mrad": 1e-3, "deg": np.pi / 180, "m": 1.0, "mm": 1e-3}
LENGTH_UNITS = ("m", "mm")
C_PLANE = ["#2a78d6", "#eb6834", "#1baf7a"]   # plane 1 / 2 / 3
C_NOMINAL = "#8a8986"
PLANE_HALF = 0.528                            # half size of a plane [m] (2 x 32 strips of 16.5 mm)
LINEAR_LIMIT = 0.05                           # rad, beyond this the linearised model is not reliable
track_angles = np.vectorize(lambda *a: tuple(TrackAngularCoordinates(*a)), otypes=[float, float])


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plane", nargs="+", required=True, help="Z1 Z2 Z3 Y1 Y2 Y3 or all")
    parser.add_argument("--param", required=True, choices=ma.PARAMS)
    parser.add_argument("--angle", type=float, required=True, help="size of the tilt (or shift)")
    parser.add_argument("--unit", default="mrad", choices=UNITS)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--runs", nargs=2, type=int, metavar=("FIRST", "LAST"))
    src.add_argument("--toy", type=int, metavar="N", help="number of toy tracks")
    parser.add_argument("--data-dir", type=Path,
                        default=Path("~/MURAVES/test_data/alignment/RECONSTRUCTED/ROSSO/v1"))
    parser.add_argument("--good-runs", type=Path, help="text file with the good runs")
    parser.add_argument("--max-events", type=int, help="per file")
    parser.add_argument("--detector", default="ROSSO", choices=("ROSSO", "NERO", "BLU"))
    parser.add_argument("--no-smear", action="store_true", help="toy: no hit resolution smearing")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "images" / "tilt_effect")
    args = parser.parse_args()
    args.plane = ma.PLANES if args.plane == ["all"] else args.plane
    bad = set(args.plane) - set(ma.PLANES)
    if bad:
        parser.error(f"unknown plane(s) {sorted(bad)}")
    if (args.param in ("du", "dw")) != (args.unit in LENGTH_UNITS):
        parser.error(f"--param {args.param} does not go with --unit {args.unit}")
    return args


# ----------------------------------------------------------------------------
# tracks
# ----------------------------------------------------------------------------
def make_sample(u):
    n = len(u["Z"])
    return ma.TrackSample(u, {L: np.zeros((n, 3)) for L in ma.LETTERS}, {p: np.zeros(n) for p in ("xz", "xy")},
                          np.zeros(n, dtype=int), np.arange(n))


def toy_geometry(detector):
    g = ma.FALLBACK_GEOMETRY[detector]
    return ma.Geometry(detector, np.array(g["x_pos"][:3], float), np.array(g["z_add"][:3], float),
                       np.zeros(3), {"Z": 0.004, "Y": 0.0035}, -0.528 + 31.5 * 0.0165)


def rotation(param, angle):
    """plane local (u, v, w) -> global (u, v, w); first order identical to ma.corrected_hits"""
    c, s = np.cos(angle), np.sin(angle)
    if param == "dalpha":   # about u: v tilts into w   (w = dalpha * v)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if param == "dbeta":    # about v: u tilts into w   (w = dbeta * u)
        return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])
    if param == "dgamma":   # about w: v tilts into -u  (u = -dgamma * v)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    return np.eye(3)


def toy_tracks(geom, n, rng):
    """lines z = z0 + bz x, y = y0 + by x (flat in slope) that cross all 6 flat planes"""
    smax = 2 * PLANE_HALF / (geom.x_pos[2] - geom.x_pos[0])
    c = {"Z": geom.u_centre["Z"].mean(), "Y": geom.u_centre["Y"].mean()}
    parts, n_have = [], 0
    while n_have < n:
        m = 4 * n
        t = {"z0": rng.uniform(-PLANE_HALF, PLANE_HALF, m) + c["Z"], "bz": rng.uniform(-smax, smax, m),
             "y0": rng.uniform(-PLANE_HALF, PLANE_HALF, m) + c["Y"], "by": rng.uniform(-smax, smax, m)}
        ok = np.ones(m, bool)
        for L, (a, b) in (("Z", ("z0", "bz")), ("Y", ("y0", "by"))):
            for k in range(3):
                for M, (a2, b2) in (("Z", ("z0", "bz")), ("Y", ("y0", "by"))):
                    pos = t[a2] + t[b2] * geom.x_of[L][k]
                    centre = geom.u_centre[L][k] if M == L else geom.v_centre[L][k]
                    ok &= np.abs(pos - centre) < PLANE_HALF
        parts.append({key: v[ok] for key, v in t.items()})
        n_have += ok.sum()
    return {key: np.concatenate([p[key] for p in parts])[:n] for key in parts[0]}


def toy_hits(geom, t, tilted=None, param=None, value=0.0):
    """raw hit positions u (N,3) per projection; plane `tilted` rotated/shifted exactly"""
    line = {"Z": (t["z0"], t["bz"], t["y0"], t["by"]), "Y": (t["y0"], t["by"], t["z0"], t["bz"])}
    u = {}
    for L in ma.LETTERS:
        u0, bu, v0, bv = line[L]
        cols = []
        for k in range(3):
            centre = np.array([geom.u_centre[L][k], geom.v_centre[L][k], geom.x_of[L][k]])
            R = np.eye(3)
            if f"{L}{k + 1}" == tilted:
                R = rotation(param, value)
                centre += {"du": [value, 0, 0], "dw": [0, 0, value]}.get(param, [0, 0, 0])
            p0 = np.stack([u0, v0, np.zeros_like(u0)], 1)
            d = np.stack([bu, bv, np.ones_like(u0)], 1)
            normal = R[:, 2]
            s = ((centre - p0) @ normal) / (d @ normal)
            local = (p0 + s[:, None] * d - centre) @ R      # = R^T (p - centre)
            cols.append(geom.u_centre[L][k] + local[:, 0])
        u[L] = np.stack(cols, 1)
    return u


# ----------------------------------------------------------------------------
# comparison
# ----------------------------------------------------------------------------
def angles(geom, fits):
    """theta, phi [deg] exactly as MURAVES_reco.py (slope_xy = b_Y, slope_xz = b_Z)"""
    return track_angles(fits["Y"]["b"], fits["Z"]["b"], geom.x_pos[0], geom.x_pos[2])


def compare(geom, reco, true, plane, model):
    th_r, ph_r = angles(geom, reco)
    th_t, ph_t = angles(geom, true)
    L = plane[0]
    return {"plane": plane, "model": model, "L": L,
            "theta": th_t, "phi": np.degrees(np.arctan(true["Y"]["b"])),
            "theta_reco": th_r, "phi_reco": np.degrees(np.arctan(reco["Y"]["b"])),
            "dtheta": np.radians(th_r - th_t) * 1e3,
            "dphi": np.radians((ph_r - ph_t + 90) % 180 - 90) * 1e3,
            "chi2_reco": reco[L]["chi2"], "chi2_true": true[L]["chi2"],
            "res_reco": reco[L]["r"].mean(0) * 1e6, "res_true": true[L]["r"].mean(0) * 1e6}


def print_table(rows, label):
    print(f"\n{label}   (Delta = reco with flat planes - correct;  angles in mrad, residuals in um)")
    print(f"{'plane':6} {'model':7} {'<dtheta>':>9} {'rms':>8} {'<dphi>':>9} {'rms':>8} "
          f"{'<chi2> correct':>15} {'<chi2> reco':>12} {'<res> reco planes 1/2/3':>27}")
    for r in rows:
        res = " ".join(f"{x:8.1f}" for x in r["res_reco"])
        print(f"{r['plane']:6} {r['model']:7} {r['dtheta'].mean():9.4f} {r['dtheta'].std():8.4f} "
              f"{r['dphi'].mean():9.4f} {r['dphi'].std():8.4f} {r['chi2_true'].mean():15.3f} "
              f"{r['chi2_reco'].mean():12.3f} {res:>27}")


# ----------------------------------------------------------------------------
# plots
# ----------------------------------------------------------------------------
def plot_comparison(rows, L, output_file, title, nbins=20):
    rows = [r for r in rows if r["L"] == L]
    fig, ax = plt.subplots(3, 2, figsize=(12, 12))
    proj = ma.PROJ_OF[L]
    for key, a, xlab in (("dtheta", ax[0, 0], r"$\Delta\theta$ [mrad]"), ("dphi", ax[0, 1], r"$\Delta\phi$ [mrad]")):
        vals = np.concatenate([r[key] for r in rows])
        lo, hi = np.percentile(vals, [0.5, 99.5])
        edges = np.linspace(min(lo, -1e-3), max(hi, 1e-3), 81)
        for r in rows:
            k = int(r["plane"][1]) - 1
            a.hist(np.clip(r[key], edges[0], edges[-1]), edges, histtype="step", lw=2, color=C_PLANE[k],
                   label=f"{r['plane']}: mean {r[key].mean():.3f}, rms {r[key].std():.3f}")
        a.set_xlabel(xlab)
        a.set_ylabel("tracks")
        a.legend(fontsize=8, frameon=False)
    for key, var, a, xlab, ylab in (("dtheta", "theta", ax[1, 0], r"$\theta$ [deg]", r"mean $\Delta\theta$ [mrad]"),
                                    ("dphi", "phi", ax[1, 1], r"$\phi$ mod 180 [deg]", r"mean $\Delta\phi$ [mrad]")):
        edges = np.linspace(*np.percentile(rows[0][var], [0.5, 99.5]), nbins + 1)
        centres = 0.5 * (edges[1:] + edges[:-1])
        for r in rows:
            k = int(r["plane"][1]) - 1
            m, sem, n = _profile(r[var], r[key], edges)
            a.errorbar(centres, m, sem, fmt="o-", ms=4, lw=2, color=C_PLANE[k], label=r["plane"])
        a.axhline(0, color=C_NOMINAL, lw=1)
        a.set_xlabel(xlab)
        a.set_ylabel(ylab)
        a.legend(fontsize=8, frameon=False)
    edges = np.linspace(0, np.percentile(rows[0]["chi2_true"], 99), 81)
    a = ax[2, 0]
    a.hist(rows[0]["chi2_true"], edges, histtype="step", lw=2, color=C_NOMINAL, label="correct geometry")
    for r in rows:
        k = int(r["plane"][1]) - 1
        a.hist(r["chi2_reco"], edges, histtype="step", lw=2, color=C_PLANE[k],
               label=f"{r['plane']} tilted, flat reco: <chi2> {r['chi2_reco'].mean():.2f}")
    a.set_xlabel(rf"$\chi^2_{{{proj}}}$")
    a.set_ylabel("tracks")
    a.set_yscale("log")
    a.legend(fontsize=8, frameon=False)
    a = ax[2, 1]
    width = 0.8 / len(rows)
    for i, r in enumerate(rows):
        k = int(r["plane"][1]) - 1
        a.bar(np.arange(3) + (i - (len(rows) - 1) / 2) * width, r["res_reco"] - r["res_true"], width * 0.9,
              color=C_PLANE[k], label=f"{r['plane']} tilted")
    a.axhline(0, color=C_NOMINAL, lw=1)
    a.set_xticks(range(3), [f"{L}{k}" for k in (1, 2, 3)])
    a.set_ylabel(r"$\Delta$ mean residual [$\mu$m]")
    a.legend(fontsize=8, frameon=False)
    for a in ax.flat:
        a.grid(alpha=0.3)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {output_file}")


def plot_image(r, output_file, title, bin_deg=1.0, min_count=10):
    phi_e = np.arange(np.floor(np.percentile(r["phi"], 0.5)), np.ceil(np.percentile(r["phi"], 99.5)) + bin_deg, bin_deg)
    th_e = np.arange(0, np.ceil(np.percentile(r["theta"], 99.5)) + bin_deg, bin_deg)
    h_true, _, _ = np.histogram2d(r["phi"], r["theta"], [phi_e, th_e])
    h_reco, _, _ = np.histogram2d(r["phi_reco"], r["theta_reco"], [phi_e, th_e])
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(h_true >= min_count, h_reco / h_true, np.nan)
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    vmax = max(h_true.max(), h_reco.max())
    for a, h, lab in ((ax[0], h_true, "correct geometry"), (ax[1], h_reco, f"{r['plane']} tilted, flat reco")):
        im = a.pcolormesh(phi_e, th_e, h.T, cmap="Blues", vmin=0, vmax=vmax)
        fig.colorbar(im, ax=a, label="tracks")
        a.set_title(lab)
    dev = np.nanpercentile(np.abs(ratio - 1), 99) if np.isfinite(ratio).any() else 0.1
    im = ax[2].pcolormesh(phi_e, th_e, ratio.T, cmap="RdBu_r",
                          norm=matplotlib.colors.TwoSlopeNorm(1.0, 1 - dev - 1e-9, 1 + dev + 1e-9))
    fig.colorbar(im, ax=ax[2], label="reco / correct")
    ax[2].set_title(f"ratio (bins with >= {min_count} tracks)")
    for a in ax:
        a.set_xlabel(r"$\phi$ mod 180 [deg]")
        a.set_ylabel(r"$\theta$ [deg]")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {output_file}")


# ----------------------------------------------------------------------------
def main() -> None:
    args = parse_args()
    value = args.angle * UNITS[args.unit]
    tilt = args.param not in ("du", "dw")
    if tilt and abs(value) > LINEAR_LIMIT:
        print(f"[warn] |{args.angle} {args.unit}| > {LINEAR_LIMIT * 1e3:.0f} mrad: the linearised model "
              f"({'data' if args.runs else 'linear rows'}) is not reliable, use the exact toy rows")

    def params_for(plane):
        p = ma.zero_params()
        p[plane][ma.PARAMS.index(args.param)] = value
        return p

    rows = []
    if args.runs:
        first, last = sorted(args.runs)
        good = gio.parse_run_list(args.good_runs)[0] if args.good_runs else None
        files, _ = gio.select_runs(args.data_dir, first, last, good, pattern="MURAVES_AnalyzedData_run*.root")
        if not files:
            raise SystemExit("no ROOT files selected")
        geom = ma.load_geometry(files[0][1], args.detector)
        sample = ma.read_samples([p for _, p in files], geom, max_events=args.max_events)
        reco = ma.evaluate(geom, sample, ma.zero_params())
        dres, dchi = ma.check_against_stored(geom, sample, reco)
        print(f"{len(sample):,} tracks; flat refit vs stored: max |dres| {dres * 1e6:.2g} um, max |dchi2| {dchi:.2g}")
        for plane in args.plane:
            rows.append(compare(geom, reco, ma.evaluate(geom, sample, params_for(plane)), plane, "data"))
        source = f"runs{first}-{last}"
    else:
        geom = toy_geometry(args.detector)
        rng = np.random.default_rng(args.seed)
        t = toy_tracks(geom, args.toy, rng)
        noise = {L: (0 if args.no_smear else geom.sigma[L]) * rng.standard_normal((args.toy, 3)) for L in ma.LETTERS}
        flat = {L: u + noise[L] for L, u in toy_hits(geom, t).items()}
        true = ma.evaluate(geom, make_sample(flat), ma.zero_params())
        for plane in args.plane:
            hits = {L: u + noise[L] for L, u in toy_hits(geom, t, plane, args.param, value).items()}
            sample = make_sample(hits)
            reco = ma.evaluate(geom, sample, ma.zero_params())
            rows.append(compare(geom, reco, true, plane, "exact"))
            rows.append(compare(geom, reco, ma.evaluate(geom, sample, params_for(plane)), plane, "linear"))
        source = f"toy{args.toy}"

    label = f"{geom.detector} {source}: {args.param} = {args.angle:g} {args.unit}"
    print_table(rows, label)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_dir / f"tilt_{geom.detector}_{source}_{args.param}_{args.angle:g}{args.unit}_"
    main_rows = [r for r in rows if r["model"] != "linear"]
    for L in ma.LETTERS:
        if any(r["L"] == L for r in main_rows):
            plot_comparison(main_rows, L, Path(f"{prefix}{L}planes.png"), f"{label}, one {L} plane tilted at a time")
    for r in main_rows:
        plot_image(r, Path(f"{prefix}{r['plane']}_image.png"), f"{label}, {r['plane']} tilted")


if __name__ == "__main__":
    main()
