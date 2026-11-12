#!/usr/bin/env python3
"""
Origin of the tails of the plane-2 residual: cluster sizes and number of candidate clusters vs track slope.

Same track selection as the alignment (best 3p track in BOTH projections, optional --forward); the residual is
the unbiased residual of plane 2 (hit minus the line through planes 1 and 3), with the alignment parameters of
--params-in applied (default: nominal geometry).  For every track it also reads, for planes 1-3 of the view:
    size and energy of the cluster on the best track   (ClusterSize_*, ClusterEnergy_* via TrackCluster_*_index)
    number of clusters on the plane                    (Nclusters_*)  -> >1 = the track finding had a choice

Questions it answers:
  * does the tail fraction (|e| > tail cut, default 3 x robust sigma) grow with |b|?
  * do the tails come with several candidate clusters (wrong cluster picked) or also with exactly one cluster on
    every plane (then: cluster position itself, e.g. split / merged / noisy clusters)?
  * which cluster sizes go with the tails, for steep vs moderate tracks?

Outputs (in --output-dir, prefix --tag):
  <tag>tails_vs_slope.png     per view: tail fraction, cluster size planes 1-3, fraction with >1 cluster planes 1-3,
                              and the mean residual split into "1 cluster on every plane" / "ambiguous" plus the
                              contribution of the tails to the mean
  <tag>tails_hist_by_size.png plane-2 residual for steep-, moderate, steep+ tracks, split by plane-2 cluster size
  <tag>tails_summary.csv      the numbers of the printed table

Example
  python tail_origin.py ~/MURAVES/test_data/alignment/RECONSTRUCTED/ROSSO/v1/MURAVES_AnalyzedData_run40[0-4]?.root \\
         --params-in geometry/split_ordered/alignment_ROSSO_runs4000-4100_ordered_parameters.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

import muraves_alignment as ma

CLASSES = ("steep-", "moderate", "steep+")
COL = {"steep-": "#3B6FB6", "moderate": "0.35", "steep+": "#D1495B"}


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="+", help="MURAVES_AnalyzedData_run*.root")
    p.add_argument("--params-in", type=Path, help="alignment parameters json (default: nominal geometry)")
    p.add_argument("--detector", choices=["NERO", "ROSSO", "BLU"])
    p.add_argument("--forward", action="store_true", help="only 90 < Phi_3p < 270")
    p.add_argument("--max-events", type=int, default=None, help="events per file")
    p.add_argument("--tail-cut", type=float, default=None,
                   help="|unbiased residual| [mm] above which a track is a tail (default: --tail-sigma x the robust "
                        "width 1.4826*MAD of the residual, per view)")
    p.add_argument("--tail-sigma", type=float, default=3.0, help="default tail cut in units of the robust width")
    p.add_argument("--steep", type=float, default=5.0,
                   help="steep tracks: b below the STEEP-th / above the (100-STEEP)-th percentile (per view)")
    p.add_argument("--moderate", nargs=2, type=float, default=[25.0, 75.0], metavar=("LO", "HI"),
                   help="moderate tracks: b between these percentiles (per view)")
    p.add_argument("--nbins", type=int, default=20, help="equal-statistics bins in b (1-99 percentile)")
    p.add_argument("--output-dir", type=Path, default=Path("images/tail_origin"))
    p.add_argument("--tag", default="", help="file name prefix")
    return p.parse_args()


def cluster_info(path, sample, max_events):
    """per track of `sample` (entries sample.event of `path`): size, energy of the track clusters and number of
    clusters on planes 1-3, for both views -> {L: {"size": (N,3), "energy": (N,3), "ncl": (N,3)}}"""
    br = ["BestTrack_3p_xz_index", "BestTrack_3p_xy_index"]
    for L in ma.LETTERS:
        for k in (1, 2, 3):
            br += [f"TrackCluster_{L.lower()}{k}_index", f"ClusterSize_{L}{k}", f"ClusterEnergy_{L}{k}",
                   f"Nclusters_{L}{k}"]
    with uproot.open(path) as f:
        arr = f[ma.TREE_NAME].arrays(br, entry_stop=max_events, library="ak")[sample.event]
    out = {}
    for L, view in (("Z", "xz"), ("Y", "xy")):
        best = ak.to_numpy(arr[f"BestTrack_3p_{view}_index"]).astype(np.int64)
        size, energy, ncl = [], [], []
        for k in (1, 2, 3):
            cl = ma._pick(arr[f"TrackCluster_{L.lower()}{k}_index"], best).astype(np.int64)
            size.append(ma._pick(arr[f"ClusterSize_{L}{k}"], cl))
            energy.append(ma._pick(arr[f"ClusterEnergy_{L}{k}"], cl))
            ncl.append(ak.to_numpy(arr[f"Nclusters_{L}{k}"]).astype(float))
        out[L] = {"size": np.stack(size, 1), "energy": np.stack(energy, 1), "ncl": np.stack(ncl, 1)}
    return out


def load(args):
    geom = ma.load_geometry(args.input[0], args.detector)
    samples, infos = [], []
    for p in args.input:
        try:
            s = ma.read_sample(p, geom, max_events=args.max_events, forward=args.forward)
        except ma.MissingBranchesError as exc:
            print(f"  [skip] {exc.path.name}: missing branches")
            continue
        samples.append(s)
        infos.append(cluster_info(p, s, args.max_events))
        print(f"  {Path(p).name}: {len(s):,} tracks")
    if not samples:
        raise SystemExit("no usable file")
    cat = lambda f: np.concatenate([f(s) for s in samples])  # noqa: E731
    sample = ma.TrackSample({L: cat(lambda s, L=L: s.u[L]) for L in ma.LETTERS},
                            {L: cat(lambda s, L=L: s.stored_res[L]) for L in ma.LETTERS},
                            {q: cat(lambda s, q=q: s.stored_chi2[q]) for q in ("xz", "xy")},
                            cat(lambda s: s.run), cat(lambda s: s.event))
    info = {L: {key: np.concatenate([i[L][key] for i in infos]) for key in ("size", "energy", "ncl")}
            for L in ma.LETTERS}
    return geom, sample, info


def b_classes(b, args):
    q = np.percentile(b, [args.steep, args.moderate[0], args.moderate[1], 100 - args.steep])
    return {"steep-": b < q[0], "moderate": (b >= q[1]) & (b <= q[2]), "steep+": b > q[3]}, q


def summary_rows(L, b, e, inf, cls, cut):
    tail = np.abs(e) > cut
    clean = (inf["ncl"] == 1).all(1)          # exactly one cluster on planes 1, 2, 3: no choice for the tracking
    rows = []
    for c in CLASSES:
        m = cls[c]
        s2 = inf["size"][:, 1]
        row ={"view": L, "class": c, "b_range": f"{b[m].min():+.2f}..{b[m].max():+.2f}", "n": int(m.sum()),
               "tail_frac_%": 100 * tail[m].mean(), "mean_e_mm": e[m].mean(),
               "mean_e_no_tail_mm": e[m & ~tail].mean(),
               "tail_contrib_to_mean_mm": (e * tail)[m].mean(),
               "frac_1cluster_all_planes_%": 100 * clean[m].mean(),
               "tail_frac_1cluster_%": 100 * tail[m & clean].mean(),
               "tail_frac_ambiguous_%": 100 * tail[m & ~clean].mean(),
               "tail_contrib_1cluster_mm": (e * (tail & clean))[m].mean(),
               "tail_contrib_ambiguous_mm": (e * (tail & ~clean))[m].mean()}
        for k in range(3):
            row[f"size_p{k + 1}"] = np.nanmean(inf["size"][m, k])
            row[f"ncl>1_p{k + 1}_%"] = 100 * (inf["ncl"][m, k] > 1).mean()
        for sz, lab in ((1, "1"), (2, "2"), (3, "3+")):
            sel = m & ((s2 >= 3) if sz == 3 else (inf["size"][:, 1] == sz))
            row[f"p2_size{lab}_%"] = 100 * sel.sum() / max(m.sum(), 1)
            row[f"tail_frac_p2_size{lab}_%"] = 100 * tail[sel].mean() if sel.any() else np.nan
        row["tail_cut_mm"] = cut
        rows.append(row)
    return rows


def print_rows(rows):
    print(f"\n{'':4s}{'class':9s}{'b range':>14s}{'N':>9s}{'tail%':>7s}{'<e>':>7s}{'<e>core':>8s}"
          f"{'tails->':>8s}{'1cl%':>6s}{'tail%|1cl':>10s}{'tail%|amb':>10s}{'->1cl':>7s}{'->amb':>7s}"
          f"{'size p1/p2/p3':>16s}{'ncl>1% p1/p2/p3':>18s}{'tail% by p2 size 1/2/3+':>26s}")
    for r in rows:
        print(f"{r['view']:4s}{r['class']:9s}{r['b_range']:>14s}{r['n']:9d}{r['tail_frac_%']:7.2f}"
              f"{r['mean_e_mm']:7.2f}{r['mean_e_no_tail_mm']:8.2f}{r['tail_contrib_to_mean_mm']:8.2f}"
              f"{r['frac_1cluster_all_planes_%']:6.1f}{r['tail_frac_1cluster_%']:10.2f}{r['tail_frac_ambiguous_%']:10.2f}"
              f"{r['tail_contrib_1cluster_mm']:7.2f}{r['tail_contrib_ambiguous_mm']:7.2f}"
              f"{r['size_p1']:7.2f}/{r['size_p2']:.2f}/{r['size_p3']:.2f}"
              f"{r['ncl>1_p1_%']:8.1f}/{r['ncl>1_p2_%']:.1f}/{r['ncl>1_p3_%']:.1f}"
              f"{r['tail_frac_p2_size1_%']:12.2f}/{r['tail_frac_p2_size2_%']:.2f}/{r['tail_frac_p2_size3+_%']:.2f}")
    print("  <e>: mean unbiased plane-2 residual [mm]; <e>core: without the tails; tails->: contribution of the tails "
          "to <e>;\n  1cl: exactly one cluster on planes 1-3 (no choice for the track finding), amb: ambiguous; "
          "->1cl / ->amb: tail contribution to <e> from those events")


def binned(b, y, edges):
    idx = np.digitize(b, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    return np.array([np.nanmean(y[ok & (idx == i)]) if (ok & (idx == i)).any() else np.nan
                     for i in range(len(edges) - 1)])


def binned_mean_err(b, y, edges):
    idx = np.digitize(b, edges) - 1
    m, s = np.full(len(edges) - 1, np.nan), np.full(len(edges) - 1, np.nan)
    for i in range(len(edges) - 1):
        yi = y[idx == i]
        if len(yi) > 10:
            m[i], s[i] = yi.mean(), yi.std() / np.sqrt(len(yi))
    return m, s


def plot_vs_slope(data, cuts, args, out, title):
    fig, axs = plt.subplots(4, 2, figsize=(13, 14), sharex="col")
    for c, (L, (b, e, inf)) in enumerate(data.items()):
        lo, hi = np.percentile(b, [1, 99])
        edges = np.unique(np.quantile(b[(b >= lo) & (b <= hi)], np.linspace(0, 1, args.nbins + 1)))
        x = 0.5 * (edges[1:] + edges[:-1])
        tail = np.abs(e) > cuts[L]
        clean = (inf["ncl"] == 1).all(1)
        ax = axs[0, c]
        ax.plot(x, 100 * binned(b, tail.astype(float), edges), "o-", color="k", ms=4, label="all")
        ax.plot(x, 100 * binned(b[clean], tail[clean].astype(float), edges), "s--", color="#2A9D5C", ms=4,
                mfc="white", label="1 cluster on planes 1-3")
        ax.plot(x, 100 * binned(b[~clean], tail[~clean].astype(float), edges), "^--", color="#D1495B", ms=4,
                mfc="white", label=">1 cluster on some plane")
        ax.set_ylabel(f"tail fraction |e|>{cuts[L]:.1f} mm [%]")
        ax.set_title(f"{L}2 ({ma.PROJ_OF[L]} view)", loc="left", fontweight="bold")
        ax.legend(fontsize=8)
        for k, mk in zip(range(3), "osD"):
            axs[1, c].plot(x, binned(b, inf["size"][:, k], edges), mk + "-", ms=4, label=f"plane {L}{k + 1}")
            axs[2, c].plot(x, 100 * binned(b, (inf["ncl"][:, k] > 1).astype(float), edges), mk + "-", ms=4,
                           label=f"plane {L}{k + 1}")
        axs[1, c].set_ylabel("mean size of the track cluster [strips]")
        axs[2, c].set_ylabel("events with >1 cluster on the plane [%]")
        axs[1, c].legend(fontsize=8)
        axs[2, c].legend(fontsize=8)
        ax = axs[3, c]
        for sel, lab, col, mk in ((np.ones_like(clean), "all", "k", "o"), (clean, "1 cluster on planes 1-3", "#2A9D5C", "s"),
                                  (~clean, ">1 cluster on some plane", "#D1495B", "^")):
            m, s = binned_mean_err(b[sel], e[sel], edges)
            ax.errorbar(x, m, yerr=s, fmt=mk + "-", ms=4, color=col, lw=1, capsize=0, label=f"mean, {lab}")
        idx = np.digitize(b, edges) - 1
        contrib = np.array([(e * tail)[idx == i].mean() if (idx == i).any() else np.nan for i in range(len(x))])
        ax.plot(x, contrib, "x:", color="#7A4BB0", label="contribution of the tails to the mean (all)")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylabel(f"unbiased residual {L}2 [mm]")
        ax.set_xlabel(f"track slope b ({ma.PROJ_OF[L]} projection)")
        ax.legend(fontsize=7.5)
        for a in axs[:, c]:
            a.grid(True, alpha=0.3, ls="--")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=120)
    plt.close(fig)


def plot_hists(data, cls_by_view, cuts, out, title):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5), sharey=True)
    for r, (L, (b, e, inf)) in enumerate(data.items()):
        lim = 1.05 * np.percentile(np.abs(e), 99.5)
        edges = np.linspace(-lim, lim, 81)
        s2 = inf["size"][:, 1]
        for c, cl in enumerate(CLASSES):
            ax = axs[r, c]
            m = cls_by_view[L][cl]
            for sel, lab, col in ((m, "all", "k"), (m & (s2 == 1), "p2 size 1", "#3B6FB6"),
                                  (m & (s2 == 2), "p2 size 2", "#2A9D5C"), (m & (s2 >= 3), "p2 size 3+", "#D1495B")):
                if sel.sum() < 10:
                    continue
                n, _ = np.histogram(e[sel], edges)
                ax.stairs(n / sel.sum(), edges, color=col, lw=1.3,
                          label=f"{lab}: {100 * sel.sum() / m.sum():.0f}%, mean {e[sel].mean():+.2f}, "
                                f"tail {100 * (np.abs(e[sel]) > cuts[L]).mean():.1f}%")
            ax.set_yscale("log")
            ax.axvline(0, color="k", lw=0.6)
            ax.set_title(f"{L}2, {cl} (b {b[m].min():+.2f}..{b[m].max():+.2f})", fontsize=9, loc="left")
            ax.set_xlabel(f"unbiased residual {L}2 [mm]")
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3, ls="--")
        axs[r, 0].set_ylabel("fraction of tracks / bin")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=120)
    plt.close(fig)


def main():
    args = parse()
    geom, sample, info = load(args)
    params = ma.zero_params()
    if args.params_in:
        params = ma.params_from_json(json.loads(args.params_in.read_text())["parameters"])
    fits = ma.evaluate(geom, sample, params)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data, cls_by_view, rows, cuts = {}, {}, [], {}
    for L in ma.LETTERS:
        b = fits[L]["b"]
        e = fits[L]["r"][:, 1] / (1.0 - fits[L]["h"][:, 1]) * 1e3
        data[L] = (b, e, info[L])
        cls_by_view[L], _ = b_classes(b, args)
        sig = 1.4826 * np.median(np.abs(e - np.median(e)))
        cuts[L] = args.tail_cut if args.tail_cut is not None else args.tail_sigma * sig
        print(f"{L}2: robust width of the unbiased residual {sig:.2f} mm, rms {e.std():.2f} mm -> tail cut |e| > {cuts[L]:.2f} mm")
        rows += summary_rows(L, b, e, info[L], cls_by_view[L], cuts[L])
    print_rows(rows)
    runs = np.unique(sample.run)
    title = (f"{geom.detector}, runs {runs.min()}-{runs.max()} ({len(runs)}), {len(sample):,} 3D 3p "
             f"{'forward ' if args.forward else ''}tracks, "
             f"{'params ' + args.params_in.name if args.params_in else 'nominal geometry'}; "
             f"steep: outer {args.steep:g}%, moderate: {args.moderate[0]:g}-{args.moderate[1]:g}% of b")
    plot_vs_slope(data, cuts, args, args.output_dir / f"{args.tag}tails_vs_slope.png", title)
    plot_hists(data, cls_by_view, cuts, args.output_dir / f"{args.tag}tails_hist_by_size.png", title)
    with open(args.output_dir / f"{args.tag}tails_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()} for r in rows)
    print(f"\nPlots and tails_summary.csv in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
