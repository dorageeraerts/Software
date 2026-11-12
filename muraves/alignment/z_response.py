#!/usr/bin/env python3
"""
Position response of the plane-2 clusters vs the incidence angle (track slope b), Z view vs Y view.

Reference position: the track predicted at plane 2 from planes 1 and 3 only (u_pred = hit - unbiased residual), with
the alignment parameters of --params-in.  eta = position of u_pred inside the strip pitch (-0.5 .. 0.5, 0 = strip
centre).  Same track selection as the alignment (best 3p track in both views, optional --forward).
Cluster position in the reco = energy-weighted mean of the strip positions; clusters are runs of strips >= s1 (6);
single-strip clusters (>= s2 = 10, plus neighbours > s3 = 2) get a uniform +-pitch/2 random smearing.

Figures (in --output-dir, prefix --tag); every printed number says which figure / panel it belongs to:
  fig1_residual_vs_eta.png        core mean of the plane-2 residual vs eta, one curve per slope class
                                  (a flat line at 0 = linear response; an S-shape = non-linear centroid)
  fig2_size_split_energy_vs_b.png per view, vs b: (a) fraction of track-cluster size 1/2/3/4+,
                                  (b) fraction with a 2nd cluster on plane 2 within 3 strips ("split cluster"),
                                  (c) track-cluster energy / path length (energy / sqrt(1+b^2))
  fig3_strip_profile.png          energy-weighted profile of ALL plane-2 cluster strips vs (strip - u_pred) in pitch
                                  units, per slope class (the response function; its centroid = the position bias)
  fig4_residual_by_size_vs_b.png  core mean of the plane-2 residual vs b, split by track-cluster size 1/2/3+ and
                                  for "split" events
  fig5_residual_vs_b_by_parity.png core mean of the plane-2 residual vs b for a predicted hit on an even / odd strip
                                  (alternating triangular bars: opposite slopes = different depth of the two bar types)
  response_summary.csv            the printed table
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

B_EDGES = [-2.2, -1.3, -0.8, -0.3, 0.3, 0.8, 1.3, 2.2]
B_LABELS = ["b<-1.3", "-1.3..-0.8", "-0.8..-0.3", "|b|<0.3", "0.3..0.8", "0.8..1.3", "b>1.3"]
B_COLORS = ["#08306B", "#2171B5", "#6BAED6", "0.3", "#FC9272", "#EF3B2C", "#67000D"]
SIZE_COL = {"1": "#3B6FB6", "2": "#2A9D5C", "3+": "#D1495B", "split": "#7A4BB0"}


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="+", help="MURAVES_AnalyzedData_run*.root")
    p.add_argument("--params-in", type=Path, help="alignment parameters json (default: nominal geometry)")
    p.add_argument("--detector", choices=["NERO", "ROSSO", "BLU"])
    p.add_argument("--forward", action="store_true", help="only 90 < Phi_3p < 270")
    p.add_argument("--max-events", type=int, default=None, help="events per file")
    p.add_argument("--core-sigma", type=float, default=3.0,
                   help="core = |residual| < CORE_SIGMA x robust width (1.4826*MAD), per view")
    p.add_argument("--split-strips", type=int, default=3,
                   help="a 2nd plane-2 cluster closer than this many strips to the track cluster = split cluster")
    p.add_argument("--output-dir", type=Path, default=Path("images/z_response"))
    p.add_argument("--tag", default="")
    return p.parse_args()


def plane2_info(path, sample, L, max_events, pitch, split_strips):
    """per track: size, energy, strip ids/energies of the track cluster on plane 2 of view L, all strips of all
    plane-2 clusters, and whether another cluster lies within split_strips strips"""
    view = "xz" if L == "Z" else "xy"
    br = [f"BestTrack_3p_{view}_index", f"TrackCluster_{L.lower()}2_index", f"ClusterSize_{L}2",
          f"ClusterEnergy_{L}2", f"ClusterPosition_{L}2", f"StripsEnergy_{L}2", f"StripsID_{L}2"]
    with uproot.open(path) as f:
        a = f[ma.TREE_NAME].arrays(br, entry_stop=max_events, library="ak")[sample.event]
    best = ak.to_numpy(a[f"BestTrack_3p_{view}_index"]).astype(np.int64)
    cl = ma._pick(a[f"TrackCluster_{L.lower()}2_index"], best).astype(np.int64)
    size = ma._pick(a[f"ClusterSize_{L}2"], cl)
    energy = ma._pick(a[f"ClusterEnergy_{L}2"], cl)
    cpos = a[f"ClusterPosition_{L}2"]
    tpos = ma._pick(cpos, cl)
    others = cpos[ak.local_index(cpos, axis=1) != cl[:, None]]
    dmin = ak.to_numpy(ak.fill_none(ak.min(abs(others - ak.Array(tpos)[:, None]), axis=1), np.inf)).astype(float)
    split = dmin < split_strips * pitch
    return {"size": size, "energy": energy, "split": split,
            "strip_id": a[f"StripsID_{L}2"], "strip_e": a[f"StripsEnergy_{L}2"]}


def core_mean(e, sel):
    x = e[sel]
    if len(x) < 30:
        return np.nan, np.nan
    return x.mean(), x.std() / np.sqrt(len(x))


def main():
    args = parse()
    geom = ma.load_geometry(args.input[0], args.detector)
    meta = ma.read_metadata(args.input[0]) or {}
    cl_cfg = next((c["parsed_content"].get("cluster_lists", {}) for c in meta.get("detector_config_files", [])
                   if c["path"].split("/")[-1].startswith("reco_parameters")), {})
    pitch = float(cl_cfg.get("adjacent_strips_distance", 0.0165))
    first = float(cl_cfg.get("first_strip_pos", -0.528))
    params = ma.zero_params()
    if args.params_in:
        params = ma.params_from_json(json.loads(args.params_in.read_text())["parameters"])

    samples, infos = [], {L: [] for L in ma.LETTERS}
    for p in args.input:
        try:
            s = ma.read_sample(p, geom, max_events=args.max_events, forward=args.forward)
        except ma.MissingBranchesError as exc:
            print(f"  [skip] {exc.path.name}: missing branches")
            continue
        samples.append(s)
        for L in ma.LETTERS:
            infos[L].append(plane2_info(p, s, L, args.max_events, pitch, args.split_strips))
        print(f"  {Path(p).name}: {len(s):,} tracks")
    cat = lambda f: np.concatenate([f(s) for s in samples])  # noqa: E731
    sample = ma.TrackSample({L: cat(lambda s, L=L: s.u[L]) for L in ma.LETTERS},
                            {L: cat(lambda s, L=L: s.stored_res[L]) for L in ma.LETTERS},
                            {q: cat(lambda s, q=q: s.stored_chi2[q]) for q in ("xz", "xy")},
                            cat(lambda s: s.run), cat(lambda s: s.event))
    fits = ma.evaluate(geom, sample, params)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = lambda name: args.output_dir / f"{args.tag}{name}"  # noqa: E731

    res = {}
    for L in ma.LETTERS:
        inf = {k: (np.concatenate([i[k] for i in infos[L]]) if k not in ("strip_id", "strip_e")
                   else ak.concatenate([i[k] for i in infos[L]])) for k in infos[L][0]}
        b = fits[L]["b"]
        e = fits[L]["r"][:, 1] / (1.0 - fits[L]["h"][:, 1])                     # unbiased residual [m]
        u_pred = sample.u[L][:, 1] - e
        off = (geom.z_add if L == "Z" else geom.y_add)[1]
        s_pred = (u_pred - off - first) / pitch                                  # predicted strip coordinate
        eta = s_pred - np.round(s_pred)
        sig = 1.4826 * np.median(np.abs(e - np.median(e)))
        core = np.abs(e) < args.core_sigma * sig
        bcls = np.digitize(b, B_EDGES) - 1
        res[L] = dict(b=b, e=e * 1e3, eta=eta, s_pred=s_pred, core=core, bcls=bcls, sig=sig * 1e3, **inf)
        print(f"{L}2: robust width {sig * 1e3:.2f} mm, core |e| < {args.core_sigma * sig * 1e3:.1f} mm")

    rows = []
    # ---------------- fig 1: residual vs eta per slope class
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    eta_edges = np.linspace(-0.5, 0.5, 11)
    ex = 0.5 * (eta_edges[1:] + eta_edges[:-1])
    for ax, L in zip(axs, ma.LETTERS):
        r = res[L]
        for c, (lab, col) in enumerate(zip(B_LABELS, B_COLORS)):
            m = (r["bcls"] == c) & r["core"]
            if m.sum() < 500:
                continue
            idx = np.digitize(r["eta"][m], eta_edges) - 1
            prof = np.array([core_mean(r["e"][m], idx == i) for i in range(10)])
            ax.errorbar(ex, prof[:, 0], yerr=prof[:, 1], fmt="o-", ms=3.5, lw=1.2, color=col, capsize=0,
                        label=f"{lab} ({m.sum() / 1e3:.0f}k)")
            amp = np.nanmax(prof[:, 0]) - np.nanmin(prof[:, 0])
            rows.append({"figure": "fig1", "view": L, "b_class": lab, "n_core": int(m.sum()),
                         "quantity": "peak-to-peak of core residual vs eta [mm]", "value": f"{amp:.2f}"})
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"fig1 {'a' if L == 'Z' else 'b'}: {L}2 ({ma.PROJ_OF[L]} view)", loc="left", fontweight="bold")
        ax.set_xlabel("eta = predicted position inside the strip pitch (0 = strip centre)")
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend(fontsize=7.5, title="track slope b", title_fontsize=8)
    axs[0].set_ylabel("core mean of the unbiased residual [mm]")
    fig.suptitle("fig1: plane-2 residual vs position inside the strip, per incidence angle "
                 "(flat at 0 = linear position response)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out("fig1_residual_vs_eta.png"), dpi=120)
    plt.close(fig)

    # ---------------- fig 2: cluster size, split clusters, energy per path length vs b
    fig, axs = plt.subplots(3, 2, figsize=(13, 11), sharex="col")
    for c_i, L in enumerate(ma.LETTERS):
        r = res[L]
        lo, hi = np.percentile(r["b"], [1, 99])
        edges = np.unique(np.quantile(r["b"][(r["b"] >= lo) & (r["b"] <= hi)], np.linspace(0, 1, 21)))
        x = 0.5 * (edges[1:] + edges[:-1])
        idx = np.digitize(r["b"], edges) - 1
        frac = lambda sel: np.array([sel[idx == i].mean() * 100 for i in range(len(x))])  # noqa: E731
        sz = r["size"]
        for lab, sel in (("size 1", sz == 1), ("size 2", sz == 2), ("size 3", sz == 3), ("size 4+", sz >= 4)):
            axs[0, c_i].plot(x, frac(sel), "o-", ms=3.5, label=lab)
        axs[1, c_i].plot(x, frac(r["split"]), "o-", ms=3.5, color=SIZE_COL["split"],
                         label=f"2nd cluster within {args.split_strips} strips")
        epl = r["energy"] / np.sqrt(1 + r["b"] ** 2)
        axs[2, c_i].plot(x, [np.median(epl[idx == i]) for i in range(len(x))], "o-", ms=3.5, color="k",
                         label="median track-cluster energy / sqrt(1+b^2)")
        axs[2, c_i].plot(x, [np.median(r["energy"][idx == i]) for i in range(len(x))], "s--", ms=3.5, color="0.5",
                         mfc="white", label="median track-cluster energy")
        for row_i, (txt, yl) in enumerate((("a: track-cluster size", "fraction [%]"),
                                           ("b: split clusters", "fraction [%]"),
                                           ("c: cluster energy per path length", "energy [ADC-like units]"))):
            ax = axs[row_i, c_i]
            ax.set_title(f"fig2 {L}2 {txt}", loc="left", fontsize=10, fontweight="bold")
            ax.set_ylabel(yl)
            ax.grid(True, alpha=0.3, ls="--")
            ax.legend(fontsize=8)
        axs[2, c_i].set_xlabel(f"track slope b ({ma.PROJ_OF[L]} projection)")
        for lab, sel in (("size1", sz == 1), ("size3+", sz >= 3), ("split", r["split"])):
            for c, cl_lab in enumerate(B_LABELS):
                m = r["bcls"] == c
                if m.sum() >= 500:
                    rows.append({"figure": "fig2", "view": L, "b_class": cl_lab, "n_core": int(m.sum()),
                                 "quantity": f"{lab} fraction [%]", "value": f"{100 * sel[m].mean():.2f}"})
    fig.suptitle("fig2: plane-2 track cluster vs incidence angle", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out("fig2_size_split_energy_vs_b.png"), dpi=120)
    plt.close(fig)

    # ---------------- fig 3: strip energy profile around the predicted position
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
    d_edges = np.linspace(-4.5, 4.5, 37)
    dx = 0.5 * (d_edges[1:] + d_edges[:-1])
    for ax, L in zip(axs, ma.LETTERS):
        r = res[L]
        d_all = r["strip_id"] - ak.Array(r["s_pred"])[:, None]                             # strip - prediction [pitch]
        for c, (lab, col) in enumerate(zip(B_LABELS, B_COLORS)):
            m = r["bcls"] == c
            if m.sum() < 500:
                continue
            d = ak.to_numpy(ak.flatten(d_all[m]))
            w = ak.to_numpy(ak.flatten(r["strip_e"][m]))
            n, _ = np.histogram(d, d_edges, weights=w)
            near = np.abs(d) < 2.5
            centroid = (d[near] * w[near]).sum() / w[near].sum()
            rms = np.sqrt((((d[near] - centroid) ** 2) * w[near]).sum() / w[near].sum())
            ax.step(dx, n / n.sum(), where="mid", color=col, lw=1.3,
                    label=f"{lab}: centroid {centroid * pitch * 1e3:+.2f} mm, rms {rms:.2f} pitch")
            rows.append({"figure": "fig3", "view": L, "b_class": lab, "n_core": int(m.sum()),
                         "quantity": "strip-energy centroid within 2.5 strips of u_pred [mm]",
                         "value": f"{centroid * pitch * 1e3:.2f}"})
            rows.append({"figure": "fig3", "view": L, "b_class": lab, "n_core": int(m.sum()),
                         "quantity": "strip-energy rms within 2.5 strips [pitch]", "value": f"{rms:.3f}"})
        ax.set_yscale("log")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_title(f"fig3 {'a' if L == 'Z' else 'b'}: {L}2 strip energy around the predicted position",
                     loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel("strip - predicted position [strip pitch]")
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend(fontsize=7)
    axs[0].set_ylabel("fraction of the strip energy")
    fig.suptitle("fig3: response function (all plane-2 cluster strips, energy-weighted), per incidence angle",
                 fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out("fig3_strip_profile.png"), dpi=120)
    plt.close(fig)

    # ---------------- fig 4: residual vs b by cluster size
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, L in zip(axs, ma.LETTERS):
        r = res[L]
        lo, hi = np.percentile(r["b"], [1, 99])
        edges = np.unique(np.quantile(r["b"][(r["b"] >= lo) & (r["b"] <= hi)], np.linspace(0, 1, 21)))
        x = 0.5 * (edges[1:] + edges[:-1])
        idx = np.digitize(r["b"], edges) - 1
        sz = r["size"]
        for lab, sel in (("1", (sz == 1) & ~r["split"]), ("2", (sz == 2) & ~r["split"]),
                         ("3+", (sz >= 3) & ~r["split"]), ("split", r["split"])):
            prof = np.array([core_mean(r["e"], (idx == i) & sel & r["core"]) for i in range(len(x))])
            ax.errorbar(x, prof[:, 0], yerr=prof[:, 1], fmt="o-", ms=3.5, lw=1, capsize=0, color=SIZE_COL[lab],
                        label=("size " + lab + ", no 2nd cluster nearby") if lab != "split" else
                        f"2nd cluster within {args.split_strips} strips")
        allp = np.array([core_mean(r["e"], (idx == i) & r["core"]) for i in range(len(x))])
        ax.errorbar(x, allp[:, 0], yerr=allp[:, 1], fmt="s--", ms=3.5, lw=1, color="k", mfc="white", label="all")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"fig4 {'a' if L == 'Z' else 'b'}: {L}2 core residual vs b by cluster type",
                     loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel(f"track slope b ({ma.PROJ_OF[L]} projection)")
        ax.set_ylabel("core mean of the unbiased residual [mm]")
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend(fontsize=7.5)
    fig.suptitle("fig4: which cluster types carry the slope dependence of the plane-2 residual", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out("fig4_residual_by_size_vs_b.png"), dpi=120)
    plt.close(fig)

    # ---------------- fig 5: residual vs b for even / odd predicted strip (alternating triangular bars)
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, L in zip(axs, ma.LETTERS):
        r = res[L]
        parity = np.round(r["s_pred"]).astype(int) % 2
        lo, hi = np.percentile(r["b"], [1, 99])
        edges = np.unique(np.quantile(r["b"][(r["b"] >= lo) & (r["b"] <= hi)], np.linspace(0, 1, 21)))
        x = 0.5 * (edges[1:] + edges[:-1])
        idx = np.digitize(r["b"], edges) - 1
        profs = {}
        for par, lab, col in ((0, "even strips", "#3B6FB6"), (1, "odd strips", "#D1495B")):
            profs[par] = np.array([core_mean(r["e"], (idx == i) & (parity == par) & r["core"]) for i in range(len(x))])
            ax.errorbar(x, profs[par][:, 0], yerr=profs[par][:, 1], fmt="o-", ms=3.5, lw=1, capsize=0, color=col,
                        label=f"predicted hit on {lab}")
            (sl, ic), cov = np.polyfit(x, profs[par][:, 0], 1, w=1 / profs[par][:, 1], cov="unscaled")
            rows.append({"figure": "fig5", "view": L, "b_class": lab, "n_core": int(((parity == par) & r["core"]).sum()),
                         "quantity": "slope of core residual vs b [mm per unit b]", "value": f"{sl:.2f}"})
        avg = 0.5 * (profs[0][:, 0] + profs[1][:, 0])
        ax.plot(x, avg, "s--", ms=3.5, color="k", mfc="white", lw=1, label="average of even and odd")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"fig5 {'a' if L == 'Z' else 'b'}: {L}2 core residual vs b, even vs odd strip",
                     loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel(f"track slope b ({ma.PROJ_OF[L]} projection)")
        ax.set_ylabel("core mean of the unbiased residual [mm]")
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend(fontsize=8)
    fig.suptitle("fig5: strip parity (alternating triangular bars): opposite residual slopes = the two bar "
                 "orientations sit at different depths", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out("fig5_residual_vs_b_by_parity.png"), dpi=120)
    plt.close(fig)

    with open(out("response_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\n{'figure':7s}{'view':5s}{'b class':12s}{'N':>9s}  quantity = value")
    for r in rows:
        print(f"{r['figure']:7s}{r['view']:5s}{r['b_class']:12s}{r['n_core']:9d}  {r['quantity']} = {r['value']}")
    print(f"\nFigures and response_summary.csv in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
