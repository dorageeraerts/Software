#!/usr/bin/env python3
"""Plots for the MURAVES track-based alignment (style follows chi2_strip_analysis_residuals.py)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

from muraves_alignment import (LETTERS, PARAMS, PARAM_LABEL, PARAM_SCALE, PROJ_OF, gauss_core_fit)

C_BEFORE, C_AFTER = "#D1495B", "#2A9D5C"   # muted crimson / emerald green


def _stats(r_mm):
    return r_mm.mean(), r_mm.std() / np.sqrt(len(r_mm)), r_mm.std()


# ----------------------------------------------------------------------------
# 1) main plot: residuals before / after + ratio panel
# ----------------------------------------------------------------------------
def plot_residuals_before_after(fits0, fits1, output_file, title, nbins=80, min_count=10):
    """3p residuals (Residue_Track3p_*, meas - fit) before and after alignment.
    Rows: Z1-Z3 (+ratio), Y1-Y3 (+ratio).  Ratio = after / before per bin (Poisson errors; the two
    histograms contain the same tracks, so the errors are conservative)."""
    fig, axs = plt.subplots(4, 3, figsize=(13, 11.5), squeeze=False, layout="constrained",
                            gridspec_kw={"height_ratios": [3, 1, 3, 1]})
    for row, L in enumerate(LETTERS):
        for k in range(3):
            ax, axr = axs[2 * row, k], axs[2 * row + 1, k]
            r0 = fits0[L]["r"][:, k] * 1e3
            r1 = fits1[L]["r"][:, k] * 1e3
            lim = 1.05 * np.percentile(np.abs(np.concatenate([r0, r1])), 99.9)
            edges = np.linspace(-lim, lim, nbins + 1)
            n0, _ = np.histogram(r0, edges)
            n1, _ = np.histogram(r1, edges)
            ctr = 0.5 * (edges[1:] + edges[:-1])
            for n, col, lab, r in ((n0, C_BEFORE, "before", r0), (n1, C_AFTER, "after", r1)):
                m, sem, rms = _stats(r)
                ax.stairs(n, edges, color=col, lw=1.5,
                          label=f"{lab}: mean={m:+.3f}$\\pm${sem:.3f} mm, rms={rms:.3f} mm\n"
                                f"median={np.median(r):+.3f} mm")
            ax.set_title(f"{L}{k + 1}: Residue_Track3p_{L.lower()}{k + 1}", fontsize=9, pad=7)
            ax.set_xlim(-lim, lim)
            ax.set_ylim(0, 1.45 * max(n0.max(), n1.max()))
            ax.grid(True, alpha=0.3, ls="--")
            ax.legend(fontsize=7, loc="upper right")
            if k == 0:
                ax.set_ylabel("tracks / bin")
            ax.tick_params(labelbottom=False)

            ok = (n0 >= min_count) & (n1 >= min_count)
            axr.axhline(1.0, color="black", lw=1.0)
            if ok.any():
                ratio = n1[ok] / n0[ok]
                err = ratio * np.sqrt(1.0 / n0[ok] + 1.0 / n1[ok])
                axr.errorbar(ctr[ok], ratio, yerr=err, fmt=".", ms=3, color="0.25", lw=0.7, capsize=0)
                m_ = float(np.clip(np.max(np.abs(np.log(ratio))), np.log(1.5), np.log(4)))
                axr.set_ylim(np.exp(-m_ * 1.15), np.exp(m_ * 1.15))
            axr.set_xlim(-lim, lim)
            axr.set_yscale("log")
            axr.set_yticks([0.5, 1, 2])
            axr.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            axr.minorticks_off()
            axr.grid(True, alpha=0.3, ls="--")
            if k == 0:
                axr.set_ylabel("after\n/ before", fontsize=8)
            axr.set_xlabel("residual (measured - fit) [mm]")
    fig.suptitle(title, fontsize=10)
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# 2) convergence of the alignment parameters (paper Fig. 4)
# ----------------------------------------------------------------------------
def plot_parameter_convergence(history, free_planes, free_params, errors, output_file, title):
    """parameter value vs iteration, one column per aligned plane (like Fig. 4 of the paper)"""
    its = [h["iteration"] for h in history]
    pars = [p for p in PARAMS if p in free_params]
    fig, axs = plt.subplots(len(pars), len(free_planes), figsize=(4.6 * len(free_planes), 2.1 * len(pars) + 0.8),
                            sharex=True, squeeze=False)
    for c, plane in enumerate(free_planes):
        for r, name in enumerate(pars):
            i = PARAMS.index(name)
            ax = axs[r, c]
            vals = np.array([h["params"][plane][i] for h in history]) * PARAM_SCALE[name]
            ax.plot(its, vals, "o-", ms=4, color="C0")
            fin = vals[-1]
            ax.axhline(fin, color="0.4", ls="--", lw=0.8)
            if errors is not None and plane in errors and np.isfinite(errors[plane][i]):
                e = errors[plane][i] * PARAM_SCALE[name]
                ax.axhspan(fin - e, fin + e, color="0.5", alpha=0.25, lw=0)
                ax.set_title(f"final {fin:+.3f} $\\pm$ {e:.3f}", loc="right", fontsize=8)
            ax.grid(True, alpha=0.3, ls="--")
            ax.set_ylabel(PARAM_LABEL[name], fontsize=9)
            if r == 0:
                ax.set_title(f"plane {plane}", loc="left", fontsize=10, fontweight="bold")
            if r == len(pars) - 1:
                ax.set_xlabel("iteration")
    fig.suptitle(title + "\n(grey band: final value $\\pm 1\\sigma_{stat}$)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# 3) quality of the fit vs iteration
# ----------------------------------------------------------------------------
def plot_convergence_quality(history, output_file, title):
    its = [h["iteration"] for h in history]
    panels = [("res_mean", 1e6, "mean residual [$\\mu$m]"), ("res_rms", 1e3, "rms residual [mm]")]
    if "res_gauss_mu" in history[0]:
        panels.append(("res_gauss_mu", 1e6, "Gaussian core mean residual [$\\mu$m]"))
    fig, axs = plt.subplots(1, 1 + len(panels), figsize=(5 * (1 + len(panels)), 4.2))
    for proj, col in (("xz", "C0"), ("xy", "C3")):
        axs[0].plot(its, [h["chi2_mean"][proj] for h in history], "o-", ms=4, color=col,
                    label=f"3p {proj}")
    axs[0].set_ylabel(r"$\langle\chi^2\rangle$ per track (1 dof)")
    axs[0].legend()
    for i, (key, scale, lab) in enumerate(panels):
        ax = axs[i + 1]
        for j, plane in enumerate(("Z1", "Z2", "Z3", "Y1", "Y2", "Y3")):
            v = np.array([h[key][plane] for h in history]) * scale      # m -> um or mm
            ax.plot(its, v, marker="o", ms=3, color=f"C{j}", ls="-" if plane[0] == "Z" else "--", label=plane)
        ax.set_ylabel(lab)
        ax.legend(ncol=2, fontsize=8)
    for ax in axs:
        ax.set_xlabel("iteration")
        ax.grid(True, alpha=0.3, ls="--")
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# 4) residual profiles vs the quantities that enter the Jacobian
# ----------------------------------------------------------------------------
def _profile(x, y, edges):
    idx = np.clip(np.digitize(x, edges) - 1, 0, len(edges) - 2)
    n = np.bincount(idx, minlength=len(edges) - 1).astype(float)
    s1 = np.bincount(idx, weights=y, minlength=len(edges) - 1)
    s2 = np.bincount(idx, weights=y * y, minlength=len(edges) - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = s1 / n
        sem = np.sqrt(np.maximum(s2 / n - m * m, 0) / n)
    return m, sem, n


def plot_profiles(geom, fits0, fits1, free_planes, output_dir, title, nbins=14, prefix=""):
    """mean residual vs track slope b (tan psi), vs u and vs v of the track at the plane.
    These are exactly the 'handles' of the method: Delta w couples to b, Delta beta to b*u,
    Delta alpha to b*v, Delta gamma to v, Delta u to a constant."""
    names = [("b", "track slope  $b=\\tan\\psi$"), ("u", "u at plane (from plane centre) [m]"),
             ("v", "v at plane (along the strips) [m]")]
    for L in LETTERS:
        fig, axs = plt.subplots(3, 3, figsize=(13, 9.5), sharey="row", squeeze=False)
        for k in range(3):
            f0, f1 = fits0[L], fits1[L]
            feat = {"b": f0["b"],
                    "u": f0["fit"][:, k] - geom.u_centre[L][k],
                    "v": f0["v_loc"][:, k]}
            for c, (key, lab) in enumerate(names):
                ax = axs[k, c]
                edges = np.unique(np.quantile(feat[key], np.linspace(0, 1, nbins + 1)))
                ctr = 0.5 * (edges[1:] + edges[:-1])
                for f, col, lbl, dx in ((f0, C_BEFORE, "before", -0.01), (f1, C_AFTER, "after", 0.01)):
                    m, sem, n = _profile(feat[key], f["r"][:, k] * 1e3, edges)
                    ax.errorbar(ctr + dx * (edges[-1] - edges[0]), m, yerr=sem, fmt="o-", ms=3.5, lw=1,
                                color=col, label=lbl, capsize=0)
                ax.axhline(0, color="black", lw=0.8)
                ax.grid(True, alpha=0.3, ls="--")
                if k == 2:
                    ax.set_xlabel(lab)
                if c == 0:
                    tag = "aligned" if f"{L}{k + 1}" in free_planes else "reference (fixed)"
                    ax.set_ylabel(f"{L}{k + 1} ({tag})\nmean residual [mm]", fontsize=9)
                if k == 0 and c == 0:
                    ax.legend(fontsize=8)
        fig.suptitle(f"{title}\n{L} planes: mean residual vs the track quantities entering the Jacobian "
                     f"(equal-statistics bins)")
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(output_dir / f"{prefix}residual_profiles_{L}.png", dpi=130, bbox_inches="tight")
        plt.close(fig)


def _core_profile(x, y, edges, window):
    """Gaussian core mean (and its error) of y in bins of x; NaN where the fit fails / too few entries"""
    m, err = np.full(len(edges) - 1, np.nan), np.full(len(edges) - 1, np.nan)
    idx = np.digitize(x, edges) - 1
    for i in range(len(edges) - 1):
        yi = y[idx == i]
        if len(yi) < 200:
            continue
        try:
            g = gauss_core_fit(yi, window)
            m[i], err[i] = g["mu"], g["mu_err"]
        except (RuntimeError, ValueError):
            pass
    return m, err


def plot_residual_vs_slope(fits0, fits1, output_file, title, planes=("Z2", "Y2"), nbins=20, b_range=(1, 99),
                           gauss_window="quantile"):
    """Mean UNBIASED residual e = r / (1 - h_kk) (hit minus the line through the other two planes) vs the track
    slope b, before / after alignment, with a straight-line fit to each profile.  Applying a shift dw to the
    plane gives e = e0 - b * dw exactly, so the slope of the profile is the dw that is still needed (before: the
    full correction, after: what the alignment left, e.g. by a rigid-station constraint); a profile that is not
    a straight line is something other than a plane shift (e.g. a slope-dependent cluster-position bias).
    Top row: plain mean per bin (sensitive to the tails).  Bottom row: mean of a Gaussian fitted to the core of
    the residuals in each bin (gauss_core_fit, window `gauss_window`), i.e. the profile without the tails.
    fits0/fits1 must hold BIASED residuals (h is used to unbias them).  Equal-statistics bins in b between the
    b_range percentiles."""
    fig, axs = plt.subplots(2, len(planes), figsize=(6.2 * len(planes), 9.2), squeeze=False, sharey="row")
    for col_i, plane in enumerate(planes):
        L, k = plane[0], int(plane[1]) - 1
        b = fits0[L]["b"]
        lo, hi = np.percentile(b, b_range)
        edges = np.unique(np.quantile(b[(b >= lo) & (b <= hi)], np.linspace(0, 1, nbins + 1)))
        ctr = 0.5 * (edges[1:] + edges[:-1])
        for row, kind in enumerate(("mean", "core")):
            ax = axs[row, col_i]
            for f, col, lbl, dx in ((fits0, C_BEFORE, "before", -0.004), (fits1, C_AFTER, "after", 0.004)):
                bb = f[L]["b"]
                e = f[L]["r"][:, k] / (1.0 - f[L]["h"][:, k]) * 1e3
                sel = (bb >= edges[0]) & (bb <= edges[-1])
                if kind == "mean":
                    m, sem, _ = _profile(bb[sel], e[sel], edges)
                else:
                    m, sem = _core_profile(bb[sel], e[sel], edges, gauss_window)
                ok = np.isfinite(m) & np.isfinite(sem) & (sem > 0)
                lab = lbl
                if ok.sum() > 2:
                    (slope, icpt), cov = np.polyfit(ctr[ok], m[ok], 1, w=1.0 / sem[ok], cov="unscaled")
                    chi2 = float((((m[ok] - (icpt + slope * ctr[ok])) / sem[ok]) ** 2).sum())
                    lab = (f"{lbl}: slope (dw still needed) = {slope:+.2f} $\\pm$ {np.sqrt(cov[0, 0]):.2f} mm, "
                           f"offset {icpt:+.2f} mm\n$\\chi^2$/ndf (line) = {chi2:.1f}/{int(ok.sum()) - 2}")
                    xx = np.array([edges[0], edges[-1]])
                    ax.plot(xx, icpt + slope * xx, color=col, lw=1, ls="--")
                ax.errorbar(ctr + dx * (edges[-1] - edges[0]), m, yerr=sem, fmt="o" if kind == "mean" else "s",
                            ms=3.5, lw=1, color=col, capsize=0, mfc=col if kind == "mean" else "white", label=lab)
            ax.axhline(0, color="black", lw=0.8)
            ax.set_title(f"{plane}: " + ("mean residual" if kind == "mean" else
                                         f"Gaussian core mean ({gauss_window} window)"),
                         loc="left", fontsize=10, fontweight="bold")
            ax.set_xlabel(f"track slope $b$ ({PROJ_OF[L]} projection)")
            ax.grid(True, alpha=0.3, ls="--")
            ax.legend(fontsize=7.5, loc="upper left")
    axs[0, 0].set_ylabel("mean unbiased residual $r/(1-h)$ [mm]")
    axs[1, 0].set_ylabel("Gaussian core mean of $r/(1-h)$ [mm]")
    fig.suptitle(f"{title}\nresidual vs track slope: a missing plane shift dw shows up as a straight line of slope dw"
                 f"\ntop: mean (incl. tails), bottom: Gaussian core mean", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# 5) track chi2
# ----------------------------------------------------------------------------
def plot_track_chi2(fits0, fits1, output_file, title, source="recomputed"):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, L in zip(axs, LETTERS):
        c0, c1 = fits0[L]["chi2"], fits1[L]["chi2"]
        hi = max(np.percentile(c0, 99.9), np.percentile(c1, 99.9))
        edges = np.linspace(0, hi, 61)
        ax.hist(c0, edges, histtype="step", color=C_BEFORE, lw=1.5, label=f"before ($\\langle\\chi^2\\rangle$={c0.mean():.4f})")
        ax.hist(c1, edges, histtype="step", color=C_AFTER, lw=1.5, label=f"after ($\\langle\\chi^2\\rangle$={c1.mean():.4f})")
        ax.set_yscale("log")
        ax.set_xlabel(f"BestTrack_3p_ChiSquare_{PROJ_OF[L]} ({source})")
        ax.set_ylabel("tracks")
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)

# ----------------------------------------------------------------------------
# 6) Gaussian core fits of the residuals (estimator "gauss")
# ----------------------------------------------------------------------------
def plot_gauss_fits(fits, output_file, title, window="quantile"):
    """residual distribution per plane with the Gaussian fitted to its core (gauss_core_fit): fit window
    shaded, fitted curve solid inside / dashed outside the window, fit parameters and chi2/ndf in the legend"""
    fig, axs = plt.subplots(2, 3, figsize=(14, 8.5), squeeze=False)
    for row, L in enumerate(LETTERS):
        for k in range(3):
            ax = axs[row, k]
            r = fits[L]["r"][:, k] * 1e3
            gf = gauss_core_fit(r, window)
            mu, sig = gf["mu"], gf["sigma"]
            bw = gf["edges"][1] - gf["edges"][0]
            lim = (mu - 4 * sig, mu + 4 * sig)
            e0 = gf["edges"][0]   # display bins aligned with the fit bins
            edges = e0 + bw * np.arange(np.floor((lim[0] - e0) / bw), np.ceil((lim[1] - e0) / bw) + 1)
            n, _ = np.histogram(r, edges)
            ax.stairs(n, edges, color="0.35", lw=1.2, label=f"residuals (N={len(r):,})")
            ax.axvspan(gf["lo"], gf["hi"], color="C0", alpha=0.12, lw=0,
                       label=f"fit window ({gf['n_in'] / len(r) * 100:.0f}% of tracks)")
            ctr = 0.5 * (gf["edges"][1:] + gf["edges"][:-1])
            ok = gf["counts"] > 0
            ax.errorbar(ctr[ok], gf["counts"][ok], yerr=np.sqrt(gf["counts"][ok]), fmt=".", ms=3, color="black",
                        lw=0.7, capsize=0)
            x_in = np.linspace(gf["lo"], gf["hi"], 200)
            x_all = np.linspace(*lim, 400)
            g = lambda x: gf["amp"] * np.exp(-0.5 * ((x - mu) / sig) ** 2)  # noqa: E731
            ax.plot(x_all, g(x_all), color="C3", lw=1, ls="--")
            ax.plot(x_in, g(x_in), color="C3", lw=2,
                    label=f"Gauss fit: $\\mu$={mu:+.3f}$\\pm${gf['mu_err']:.3f} mm\n"
                          f"$\\sigma$={sig:.3f}$\\pm${gf['sigma_err']:.3f} mm\n"
                          f"A={gf['amp']:.0f}$\\pm${gf['amp_err']:.0f} / bin\n"
                          f"$\\chi^2$/ndf={gf['chi2']:.1f}/{gf['ndf']} = {gf['chi2'] / max(gf['ndf'], 1):.2f}")
            ax.axvline(mu, color="C3", lw=0.8)
            ax.axvline(r.mean(), color="C2", lw=0.8, ls=":", label=f"sample mean={r.mean():+.3f} mm")
            ax.axvline(np.median(r), color="C1", lw=0.8, ls="-.", label=f"median={np.median(r):+.3f} mm")
            ax.set_xlim(*lim)
            ax.set_ylim(0, 1.6 * n.max())
            ax.set_title(f"{L}{k + 1}: Residue_Track3p_{L.lower()}{k + 1}", fontsize=9)
            ax.set_xlabel("residual (measured - fit) [mm]")
            if k == 0:
                ax.set_ylabel(f"tracks / {bw:.3f} mm")
            ax.grid(True, alpha=0.3, ls="--")
            ax.legend(fontsize=6.5, loc="upper left")
    wtxt = "central 68% of the residuals" if window == "quantile" else "iterative $\\mu\\pm\\sigma$ of the fit"
    fig.suptitle(f"{title}\nGaussian core fit, window: {wtxt}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output_file, dpi=130, bbox_inches="tight")
    plt.close(fig)
