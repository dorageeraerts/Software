#!/usr/bin/env python3
"""
Track-based sensor alignment for MURAVES (planes 1-3 only, "3p" tracks).

Implements the method of V. Karimaki et al., "Sensor Alignment by Tracks",
CHEP03 (arXiv:physics/0306034) for the single-sided strip planes of MURAVES.

------------------------------------------------------------------------------
Conventions (all verified against the ROOT file, see check_against_stored())
------------------------------------------------------------------------------
* x is the telescope axis (plane-to-plane direction, x_pos = [-0.26, 0, 0.262] m).
  Z planes measure z, Y planes measure y.  Each projection has its own straight
  line  u(x) = a + b*x  fitted through its 3 hits (equal weights).
* hit position used by the reconstruction = ClusterPosition + offset
  (Z planes: z_add, Y planes: y_add).  Offsets are read from the METADATA of the
  ROOT file (reco_parameters.json), so the alignment parameters found here are
  *corrections on top of the reconstruction geometry*.
* Residual convention = Residue_Track3p_* of the ROOT file:   r = measured - fit.
* chi2 = sum_planes r^2 / sigma^2  (sigma_z = 4.0 mm, sigma_y = 3.5 mm, from METADATA)
  -> reproduces BestTrack_3p_ChiSquare_{xz,xy} exactly (1 degree of freedom).

Paper -> MURAVES mapping (local frame of a plane = (u, v, w)):
    u = precise (strip) coordinate   -> z for Z planes, y for Y planes
    v = coarse coordinate (along the strips, NOT measured)  -> taken from the
        *other* projection of the same 3D track (this is what the 3D track gives us)
    w = plane normal                 -> x
Per plane 5 parameters (Delta v is redundant for single-sided strips, paper IV.B):
    du, dw  (shifts, metres),  dalpha, dbeta, dgamma (tilts about u, v, w; radians)
Linearised model (paper Eqs. 19/20, "tilt formalism"):
    corrected hit  :  x'  = x_plane + dw + dbeta*u + dalpha*v
                      u'  = u_meas + du - dgamma*v
    residual model :  r(p) = r0 + du - b*dw - b*v*dalpha - b*u*dbeta - v*dgamma
    Jacobian (Eq. 17, in the r = meas - fit convention): J = [1, -b, -b*v, -b*u, -v]
    with b = tan(psi) = track slope in that projection, u, v measured from the
    centre of the plane.

NOTE on a sign in the paper: Eq. (15) is printed as  dp = +(J^T V^-1 J)^-1 J^T V^-1 eps.
For eps = fit - meas and J = d eps/dp the Gauss-Newton step is the NEGATIVE of that,
which is what is implemented here (checked with the closure test).

------------------------------------------------------------------------------
What can and cannot be determined with only 3 planes
------------------------------------------------------------------------------
With three planes per projection a straight-line fit absorbs a common shift and a
common tilt of the three planes, so only ONE combination of the plane positions
("bow", roughly u1 - 2*u2 + u3) is observable per feature.  Reference planes must
therefore be fixed (exactly like the reference sensors of the paper).  The default
is: Z1, Z3, Y1, Y3 fixed, Z2, Y2 aligned.  Other choices are possible through
`free_planes`, but the degenerate directions are then NOT constrained by the data.
The 4th plane (far lever arm) is what would lift this degeneracy later.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

try:
    import awkward as ak
    import uproot
except ImportError:  # only needed for reading ROOT files, not for the toy test
    ak = None
    uproot = None

# ----------------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------------
TREE_NAME = "AnalyzedData"
PROJ_OF = {"Z": "xz", "Y": "xy"}           # plane letter -> track projection
OTHER = {"Z": "Y", "Y": "Z"}
LETTERS = ("Z", "Y")
PLANES = [f"{l}{k}" for l in LETTERS for k in (1, 2, 3)]
PARAMS = ["du", "dw", "dalpha", "dbeta", "dgamma"]
PARAM_LABEL = {"du": r"$\Delta u$ [mm]", "dw": r"$\Delta w$ [mm]",
               "dalpha": r"$\Delta\alpha$ [mrad]", "dbeta": r"$\Delta\beta$ [mrad]",
               "dgamma": r"$\Delta\gamma$ [mrad]"}
PARAM_SCALE = {"du": 1e3, "dw": 1e3, "dalpha": 1e3, "dbeta": 1e3, "dgamma": 1e3}  # m->mm, rad->mrad

# fallback geometry if METADATA is missing (copy of reco_parameters.json)
FALLBACK_GEOMETRY = {
    "ROSSO": {"z_add": [0.292, 0.251, 0.207, 0.0], "x_pos": [-0.265, 0.0, 0.262, 1.475]},
    "NERO": {"z_add": [0.293, 0.251, 0.21, 0.0], "x_pos": [-0.26, 0.0, 0.262, 1.492]},
    "BLU": {"z_add": [0.2712, 0.2312, 0.1892, 0.0], "x_pos": [-0.26, 0.0, 0.262, 1.492]},
}
# x(Y layer) - x(Z layer) within a station [m] (approx. 5.5 cm, same for all stations)
YZ_SPACING = 0.055


# ----------------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------------
@dataclass
class Geometry:
    detector: str
    x_pos: np.ndarray                 # (3,) x of planes 1..3 [m] (Z planes)
    z_add: np.ndarray                 # (3,) offset added to Z cluster positions [m]
    y_add: np.ndarray                 # (3,) offset added to Y cluster positions [m]
    sigma: dict                       # {"Z": sigma_z, "Y": sigma_y} [m]
    plane_centre_raw: float           # centre of the strip array in raw coordinates [m]
    x_pos_y: np.ndarray = None        # (3,) x of the Y planes [m] (reco "x_pos_y"; default: x_pos)
    yz_spacing: float = 0.0           # x(Y layer) - x(Z layer) of a station [m], added on top of x_pos_y in the model only
    u_centre: dict = field(default_factory=dict)   # {"Z"|"Y": (3,)} centre of plane in its measured coordinate
    v_centre: dict = field(default_factory=dict)   # {"Z"|"Y": (3,)} centre in the coordinate along the strips
    x_of: dict = field(default_factory=dict)       # {"Z"|"Y": (3,)} x of the planes of that projection (model)
    x_reco: dict = field(default_factory=dict)     # {"Z"|"Y": (3,)} x of the planes as in the reco (x_pos, x_pos_y)

    def __post_init__(self):
        if self.x_pos_y is None:
            self.x_pos_y = self.x_pos.copy()
        # A spacing common to all stations is a constant x shift of the whole Y projection: the residuals do
        # not change (so the check against the ROOT file still holds), only v (the coordinate along the strips,
        # taken from the other projection at the plane's x) does.  Outputs for the reco use x_reco.
        self.x_reco = {"Z": self.x_pos, "Y": self.x_pos_y}
        self.x_of = {"Z": self.x_pos, "Y": self.x_pos_y + self.yz_spacing}
        c0 = self.plane_centre_raw
        # Z plane k: measures z (offset z_add), strips run along y
        # Y plane k: measures y (offset y_add), strips run along z (same station => z_add of that station)
        self.u_centre = {"Z": c0 + self.z_add, "Y": c0 + self.y_add}
        self.v_centre = {"Z": c0 + self.y_add, "Y": c0 + self.z_add}

    def describe(self):
        return (f"{self.detector}: x_pos={self.x_pos.tolist()} x_pos_y={self.x_pos_y.tolist()} "
                f"yz_spacing={self.yz_spacing * 1e3:.1f} mm z_add={self.z_add.tolist()} "
                f"y_add={self.y_add.tolist()} sigma_z={self.sigma['Z']*1e3:.2f} mm "
                f"sigma_y={self.sigma['Y']*1e3:.2f} mm")


def read_metadata(path):
    """METADATA TObjString of the ROOT file -> dict (or None)."""
    if uproot is None:
        return None
    try:
        with uproot.open(path) as f:
            if "METADATA" not in f:
                return None
            return json.loads(str(f["METADATA"]))
    except Exception as exc:  # noqa: BLE001
        print(f"  [warn] could not parse METADATA of {path}: {exc}")
        return None


def load_geometry(path, detector=None, yz_spacing=YZ_SPACING) -> Geometry:
    """Geometry from the METADATA of a ROOT file (fallback: built-in table).
    yz_spacing: x(Y layer) - x(Z layer) within a station [m], added to the reco x_pos_y in the alignment model."""
    meta = read_metadata(path)
    reco, det_auto = None, None
    if meta is not None:
        for cfg in meta.get("detector_config_files", []):
            name = cfg["path"].split("/")[-1]
            # telescopeROSSO.cfg or an aligned copy such as telescopeROSSO_aligned.cfg
            m = re.match(r"telescope([A-Za-z]+)\w*\.cfg", name)
            if m:
                det_auto = m.group(1).upper()
            # reco_parameters.json or an aligned copy such as reco_parameters_aligned_runs4000-4100.json
            if name.startswith("reco_parameters") and name.endswith(".json"):
                reco = cfg["parsed_content"]
    det = (detector or det_auto)
    if det is None:
        raise SystemExit("Cannot determine the detector from METADATA; use --detector NERO|ROSSO|BLU "
                         "(run numbers alone are ambiguous between detectors).")
    det = det.upper()
    if detector and det_auto and det != det_auto:
        print(f"  [warn] --detector {det} differs from the detector in METADATA ({det_auto})")
    if reco is not None and det in reco["detector_geometry"]:
        g = reco["detector_geometry"][det]
        rc = reco["reconstruction"]
        sig = {"Z": float(rc["sigma_z"]), "Y": float(rc["sigma_y"])}
        y_add = np.array(rc.get("y_add", [0, 0, 0])[:3], float)
        cl = reco.get("cluster_lists", {})
        c0 = float(cl.get("first_strip_pos", -0.528)) + 0.5 * (int(cl.get("n_strips", 64)) - 1) * float(
            cl.get("adjacent_strips_distance", 0.0165))
    else:
        print("  [warn] no reco_parameters in METADATA, using built-in geometry table")
        g = FALLBACK_GEOMETRY[det]
        sig = {"Z": 0.004, "Y": 0.0035}
        y_add = np.zeros(3)
        c0 = -0.528 + 31.5 * 0.0165
    return Geometry(det, np.array(g["x_pos"][:3], float), np.array(g["z_add"][:3], float),
                    y_add, sig, c0, np.array(g.get("x_pos_y", g["x_pos"])[:3], float), float(yz_spacing))


# ----------------------------------------------------------------------------
# reading the tracks
# ----------------------------------------------------------------------------
@dataclass
class TrackSample:
    u: dict            # {"Z": (N,3), "Y": (N,3)} hit positions as used by the reco (cluster pos + offset) [m]
    stored_res: dict   # {"Z": (N,3), "Y": (N,3)} Residue_Track3p_* from the file (meas - fit) [m]
    stored_chi2: dict  # {"xz": (N,), "xy": (N,)} BestTrack_3p_ChiSquare_*
    run: np.ndarray    # (N,)
    event: np.ndarray  # (N,) entry number in the file

    def __len__(self):
        return len(self.run)

    def subset(self, mask):
        return TrackSample({k: v[mask] for k, v in self.u.items()},
                           {k: v[mask] for k, v in self.stored_res.items()},
                           {k: v[mask] for k, v in self.stored_chi2.items()},
                           self.run[mask], self.event[mask])


class MissingBranchesError(RuntimeError):
    """A ROOT file does not contain the branches required for alignment."""

    def __init__(self, path, missing):
        self.path = Path(path)
        self.missing = tuple(sorted(missing))
        super().__init__(f"{self.path}: missing branches {list(self.missing)}")


def _pick(jag, idx):
    """jag (n, var), idx (n,) -> float (n,), NaN where not available."""
    mask = ak.local_index(jag, axis=1) == idx[:, None]
    return ak.to_numpy(ak.fill_none(ak.firsts(jag[mask], axis=1), np.nan)).astype(float)


def _needed_branches(corr=False):
    """corr=True: the shoulder-corrected *_corr branches (per-track hit positions TrackClusterPosition_*_corr)"""
    c = "_corr" if corr else ""
    b = ["Run", f"BestTrack_3p_ChiSquare_xz{c}", f"BestTrack_3p_ChiSquare_xy{c}",
         f"BestTrack_3p_xz_index{c}", f"BestTrack_3p_xy_index{c}"]
    for k in (1, 2, 3):
        if corr:
            b += [f"TrackClusterPosition_z{k}_corr", f"TrackClusterPosition_y{k}_corr"]
        else:
            b += [f"TrackCluster_z{k}_index", f"TrackCluster_y{k}_index",
                  f"ClusterPosition_Z{k}", f"ClusterPosition_Y{k}"]
        b += [f"Residue_Track3p_z{k}{c}", f"Residue_Track3p_y{k}{c}"]
    return b


def read_sample(path, geom: Geometry, max_events=None, first_event=0, corr=False, forward=False) -> TrackSample:
    """3D 3p tracks (chi2_xz > 0 AND chi2_xy > 0, as in chi2_strip_analysis_residuals.py)
    of one ROOT file.  Only the best 3p track of each event (BestTrack_3p_*_index) is used.
    corr=True: use the shoulder-corrected *_corr branches instead.
    forward=True: only tracks with 90 < Phi_3p < 270 deg (Phi_3p is computed in the reco from the same two
    best 3p tracks; it is in (90, 270) exactly when the xz slope is negative)."""
    c = "_corr" if corr else ""
    stop = None if max_events is None else first_event + max_events
    branches = _needed_branches(corr) + (["Phi_3p"] if forward else [])
    with uproot.open(path) as f:
        tree = f[TREE_NAME]
        missing = set(branches) - set(tree.keys())
        if missing:
            raise MissingBranchesError(path, missing)
        arr = tree.arrays(branches, entry_start=first_event, entry_stop=stop, library="ak")
    entry = np.arange(len(arr)) + first_event
    cxz = ak.to_numpy(arr[f"BestTrack_3p_ChiSquare_xz{c}"]).astype(float)
    cxy = ak.to_numpy(arr[f"BestTrack_3p_ChiSquare_xy{c}"]).astype(float)
    sel = np.isfinite(cxz) & np.isfinite(cxy) & (cxz > 0) & (cxy > 0)
    if forward:
        phi = ak.to_numpy(arr["Phi_3p"]).astype(float)
        sel &= (phi > 90) & (phi < 270)
    arr, entry, cxz, cxy = arr[sel], entry[sel], cxz[sel], cxy[sel]
    idx = {"Z": ak.to_numpy(arr[f"BestTrack_3p_xz_index{c}"]).astype(np.int64),
           "Y": ak.to_numpy(arr[f"BestTrack_3p_xy_index{c}"]).astype(np.int64)}
    offs = {"Z": geom.z_add, "Y": geom.y_add}
    u, res = {}, {}
    for L in LETTERS:
        cols_u, cols_r = [], []
        for k in (1, 2, 3):
            if corr:
                pos = _pick(arr[f"TrackClusterPosition_{L.lower()}{k}_corr"], idx[L])
            else:
                cl = _pick(arr[f"TrackCluster_{L.lower()}{k}_index"], idx[L])
                ok = np.isfinite(cl)
                pos = np.full(len(arr), np.nan)
                if ok.any():
                    pos[ok] = _pick(arr[f"ClusterPosition_{L}{k}"][ok], cl[ok].astype(np.int64))
            cols_u.append(pos + offs[L][k - 1])
            cols_r.append(_pick(arr[f"Residue_Track3p_{L.lower()}{k}{c}"], idx[L]))
        u[L], res[L] = np.stack(cols_u, 1), np.stack(cols_r, 1)
    good = np.isfinite(u["Z"]).all(1) & np.isfinite(u["Y"]).all(1)
    run = ak.to_numpy(arr["Run"]).astype(int)
    out = TrackSample(u, res, {"xz": cxz, "xy": cxy}, run, entry)
    n_bad = int((~good).sum())
    if n_bad:
        print(f"  [info] {Path(path).name}: dropped {n_bad} tracks with an unavailable cluster")
    return out.subset(good)


def read_samples(paths, geom, max_events=None, chi2_max=None, first_event=0, corr=False, forward=False):
    parts = []
    for p in paths:
        try:
            s = read_sample(p, geom, max_events=max_events, first_event=first_event, corr=corr, forward=forward)
        except MissingBranchesError as exc:
            print(f"  [skip] {exc.path.name}: missing branches; no usable alignment tree")
            continue
        print(f"  {Path(p).name}: {len(s):,} 3D 3p tracks")
        parts.append(s)
    if not parts:
        return TrackSample({L: np.empty((0, 3)) for L in LETTERS},
                           {L: np.empty((0, 3)) for L in LETTERS},
                           {p: np.empty(0) for p in ("xz", "xy")},
                           np.empty(0, dtype=int), np.empty(0, dtype=int))
    cat = lambda f: np.concatenate([f(s) for s in parts])  # noqa: E731
    sample = TrackSample({L: cat(lambda s, L=L: s.u[L]) for L in LETTERS},
                         {L: cat(lambda s, L=L: s.stored_res[L]) for L in LETTERS},
                         {p: cat(lambda s, p=p: s.stored_chi2[p]) for p in ("xz", "xy")},
                         cat(lambda s: s.run), cat(lambda s: s.event))
    if chi2_max is not None:
        keep = (sample.stored_chi2["xz"] < chi2_max) & (sample.stored_chi2["xy"] < chi2_max)
        print(f"  chi2 < {chi2_max}: keeping {keep.sum():,} / {len(sample):,} tracks")
        sample = sample.subset(keep)
    return sample


# ----------------------------------------------------------------------------
# alignment parameters
# ----------------------------------------------------------------------------
def zero_params():
    return {p: np.zeros(len(PARAMS)) for p in PLANES}


def params_to_json(params, errors=None):
    out = {}
    for p, v in params.items():
        out[p] = {n: float(v[i]) for i, n in enumerate(PARAMS)}
        if errors is not None and p in errors:
            out[p].update({f"{n}_err": float(errors[p][i]) for i, n in enumerate(PARAMS)})
    return out


def params_from_json(d):
    out = zero_params()
    for p in PLANES:
        if p in d:
            out[p] = np.array([d[p][n] for n in PARAMS], float)
    return out


# ----------------------------------------------------------------------------
# track model with the current alignment
# ----------------------------------------------------------------------------
def corrected_hits(geom, L, u_raw, v_loc, params):
    """Hits moved to their corrected positions (see module docstring).
    u_raw (N,3) raw hit positions; v_loc (N,3) coordinate along the strips relative to the
    plane centre (from the other projection of the track).  Returns x', u' (N,3)."""
    xs = np.empty_like(u_raw)
    us = np.empty_like(u_raw)
    for k in range(3):
        du, dw, da, db, dg = params[f"{L}{k + 1}"]
        ul = u_raw[:, k] - geom.u_centre[L][k]
        vl = v_loc[:, k]
        xs[:, k] = geom.x_of[L][k] + dw + db * ul + da * vl
        us[:, k] = u_raw[:, k] + du - dg * vl
    return xs, us


def fit_lines(xs, us):
    """Equal-weight straight-line fit u = a + b x through 3 points per track.
    Returns a, b, fit (N,3), residuals r = u - fit, leverage h (N,3)."""
    xm = xs.mean(1, keepdims=True)
    um = us.mean(1, keepdims=True)
    sxx = ((xs - xm) ** 2).sum(1, keepdims=True)
    b = ((xs - xm) * (us - um)).sum(1, keepdims=True) / sxx
    a = um - b * xm
    fit = a + b * xs
    h = 1.0 / xs.shape[1] + (xs - xm) ** 2 / sxx
    return a[:, 0], b[:, 0], fit, us - fit, h


def fit_projection(geom, sample, L, v_loc, params):
    xs, us = corrected_hits(geom, L, sample.u[L], v_loc, params)
    a, b, fit, r, h = fit_lines(xs, us)
    return {"a": a, "b": b, "xs": xs, "us": us, "fit": fit, "r": r, "h": h,
            "chi2": (r ** 2).sum(1) / geom.sigma[L] ** 2}


def fit_all(geom, sample, params, v_est):
    """both projections; v_est[L] = coordinate along the strips of the Z (Y) planes, i.e. the
    *global* coordinate of the other projection's track at the plane x positions (N,3)."""
    fits = {}
    for L in LETTERS:
        v_loc = v_est[L] - geom.v_centre[L][None, :]
        fits[L] = fit_projection(geom, sample, L, v_loc, params)
        fits[L]["v_loc"] = v_loc
    return fits


def update_v_est(geom, fits):
    """v of the Z planes = y of the Y track at the Z-plane x, v of the Y planes = z of the Z track at the Y-plane x."""
    out = {}
    for L in LETTERS:
        o = fits[OTHER[L]]
        out[L] = o["a"][:, None] + o["b"][:, None] * geom.x_of[L][None, :]
    return out


def initial_v_est(n):
    return {L: np.zeros((n, 3)) for L in LETTERS}


def settle(geom, sample, params, n_pass=3):
    """tracks + v estimates consistent with `params`"""
    v_est = initial_v_est(len(sample))
    # raw v is irrelevant for zero tilts; a few passes make it consistent otherwise
    for _ in range(n_pass):
        fits = fit_all(geom, sample, params, v_est)
        v_est = update_v_est(geom, fits)
    fits = fit_all(geom, sample, params, v_est)
    return fits, v_est


def jacobian(geom, fit, L, k):
    """J (N,5) for plane k of projection L at the current track (paper Eq. 17, meas-fit convention)"""
    b = fit["b"]
    u = fit["fit"][:, k] - geom.u_centre[L][k]
    v = fit["v_loc"][:, k]
    return np.stack([np.ones_like(b), -b, -b * v, -b * u, -v], axis=1)


# ----------------------------------------------------------------------------
# the iterative alignment (paper IV.B, Eq. 15)
# ----------------------------------------------------------------------------
def unbiased_residuals(fits):
    """copy of `fits` with r replaced by the unbiased residual e = r / (1 - h_kk): hit minus the line
    through the OTHER two planes (exact identity of the least-squares fit, equal weights)"""
    return {L: {**f, "r": f["r"] / (1.0 - f["h"])} for L, f in fits.items()}


def plane_normal_eq(geom, fits, L, k, free_mask, leverage, unbiased=False):
    """returns JTWJ, JTWr for one plane (sensor).  Independent 5x5 problem per sensor (paper).
    unbiased=True: the residual is e = r / (1 - h_kk) (hit not in the track fit), unit weights."""
    fit = fits[L]
    J = jacobian(geom, fit, L, k)[:, free_mask]
    r = fit["r"][:, k]
    if unbiased:
        r = r / (1.0 - fit["h"][:, k])
    s2 = geom.sigma[L] ** 2
    w = (1.0 - fit["h"][:, k]) if leverage else np.ones(len(r))
    return (J * w[:, None]).T @ J / s2, J.T @ r / s2, J, r, w


def solve_step(A, g):
    try:
        return -np.linalg.solve(A, g)
    except np.linalg.LinAlgError:
        print("  [warn] singular normal matrix: using pseudo-inverse")
        return -np.linalg.pinv(A) @ g


def median_step(J, r, w, n_iter=100, eps=1e-7):
    """median analogue of solve_step: the step d minimising sum_i |r_i + w_i J_i d| (least absolute
    deviations, by iteratively reweighted least squares).  For du alone this is -median(r) / w, i.e. the
    alignment converges to zero MEDIAN residual instead of zero mean.  eps [m] regularises |residual| ~ 0."""
    Jw = J * w[:, None]
    d = solve_step(Jw.T @ Jw, Jw.T @ r)
    for _ in range(n_iter):
        u = 1.0 / np.maximum(np.abs(r + Jw @ d), eps)
        d_new = solve_step((Jw * u[:, None]).T @ Jw, (Jw * u[:, None]).T @ r)
        done = np.abs(d_new - d).max() < 1e-9
        d = d_new
        if done:
            break
    return d


def gauss_core_fit(r, window="quantile", nbins=40, max_iter=20):
    """Gaussian fit A*exp(-(x-mu)^2 / (2 sigma^2)) to the core of the residual distribution r (histogram,
    Poisson errors).  Fit window:
      "quantile"  : the central 68% of the residuals (16th - 84th percentile)
      "iterative" : start from the quantile window, then refit in [mu - sigma, mu + sigma] of the previous fit
                    until the window is stable.
    Returns dict(mu, sigma, amp, *_err, chi2, ndf, lo, hi, edges, counts, n_in)."""
    from scipy.optimize import curve_fit

    def gauss(x, amp, mu, sigma):
        return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    lo, hi = np.percentile(r, [16, 84])
    for i in range(max_iter if window == "iterative" else 1):
        new_edges = np.linspace(lo, hi, nbins + 1)
        new_n, _ = np.histogram(r, new_edges)
        new_ok = new_n > 0
        ctr = 0.5 * (new_edges[1:] + new_edges[:-1])
        p0 = [new_n.max(), r[(r >= lo) & (r <= hi)].mean(), 0.5 * (hi - lo)]
        try:
            popt_i, pcov_i = curve_fit(gauss, ctr[new_ok], new_n[new_ok], p0=p0, sigma=np.sqrt(new_n[new_ok]),
                                       absolute_sigma=True)
        except (RuntimeError, ValueError):   # no convergence, or empty window
            if i == 0:
                raise
            # iterative refit failed (low statistics): keep the previous window's fit
            lo, hi = edges[0], edges[-1]
            break
        edges, n, ok, popt, pcov = new_edges, new_n, new_ok, popt_i, pcov_i
        amp, mu, sigma = popt[0], popt[1], abs(popt[2])
        new_lo, new_hi = mu - sigma, mu + sigma
        done = max(abs(new_lo - lo), abs(new_hi - hi)) < 1e-3 * (hi - lo)
        if window == "iterative" and not done:
            lo, hi = new_lo, new_hi
        else:
            break
    ctr = 0.5 * (edges[1:] + edges[:-1])
    err = np.sqrt(np.diag(pcov))
    chi2 = float((((n[ok] - gauss(ctr[ok], *popt)) / np.sqrt(n[ok])) ** 2).sum())
    return {"mu": mu, "sigma": sigma, "amp": amp, "mu_err": err[1], "sigma_err": err[2], "amp_err": err[0],
            "chi2": chi2, "ndf": int(ok.sum()) - 3, "lo": lo, "hi": hi, "edges": edges, "counts": n,
            "n_in": int(n.sum())}


def gauss_step(A, J, r, s2, window="quantile"):
    """solve_step with the mean of the residuals replaced by the mean of a Gaussian fitted to their core
    (gauss_core_fit).  Shifting all residuals by a constant only changes the du component of the step, so du
    converges to zero GAUSSIAN-CORE mean residual while the other parameters are the least-squares ones."""
    mu = gauss_core_fit(r, window)["mu"]
    return solve_step(A, J.T @ (r - r.mean() + mu) / s2)


def snapshot(it, params, fits, geom, delta=None, gauss_window=None):
    st = {"iteration": it,
          "params": {p: v.copy() for p, v in params.items()},
          "chi2_mean": {PROJ_OF[L]: float(fits[L]["chi2"].mean()) for L in LETTERS},
          "res_mean": {f"{L}{k + 1}": float(fits[L]["r"][:, k].mean()) for L in LETTERS for k in range(3)},
          "res_rms": {f"{L}{k + 1}": float(fits[L]["r"][:, k].std()) for L in LETTERS for k in range(3)}}
    if gauss_window is not None:
        st["res_gauss_mu"] = {f"{L}{k + 1}": float(gauss_core_fit(fits[L]["r"][:, k], gauss_window)["mu"])
                              for L in LETTERS for k in range(3)}
    if delta is not None:
        st["delta"] = {p: d.copy() for p, d in delta.items()}
    return st


def check_free_planes(free_planes, leverage, unbiased=False):
    """True if some projection has more than one free plane (degenerate for 3-plane tracks).
    Raises ValueError for the combinations that are not valid (leverage weights or unbiased residuals
    + several free planes; leverage + unbiased)."""
    if leverage and unbiased:
        raise ValueError("leverage and unbiased are alternatives (both correct for the (1-h_kk) of the "
                         "biased residuals); use only one of them.")
    multi = max(sum(p[0] == L for p in free_planes) for L in LETTERS) > 1
    if multi and leverage:
        raise ValueError("leverage=True is only valid with ONE free plane per projection (the (1-h_kk) "
                         "weights ignore the coupling between free planes and the iteration diverges).")
    if multi and unbiased:
        raise ValueError("unbiased=True is only valid with ONE free plane per projection (the unbiased "
                         "residual of a plane is measured w.r.t. the other two planes, which must be fixed).")
    return multi


def station_maps(free_mask):
    """Rigid station: the Z and Y layer share the shifts (dx, dy, dz); the tilts stay per layer.
    Station parameters q = (dx, dy, dz, dalpha_Z, dbeta_Z, dgamma_Z, dalpha_Y, dbeta_Y, dgamma_Y).
    Returns {L: C_L} with p_L[free_mask] = C_L @ q[qmask]   (Z: du = dz, dw = dx;  Y: du = dy, dw = dx)."""
    C = {L: np.zeros((len(PARAMS), 9)) for L in LETTERS}
    for L, i_u, i_t in (("Z", 2, 3), ("Y", 1, 6)):
        C[L][0, i_u] = 1.0
        C[L][1, 0] = 1.0
        C[L][2:, i_t:i_t + 3] = np.eye(3)
        C[L] = C[L][free_mask]
    qmask = (C["Z"] != 0).any(0) | (C["Y"] != 0).any(0)
    return {L: C[L][:, qmask] for L in LETTERS}


def check_rigid_stations(rigid_stations, free_planes, estimator="mean"):
    for s in rigid_stations:
        if f"Z{s}" not in free_planes or f"Y{s}" not in free_planes:
            raise ValueError(f"rigid station {s}: both Z{s} and Y{s} must be free planes")
    if rigid_stations and estimator != "mean":
        raise ValueError("rigid stations are only implemented for estimator='mean'")


def run_alignment(geom, sample, free_planes=("Z2", "Y2"), free_params=tuple(PARAMS),
                  params0=None, max_iter=30, tol_len=2e-6, tol_ang=2e-6,
                  leverage=False, unbiased=False, estimator="mean", gauss_window="quantile",
                  rigid_stations=(), verbose=True):
    """Iterate:  tracks refitted with the current alignment -> residuals -> dp per sensor.

    leverage=False : the update exactly as in the paper (residuals of the refitted tracks,
                     tracks held fixed in the Jacobian).  Converges geometrically (factor
                     ~1/3 per iteration for the middle plane of 3 equidistant planes).
    leverage=True  : same minimum, faster: the normal matrix is weighted with (1 - h_kk), the
                     part of the hit that the track fit cannot absorb.
    unbiased=True  : the update uses the unbiased residual e = r / (1 - h_kk) (hit of the aligned
                     plane not in the track fit).  Same minimum as leverage=True (h_kk is the same
                     for all tracks), and the history holds unbiased residuals.
    estimator      : "mean" = least-squares step (paper, zero mean residual at convergence);
                     "median" = least-absolute-deviations step (zero median residual, see median_step);
                     "gauss" = zero mean of a Gaussian fitted to the residual core (see gauss_step),
                     fit window gauss_window = "quantile" (central 68%) or "iterative" (mu +- sigma).
    rigid_stations : stations (1-3) whose Z and Y layer move together: ONE shift (dx, dy, dz) for the station
                     (Z: du = dz, dw = dx; Y: du = dy, dw = dx), solved from the summed normal equations of
                     both projections (station_maps); the tilts, if free, stay per layer.
    Returns dict(params, errors, history, fits, free_planes, ...)"""
    free_planes = [p.upper() for p in free_planes]
    free_mask = np.array([n in free_params for n in PARAMS])
    rigid_stations = [int(s) for s in rigid_stations]
    check_rigid_stations(rigid_stations, free_planes, estimator)
    C = station_maps(free_mask)
    if check_free_planes(free_planes, leverage, unbiased) and verbose:
        print("  [warn] more than one free plane in a projection: 3-plane tracks only constrain ONE "
              "combination of the plane positions per feature.\n"
              "         The result is one of many equivalent solutions (the degenerate directions keep "
              "drifting, convergence is slow) and the printed errors assume all other planes fixed.\n"
              "         Prefer fixing the outer planes (default: --free Z2 Y2).")
    params = zero_params() if params0 is None else {p: v.copy() for p, v in params0.items()}
    fits, v_est = settle(geom, sample, params)
    history = []
    converged = False
    for it in range(max_iter + 1):
        fits = fit_all(geom, sample, params, v_est)
        v_est = update_v_est(geom, fits)
        delta = {}
        for s in rigid_stations:
            Aq, gq = 0.0, 0.0
            for L in LETTERS:
                A, g, *_ = plane_normal_eq(geom, fits, L, s - 1, free_mask, leverage, unbiased)
                Aq, gq = Aq + C[L].T @ A @ C[L], gq + C[L].T @ g
            dq = solve_step(Aq, gq)
            for L in LETTERS:
                delta[f"{L}{s}"] = np.zeros(len(PARAMS))
                delta[f"{L}{s}"][free_mask] = C[L] @ dq
        for plane in free_planes:
            if plane in delta:
                continue
            L, k = plane[0], int(plane[1]) - 1
            A, g, J, r, w = plane_normal_eq(geom, fits, L, k, free_mask, leverage, unbiased)
            d = np.zeros(len(PARAMS))
            if estimator == "mean":
                d[free_mask] = solve_step(A, g)
            elif estimator == "gauss":
                d[free_mask] = gauss_step(A, J, r, geom.sigma[L] ** 2, gauss_window)
            else:
                d[free_mask] = median_step(J, r, w)
            delta[plane] = d
        history.append(snapshot(it, params, unbiased_residuals(fits) if unbiased else fits, geom, delta,
                                gauss_window if estimator == "gauss" else None))
        step_len = max([np.abs(d[:2]).max() for d in delta.values()] + [0.0])
        step_ang = max([np.abs(d[2:]).max() for d in delta.values()] + [0.0])
        if verbose:
            print(f"  iter {it:2d}: <chi2> xz={history[-1]['chi2_mean']['xz']:.4f} "
                  f"xy={history[-1]['chi2_mean']['xy']:.4f}  max step: "
                  f"{step_len * 1e6:9.2f} um / {step_ang * 1e6:9.2f} urad")
        if step_len < tol_len and step_ang < tol_ang:
            converged = True
            break
        if it == max_iter:
            break
        for plane, d in delta.items():
            params[plane] = params[plane] + d
    if verbose:
        print("  converged." if converged else f"  [warn] not converged after {max_iter} iterations")
    errors = parameter_errors(geom, fits, free_planes, free_mask, rigid_stations)
    return {"params": params, "errors": errors, "history": history, "fits": fits,
            "free_planes": free_planes, "free_params": [n for n in PARAMS if n in free_params],
            "rigid_stations": rigid_stations, "converged": converged, "v_est": v_est}


def _plane_hessian(geom, fits, L, k, free_mask):
    """sum (1 - h_kk) J^T J / s^2 with s^2 = < r^2 / (1 - h_kk) > estimated from the data"""
    fit = fits[L]
    J = jacobian(geom, fit, L, k)[:, free_mask]
    lev = 1.0 - fit["h"][:, k]
    s2 = np.mean(fit["r"][:, k] ** 2 / lev)
    return (J * lev[:, None]).T @ J / s2


def parameter_errors(geom, fits, free_planes, free_mask, rigid_stations=()):
    """Statistical errors of the parameters of the free planes.
    The hessian of the chi2 with the track parameters profiled out is
        sum_tracks (1 - h_kk) J^T J / sigma^2
    and the hit resolution is estimated from the data: s^2 = < r^2 / (1 - h_kk) >,
    so  Cov = s^2 [ sum (1 - h_kk) J^T J ]^-1.   (stat. only; independent tracks assumed)
    Rigid station: Cov_q = [ sum_L C_L^T H_L C_L ]^-1, propagated to each layer as C_L Cov_q C_L^T."""
    err = {}
    C = station_maps(free_mask)
    for s in rigid_stations:
        H = {L: _plane_hessian(geom, fits, L, s - 1, free_mask) for L in LETTERS}
        cov_q = np.linalg.pinv(sum(C[L].T @ H[L] @ C[L] for L in LETTERS))
        for L in LETTERS:
            e = np.full(len(PARAMS), np.nan)
            e[free_mask] = np.sqrt(np.diag(C[L] @ cov_q @ C[L].T))
            err[f"{L}{s}"] = e
    for plane in free_planes:
        if plane in err:
            continue
        L, k = plane[0], int(plane[1]) - 1
        e = np.full(len(PARAMS), np.nan)
        e[free_mask] = np.sqrt(np.diag(np.linalg.pinv(_plane_hessian(geom, fits, L, k, free_mask))))
        err[plane] = e
    return err


def evaluate(geom, sample, params):
    """fits for a given set of parameters (no fitting of the parameters)"""
    fits, v_est = settle(geom, sample, params)
    return fits


# ----------------------------------------------------------------------------
# sanity check of the track model against the numbers stored in the ROOT file
# ----------------------------------------------------------------------------
def check_against_stored(geom, sample, fits0):
    """largest |my residual - Residue_Track3p_*| and |my chi2 - BestTrack_3p_ChiSquare_*|"""
    dres = max(np.abs(fits0[L]["r"] - sample.stored_res[L]).max() for L in LETTERS)
    dchi = max(np.abs(fits0[L]["chi2"] - sample.stored_chi2[PROJ_OF[L]]).max() for L in LETTERS)
    return dres, dchi


def reco_equivalent_offsets(geom, params):
    """The parts of the alignment that can be expressed as the offsets of the reconstruction geometry
    (reco_parameters.json): z_add / y_add  (+= du) and the plane x position (+= dw).
    The tilts dalpha, dbeta, dgamma have no counterpart there.  The fit gives dw separately for the Z and
    the Y plane of a station: reco x_pos (Z view) and x_pos_y (Y view)."""
    k3 = range(3)
    return {
        "z_add": [float(geom.z_add[k] + params[f"Z{k + 1}"][0]) for k in k3],
        "y_add": [float(geom.y_add[k] + params[f"Y{k + 1}"][0]) for k in k3],
        "x_pos_for_Z_planes": [float(geom.x_reco["Z"][k] + params[f"Z{k + 1}"][1]) for k in k3],
        "x_pos_for_Y_planes": [float(geom.x_reco["Y"][k] + params[f"Y{k + 1}"][1]) for k in k3],
        "note": "tilts (dalpha, dbeta, dgamma) cannot be expressed as offsets; metres",
    }


def same_geometry(g1, g2):
    return (np.allclose(g1.x_pos, g2.x_pos) and np.allclose(g1.x_pos_y, g2.x_pos_y)
            and np.allclose(g1.z_add, g2.z_add)
            and np.allclose(g1.y_add, g2.y_add) and g1.detector == g2.detector
            and np.isclose(g1.yz_spacing, g2.yz_spacing)
            and all(np.isclose(g1.sigma[L], g2.sigma[L]) for L in LETTERS))