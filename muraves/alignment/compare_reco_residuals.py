#!/usr/bin/env python3
"""Compare stored 3p residuals and track chi2 before and after reconstruction alignment."""
# python compare_reco_residuals.py --runs 4000 4100

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

import alignment_plots as ap
import geometry_io as gio
import muraves_alignment as ma


RUN_RE = re.compile(r"MURAVES_AnalyzedData_run(\d+)\.root$", re.IGNORECASE)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--before-dir",
        type=Path,
        default=Path("/user/dgeeraer/MURAVES/test_data/alignment/RECONSTRUCTED/ROSSO/v1"),
        help="directory containing the pre-alignment ROOT files",
    )
    parser.add_argument(
        "--after-dir",
        type=Path,
        default=Path("/user/dgeeraer/MURAVES/test_data/alignment/AFTER_ALIGNMENT/RECONSTRUCTED/ROSSO/v1"),
        help="directory containing the post-alignment ROOT files",
    )
    parser.add_argument("--runs", nargs=2, type=int, metavar=("FIRST", "LAST"), default=(4000, 4100))
    parser.add_argument("--good-runs", type=Path,
                        help="text file with the good runs; without it ALL runs in the range are used")
    parser.add_argument("--detector", default="ROSSO", choices=("ROSSO", "NERO", "BLU"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/user/dgeeraer/MURAVES/Software/muraves/alignment/images/AFTER_RECO"),
    )
    return parser.parse_args()


def files_by_run(directory: Path, first: int, last: int) -> dict[int, Path]:
    if not directory.is_dir():
        raise SystemExit(f"directory not found: {directory}")
    found = {}
    for path in directory.rglob("MURAVES_AnalyzedData_run*.root"):
        match = RUN_RE.fullmatch(path.name)
        if match is None:
            continue
        run = int(match.group(1))
        if first <= run <= last:
            if run in found:
                raise SystemExit(f"two analyzed ROOT files for run {run} under {directory}")
            found[run] = path
    return found


def read_tracks(path: Path, detector: str) -> dict[str, dict[str, np.ndarray]] | None:
    """stored residuals r (N,3) [m] and chi2 (N,) per projection, in the format of the alignment fits"""
    try:
        geom = ma.load_geometry(path, detector)
        sample = ma.read_sample(path, geom)
    except ma.MissingBranchesError as exc:
        print(f"[skip] {exc.path.name}: missing branches")
        return None
    if len(sample) == 0:
        print(f"[skip] {path.name}: no 3D 3p tracks")
        return None
    ok = np.isfinite(sample.stored_res["Z"]).all(1) & np.isfinite(sample.stored_res["Y"]).all(1)
    return {L: {"r": sample.stored_res[L][ok], "chi2": sample.stored_chi2[ma.PROJ_OF[L]][ok]}
            for L in ma.LETTERS}


def concat(parts: list[dict]) -> dict[str, dict[str, np.ndarray]]:
    return {L: {key: np.concatenate([p[L][key] for p in parts]) for key in ("r", "chi2")} for L in ma.LETTERS}


def main() -> None:
    args = parse_args()
    first, last = sorted(args.runs)
    before_files = files_by_run(args.before_dir, first, last)
    after_files = files_by_run(args.after_dir, first, last)
    common = sorted(set(before_files) & set(after_files))
    if args.good_runs:
        good, _ = gio.parse_run_list(args.good_runs)
        print(f"good-run list {args.good_runs}: keeping {len(set(common) & good)} of {len(common)} runs")
        common = [run for run in common if run in good]
    if not common:
        raise SystemExit("no matching before/after analyzed ROOT files found")

    before_parts, after_parts = [], []
    for run in common:
        print(f"run {run}:")
        old = read_tracks(before_files[run], args.detector)
        new = read_tracks(after_files[run], args.detector)
        if old is None or new is None:
            continue
        before_parts.append(old)
        after_parts.append(new)

    if not before_parts:
        raise SystemExit("no complete before/after residual sample found")
    before, after = concat(before_parts), concat(after_parts)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_dir / f"alignment_{args.detector}_runs{first}-{last}_"
    n_tracks = f"{len(before['Z']['r']):,} / {len(after['Z']['r']):,} 3D 3p tracks (before / after)"
    title = f"{args.detector}: stored reco {{}} before/after alignment, runs {first}-{last}\n" \
            f"{len(before_parts)} runs, {n_tracks}"

    output = Path(f"{prefix}residuals_3p_before_after.png")
    ap.plot_residuals_before_after(before, after, output, title.format("3p residuals"))
    print(f"saved {output}")
    output = Path(f"{prefix}track_chi2_before_after.png")
    ap.plot_track_chi2(before, after, output, title.format("track chi2"), source="stored")
    print(f"saved {output}")


if __name__ == "__main__":
    main()
