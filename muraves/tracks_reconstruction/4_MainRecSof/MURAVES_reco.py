from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Sequence
import numpy as np

from SearchFileName import Search_File
from ClusterLists import CreateClusterList, DeterministicSmearingRNG
from Evaluate_Angular_Coordinates import TrackAngularCoordinates
from ReadEvent import ReadEvent
from Tracking import MakeTracks
from reco_config import get_reco_config, resolve_first_existing, set_runtime_config_path

from muraves_lib import file_handler, root_wrapper
from multiprocessing import Pool
from functools import partial
import argparse as argp
import logging
logger = logging.getLogger(__name__)


def _mean_rms(values: Sequence[float]) -> tuple[float, float]:
    """Return arithmetic mean and population RMS for a sequence."""
    if not values:
        return 0.0, 0.0
    mean = sum(values) / len(values)
    rms = math.sqrt(sum((value - mean) ** 2 for value in values) / len(values))
    return mean, rms


def _mirror_azimuth(phi: float, axis: float = 180.0) -> float:
    """Reflect an azimuth (degrees, in [0, 360)) around phi == axis."""
    if phi < 0.0:
        # Sentinel value (e.g. -1.0) meaning "no track found": leave untouched.
        return phi
    return (2.0 * axis - phi) % 360.0


def _read_slow_control(
    slow_control_path: Path, run: int, config: dict
) -> tuple[float, float, float]:
    trigger_rate = 0.0
    temperature = 0.0
    working_point = 0.0

    if not slow_control_path.exists():
        return trigger_rate, temperature, working_point

    slow_control_cfg = config["slow_control"]
    delimiter = str(slow_control_cfg["delimiter"])
    minimum_fields = int(slow_control_cfg["minimum_fields"])
    run_index = int(slow_control_cfg["run_index"])
    trigger_rate_index = int(slow_control_cfg["trigger_rate_index"])
    temperature_index = int(slow_control_cfg["temperature_index"])
    working_point_index = int(slow_control_cfg["working_point_index"])

    with slow_control_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split(delimiter)
            if len(fields) < minimum_fields:
                continue

            try:
                sc_run = int(float(fields[run_index]))
            except ValueError:
                continue

            if sc_run == run:
                try:
                    trigger_rate = float(fields[trigger_rate_index])
                except ValueError:
                    trigger_rate = 0.0
                try:
                    temperature = float(fields[temperature_index])
                except ValueError:
                    temperature = 0.0
                try:
                    working_point = float(fields[working_point_index])
                except ValueError:
                    working_point = 0.0
                break

    return trigger_rate, temperature, working_point


def _load_spiroc_mapping(spiroc_cfg: Path) -> list[int]:
    strips: list[int] = []
    channels: list[int] = []

    with spiroc_cfg.open("r", encoding="utf-8", errors="ignore") as handle:
        for idx, line in enumerate(handle):
            if idx == 0:
                continue
            fields = line.strip().split()
            if len(fields) < 2:
                continue
            try:
                strip = int(fields[0])
                channel = int(fields[1])
            except ValueError:
                continue
            strips.append(strip)
            channels.append(channel)

    # Reorder channel ids by strip number to match the strip-wise ordering used downstream.
    ordered = sorted(range(len(strips)), key=lambda i: strips[i])
    return [channels[i] for i in ordered]


def _load_telescope_config(telescope_cfg: Path) -> tuple[list[int], list[str]]:
    n_stations: list[int] = []
    views: list[str] = []

    with telescope_cfg.open("r", encoding="utf-8", errors="ignore") as handle:
        for idx, line in enumerate(handle):
            if idx == 0:
                continue
            fields = line.split()
            # Expected columns: n_mod n_scheda n_station n_plane type ...
            if len(fields) < 5:
                continue
            try:
                n_station = int(fields[2])
            except ValueError:
                continue
            n_stations.append(n_station)
            views.append(fields[4])

    if not n_stations or not views:
        raise ValueError(
            f"No valid board configuration could be parsed from telescope file: {telescope_cfg}"
        )

    if len(n_stations) != len(views):
        raise ValueError(
            f"Malformed telescope configuration in {telescope_cfg}: "
            f"n_stations={len(n_stations)} differs from views={len(views)}"
        )

    return n_stations, views


def _load_pedestals(
    pedestal_folder: Path,
    n_boards: int,
    sorted_channels: Sequence[int],
) -> tuple[list[list[float]], list[list[float]], list[list[float]]]:
    boards_peds: list[list[float]] = []
    boards_onephes: list[list[float]] = []
    boards_is_copy: list[list[float]] = []

    for board_n in range(n_boards):
        ped_file = pedestal_folder / f"pedestal_{board_n}"

        peds: list[float] = []
        onephes: list[float] = []
        is_copy: list[float] = []

        if ped_file.exists():
            with ped_file.open("r", encoding="utf-8", errors="ignore") as handle:
                for idx, line in enumerate(handle):
                    if idx == 0:
                        continue
                    fields = line.strip().split()
                    if len(fields) < 3:
                        continue
                    try:
                        peds.append(float(fields[1]))
                    except ValueError:
                        peds.append(0.0)
                    try:
                        onephes.append(float(fields[2]))
                    except ValueError:
                        onephes.append(1.0)

                    if len(fields) > 3:
                        try:
                            is_copy.append(float(fields[3]))
                        except ValueError:
                            is_copy.append(0.0)
                    else:
                        is_copy.append(0.0)

        # Keep the same strip ordering used by C++ via sorted channel map.
        sorted_peds = [peds[ch] if 0 <= ch < len(peds) else 0.0 for ch in sorted_channels]
        sorted_onephes = [onephes[ch] if 0 <= ch < len(onephes) else 1.0 for ch in sorted_channels]
        sorted_is_copy = [is_copy[ch] if 0 <= ch < len(is_copy) else 0.0 for ch in sorted_channels]

        boards_peds.append(sorted_peds)
        boards_onephes.append(sorted_onephes)
        boards_is_copy.append(sorted_is_copy)

    return boards_peds, boards_onephes, boards_is_copy


def _safe_get(container: Sequence[Sequence[float]], idx: int) -> Sequence[float]:
    return container[idx] if 0 <= idx < len(container) else []


# ---------------------------------------------------------------------------
# AnalyzedData branches whose per-event array length (cardinality) is
# identical to one or more sibling branches. uproot.mktree() gives every
# jagged/variable-length branch its own private "n<branch>" counter branch by
# default; when several branches share the exact same per-event length this
# duplicates the same length information many times over in the output ROOT
# file (e.g. nClusterPosition_Z1 == nClusterEnergy_Z1 == ... == Nclusters_Z1
# for every event). Grouping branches here routes them to one shared counter
# branch instead.
#
# Every group below was verified empirically against a full reconstruction of
# run 2500 (40000 real events, test_data/RAW_GZ/NERO) by comparing per-event
# list lengths across the whole run, not just assumed from the source code.
# Only groups that matched for *all* events were kept; anything that did not
# match on every single event was left ungrouped (falls back to its own
# default "n<branch>" counter, i.e. today's behaviour), out of caution with
# real physics data.
_SHARED_COUNTER_GROUPS: dict[str, list[str]] = {
    # Per-plane cluster branches: length == Nclusters_<plane> (verified, 0 mismatches/40000).
    # isInTrack_3p_cl<plane> / isInTrack_4p_cl<plane> are only added for planes 1 (Z1/Y1):
    # for planes 2/3 the same flags come out length-0 whenever the *neighbouring* plane has
    # 0 clusters (a pre-existing quirk of Tracking.MakeTracks, unrelated to this branch), so
    # their length does not reliably equal Nclusters_<plane> and grouping them would risk a
    # ValueError (or worse, silently wrong data) at write time -- verified mismatch on
    # 1483-2036 out of 40000 events depending on plane/projection.
    "ClusterGroup_Z1": [
        "ClusterSize_Z1", "ClusterZ1_Texp", "ClusterEnergy_Z1", "ClusterPosition_Z1",
        "isInTrack_3p_clZ1", "isInTrack_4p_clZ1",
    ],
    "ClusterGroup_Z2": ["ClusterSize_Z2", "ClusterZ2_Texp", "ClusterEnergy_Z2", "ClusterPosition_Z2"],
    "ClusterGroup_Z3": ["ClusterSize_Z3", "ClusterZ3_Texp", "ClusterEnergy_Z3", "ClusterPosition_Z3"],
    "ClusterGroup_Z4": ["ClusterSize_Z4", "ClusterZ4_Texp", "ClusterEnergy_Z4", "ClusterPosition_Z4"],
    "ClusterGroup_Y1": [
        "ClusterSize_Y1", "ClusterY1_Texp", "ClusterEnergy_Y1", "ClusterPosition_Y1",
        "isInTrack_3p_clY1", "isInTrack_4p_clY1",
    ],
    "ClusterGroup_Y2": ["ClusterSize_Y2", "ClusterY2_Texp", "ClusterEnergy_Y2", "ClusterPosition_Y2"],
    "ClusterGroup_Y3": ["ClusterSize_Y3", "ClusterY3_Texp", "ClusterEnergy_Y3", "ClusterPosition_Y3"],
    "ClusterGroup_Y4": ["ClusterSize_Y4", "ClusterY4_Texp", "ClusterEnergy_Y4", "ClusterPosition_Y4"],

    # Per-plane strip branches (flattened nested lists, see encode_nested()). There is no
    # pre-existing explicit "total strip count" field for these, but the 3 branches always
    # share the same per-event flattened length within a plane (verified, 0 mismatches/40000).
    # Note: the per-cluster grouping of strips was already lost before this change (the
    # encode_nested __counts breakdown is intentionally left disabled, "Store only flat
    # values like in c++ code") -- that is pre-existing behaviour, unrelated to this fix.
    "StripsGroup_Z1": ["StripsEnergy_Z1", "StripsPosition_Z1", "StripsID_Z1"],
    "StripsGroup_Z2": ["StripsEnergy_Z2", "StripsPosition_Z2", "StripsID_Z2"],
    "StripsGroup_Z3": ["StripsEnergy_Z3", "StripsPosition_Z3", "StripsID_Z3"],
    "StripsGroup_Z4": ["StripsEnergy_Z4", "StripsPosition_Z4", "StripsID_Z4"],
    "StripsGroup_Y1": ["StripsEnergy_Y1", "StripsPosition_Y1", "StripsID_Y1"],
    "StripsGroup_Y2": ["StripsEnergy_Y2", "StripsPosition_Y2", "StripsID_Y2"],
    "StripsGroup_Y3": ["StripsEnergy_Y3", "StripsPosition_Y3", "StripsID_Y3"],
    "StripsGroup_Y4": ["StripsEnergy_Y4", "StripsPosition_Y4", "StripsID_Y4"],

    # Per-3-point-track branches: length == Ntracks_3p_xz / Ntracks_3p_xy (verified, 0
    # mismatches/40000). ExpectedPosition_OnPlane4th_xz/xy is intentionally excluded: it is
    # only appended when the 4th-plane extrapolation lands within the detector acceptance, so
    # its length is <= Ntracks_3p_* and does NOT match it for every event (verified: matches
    # only 9510/40000 xz events and 21566/40000 xy events) -- it keeps its own private counter.
    "Track3pGroup_xz": [
        "Intercept_3p_xz", "Slope_3p_xz", "chiSquare_3p_xz",
        "TrackCluster_z1_index", "TrackCluster_z2_index", "TrackCluster_z3_index",
        "Residue_Track3p_z1", "Residue_Track3p_z2", "Residue_Track3p_z3",
        "TrackEnergy_3p_xz", "Plane4th_isIntercepted_xz",
    ],
    "Track3pGroup_xy": [
        "Intercept_3p_xy", "Slope_3p_xy", "chiSquare_3p_xy",
        "TrackCluster_y1_index", "TrackCluster_y2_index", "TrackCluster_y3_index",
        "Residue_Track3p_y1", "Residue_Track3p_y2", "Residue_Track3p_y3",
        "TrackEnergy_3p_xy", "Plane4th_isIntercepted_xy",
    ],

    # Per-4-point-candidate branches (flattened nested lists), no pre-existing explicit count
    # field, but always share the same per-event flattened length within a side (verified, 0
    # mismatches/40000).
    "Track4pGroup_xz": [
        "Intercept_4p_xz", "Slope_4p_xz", "chiSquare_4p_xz",
        "displacement_p4_xz", "cluster_c4_index_xz", "ScatteringAngle_xz",
    ],
    "Track4pGroup_xy": [
        "Intercept_4p_xy", "Slope_4p_xy", "chiSquare_4p_xy",
        "displacement_p4_xy", "cluster_c4_index_xy", "ScatteringAngle_xy",
    ],
}


def _apply_shared_counter_groups(ak_module, payload: dict[str, object]) -> dict[str, object]:
    """Merge grouped branches into single record-jagged branches sharing one counter.

    uproot.mktree() gives every jagged/variable-length branch its own private
    "n<branch>" counter by default. Passing a plain ``counter_name`` callable
    that maps several *different* branch names to the same counter name is
    NOT safe with the installed uproot version (5.6.9): it corrupts internal
    TTree bookkeeping and raises ``struct.error: required argument is not an
    integer`` at write time (verified with a minimal 3-branch reproduction,
    independent of this pipeline's data).

    The supported way to make several jagged branches share exactly one
    counter branch is to zip them into a single Awkward record array (one
    List[Record] branch instead of several parallel List[float] branches);
    uproot's writer natively splits a List[Record] branch back out into one
    physical TBranch per record field, all driven by one shared counter
    branch (verified with the same minimal reproduction: fields keep their
    original names and round-trip identically, and only one extra counter
    branch is created for the whole group).

    Only branches present in ``payload`` are merged; groups with fewer than
    two members actually present are left untouched.
    """
    merged = dict(payload)
    for group_name, members in _SHARED_COUNTER_GROUPS.items():
        available = [m for m in members if m in merged]
        if len(available) < 2:
            continue
        fields = {m: merged.pop(m) for m in available}
        merged[group_name] = ak_module.zip(fields, depth_limit=None)
    return merged


def write_root_outputs(
    analysis_jsonl: Path,
    mini_summary_json: Path,
    output_dir: Path | None = None,
    config_files: list[Path] | None = None,
    reco_config_file: Path | None = None,
) -> tuple[Path, Path]:
    """Write ROOT TTrees from JSON outputs while preserving jagged event content."""
    import importlib

    try:
        ak = importlib.import_module("awkward")
        np = importlib.import_module("numpy")
        uproot = importlib.import_module("uproot")
    except ModuleNotFoundError as exc:
        missing = getattr(exc, "name", "required dependency")
        raise RuntimeError(
            f"ROOT export requires optional dependency '{missing}'. "
            "Install uproot and awkward in your environment."
        ) from exc

    analysis_jsonl = Path(analysis_jsonl)
    mini_summary_json = Path(mini_summary_json)

    _root_export_t0 = time.time()
    print(f"[progress] ROOT export started for {analysis_jsonl.name}", flush=True)

    if not analysis_jsonl.exists():
        raise FileNotFoundError(f"Missing analysis JSONL file: {analysis_jsonl}")
    if not mini_summary_json.exists():
        raise FileNotFoundError(f"Missing summary JSON file: {mini_summary_json}")

    with mini_summary_json.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)

    event_records: list[dict] = []
    with analysis_jsonl.open("r", encoding="utf-8", errors="ignore") as handle:
        for line_idx, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            try:
                event_records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Invalid JSON on line {line_idx} of {analysis_jsonl}: {exc}"
                ) from exc
            if line_idx % 50000 == 0:
                print(f"[progress] ROOT export parsed {line_idx} JSON lines", flush=True)

    _t_json_read = time.time()
    print(f"[progress] ROOT export: JSON read in {_t_json_read - _root_export_t0:.2f} s ({len(event_records)} events)", flush=True)

    run = int(summary.get("Run", 0))
    export_dir = Path(output_dir) if output_dir is not None else analysis_jsonl.parent
    export_dir.mkdir(parents=True, exist_ok=True)

    def to_root_branch(values: Sequence[object]) -> object:
        if all(isinstance(v, bool) for v in values):
            return np.asarray(values, dtype=np.bool_)
        if all(isinstance(v, int) and not isinstance(v, bool) for v in values):
            return np.asarray(values, dtype=np.int64)
        if all(isinstance(v, (int, float, bool)) for v in values):
            return np.asarray(values, dtype=np.float64)
        return ak.Array(values)

    def is_sequence(value: object) -> bool:
        return isinstance(value, (list, tuple))

    def needs_nested_encoding(values: Sequence[object]) -> bool:
        for value in values:
            if not is_sequence(value):
                continue
            if any(is_sequence(item) for item in value):
                return True
        return False

    def encode_nested(values: Sequence[object]) -> tuple[list[list[object]], list[list[int]]]:
        # uproot cannot directly write list[list[x]] as a single branch in this schema.
        # Flatten inner lists and store their sizes in a companion __counts branch.
        flat_values: list[list[object]] = []
        item_counts: list[list[int]] = []

        for value in values:
            if not is_sequence(value):
                flat_values.append([])
                item_counts.append([])
                continue

            flat_event: list[object] = []
            count_event: list[int] = []

            for item in value:
                if is_sequence(item):
                    count_event.append(len(item))
                    flat_event.extend(item)
                elif item is None:
                    count_event.append(0)
                else:
                    count_event.append(1)
                    flat_event.append(item)

            flat_values.append(flat_event)
            item_counts.append(count_event)

        return flat_values, item_counts

    def add_payload_branches(
        payload: dict[str, object],
        key: str,
        values: Sequence[object],
        nested_keys: list[str],
    ) -> None:
        if needs_nested_encoding(values):
            #print(f"Branch {key}, needs nested encoding!")
            flat_values, item_counts = encode_nested(values)
            payload[key] = to_root_branch(flat_values)
            # Store only flat values like in c++ code
            #payload[f"{key}__counts"] = to_root_branch(item_counts)
            nested_keys.append(key)
            return
        payload[key] = to_root_branch(values)

    analyzed_payload: dict[str, object] = {}
    analyzed_nested_keys: list[str] = []
    if event_records:
        # Build one ROOT branch per JSON key.
        for key in event_records[0].keys():
            add_payload_branches(
                analyzed_payload,
                key,
                [record.get(key) for record in event_records],
                analyzed_nested_keys,
            )
    else:
        analyzed_payload["Run"] = np.asarray([], dtype=np.int64)
        analyzed_payload["Event"] = np.asarray([], dtype=np.int64)

    if analyzed_nested_keys:
        preview = ", ".join(analyzed_nested_keys[:8])
        if len(analyzed_nested_keys) > 8:
            preview += ", ..."
        print(
            "[progress] Encoded nested event branches for ROOT compatibility: "
            f"{preview}",
            flush=True,
        )

    run_info_payload: dict[str, object] = {}
    run_info_nested_keys: list[str] = []
    for key, value in summary.items():
        if isinstance(value, str):
            # Keep Run_info numeric/object payload ROOT-friendly.
            continue
        add_payload_branches(run_info_payload, key, [value], run_info_nested_keys)

    if not run_info_payload:
        run_info_payload["Run"] = np.asarray([run], dtype=np.int64)

    if run_info_nested_keys:
        print(
            "[progress] Encoded nested Run_info branches for ROOT compatibility: "
            + ", ".join(run_info_nested_keys),
            flush=True,
        )

    _t_arrays = time.time()
    print(f"[progress] ROOT export: array building in {_t_arrays - _t_json_read:.2f} s", flush=True)

    analyzed_root = export_dir / f"MURAVES_AnalyzedData_run{run}.root"
    run_info_root = export_dir / f"MURAVES_miniRunTree_run{run}.root"

    all_config_files = list(config_files) if config_files else []
    if reco_config_file is not None:
        all_config_files.append(reco_config_file)

    # Metadata must be appended to the scratch tmp_path (local disk) while it is
    # still open under temp_to_output, before that context manager copies the
    # finished file out to its final destination (which may be on pnfs).
    # add_metadata_to_root reopens the file with ROOT.TFile(..., "UPDATE") for a
    # random-access write, which pnfs does not support reliably once the file has
    # already been copied there.
    meta = None
    # Merge branches with an identical per-event cardinality (see
    # _SHARED_COUNTER_GROUPS above) into single record-jagged branches so uproot writes
    # one shared counter per group instead of a private counter per branch, without
    # renaming any of the original branches (field_name below strips the synthetic
    # group key and keeps each field's original name).
    analyzed_payload_for_root = _apply_shared_counter_groups(ak, analyzed_payload)
    with file_handler.temp_to_output(analyzed_root) as tmp_path:
        with uproot.recreate(tmp_path) as root_file:
            # mktree(payload) writes once; avoid extend() here to prevent duplicated entries.
            tree = root_file.mktree(
                "AnalyzedData",
                analyzed_payload_for_root,
                field_name=lambda outer, inner: inner,
            )
        if all_config_files:
            # Build git + config metadata once and reuse for the second ROOT file
            # to avoid running git subprocess calls and SHA256 twice.
            meta = root_wrapper.add_metadata_to_root(tmp_path, all_config_files)
    _t_analyzed_written = time.time()
    print(f"[progress] ROOT export: AnalyzedData written in {_t_analyzed_written - _t_arrays:.2f} s", flush=True)

    with file_handler.temp_to_output(run_info_root) as tmp_path:
        with uproot.recreate(tmp_path) as root_file:
            # Same pattern for run-level tree.
            tree = root_file.mktree("Run_info", run_info_payload)
            #tree.extend(run_info_payload)
        if all_config_files:
            root_wrapper.add_metadata_to_root(tmp_path, all_config_files, prebuilt_metadata=meta)

    _t_runinfo_written = time.time()
    print(f"[progress] ROOT export: Run_info written in {_t_runinfo_written - _t_analyzed_written:.2f} s", flush=True)

    _t_end = time.time()
    print(
        f"[progress] ROOT export completed in {_t_end - _root_export_t0:.2f} s total: "
        f"{analyzed_root.name}, {run_info_root.name}",
        flush=True,
    )

    return analyzed_root, run_info_root


def run_reconstruction(
    color: str,
    run: int,
    adc_file: Path,
    pedestal_folder: Path,
    reconstructed_path: Path,
    slow_control_file: Path,
    spiroc_mapping_file: Path,
    telescope_config_file: Path,
    write_root: bool = False,
    progress_every: int = 1000,
    cluster_smearing_seed: int | None = None,
    config: dict | None = None,
    reco_config_file: Path | None = None,
) -> tuple[Path, Path]:
    """Run full event reconstruction for one run and emit JSON (and optional ROOT)."""
    print(" ~~~~~~~  Welcome to the MURAVES reconstruction (Python) ~~~~~~~~")
    total_start_wall_time = time.time()
    json_event_output_elapsed = 0.0

    cfg = config or get_reco_config()
    geometry_cfg = cfg["detector_geometry"]
    reco_cfg = cfg["reconstruction"]

    # Geometry
    if color not in geometry_cfg:
        raise ValueError(f"Unsupported detector color: {color}")
    z_add = [float(v) for v in geometry_cfg[color]["z_add"]]
    x_pos = [float(v) for v in geometry_cfg[color]["x_pos"]]

    y_add = [float(v) for v in reco_cfg["y_add"]]

    # SPACIAL RESOLUTION PARAMETERS (C++ values: sigma_z=0.0040, sigma_y=0.0035)
    sigma_z = float(reco_cfg["sigma_z"])
    sigma_y = float(reco_cfg["sigma_y"])

    # CLUSTERING PARAMETERS (C++ values: s1=6, s2=10, s3=2)
    cluster_cfg = reco_cfg["cluster_thresholds"]
    s1 = float(cluster_cfg["s1"]) # cluster strips must have at least this energy to be seed
    s2 = float(cluster_cfg["s2"]) # single strip clusters must have at least this energy to be accepted as 1-strip cluster
    s3 = float(cluster_cfg["s3"]) # adiacent strips must have at least this energy to be merged into cluster

    # TRACKING PARAMETERS
    # C++ code uses 5*sigma_z and 5*sigma_y as proximity cuts for track building, we keep the same values here.
    proximity_multiplier = float(reco_cfg["proximity_sigma_multiplier"])
    proximity_cut_xz = proximity_multiplier * sigma_z # from C++ comment: "point-trackcandidate distance requirement"
    proximity_cut_xy = proximity_multiplier * sigma_y # from C++ comment: "point-trackcandidate distance requirement"

    # CONSTANT PARAMETERS
    n_boards = int(reco_cfg["n_boards"])
    n_channels = int(reco_cfg["n_channels"])
    cluster_lists_cfg = cfg["cluster_lists"]
    read_events_cfg = cfg["read_event"]
    tracking_cfg = cfg["tracking"]
    # additional parameter available in C++ but not used: 
    #const int nInfo = 168;
    #const int nChInfo = 5;
    #const double FirstStripPos = -0.528;
    #const double AdiacentStripsDistance = 0.0165;


    # PATHS
    reconstructed_path.mkdir(parents=True, exist_ok=True)

    analysis_jsonl = reconstructed_path / f"MURAVES_AnalyzedData_run{run}.jsonl"
    mini_summary_json = reconstructed_path / f"MURAVES_miniRunTree_run{run}.json"
    smearing_rng = DeterministicSmearingRNG(cluster_smearing_seed) if cluster_smearing_seed is not None else None

    if cluster_smearing_seed is not None:
        print(f"Using deterministic single-strip smearing seed: {cluster_smearing_seed}")

    adc_file = Path(adc_file)
    if not adc_file.exists():
        raise FileNotFoundError(f"ADC input file not found: {adc_file}")

    print(f"[progress] Input ADC file: {adc_file}", flush=True)

    trigger_rate, temperature, working_point = _read_slow_control(slow_control_file, run, cfg)

    spiroc_cfg = Path(spiroc_mapping_file)
    # The SPiROC mapping file defines the strip-to-channel mapping for each board, which is crucial for correctly interpreting the ADC data and applying pedestals. The C++ code relies on this mapping to reorder channels and access pedestals in the correct order, so we must load it before processing events.
    sorted_channels = _load_spiroc_mapping(spiroc_cfg)

    # Pedestal files contain the baseline ADC values (pedestals) and the conversion factors to photoelectrons (onephes) for each channel. These are essential for calibrating the raw ADC counts into physical energy deposits. The C++ code loads these pedestals before event reconstruction, so we do the same here to ensure that we can apply the correct calibration to each channel's ADC counts during event processing.
    boards_peds, boards_onephes, boards_is1phe_copy = _load_pedestals(
        pedestal_folder, n_boards, sorted_channels
    )

    # Build fixed-size NumPy views once; reused for every event in the run.
    boards_peds_np: list[np.ndarray] = []
    boards_onephes_np: list[np.ndarray] = []
    for board_n in range(n_boards):
        peds_board = boards_peds[board_n] if board_n < len(boards_peds) else []
        onephe_board = boards_onephes[board_n] if board_n < len(boards_onephes) else []

        peds_np = np.zeros(n_channels, dtype=np.float64)
        onephes_np = np.ones(n_channels, dtype=np.float64)

        if peds_board:
            peds_len = min(n_channels, len(peds_board))
            peds_np[:peds_len] = peds_board[:peds_len]
        if onephe_board:
            onephe_len = min(n_channels, len(onephe_board))
            onephes_np[:onephe_len] = onephe_board[:onephe_len]

        boards_peds_np.append(peds_np)
        boards_onephes_np.append(onephes_np)

    telescope_cfg = Path(telescope_config_file)
    n_stations, views = _load_telescope_config(telescope_cfg)
    if len(n_stations) < n_boards or len(views) < n_boards:
        raise ValueError(
            f"Incomplete telescope configuration in {telescope_cfg}: "
            f"expected at least {n_boards} board rows, got {len(n_stations)}"
        )

    # Prepare variables for ADC file reading loop and event reconstruction.
    ev = 0
    n_events_with_3p = 0
    n_events_with_4p = 0
    n_events_with_3p_xy = 0
    n_events_with_3p_xz = 0
    
    # Cluster counts per plane (?) x1, x2, x3, x4: ? - y1, y2, y3, y4: plane with scintillating bars directed along y-axis, i.e. measuring x coordinate ?
    n_events_with_x1_cl = 0
    n_events_with_x2_cl = 0
    n_events_with_x3_cl = 0
    n_events_with_x4_cl = 0
    n_events_with_y1_cl = 0
    n_events_with_y2_cl = 0
    n_events_with_y3_cl = 0
    n_events_with_y4_cl = 0

    n_tracks_xy_3p = 0
    n_tracks_xz_3p = 0
    n_tracks_3p = 0
    n_tracks_xy_4p = 0
    n_tracks_xz_4p = 0
    n_tracks_4p = 0

    n_good_tracks_3p = 0
    n_good_tracks_4p = 0

    first_timestamp = None
    last_timestamp = None

    # (?)
    cluster_energy_acc = {
        "Z1": [],
        "Z2": [],
        "Z3": [],
        "Z4": [],
        "Y1": [],
        "Y2": [],
        "Y3": [],
        "Y4": [],
    }
    ncluster_acc = {
        "Z1": [],
        "Z2": [],
        "Z3": [],
        "Z4": [],
        "Y1": [],
        "Y2": [],
        "Y3": [],
        "Y4": [],
    }


    # READ ADC FILE AND RECONSTRUCT EVENTS
    with adc_file.open("r", encoding="utf-8", errors="ignore") as adc_handle, analysis_jsonl.open(
        "w", encoding="utf-8"
    ) as out_events:
        start_wall_time = time.time()
        for event in adc_handle:
            ev += 1

            event_info = ReadEvent(event, sorted_channels, events_cfg=read_events_cfg)
            timestamp = event_info.timeStamp

            if first_timestamp is None:
                first_timestamp = timestamp
            last_timestamp = timestamp

            trigger_mask_strips = event_info.TrMask_strips
            trigger_mask_channels = event_info.TrMask_channels
            trigger_mask_size = event_info.TrMask_size

            # Deposits per module 
            deposits_p1x1: list[float] = []
            deposits_p1x2: list[float] = []
            deposits_p2x1: list[float] = []
            deposits_p2x2: list[float] = []
            deposits_p3x1: list[float] = []
            deposits_p3x2: list[float] = []
            deposits_p4x1: list[float] = []
            deposits_p4x2: list[float] = []
            deposits_p1y1: list[float] = []
            deposits_p1y2: list[float] = []
            deposits_p2y1: list[float] = []
            deposits_p2y2: list[float] = []
            deposits_p3y1: list[float] = []
            deposits_p3y2: list[float] = []
            deposits_p4y1: list[float] = []
            deposits_p4y2: list[float] = []

            # Time expansions per module
            texp_1x1 = texp_1x2 = texp_2x1 = texp_2x2 = 0.0
            texp_3x1 = texp_3x2 = texp_4x1 = texp_4x2 = 0.0
            texp_1y1 = texp_1y2 = texp_2y1 = texp_2y2 = 0.0
            texp_3y1 = texp_3y2 = texp_4y1 = texp_4y2 = 0.0

            for b in range(n_boards):
                if b >= len(views) or b >= len(n_stations):
                    continue

                view = views[b]
                n_station = n_stations[b]

                adc_counts = event_info.boards[b] if b < len(event_info.boards) else []
                peds_board_np = boards_peds_np[b] if b < len(boards_peds_np) else np.zeros(n_channels, dtype=np.float64)
                onephe_board_np = boards_onephes_np[b] if b < len(boards_onephes_np) else np.ones(n_channels, dtype=np.float64)

                time_exp_value = event_info.timeExp[b] if b < len(event_info.timeExp) else 0.0

                if n_station == 1:
                    if view == "x1":
                        texp_1x1 = time_exp_value
                    if view == "x2":
                        texp_1x2 = time_exp_value
                    if view == "y1":
                        texp_1y1 = time_exp_value
                    if view == "y2":
                        texp_1y2 = time_exp_value
                if n_station == 2:
                    if view == "x1":
                        texp_2x1 = time_exp_value
                    if view == "x2":
                        texp_2x2 = time_exp_value
                    if view == "y1":
                        texp_2y1 = time_exp_value
                    if view == "y2":
                        texp_2y2 = time_exp_value
                if n_station == 3:
                    if view == "x1":
                        texp_3x1 = time_exp_value
                    if view == "x2":
                        texp_3x2 = time_exp_value
                    if view == "y1":
                        texp_3y1 = time_exp_value
                    if view == "y2":
                        texp_3y2 = time_exp_value
                if n_station == 4:
                    if view == "x1":
                        texp_4x1 = time_exp_value
                    if view == "x2":
                        texp_4x2 = time_exp_value
                    if view == "y1":
                        texp_4y1 = time_exp_value
                    if view == "y2":
                        texp_4y2 = time_exp_value

                adc_np = np.zeros(n_channels, dtype=np.float64)
                if adc_counts:
                    adc_len = min(n_channels, len(adc_counts))
                    adc_np[:adc_len] = adc_counts[:adc_len]

                # Convert ADC counts to photoelectron-equivalent deposits in one vectorized step.
                deposits_np = np.divide(
                    adc_np - peds_board_np,
                    onephe_board_np,
                    out=np.zeros_like(adc_np),
                    where=onephe_board_np != 0.0,
                )
                deposits_board = deposits_np.tolist()

                if n_station == 1:
                    if view == "x1":
                        deposits_p1x1.extend(deposits_board)
                    if view == "x2":
                        deposits_p1x2.extend(deposits_board)
                    if view == "y1":
                        deposits_p1y1.extend(deposits_board)
                    if view == "y2":
                        deposits_p1y2.extend(deposits_board)
                if n_station == 2:
                    if view == "x1":
                        deposits_p2x1.extend(deposits_board)
                    if view == "x2":
                        deposits_p2x2.extend(deposits_board)
                    if view == "y1":
                        deposits_p2y1.extend(deposits_board)
                    if view == "y2":
                        deposits_p2y2.extend(deposits_board)
                if n_station == 3:
                    if view == "x1":
                        deposits_p3x1.extend(deposits_board)
                    if view == "x2":
                        deposits_p3x2.extend(deposits_board)
                    if view == "y1":
                        deposits_p3y1.extend(deposits_board)
                    if view == "y2":
                        deposits_p3y2.extend(deposits_board)
                if n_station == 4:
                    if view == "x1":
                        deposits_p4x1.extend(deposits_board)
                    if view == "x2":
                        deposits_p4x2.extend(deposits_board)
                    if view == "y1":
                        deposits_p4y1.extend(deposits_board)
                    if view == "y2":
                        deposits_p4y2.extend(deposits_board)

            # Concatenate sub-plane energy deposits to form a full plane deposit. 
            deposits_p1x = deposits_p1x1 + deposits_p1x2
            deposits_p2x = deposits_p2x1 + deposits_p2x2
            deposits_p3x = deposits_p3x1 + deposits_p3x2
            deposits_p4x = deposits_p4x1 + deposits_p4x2

            # Y view uses reversed sub-plane concatenation, same was done in the C++ code. 
            # !!This is a concatenation of lists so it's not commutative!!
            deposits_p1y = deposits_p1y2 + deposits_p1y1
            deposits_p2y = deposits_p2y2 + deposits_p2y1
            deposits_p3y = deposits_p3y2 + deposits_p3y1
            deposits_p4y = deposits_p4y2 + deposits_p4y1

	        #/////////////////////////////////////////////////////////////////////////////////////
	        #//////////////////////////////////            ///////////////////////////////////////
	        #///////////////////////////////// CLUSTERING ///////////////////////////////////////
	        #////////////////////////////////            ///////////////////////////////////////
	        #///////////////////////////////////////////////////////////////////////////////////

            results_p1x = CreateClusterList(
                deposits_p1x,
                s1,
                s2,
                s3,
                texp_1x1,
                texp_1x2,
                # board 0 and 1 are x1 sub-planes
                _safe_get(trigger_mask_strips, 0),
                _safe_get(trigger_mask_strips, 1),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )
            results_p2x = CreateClusterList(
                deposits_p2x,
                s1,
                s2,
                s3,
                texp_2x1,
                texp_2x2,
                # board 4 and 5 are x2 sub-planes
                _safe_get(trigger_mask_strips, 4),
                _safe_get(trigger_mask_strips, 5),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )
            results_p3x = CreateClusterList(
                deposits_p3x,
                s1,
                s2,
                s3,
                texp_3x1,
                texp_3x2,
                # board 8 and 9 are x3 sub-planes
                _safe_get(trigger_mask_strips, 8),
                _safe_get(trigger_mask_strips, 9),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )
            results_p4x = CreateClusterList(
                deposits_p4x,
                s1,
                s2,
                s3,
                texp_4x1,
                texp_4x2,
                # board 12 and 13 are x4 sub-planes
                _safe_get(trigger_mask_strips, 12),
                _safe_get(trigger_mask_strips, 13),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )

            results_p1y = CreateClusterList(
                deposits_p1y,
                s1,
                s2,
                s3,
                texp_1y2,
                texp_1y1,
                # board 2 and 3 are y1 sub-planes, note the reversed order for y views
                _safe_get(trigger_mask_strips, 3),
                _safe_get(trigger_mask_strips, 2),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )
            results_p2y = CreateClusterList(
                deposits_p2y,
                s1,
                s2,
                s3,
                texp_2y2,
                texp_2y1,
                # board 6 and 7 are y2 sub-planes, note the reversed order for y views
                _safe_get(trigger_mask_strips, 7),
                _safe_get(trigger_mask_strips, 6),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )
            results_p3y = CreateClusterList(
                deposits_p3y,
                s1,
                s2,
                s3,
                texp_3y2,
                texp_3y1,
                # board 10 and 11 are y3 sub-planes, note the reversed order for y views
                _safe_get(trigger_mask_strips, 11),
                _safe_get(trigger_mask_strips, 10),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )
            results_p4y = CreateClusterList(
                deposits_p4y,
                s1,
                s2,
                s3,
                texp_4y2,
                texp_4y1,
                # board 14 and 15 are y4 sub-planes, note the reversed order for y views
                _safe_get(trigger_mask_strips, 15),
                _safe_get(trigger_mask_strips, 14),
                smearing_rng=smearing_rng,
                cluster_cfg=cluster_lists_cfg,
            )

            nclusters_z1 = len(results_p1x.ClustersEnergy)
            nclusters_z2 = len(results_p2x.ClustersEnergy)
            nclusters_z3 = len(results_p3x.ClustersEnergy)
            nclusters_z4 = len(results_p4x.ClustersEnergy)

            nclusters_y1 = len(results_p1y.ClustersEnergy)
            nclusters_y2 = len(results_p2y.ClustersEnergy)
            nclusters_y3 = len(results_p3y.ClustersEnergy)
            nclusters_y4 = len(results_p4y.ClustersEnergy)

            # (?) Is this really necessary this change of name from x to z ? To me it's just confusing. 
            # Keeping like this for now to match the C++ code but we should consider refactoring this in the future for clarity.
            if nclusters_z1 > 0:
                n_events_with_x1_cl += 1
            if nclusters_z2 > 0:
                n_events_with_x2_cl += 1
            if nclusters_z3 > 0:
                n_events_with_x3_cl += 1
            if nclusters_z4 > 0:
                n_events_with_x4_cl += 1
            if nclusters_y1 > 0:
                n_events_with_y1_cl += 1
            if nclusters_y2 > 0:
                n_events_with_y2_cl += 1
            if nclusters_y3 > 0:
                n_events_with_y3_cl += 1
            if nclusters_y4 > 0:
                n_events_with_y4_cl += 1

            cluster_energy_acc["Z1"].extend(results_p1x.ClustersEnergy)
            cluster_energy_acc["Z2"].extend(results_p2x.ClustersEnergy)
            cluster_energy_acc["Z3"].extend(results_p3x.ClustersEnergy)
            cluster_energy_acc["Z4"].extend(results_p4x.ClustersEnergy)
            cluster_energy_acc["Y1"].extend(results_p1y.ClustersEnergy)
            cluster_energy_acc["Y2"].extend(results_p2y.ClustersEnergy)
            cluster_energy_acc["Y3"].extend(results_p3y.ClustersEnergy)
            cluster_energy_acc["Y4"].extend(results_p4y.ClustersEnergy)

            ncluster_acc["Z1"].append(nclusters_z1)
            ncluster_acc["Z2"].append(nclusters_z2)
            ncluster_acc["Z3"].append(nclusters_z3)
            ncluster_acc["Z4"].append(nclusters_z4)
            ncluster_acc["Y1"].append(nclusters_y1)
            ncluster_acc["Y2"].append(nclusters_y2)
            ncluster_acc["Y3"].append(nclusters_y3)
            ncluster_acc["Y4"].append(nclusters_y4)

####################################################################################################################################
############################################ DEEPLEY CHECKED UP TO HERE ############################################################
####################################################################################################################################
            tracks_xz = MakeTracks(
                results_p1x.ClustersPositions,
                results_p2x.ClustersPositions,
                results_p3x.ClustersPositions,
                results_p4x.ClustersPositions,
                results_p1x.ClustersEnergy,
                results_p2x.ClustersEnergy,
                results_p3x.ClustersEnergy,
                results_p4x.ClustersEnergy,
                results_p1x.TimeExpansions,
                results_p2x.TimeExpansions,
                results_p3x.TimeExpansions,
                results_p4x.TimeExpansions,
                proximity_cut_xz,
                x_pos,
                z_add,
                sigma_z,
                tracking_cfg=tracking_cfg
            )
            tracks_xy = MakeTracks(
                results_p1y.ClustersPositions,
                results_p2y.ClustersPositions,
                results_p3y.ClustersPositions,
                results_p4y.ClustersPositions,
                results_p1y.ClustersEnergy,
                results_p2y.ClustersEnergy,
                results_p3y.ClustersEnergy,
                results_p4y.ClustersEnergy,
                results_p1y.TimeExpansions,
                results_p2y.TimeExpansions,
                results_p3y.TimeExpansions,
                results_p4y.TimeExpansions,
                proximity_cut_xy,
                x_pos,
                y_add,
                sigma_y,
                tracking_cfg=tracking_cfg
            )

            ntracks_3p_xz = len(tracks_xz.intercepts_3p)
            ntracks_3p_xy = len(tracks_xy.intercepts_3p)

            # CALCULATE THE TOTAL NUMBER OF TRACKS IN THE RUN
            ##### we are still inside the event loop, so this is an accumulation of tracks across events.
            n_tracks_xy_3p += ntracks_3p_xy
            n_tracks_xz_3p += ntracks_3p_xz
            n_tracks_3p += ntracks_3p_xy * ntracks_3p_xz
            ###############################################

            if ntracks_3p_xz > 0 and ntracks_3p_xy > 0:
                n_events_with_3p += 1
            if ntracks_3p_xz > 0:
                n_events_with_3p_xz += 1
            if ntracks_3p_xy > 0:
                n_events_with_3p_xy += 1

            best_track_xy_ind = -1
            best_track_xz_ind = -1
            theta_3p = -1.0
            phi_3p = -1.0

            best_track_3p_chi_xy = -1.0
            best_track_3p_chi_xz = -1.0

            # /////////////////////////// 3D BEST TRACK ////////////////////////////////////////
            if ntracks_3p_xz > 0 and ntracks_3p_xy > 0:
                if tracks_xy.BestChi_index == tracks_xy.BestEnergy_index:
                    best_track_xy_ind = tracks_xy.BestEnergy_index
                else:
                    #//// CRITERIA TO CHOOSE THE BEST TRACK BETWEEN BEST CHI AND BEST ENERGY ///////
                    best_track_xy_ind = tracks_xy.BestChi_index

                if tracks_xz.BestChi_index == tracks_xz.BestEnergy_index:
                    best_track_xz_ind = tracks_xz.BestEnergy_index
                else:
                    #//// CRITERIA TO CHOOSE THE BEST TRACK BETWEEN BEST CHI AND BEST ENERGY ///////
                    best_track_xz_ind = tracks_xz.BestChi_index

                if best_track_xz_ind >= 0 and best_track_xy_ind >= 0:
                    coordinates_3p = TrackAngularCoordinates(
                        tracks_xy.slopes_3p[best_track_xy_ind],
                        tracks_xz.slopes_3p[best_track_xz_ind],
                        x_pos[0],
                        x_pos[2],
                    )
                    theta_3p, phi_3p = coordinates_3p
                    n_good_tracks_3p += 1

                    best_track_3p_chi_xy = tracks_xy.chiSquares_3p[best_track_xy_ind]
                    best_track_3p_chi_xz = tracks_xz.chiSquares_3p[best_track_xz_ind]

            n4p_xy = sum(1 for n in tracks_xy.Ntracks_4p if n > 0)
            n4p_xz = sum(1 for n in tracks_xz.Ntracks_4p if n > 0)
            if n4p_xz > 0 and n4p_xy > 0:
                n_events_with_4p += 1
                n_tracks_xy_4p += n4p_xy
                n_tracks_xz_4p += n4p_xz
                n_tracks_4p += n4p_xz * n4p_xy

            # Best 4p track estimation equivalent to C++ selection.
            theta_4p = -1.0
            phi_4p = -1.0
            best_track_4p_chi_xy = -1.0
            best_track_4p_chi_xz = -1.0
            best_scattering_angle_xy = -1.0
            best_scattering_angle_xz = -1.0

            track_3p_of_4p_idx_xy = -1
            track_3p_of_4p_idx_xz = -1
            track_4p_idx_xy = -1
            track_4p_idx_xz = -1
            best_track_4p_is_texp_null_xy = -1
            best_track_4p_is_texp_null_xz = -1

            cluster_c4_texp_xz: list[list[float]] = []
            for idx_row, row in enumerate(tracks_xz.cluster_indices_4):
                texp_row = []
                for idx_val in row:
                    idx_int = int(idx_val)
                    if 0 <= idx_int < len(results_p4x.TimeExpansions):
                        texp_row.append(results_p4x.TimeExpansions[idx_int])
                    else:
                        texp_row.append(-1.0)
                cluster_c4_texp_xz.append(texp_row)

            cluster_c4_texp_xy: list[list[float]] = []
            for idx_row, row in enumerate(tracks_xy.cluster_indices_4):
                texp_row = []
                for idx_val in row:
                    idx_int = int(idx_val)
                    if 0 <= idx_int < len(results_p4y.TimeExpansions):
                        texp_row.append(results_p4y.TimeExpansions[idx_int])
                    else:
                        texp_row.append(-1.0)
                cluster_c4_texp_xy.append(texp_row)
                
	        #////////////////////////////////////////////////////////////////////////////////
	        #////////////////////////////// BEST TRACK 4 PLANES /////////////////////////////
            if tracks_xz.Track_3p_to_4p_index and tracks_xy.Track_3p_to_4p_index:
                # Only promote 4p candidates attached to the selected best 3p tracks.
                if best_track_xy_ind in tracks_xy.Track_3p_to_4p_index and best_track_xz_ind in tracks_xz.Track_3p_to_4p_index:
                    track_3p_of_4p_idx_xy = tracks_xy.Track_3p_to_4p_index.index(best_track_xy_ind)
                    track_3p_of_4p_idx_xz = tracks_xz.Track_3p_to_4p_index.index(best_track_xz_ind)

                    best_idx_xy = 0
                    best_track_4p_is_texp_null_xy = 1
                    if track_3p_of_4p_idx_xy < len(cluster_c4_texp_xy):
                        for k_idx, texp in enumerate(cluster_c4_texp_xy[track_3p_of_4p_idx_xy]):
                            if texp > 0:
                                best_track_4p_is_texp_null_xy = 0
                                best_idx_xy = k_idx
                                break
                    track_4p_idx_xy = best_idx_xy

                    best_idx_xz = 0
                    best_track_4p_is_texp_null_xz = 1
                    if track_3p_of_4p_idx_xz < len(cluster_c4_texp_xz):
                        for k_idx, texp in enumerate(cluster_c4_texp_xz[track_3p_of_4p_idx_xz]):
                            if texp > 0:
                                best_track_4p_is_texp_null_xz = 0
                                best_idx_xz = k_idx
                                break
                    track_4p_idx_xz = best_idx_xz

                    if (
                        track_3p_of_4p_idx_xy < len(tracks_xy.chiSquares_4p)
                        and track_3p_of_4p_idx_xz < len(tracks_xz.chiSquares_4p)
                        and track_4p_idx_xy < len(tracks_xy.chiSquares_4p[track_3p_of_4p_idx_xy])
                        and track_4p_idx_xz < len(tracks_xz.chiSquares_4p[track_3p_of_4p_idx_xz])
                    ):
                        best_track_4p_chi_xy = tracks_xy.chiSquares_4p[track_3p_of_4p_idx_xy][track_4p_idx_xy]
                        best_track_4p_chi_xz = tracks_xz.chiSquares_4p[track_3p_of_4p_idx_xz][track_4p_idx_xz]

                    if (
                        track_3p_of_4p_idx_xy < len(tracks_xy.ScatteringAngles)
                        and track_3p_of_4p_idx_xz < len(tracks_xz.ScatteringAngles)
                        and track_4p_idx_xy < len(tracks_xy.ScatteringAngles[track_3p_of_4p_idx_xy])
                        and track_4p_idx_xz < len(tracks_xz.ScatteringAngles[track_3p_of_4p_idx_xz])
                    ):
                        best_scattering_angle_xy = tracks_xy.ScatteringAngles[track_3p_of_4p_idx_xy][track_4p_idx_xy]
                        best_scattering_angle_xz = tracks_xz.ScatteringAngles[track_3p_of_4p_idx_xz][track_4p_idx_xz]

                    if (
                        track_3p_of_4p_idx_xy < len(tracks_xy.slope_4p)
                        and track_3p_of_4p_idx_xz < len(tracks_xz.slope_4p)
                        and track_4p_idx_xy < len(tracks_xy.slope_4p[track_3p_of_4p_idx_xy])
                        and track_4p_idx_xz < len(tracks_xz.slope_4p[track_3p_of_4p_idx_xz])
                    ):
                        coordinates_4p = TrackAngularCoordinates(
                            tracks_xy.slope_4p[track_3p_of_4p_idx_xy][track_4p_idx_xy],
                            tracks_xz.slope_4p[track_3p_of_4p_idx_xz][track_4p_idx_xz],
                            x_pos[0],
                            x_pos[3],
                        )
                        theta_4p, phi_4p = coordinates_4p
                        n_good_tracks_4p += 1

            event_record = {
                # Event-level payload kept close to original C++ tree naming.
                "Run": run,

                "Nclusters_Z1": nclusters_z1,
                "Nclusters_Z2": nclusters_z2,
                "Nclusters_Z3": nclusters_z3,
                "Nclusters_Z4": nclusters_z4,
                "Nclusters_Y1": nclusters_y1,
                "Nclusters_Y2": nclusters_y2,
                "Nclusters_Y3": nclusters_y3,
                "Nclusters_Y4": nclusters_y4,

                "ClusterSize_Z1": results_p1x.ClustersSize,
                "ClusterSize_Z2": results_p2x.ClustersSize,
                "ClusterSize_Z3": results_p3x.ClustersSize,
                "ClusterSize_Z4": results_p4x.ClustersSize,
                "ClusterZ1_Texp": results_p1x.TimeExpansions,
                "ClusterZ2_Texp": results_p2x.TimeExpansions,
                "ClusterZ3_Texp": results_p3x.TimeExpansions,
                "ClusterZ4_Texp": results_p4x.TimeExpansions,
                "ClusterEnergy_Z1": results_p1x.ClustersEnergy,
                "ClusterEnergy_Z2": results_p2x.ClustersEnergy,
                "ClusterEnergy_Z3": results_p3x.ClustersEnergy,
                "ClusterEnergy_Z4": results_p4x.ClustersEnergy,
                "ClusterPosition_Z1": results_p1x.ClustersPositions,
                "ClusterPosition_Z2": results_p2x.ClustersPositions,
                "ClusterPosition_Z3": results_p3x.ClustersPositions,
                "ClusterPosition_Z4": results_p4x.ClustersPositions,


                "ClusterSize_Y1": results_p1y.ClustersSize,
                "ClusterSize_Y2": results_p2y.ClustersSize,
                "ClusterSize_Y3": results_p3y.ClustersSize,
                "ClusterSize_Y4": results_p4y.ClustersSize,
                "ClusterY1_Texp": results_p1y.TimeExpansions,
                "ClusterY2_Texp": results_p2y.TimeExpansions,
                "ClusterY3_Texp": results_p3y.TimeExpansions,
                "ClusterY4_Texp": results_p4y.TimeExpansions,
                "ClusterEnergy_Y1": results_p1y.ClustersEnergy,
                "ClusterEnergy_Y2": results_p2y.ClustersEnergy,
                "ClusterEnergy_Y3": results_p3y.ClustersEnergy,
                "ClusterEnergy_Y4": results_p4y.ClustersEnergy,
                "ClusterPosition_Y1": results_p1y.ClustersPositions,
                "ClusterPosition_Y2": results_p2y.ClustersPositions,
                "ClusterPosition_Y3": results_p3y.ClustersPositions,
                "ClusterPosition_Y4": results_p4y.ClustersPositions,

                "StripsEnergy_Z1": results_p1x.StripsEnergy,
                "StripsEnergy_Z2": results_p2x.StripsEnergy,
                "StripsEnergy_Z3": results_p3x.StripsEnergy,
                "StripsEnergy_Z4": results_p4x.StripsEnergy,
                "StripsEnergy_Y1": results_p1y.StripsEnergy,
                "StripsEnergy_Y2": results_p2y.StripsEnergy,
                "StripsEnergy_Y3": results_p3y.StripsEnergy,
                "StripsEnergy_Y4": results_p4y.StripsEnergy,
                "StripsPosition_Z1": results_p1x.StripsPositions,
                "StripsPosition_Z2": results_p2x.StripsPositions,
                "StripsPosition_Z3": results_p3x.StripsPositions,
                "StripsPosition_Z4": results_p4x.StripsPositions,
                "StripsPosition_Y1": results_p1y.StripsPositions,
                "StripsPosition_Y2": results_p2y.StripsPositions,
                "StripsPosition_Y3": results_p3y.StripsPositions,
                "StripsPosition_Y4": results_p4y.StripsPositions,
                "StripsID_Z1": results_p1x.StripsID,
                "StripsID_Z2": results_p2x.StripsID,
                "StripsID_Z3": results_p3x.StripsID,
                "StripsID_Z4": results_p4x.StripsID,
                "StripsID_Y1": results_p1y.StripsID,
                "StripsID_Y2": results_p2y.StripsID,
                "StripsID_Y3": results_p3y.StripsID,
                "StripsID_Y4": results_p4y.StripsID,

                # Tracking branches
                "Ntracks_3p_xz": ntracks_3p_xz,
                "Ntracks_3p_xy": ntracks_3p_xy,
                "Intercept_3p_xz": tracks_xz.intercepts_3p,
                "Slope_3p_xz": tracks_xz.slopes_3p,
                "chiSquare_3p_xz": tracks_xz.chiSquares_3p,
                "Intercept_3p_xy": tracks_xy.intercepts_3p,
                "Slope_3p_xy": tracks_xy.slopes_3p,
                "chiSquare_3p_xy": tracks_xy.chiSquares_3p,
                "TrackCluster_z1_index": tracks_xz.cluster_index_1,
                "TrackCluster_z2_index": tracks_xz.cluster_index_2,
                "TrackCluster_z3_index": tracks_xz.cluster_index_3,
                "TrackCluster_y1_index": tracks_xy.cluster_index_1,
                "TrackCluster_y2_index": tracks_xy.cluster_index_2,
                "TrackCluster_y3_index": tracks_xy.cluster_index_3,

                "Residue_Track3p_z1": tracks_xz.residue_c1,
                "Residue_Track3p_z2": tracks_xz.residue_c2,
                "Residue_Track3p_z3": tracks_xz.residue_c3,
                "Residue_Track3p_y1": tracks_xy.residue_c1,
                "Residue_Track3p_y2": tracks_xy.residue_c2,
                "Residue_Track3p_y3": tracks_xy.residue_c3,
                "TrackEnergy_3p_xy": tracks_xy.TrackEnergy_3p,
                "TrackEnergy_3p_xz": tracks_xz.TrackEnergy_3p,
                "Plane4th_isIntercepted_xz": tracks_xz.Plane4th_isIntercepted,
                "Plane4th_isIntercepted_xy": tracks_xy.Plane4th_isIntercepted,
                "BestChi_xy_index": tracks_xy.BestChi_index,
                "BestEnergy_xy_index": tracks_xy.BestEnergy_index,
                "BestChi_xz_index": tracks_xz.BestChi_index,
                "BestEnergy_xz_index": tracks_xz.BestEnergy_index,
                "BestChi_xy": tracks_xy.BestChi,
                "BestEnergy_xy": tracks_xy.BestEnergy,
                "BestChi_xz": tracks_xz.BestChi,
                "BestEnergy_xz": tracks_xz.BestEnergy,
                "ExpectedPosition_OnPlane4th_xy": tracks_xy.ExpectedPosition_OnPlane4th,
                "ExpectedPosition_OnPlane4th_xz": tracks_xz.ExpectedPosition_OnPlane4th,

                "Intercept_4p_xz": tracks_xz.intercept_4p,
                "Slope_4p_xz": tracks_xz.slope_4p,
                "chiSquare_4p_xz": tracks_xz.chiSquares_4p,
                "displacement_p4_xz": tracks_xz.displacement_p4,
                "Intercept_4p_xy": tracks_xy.intercept_4p,
                "Slope_4p_xy": tracks_xy.slope_4p,
                "chiSquare_4p_xy": tracks_xy.chiSquares_4p,
                "displacement_p4_xy": tracks_xy.displacement_p4,
                "cluster_c4_index_xz": tracks_xz.cluster_indices_4,
                "cluster_c4_index_xy": tracks_xy.cluster_indices_4,
                "ScatteringAngle_xz": tracks_xz.ScatteringAngles,
                "ScatteringAngle_xy": tracks_xy.ScatteringAngles,

                "isInTrack_3p_clZ1": tracks_xz.IsInTrack_clusters1,
                "isInTrack_3p_clZ2": tracks_xz.IsInTrack_clusters2,
                "isInTrack_3p_clZ3": tracks_xz.IsInTrack_clusters3,
                "isInTrack_3p_clY1": tracks_xy.IsInTrack_clusters1,
                "isInTrack_3p_clY2": tracks_xy.IsInTrack_clusters2,   
                "isInTrack_3p_clY3": tracks_xy.IsInTrack_clusters3,
	            "isInTrack_4p_clZ1": tracks_xz.IsInTrack_clusters1_4p,
	            "isInTrack_4p_clZ2": tracks_xz.IsInTrack_clusters2_4p,
	            "isInTrack_4p_clZ3": tracks_xz.IsInTrack_clusters3_4p,
	            "isInTrack_4p_clY1": tracks_xy.IsInTrack_clusters1_4p,
	            "isInTrack_4p_clY2": tracks_xy.IsInTrack_clusters2_4p,
	            "isInTrack_4p_clY3": tracks_xy.IsInTrack_clusters3_4p,

                # BEST TRACKS
                "Theta_3p": theta_3p,
                "Theta_4p": theta_4p,
                "Phi_3p": phi_3p,
                "Phi_4p": phi_4p,
                "Phi_3p_Flipped": _mirror_azimuth(phi_3p),
                "Phi_4p_Flipped": _mirror_azimuth(phi_4p),
    
                "BestTrack_3p_xy_index": best_track_xy_ind,
                "BestTrack_3p_xz_index": best_track_xz_ind,
                "Track_3p_of_4p_index_xy": track_3p_of_4p_idx_xy,
                "Track_3p_of_4p_index_xz": track_3p_of_4p_idx_xz,
                "Track_4p_index_xy": track_4p_idx_xy,
                "Track_4p_index_xz": track_4p_idx_xz,

                "BestTrack_3p_ChiSquare_xy": best_track_3p_chi_xy,
                "BestTrack_3p_ChiSquare_xz": best_track_3p_chi_xz,
                "BestTrack_4p_ChiSquare_xy": best_track_4p_chi_xy,
                "BestTrack_4p_ChiSquare_xz": best_track_4p_chi_xz,
                "BestTracks_ScatteringAngle_xy": best_scattering_angle_xy,
                "BestTracks_ScatteringAngle_xz": best_scattering_angle_xz,
                "Best_track_4p_isTexpNULL_xy": best_track_4p_is_texp_null_xy,
                "Best_track_4p_isTexpNULL_xz": best_track_4p_is_texp_null_xz,

                "timestamp": timestamp,
                "WorkingPoint": int(working_point),
                "Temperature": temperature,
                "TriggerRate": trigger_rate,

                "TriggerMaskChannels": trigger_mask_channels,
                "TriggerMaskStrips": trigger_mask_strips,
                "TriggerMaskSize": trigger_mask_size,

                #only missing branches with respect to the C++ code are the "isBestTrack_3p_xy" and "isBestTrack_3p_xz" flags, but they are redundant since we have the indices of the best tracks, so we can easily derive them in the analysis if needed.

            }

            json_output_start = time.time()
            out_events.write(json.dumps(event_record) + "\n")
            json_event_output_elapsed += time.time() - json_output_start

            if progress_every > 0 and ev % progress_every == 0:
                elapsed = max(1e-9, time.time() - start_wall_time)
                rate = ev / elapsed
                print(
                    f"[progress] run={run} events={ev} elapsed={elapsed/60.0:.1f}m "
                    f"rate={rate:.1f} ev/s tracks3p_events={n_events_with_3p} "
                    f"tracks4p_events={n_events_with_4p}",
                    flush=True,
                )

    run_duration_min = 0.0
    if first_timestamp is not None and last_timestamp is not None:
        run_duration_min = max(0.0, (last_timestamp - first_timestamp) / 60.0)

    mean_cluster_energy_z1, rms_cluster_energy_z1 = _mean_rms(cluster_energy_acc["Z1"])
    mean_cluster_energy_z2, rms_cluster_energy_z2 = _mean_rms(cluster_energy_acc["Z2"])
    mean_cluster_energy_z3, rms_cluster_energy_z3 = _mean_rms(cluster_energy_acc["Z3"])
    mean_cluster_energy_z4, rms_cluster_energy_z4 = _mean_rms(cluster_energy_acc["Z4"])

    mean_cluster_energy_y1, rms_cluster_energy_y1 = _mean_rms(cluster_energy_acc["Y1"])
    mean_cluster_energy_y2, rms_cluster_energy_y2 = _mean_rms(cluster_energy_acc["Y2"])
    mean_cluster_energy_y3, rms_cluster_energy_y3 = _mean_rms(cluster_energy_acc["Y3"])
    mean_cluster_energy_y4, rms_cluster_energy_y4 = _mean_rms(cluster_energy_acc["Y4"])

    mean_nclusters_z1, rms_nclusters_z1 = _mean_rms(ncluster_acc["Z1"])
    mean_nclusters_z2, rms_nclusters_z2 = _mean_rms(ncluster_acc["Z2"])
    mean_nclusters_z3, rms_nclusters_z3 = _mean_rms(ncluster_acc["Z3"])
    mean_nclusters_z4, rms_nclusters_z4 = _mean_rms(ncluster_acc["Z4"])

    mean_nclusters_y1, rms_nclusters_y1 = _mean_rms(ncluster_acc["Y1"])
    mean_nclusters_y2, rms_nclusters_y2 = _mean_rms(ncluster_acc["Y2"])
    mean_nclusters_y3, rms_nclusters_y3 = _mean_rms(ncluster_acc["Y3"])
    mean_nclusters_y4, rms_nclusters_y4 = _mean_rms(ncluster_acc["Y4"])

    summary = {
        # Run-level aggregates: monitoring counters, mean/RMS observables, and config cuts.
        "Run": run,
        "Nevents": ev,
        "boards_OnePhes": boards_onephes,
        "IsOnePhe_aCopy": boards_is1phe_copy,
        "WorkingPoint": int(working_point),
        "Temperature": temperature,
        "TriggerRate": trigger_rate,
        "Nev_withAtrack3p": float(n_events_with_3p),
        "Nev_withAtrack3p_xy": float(n_events_with_3p_xy),
        "Nev_withAtrack3p_xz": float(n_events_with_3p_xz),
        "Nev_withAtrack4p": float(n_events_with_4p),
        "Ntracks_xy_3p": float(n_tracks_xy_3p),
        "Ntracks_xz_3p": float(n_tracks_xz_3p),
        "Ntracks_3p": float(n_tracks_3p),
        "Ntracks_xy_4p": float(n_tracks_xy_4p),
        "Ntracks_xz_4p": float(n_tracks_xz_4p),
        "Ntracks_4p": float(n_tracks_4p),
        "Nev_withAcluster_Z1": float(n_events_with_x1_cl),
        "Nev_withAcluster_Z2": float(n_events_with_x2_cl),
        "Nev_withAcluster_Z3": float(n_events_with_x3_cl),
        "Nev_withAcluster_Z4": float(n_events_with_x4_cl),
        "Nev_withAcluster_Y1": float(n_events_with_y1_cl),
        "Nev_withAcluster_Y2": float(n_events_with_y2_cl),
        "Nev_withAcluster_Y3": float(n_events_with_y3_cl),
        "Nev_withAcluster_Y4": float(n_events_with_y4_cl),
        "EnergyCut_clusterStrip": s1,
        "EnergyCut_singleStrip": s2,
        "EnergyCut_additionalStrip": s3,
        "proximity_cut_xz": proximity_cut_xz,
        "proximity_cut_xy": proximity_cut_xy,
        "RunDuration": run_duration_min,
        "Mean_clusterEnergy_Z1": mean_cluster_energy_z1,
        "Mean_clusterEnergy_Z2": mean_cluster_energy_z2,
        "Mean_clusterEnergy_Z3": mean_cluster_energy_z3,
        "Mean_clusterEnergy_Z4": mean_cluster_energy_z4,
        "RMS_clusterEnergy_Z1": rms_cluster_energy_z1,
        "RMS_clusterEnergy_Z2": rms_cluster_energy_z2,
        "RMS_clusterEnergy_Z3": rms_cluster_energy_z3,
        "RMS_clusterEnergy_Z4": rms_cluster_energy_z4,
        "Mean_clusterEnergy_Y1": mean_cluster_energy_y1,
        "Mean_clusterEnergy_Y2": mean_cluster_energy_y2,
        "Mean_clusterEnergy_Y3": mean_cluster_energy_y3,
        "Mean_clusterEnergy_Y4": mean_cluster_energy_y4,
        "RMS_clusterEnergy_Y1": rms_cluster_energy_y1,
        "RMS_clusterEnergy_Y2": rms_cluster_energy_y2,
        "RMS_clusterEnergy_Y3": rms_cluster_energy_y3,
        "RMS_clusterEnergy_Y4": rms_cluster_energy_y4,
        "Mean_Nclusters_Z1": mean_nclusters_z1,
        "RMS_Nclusters_Z1": rms_nclusters_z1,
        "Mean_Nclusters_Z2": mean_nclusters_z2,
        "RMS_Nclusters_Z2": rms_nclusters_z2,
        "Mean_Nclusters_Z3": mean_nclusters_z3,
        "RMS_Nclusters_Z3": rms_nclusters_z3,
        "Mean_Nclusters_Z4": mean_nclusters_z4,
        "RMS_Nclusters_Z4": rms_nclusters_z4,
        "Mean_Nclusters_Y1": mean_nclusters_y1,
        "RMS_Nclusters_Y1": rms_nclusters_y1,
        "Mean_Nclusters_Y2": mean_nclusters_y2,
        "RMS_Nclusters_Y2": rms_nclusters_y2,
        "Mean_Nclusters_Y3": mean_nclusters_y3,
        "RMS_Nclusters_Y3": rms_nclusters_y3,
        "Mean_Nclusters_Y4": mean_nclusters_y4,
        "RMS_Nclusters_Y4": rms_nclusters_y4,
        "NGoodTracks3p": n_good_tracks_3p,
        "NGoodTracks4p": n_good_tracks_4p,
        "analysis_file": str(analysis_jsonl),
    }

    event_loop_elapsed = max(1e-9, time.time() - start_wall_time)

    summary_output_start = time.time()
    with mini_summary_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    summary_output_elapsed = time.time() - summary_output_start

    analyzed_root = None
    run_info_root = None
    root_fill_elapsed = 0.0
    if write_root:
        root_start_wall_time = time.time()
        analyzed_root, run_info_root = write_root_outputs(
            analysis_jsonl=analysis_jsonl,
            mini_summary_json=mini_summary_json,
            config_files=[spiroc_cfg, telescope_cfg],
            reco_config_file=reco_config_file,
        )
        root_fill_elapsed = max(0.0, time.time() - root_start_wall_time)

    elapsed_total = max(1e-9, time.time() - total_start_wall_time)
    json_output_elapsed = json_event_output_elapsed + summary_output_elapsed
    output_production_elapsed = json_output_elapsed + root_fill_elapsed
    reconstruction_core_elapsed = max(0.0, elapsed_total - output_production_elapsed)

    summary["EventLoopWallTimeSeconds"] = event_loop_elapsed
    summary["JsonOutputWallTimeSeconds"] = json_output_elapsed
    summary["RootFillWallTimeSeconds"] = root_fill_elapsed
    summary["OutputProductionWallTimeSeconds"] = output_production_elapsed
    summary["ReconstructionCoreWallTimeSeconds"] = reconstruction_core_elapsed
    summary["ProcessingWallTimeSeconds"] = elapsed_total
    summary["ReconstructionRateEventsPerSecond"] = (
        ev / reconstruction_core_elapsed if ev > 0 and reconstruction_core_elapsed > 0 else 0.0
    )
    summary["ProcessingRateEventsPerSecond"] = ev / elapsed_total if ev > 0 else 0.0

    with file_handler.temp_to_output(mini_summary_json) as tmp_path:
        with tmp_path.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)

    print(f"Number of events: {ev}")
    print(f"[progress] Reconstruction core wall time: {reconstruction_core_elapsed:.2f} s", flush=True)
    print(f"[progress] JSON output wall time: {json_output_elapsed:.2f} s", flush=True)
    if write_root:
        print(f"[progress] ROOT fill wall time: {root_fill_elapsed:.2f} s", flush=True)
    print(f"[progress] Output production wall time: {output_production_elapsed:.2f} s", flush=True)
    print(f"[progress] Total wall time: {elapsed_total:.2f} s", flush=True)
    if ev > 0:
        print(f"[progress] Reconstruction core rate: {summary['ReconstructionRateEventsPerSecond']:.1f} ev/s", flush=True)
        print(f"[progress] Average processing rate: {ev / elapsed_total:.1f} ev/s", flush=True)
    print(f"RUN {run} COMPLETED. Find your data in ---> {analysis_jsonl}")
    print(f"Find your minitree in ---> {mini_summary_json}")
    if analyzed_root is not None and run_info_root is not None:
        print(f"Find your ROOT analyzed data in ---> {analyzed_root}")
        print(f"Find your ROOT minitree in ---> {run_info_root}")

    return analysis_jsonl, mini_summary_json


def _parse_bool_flag(value: str) -> bool:
    return str(value).strip().lower() == "true"


def _extract_run_number(path: Path) -> int:
    stem = path.stem
    if "run" not in stem:
        raise ValueError(f"Unable to extract run number from path: {path}")
    return int(stem.split("run")[-1])


def _configure_batch_logging(log_on_console: str, verbose: str, output_filename: Path) -> None:
    if log_on_console == "True":
        logging.basicConfig(
            level=getattr(logging, verbose.upper()),
            format="%(asctime)s [%(levelname)s] %(message)s",
        )
        return

    log_file = Path("logs/RECONSTRUCTED") / output_filename.with_suffix(".log").name
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text("", encoding="utf-8")
    logging.basicConfig(
        filename=log_file,
        level=getattr(logging, verbose.upper()),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def _build_pedestal_folder_map(pedestal_batch_file: Path) -> dict[int, Path]:
    with pedestal_batch_file.open("r", encoding="utf-8") as handle:
        pedestal_outputs = [Path(line.strip()) for line in handle if line.strip()]

    return {
        int(Path(pedestal_output).parent.stem): pedestal_output.parent
        for pedestal_output in pedestal_outputs
    }


def _process_batch_run(
    adc_file_str: str,
    color: str,
    pedestal_folders_by_run: dict[int, Path],
    fixed_pedestal_folder: Path | None,
    reconstructed_base_dir: Path,
    raw_base: Path,
    spiroc_mapping_file: Path,
    telescope_config_file: Path,
    slow_control_prefix: str,
    overwrite_outputs: bool,
    write_root: bool,
    progress_every: int,
    cluster_smearing_seed: int | None,
    config: dict,
    reco_config_file: Path | None = None,
) -> tuple[Path, int]:
    adc_file = Path(adc_file_str)
    run = _extract_run_number(adc_file)

    # A telescope-level pedestal folder (produced by 3_PedAna/PedestalReader.py) applies to
    # every run of that telescope, so it takes precedence over the per-run pedestal map.
    if fixed_pedestal_folder is not None:
        pedestal_folder = fixed_pedestal_folder
    else:
        pedestal_folder = pedestal_folders_by_run.get(run)
        if pedestal_folder is None:
            raise FileNotFoundError(
                f"No pedestal outputs found for run {run} in the provided pedestal batch file"
            )

    reconstructed_path = reconstructed_base_dir
    analysis_jsonl = reconstructed_path / f"MURAVES_AnalyzedData_run{run}.jsonl"
    analyzed_root = reconstructed_path / f"MURAVES_AnalyzedData_run{run}.root"
    primary_output = analyzed_root if write_root else analysis_jsonl

    should_process = overwrite_outputs or not primary_output.exists()
    if not should_process:
        logger.info("Output file %s already exists. Skipping run %s.", primary_output, run)
        return primary_output, run

    logger.info("Processing run %s from %s", run, adc_file)
    run_reconstruction(
        color=color,
        run=run,
        adc_file=adc_file,
        pedestal_folder=pedestal_folder,
        reconstructed_path=reconstructed_path,
        slow_control_file=raw_base / color / f"{slow_control_prefix}{run}",
        spiroc_mapping_file=spiroc_mapping_file,
        telescope_config_file=telescope_config_file,
        write_root=write_root,
        progress_every=progress_every,
        cluster_smearing_seed=cluster_smearing_seed,
        config=config,
        reco_config_file=reco_config_file,
    )
    return primary_output, run


def parse_args(config_defaults: dict | None = None) -> argparse.Namespace:
    if config_defaults is None:
        config_defaults = get_reco_config()
    parser = argparse.ArgumentParser(
        description=(
            "Batch MURAVES reconstruction entrypoint for Snakemake. "
            "Consumes batch stamp files and reconstructs all runs listed in the batch."
        )
    )
    parser.add_argument("color", type=str, help="Detector color: ROSSO, NERO, or BLU")
    parser.add_argument(
        "-i",
        "--input-filename",
        dest="input_filename",
        type=Path,
        required=True,
        help="Batch stamp file listing parsed ADC input files.",
    )
    parser.add_argument(
        "-p",
        "--pedestal-input-filename",
        dest="pedestal_input_filename",
        type=Path,
        default=None,
        help="Batch stamp file listing per-run pedestal outputs for the same runs. "
             "Mutually exclusive with --pedestal-folder.",
    )
    parser.add_argument(
        "--pedestal-folder",
        dest="pedestal_folder",
        type=Path,
        default=None,
        help="Single telescope-level pedestal folder (containing pedestal_N files) applied to "
             "every run in the batch, as produced by 3_PedAna/PedestalReader.py. "
             "Mutually exclusive with -p/--pedestal-input-filename.",
    )
    parser.add_argument(
        "-o",
        "--output-filename",
        dest="output_filename",
        type=Path,
        required=True,
        help="Batch stamp file to write with reconstructed outputs.",
    )
    parser.add_argument(
        "--raw-base",
        type=Path,
        required=True,
        help="Base RAW folder that contains <COLOR>/SLOWCONTROL_run<run>",
    )
    parser.add_argument(
        "--tracks-base",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="tracks_reconstruction base path (contains AncillaryFiles)",
    )
    parser.add_argument(
        "--spiroc-mapping-file",
        type=Path,
        required=True,
        help="Path to SPiROC strip-channel mapping file.",
    )
    parser.add_argument(
        "--telescope-config-file",
        type=Path,
        required=True,
        help="Path to telescope board configuration file for the selected color.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional JSON override for reconstruction parameters.",
    )
    parser.add_argument(
        "--base-config",
        type=Path,
        default=None,
        help="Optional base JSON configuration for reconstruction parameters.",
    )
    parser.add_argument(
        "--write-root",
        action="store_true",
        help="Also export ROOT files from the generated JSON outputs.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=int(config_defaults["reconstruction"]["progress_every_default"]),
        help="Print progress every N events within each run (set 0 to disable).",
    )
    parser.add_argument(
        "--slow-control-prefix",
        dest="slow_control_prefix",
        type=str,
        default="SLOWCONTROL_run",
        help="Filename prefix for slow-control files inside <raw-base>/<color>/. Default: SLOWCONTROL_run",
    )
    parser.add_argument(
        "--cluster-smearing-seed",
        type=int,
        default=None,
        help=(
            "Optional deterministic seed for single-strip cluster smearing. "
            "Use the same value in C++ and Python to align the smearing sequence."
        ),
    )
    parser.add_argument(
        "-l",
        "--log_on_console",
        dest="log_on_console",
        required=True,
        help="If true logs are printed on terminal, if False they are printed on a file.",
    )
    parser.add_argument(
        "-ow",
        "--overwrite_outputs",
        dest="overwrite_outputs",
        required=True,
        help="If true, existing output files will be overwritten. If False, existing output files will be kept and the corresponding runs will be skipped.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        required=False,
        default="info",
        help="Logging level: debug/info/warning/error/critical (default = info)",
    )
    parser.add_argument(
        "-th",
        "--num_threads",
        dest="num_threads",
        type=int,
        default=1,
        help="Number of threads/cores to use for processing (default = 1)",
    )
    parsed = parser.parse_args()
    if (parsed.pedestal_input_filename is None) == (parsed.pedestal_folder is None):
        parser.error(
            "exactly one of -p/--pedestal-input-filename or --pedestal-folder is required"
        )
    return parsed


def main() -> None:
    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument("--config", type=Path, default=None)
    pre_parser.add_argument("--base-config", type=Path, default=None)
    pre_args, _ = pre_parser.parse_known_args()

    preloaded_config = get_reco_config(pre_args.config, pre_args.base_config)
    args = parse_args(preloaded_config)

    if args.config is not None or args.base_config is not None:
        set_runtime_config_path(args.config, args.base_config)
    config = get_reco_config(args.config, args.base_config)
    color = args.color.upper()
    overwrite_outputs = _parse_bool_flag(args.overwrite_outputs)

    _configure_batch_logging(args.log_on_console, args.verbose, args.output_filename)

    with args.input_filename.open("r", encoding="utf-8") as handle:
        adc_files = [line.strip() for line in handle if line.strip()]

    if args.pedestal_folder is not None:
        pedestal_folders_by_run: dict[int, Path] = {}
        logger.info("Using telescope-level pedestal folder %s for all runs", args.pedestal_folder)
    else:
        pedestal_folders_by_run = _build_pedestal_folder_map(args.pedestal_input_filename)

    if not adc_files:
        with file_handler.temp_to_output(args.output_filename) as tmp:
            Path(tmp).write_text("", encoding="utf-8")
        logger.info("No ADC files found in batch input %s", args.input_filename)
        return

    worker = partial(
        _process_batch_run,
        color=color,
        pedestal_folders_by_run=pedestal_folders_by_run,
        fixed_pedestal_folder=args.pedestal_folder,
        reconstructed_base_dir=args.output_filename.parent,
        raw_base=args.raw_base,
        spiroc_mapping_file=args.spiroc_mapping_file,
        telescope_config_file=args.telescope_config_file,
        slow_control_prefix=args.slow_control_prefix,
        overwrite_outputs=overwrite_outputs,
        write_root=args.write_root,
        progress_every=args.progress_every,
        cluster_smearing_seed=args.cluster_smearing_seed,
        config=config,
        reco_config_file=args.base_config or args.config,
    )

    with Pool(args.num_threads) as pool:
        results = pool.map(worker, adc_files)

    output_files, run_list = zip(*results)
    with file_handler.temp_to_output(args.output_filename) as tmp:
        with Path(tmp).open("w", encoding="utf-8") as handle:
            for output_file in output_files:
                handle.write(f"{output_file}\n")

    print(
        f"Batch {args.output_filename.stem} completed with {args.num_threads} thread(s). Runs: {run_list}",
        flush=True,
    )


if __name__ == "__main__":
    main()
