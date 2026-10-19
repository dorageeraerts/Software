__all__ = ["save_root_with_metadata", "load_root_metadata", "add_metadata_to_root"]

import ROOT
import hashlib
import json
import os
import subprocess
from datetime import datetime
import re


def git(cmd):
    return subprocess.check_output(cmd).decode().strip()

def get_git_metadata():
    """Collect git provenance in as few subprocess calls as possible."""

    # One call: get all tags with their tagger/commit date, sorted newest first.
    # Format: "<refname:short> <version:refname:short>"
    try:
        raw = subprocess.check_output(
            ["git", "tag", "--list", "--sort=-version:refname"],
            stderr=subprocess.DEVNULL,
        ).decode().splitlines()
        all_tags = [t.strip() for t in raw if t.strip()]
    except subprocess.CalledProcessError:
        all_tags = []

    # Split into algorithm tags (contain "-") and semver-only tags (no "-").
    algo_tags = [t for t in all_tags if "-" in t]

    # Build algorithms dict: algo_id -> latest tag for that algo.
    algorithms: dict = {}
    for tag in algo_tags:
        algo_id = tag.split("-")[0]
        if algo_id not in algorithms:          # list is sorted newest-first
            algorithms[algo_id] = tag

    # Latest release = first tag in the full sorted list (semver-aware sort).
    latest_release = all_tags[0] if all_tags else None

    # One more call: current short commit hash.
    try:
        commit = git(["git", "rev-parse", "--short", "HEAD"])
    except subprocess.CalledProcessError:
        commit = None

    return {
        "algorithms": algorithms,
        "latest_release": latest_release,
        "generated_at": datetime.utcnow().isoformat(),
        "commit": commit,
    }


def _sha256sum(path, chunk_size=1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _parse_cfg_table(path):
    """Parse a whitespace-separated table with a header row into a list of dicts."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = [l.rstrip("\n") for l in handle]
    if not lines:
        return []
    headers = lines[0].split()
    rows = []
    for line in lines[1:]:
        fields = line.split()
        if not fields:
            continue
        rows.append(dict(zip(headers, fields)))
    return rows


def _build_config_file_metadata(
    config_files,
    include_contents=False,
    max_inline_bytes=100_000,
):
    config_metadata = []
    for cfg_file in config_files:
        abs_path = os.path.abspath(cfg_file)
        entry = {
            "path": abs_path,
            "exists": os.path.isfile(abs_path),
        }

        if not entry["exists"]:
            config_metadata.append(entry)
            continue

        size_bytes = os.path.getsize(abs_path)
        mtime = datetime.utcfromtimestamp(os.path.getmtime(abs_path)).isoformat()
        entry.update(
            {
                "size_bytes": size_bytes,
                "modified_at": mtime,
                "sha256": _sha256sum(abs_path),
            }
        )

        ext = os.path.splitext(abs_path)[1].lower()
        if ext == ".json":
            try:
                with open(abs_path, "r", encoding="utf-8") as handle:
                    entry["parsed_content"] = json.load(handle)
            except Exception as exc:
                entry["parse_error"] = str(exc)
        elif ext == ".cfg":
            try:
                entry["parsed_content"] = _parse_cfg_table(abs_path)
            except Exception as exc:
                entry["parse_error"] = str(exc)

        config_metadata.append(entry)

    return config_metadata


def save_root_with_metadata(
    path,
    *root_objects,
    config_files=None,
    include_config_contents=False,
    max_inline_config_bytes=100_000,
):
    metadata = get_git_metadata()

    if not config_files:
        raise ValueError("config_files is required and must contain at least one file path")

    metadata["detector_config_files"] = _build_config_file_metadata(
        config_files,
        include_contents=include_config_contents,
        max_inline_bytes=max_inline_config_bytes,
    )

    f = ROOT.TFile(path, "RECREATE")

    # salva gli oggetti ROOT (istogrammi, tree, ecc.)
    for obj in root_objects:
        obj.Write()

    # salva i metadati come TObjString
    meta_str = ROOT.TObjString(json.dumps(metadata))
    meta_str.Write("METADATA")

    f.Close()
    return metadata

def add_metadata_to_root(
    path,
    config_files,
    include_config_contents=False,
    max_inline_config_bytes=100_000,
    prebuilt_metadata=None,
):
    """Append METADATA to an existing ROOT file without overwriting its contents.

    If *prebuilt_metadata* is provided it is written as-is (no git/config
    queries), which avoids redundant work when writing multiple ROOT files
    in the same run.
    """
    if prebuilt_metadata is not None:
        metadata = prebuilt_metadata
    else:
        if not config_files:
            raise ValueError("config_files is required and must contain at least one file path")
        metadata = get_git_metadata()
        metadata["detector_config_files"] = _build_config_file_metadata(
            config_files,
            include_contents=include_config_contents,
            max_inline_bytes=max_inline_config_bytes,
        )

    f = ROOT.TFile(str(path), "UPDATE")
    meta_str = ROOT.TObjString(json.dumps(metadata))
    meta_str.Write("METADATA")
    f.Close()
    return metadata


def load_root_metadata(path):
    f = ROOT.TFile(str(path))
    meta_obj = f.Get("METADATA")
    if not meta_obj:
        f.Close()
        return None
    metadata = json.loads(meta_obj.GetString().Data())
    f.Close()
    return metadata


"""
# EXAMPLE: Save ROOT with metadata
import ROOT
from root_wrapper import save_root_with_metadata

# creo un istogramma dummy
h = ROOT.TH1F("h", "Example", 100, 0, 1)
h.FillRandom("gaus", 1000)

# genero sample con metadati
meta = save_root_with_metadata("sample.root", "algoCLUSRECO", h)
print("Metadati salvati:", meta)
"""

"""
# EXEMPLE: Read ROOT metadata
f = ROOT.TFile("sample.root")
meta_obj = f.Get("METADATA")
metadata = json.loads(meta_obj.GetString().Data())
print(metadata["algorithms"])
"""