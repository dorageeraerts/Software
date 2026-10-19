"""Print METADATA stored inside one or more MURAVES ROOT output files.

Usage
-----
    python print_metadata.py path/to/file.root [path/to/other.root ...]
    python print_metadata.py path/to/file.root --json          # raw JSON dump
    python print_metadata.py path/to/file.root --config-only   # only detector config files
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from muraves_lib import root_wrapper


def _print_table(rows: list[dict], indent: str = "      ") -> None:
    if not rows:
        return
    headers = list(rows[0].keys())
    col_widths = {h: max(len(h), max(len(str(r.get(h, ""))) for r in rows)) for h in headers}
    header_line = "  ".join(h.ljust(col_widths[h]) for h in headers)
    separator = "  ".join("-" * col_widths[h] for h in headers)
    print(indent + header_line)
    print(indent + separator)
    for row in rows:
        print(indent + "  ".join(str(row.get(h, "")).ljust(col_widths[h]) for h in headers))


def _print_metadata(meta: dict, json_mode: bool, config_only: bool) -> None:
    if json_mode:
        print(json.dumps(meta, indent=2))
        return

    if not config_only:
        print(f"  commit        : {meta.get('commit', 'N/A')}")
        print(f"  latest_release: {meta.get('latest_release', 'N/A')}")
        print(f"  generated_at  : {meta.get('generated_at', 'N/A')}")
        algorithms = meta.get("algorithms", {})
        if algorithms:
            print(f"  algorithms    :")
            for algo, tag in algorithms.items():
                print(f"    {algo:>6} -> {tag}")
        else:
            print("  algorithms    : (no tags found)")

    cfg_files = meta.get("detector_config_files", [])
    print(f"  detector_config_files ({len(cfg_files)} file{'s' if len(cfg_files) != 1 else ''}):")
    for cfg in cfg_files:
        status = "OK" if cfg.get("exists", True) else "MISSING"
        sha = cfg.get("sha256", "N/A")
        sha_short = sha[:16] + "..." if len(sha) > 16 else sha
        size = cfg.get("size_bytes")
        size_str = f"  {size} B" if size is not None else ""
        modified = cfg.get("modified_at", "")
        modified_str = f"  modified={modified}" if modified else ""
        print(f"    [{status}] {cfg['path']}")
        print(f"           sha256={sha_short}{size_str}{modified_str}")
        if "parse_error" in cfg:
            print(f"           parse_error: {cfg['parse_error']}")
        elif "parsed_content" in cfg:
            content = cfg["parsed_content"]
            if isinstance(content, list) and content:
                _print_table(content, indent="           ")
            elif isinstance(content, dict):
                print(json.dumps(content, indent=4).replace("\n", "\n           "))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print METADATA from MURAVES ROOT output files."
    )
    parser.add_argument(
        "root_files",
        nargs="+",
        type=Path,
        metavar="FILE.root",
        help="One or more ROOT files to inspect.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Dump raw metadata as JSON instead of formatted output.",
    )
    parser.add_argument(
        "--config-only",
        action="store_true",
        help="Print only the detector_config_files section.",
    )
    args = parser.parse_args()

    any_error = False
    for path in args.root_files:
        print(f"\n=== {path} ===")
        if not path.exists():
            print("  ERROR: file not found")
            any_error = True
            continue

        meta = root_wrapper.load_root_metadata(path)
        if meta is None:
            print("  WARNING: no METADATA key found in this file")
            continue

        _print_metadata(meta, json_mode=args.json, config_only=args.config_only)

    sys.exit(1 if any_error else 0)


if __name__ == "__main__":
    main()
