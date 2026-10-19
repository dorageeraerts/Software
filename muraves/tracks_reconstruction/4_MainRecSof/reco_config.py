from __future__ import annotations

import functools
import json
import os
from copy import deepcopy
from pathlib import Path
from typing import Any

RUNTIME_CONFIG_ENV = "MURAVES_RECO_CONFIG"
RUNTIME_BASE_CONFIG_ENV = "MURAVES_RECO_BASE_CONFIG"


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Configuration file {path} must contain a JSON object")
    return data


def set_runtime_config_path(
    config_path: Path | str | None,
    base_config_path: Path | str | None = None,
) -> None:
    if config_path is None:
        os.environ.pop(RUNTIME_CONFIG_ENV, None)
    else:
        os.environ[RUNTIME_CONFIG_ENV] = str(Path(config_path).expanduser().resolve())

    if base_config_path is None:
        os.environ.pop(RUNTIME_BASE_CONFIG_ENV, None)
    else:
        os.environ[RUNTIME_BASE_CONFIG_ENV] = str(Path(base_config_path).expanduser().resolve())

    # Clear the cache when config paths change
    _load_reco_config_cached.cache_clear()


@functools.lru_cache(maxsize=1)
def _load_reco_config_cached(
    config_path: str | None,
    base_config_path: str | None,
) -> dict[str, Any]:
    """Internal cached version of config loading. Paths must be strings for hashability."""
    base_candidate: Path | None
    if base_config_path is not None:
        base_candidate = Path(base_config_path).expanduser().resolve()
    else:
        configured_base = os.environ.get(RUNTIME_BASE_CONFIG_ENV)
        base_candidate = Path(configured_base).expanduser().resolve() if configured_base else None

    candidate: Path | None
    if config_path is not None:
        candidate = Path(config_path).expanduser().resolve()
    else:
        configured = os.environ.get(RUNTIME_CONFIG_ENV)
        candidate = Path(configured).expanduser().resolve() if configured else None

    if base_candidate is None and candidate is None:
        raise ValueError(
            "No reconstruction configuration file provided. "
            "Pass --config in the entrypoint or set runtime config environment variables."
        )

    if base_candidate is None:
        if candidate is None:
            raise ValueError("No reconstruction configuration file provided")
        if not candidate.exists():
            raise FileNotFoundError(f"Configuration file not found: {candidate}")
        return _load_json(candidate)

    if not base_candidate.exists():
        raise FileNotFoundError(f"Base configuration file not found: {base_candidate}")

    base = _load_json(base_candidate)

    if candidate is None or candidate == base_candidate:
        return base

    if not candidate.exists():
        raise FileNotFoundError(f"Configuration override file not found: {candidate}")

    override = _load_json(candidate)
    return _deep_merge(base, override)


def get_reco_config(
    config_path: Path | str | None = None,
    base_config_path: Path | str | None = None,
) -> dict[str, Any]:
    """Load reconstruction configuration with caching.
    
    Configuration is cached after first load to avoid repeated file I/O.
    Call set_runtime_config_path() with new paths to invalidate the cache.
    """
    # Convert to strings for hashability in lru_cache
    config_str = str(Path(config_path).expanduser().resolve()) if config_path else None
    base_str = str(Path(base_config_path).expanduser().resolve()) if base_config_path else None
    
    return _load_reco_config_cached(config_str, base_str)


def resolve_first_existing(paths: list[str]) -> Path:
    for item in paths:
        candidate = Path(item)
        if candidate.exists():
            return candidate
    return Path(paths[-1])
