# 4_MainRecSof Pythonized Modules

This folder now includes Python ports of the core C++ reconstruction modules.

## Added Python Files

- `cluster_lists.py` (`CreateClusterList`, `ClusterPosition`, `ClusterEnergy`, `SortIndices`)
- `evaluate_angular_coordinates.py` (`TrackAngularCoordinates`)
- `read_event.py` (`ReadEvent`)
- `tracking.py` (`MakeTracks`)
- `MURAVES_reco_v2.py` (Python entrypoint)

Compatibility wrappers with legacy names are also provided:

- `ClusterLists.py`
- `EvaluateAngularCoordinates.py`
- `ReadEvent.py`
- `Tracking.py`

## Entrypoint

```bash
python 4_MainRecSof/MURAVES_reco_v2.py <COLOR> <RUN>
```

Example:

```bash
python 4_MainRecSof/MURAVES_reco_v2.py BLU 2512
```

Optional path overrides:

```bash
python 4_MainRecSof/MURAVES_reco_v2.py BLU 2512 \
  --output-base /workspace/muraves_outputs \
  --raw-base /data/RAW_GZ \
  --tracks-base /workspace/Software/muraves/tracks_reconstruction
```

Optional ROOT export (on demand):

```bash
python 4_MainRecSof/MURAVES_reco_v2.py BLU 2512 --write-root
```

Progress logging (prints every 500 events):

```bash
python 4_MainRecSof/MURAVES_reco_v2.py NERO 2546 --write-root --progress-every 500
```

Disable periodic progress logs:

```bash
python 4_MainRecSof/MURAVES_reco_v2.py NERO 2546 --progress-every 0
```

Programmatic on-demand export function:

```python
from pathlib import Path
from MURAVES_reco_v2 import write_root_outputs

write_root_outputs(
    analysis_jsonl=Path("/path/to/MURAVES_AnalyzedData_run2512.jsonl"),
    mini_summary_json=Path("/path/to/MURAVES_miniRunTree_run2512.json"),
)
```

## Outputs

The Python entrypoint writes:

- Event-level file: `MURAVES_AnalyzedData_run<run>.jsonl`
- Run summary file: `MURAVES_miniRunTree_run<run>.json`

If `--write-root` is used (or `write_root_outputs(...)` is called), it also writes:

- Event ROOT tree: `MURAVES_AnalyzedData_run<run>.root` (`AnalyzedData`)
- Run ROOT tree: `MURAVES_miniRunTree_run<run>.root` (`Run_info`)

Note on nested list branches:

- Some fields that are nested lists (for example `list[list[float]]`) are stored as flattened values in the original branch name, plus a companion branch `<name>__counts` to reconstruct the inner grouping.

Both files are created in:

- `<output-base>/RECONSTRUCTED/<COLOR>/`

## Notes

- The C++ sources and executable are untouched.
- This allows incremental migration: existing C++ workflows can continue while Python modules are validated.
