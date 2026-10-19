from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence
from pathlib import Path


@dataclass
class PedestalStruct:
    """Store pedestal and 1phe data for a board, reordered by channel mapping."""
    channel_numbers: list[int] = field(default_factory=list)
    pedestal_values: list[float] = field(default_factory=list)
    onephe_values: list[float] = field(default_factory=list)
    flags: list[int] = field(default_factory=list)


def _safe_float(value: str, default: float = 0.0) -> float:
    """Safely convert string to float, returning default on error."""
    try:
        val = float(value)
        # Check for NaN
        if math.isnan(val):
            return default
        return val
    except (TypeError, ValueError):
        return default


def _parse_value(value: str) -> float | None:
    """Parse a string as float; return None if NaN or invalid."""
    try:
        val = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(val) else val


def read_pedestal_file(
    filepath: str | Path,
    stripIndices: Sequence[int],
    board_idx: int = 0,
) -> PedestalStruct:
    """
    Read a pedestal/1phe file, keeping channels in the same order as in the input file.

    No channel remapping is applied here: OnePhE_evaluator.py's write_results (which
    produces the reference pedestal_N files) writes channels in plain raw ADC order
    (0..n_channels-1), so this function preserves that same order rather than
    reordering through the spiroc-hybrid-map.cfg strip mapping. This beacuse the re-ordering
    is done at reconstruction level.

    Args:
        filepath:     Path to the txt file with 3 columns: channel_number, pedestal, 1phe
        stripIndices: Only its length is used, as the number of channels per board.
        board_idx:    Board index (0-based), used to compute the per-board channel offset.
                     Input files use global channel numbering (board 0: 1-32, board 1: 33-64, ...),
                     so the offset board_idx * n_channels is subtracted before indexing.

    Returns:
        PedestalStruct containing channel data in raw file order:
        - channel_numbers: local 0-based channel indices
        - pedestal_values: pedestal values per channel
        - onephe_values: 1phe values per channel
        - flags: 1 if pedestal or 1phe was NaN, or if 1phe was 0, else 0

    Notes:
        - NaN pedestal/1phe values are replaced with the average of the other
          (non-NaN) values for that column found in the same input file
        - A 1phe value of 0 is replaced with 1000 and flagged (flag=1), since it
          usually indicates a failed calibration rather than a genuine reading
        - Channels not present in the file are padded with 0.0
    """
    filepath = Path(filepath)

    # Initialize raw data arrays matching stripIndices length
    raw_channels = [0] * len(stripIndices)
    raw_pedestals = [0.0] * len(stripIndices)
    raw_onephes = [0.0] * len(stripIndices)
    raw_flags = [0] * len(stripIndices)
    present = [False] * len(stripIndices)
    pedestal_is_nan = [False] * len(stripIndices)
    onephe_is_nan = [False] * len(stripIndices)

    # Read file and populate raw data (channels are 1-indexed in input)
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            # Input channels are globally numbered: board 0 → 1-32, board 1 → 33-64, etc.
            # Convert to local 0-based index by subtracting the board offset.
            n_channels = len(stripIndices)
            channel = int(_safe_float(parts[0], default=-1.0)) - 1 - board_idx * n_channels
            if channel < 0 or channel >= len(stripIndices):
                continue

            pedestal = _parse_value(parts[1])
            onephe = _parse_value(parts[2])

            raw_channels[channel] = channel
            present[channel] = True
            if pedestal is None:
                pedestal_is_nan[channel] = True
            else:
                raw_pedestals[channel] = pedestal
            if onephe is None:
                onephe_is_nan[channel] = True
            elif onephe == 0.0:
                raw_onephes[channel] = 1000.0
            else:
                raw_onephes[channel] = onephe
            raw_flags[channel] = 1 if (pedestal is None or onephe is None or onephe == 0.0) else 0

    # Impute NaN entries with the average of the other valid values in this file
    valid_pedestals = [
        v for v, is_present, is_nan in zip(raw_pedestals, present, pedestal_is_nan)
        if is_present and not is_nan
    ]
    valid_onephes = [
        v for v, is_present, is_nan in zip(raw_onephes, present, onephe_is_nan)
        if is_present and not is_nan
    ]
    pedestal_avg = sum(valid_pedestals) / len(valid_pedestals) if valid_pedestals else 0.0
    onephe_avg = sum(valid_onephes) / len(valid_onephes) if valid_onephes else 0.0
    for channel in range(len(stripIndices)):
        if pedestal_is_nan[channel]:
            raw_pedestals[channel] = pedestal_avg
        if onephe_is_nan[channel]:
            raw_onephes[channel] = onephe_avg

    # No reordering: keep raw file order, consistent with OnePhE_evaluator.py's write_results
    return PedestalStruct(
        channel_numbers=raw_channels,
        pedestal_values=raw_pedestals,
        onephe_values=raw_onephes,
        flags=raw_flags,
    )


def write_pedestal_file(
    output_filepath: str | Path,
    pedestal_data: PedestalStruct,
    board_idx: int = 0,
) -> None:
    """
    Write pedestal data to file in the format produced by OnePhE_evaluator.write_results.

    Values that were NaN in the source file (flag==1) have already been replaced
    with the average of the other valid values in that file by read_pedestal_file;
    they are written here like any other value, with flag=1 preserved for traceability.

    Args:
        output_filepath: Path to output file
        pedestal_data:   PedestalStruct to write
        board_idx:       Board index (0-based, reserved for future use)
    """
    output_filepath = Path(output_filepath)
    output_filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(output_filepath, 'w') as f:
        f.write("ch \t ped \t 1pe \n")

        for ch_idx, (ped, onephe, flag) in enumerate(zip(
            pedestal_data.pedestal_values,
            pedestal_data.onephe_values,
            pedestal_data.flags,
        )):
            f.write(f"{ch_idx}\t{round(ped)}\t{round(onephe)}\t {flag}\n")


# Commented out for now: this plotting workflow assumed a run-numbered output layout.
# It will be reworked for the telescope-based workflow later.
#
# def _read_output_pedestal_file(filepath: str | Path) -> tuple[list[int], list[float], list[int]]:
#     """Read a pedestal_N file (as written by write_pedestal_file) and return (channels, onephe_values, flags)."""
#     filepath = Path(filepath)
#     channels: list[int] = []
#     onephes: list[float] = []
#     flags: list[int] = []
#
#     with filepath.open('r') as f:
#         next(f, None)  # skip header line
#         for line in f:
#             parts = line.split()
#             if len(parts) < 3:
#                 continue
#             channels.append(int(parts[0]))
#             onephes.append(float(parts[2]))
#             flags.append(int(parts[3]) if len(parts) > 3 else 0)
#
#     return channels, onephes, flags
#
#
# def plot_onephe_comparison(
#     generated_dir: str | Path,
#     reference_dir: str | Path,
#     output_path: str | Path | None = None,
#     num_boards: int = 16,
#     channels_per_board: int = 32,
#     file_pattern: str = "pedestal_{}",
#     show: bool = True,
# ) -> None:
#     """
#     Plot a comparison of 1phe values between two sets of pedestal_N output files
#     (in the format written by write_pedestal_file), laid out along a shared x axis
#     of board + channel (global index = board_idx * channels_per_board + channel).
#
#     Args:
#         generated_dir:      Directory with this script's pedestal_N output files
#         reference_dir:      Directory with the reference pedestal_N files to compare against
#         output_path:        If given, save the figure to this path
#         num_boards:         Number of boards to compare (default 16)
#         channels_per_board: Channels per board, used to lay out the x axis (default 32)
#         file_pattern:       Filename pattern for both directories (use {} for board number, 0-indexed)
#         show:                Whether to display the plot interactively
#     """
#     import matplotlib.patches as mpatches
#     import matplotlib.pyplot as plt
#     from matplotlib.lines import Line2D
#
#     generated_dir = Path(generated_dir)
#     reference_dir = Path(reference_dir)
#
#     gen_x: list[int] = []
#     gen_y: list[float] = []
#     ref_x: list[int] = []
#     ref_y: list[float] = []
#     board_starts: list[float] = []
#     ref_flagged_x: set[int] = set()
#     gen_flagged_x: set[int] = set()
#
#     for board_idx in range(num_boards):
#         gen_file = generated_dir / file_pattern.format(board_idx)
#         ref_file = reference_dir / file_pattern.format(board_idx)
#
#         if gen_file.exists():
#             channels, onephes, flags = _read_output_pedestal_file(gen_file)
#             xs = [board_idx * channels_per_board + ch for ch in channels]
#             gen_x.extend(xs)
#             gen_y.extend(onephes)
#             gen_flagged_x.update(x for x, flag in zip(xs, flags) if flag == 1)
#         else:
#             print(f"Warning: {gen_file} not found, skipping board {board_idx} (generated)")
#
#         if ref_file.exists():
#             channels, onephes, flags = _read_output_pedestal_file(ref_file)
#             xs = [board_idx * channels_per_board + ch for ch in channels]
#             ref_x.extend(xs)
#             ref_y.extend(onephes)
#             ref_flagged_x.update(x for x, flag in zip(xs, flags) if flag == 1)
#         else:
#             print(f"Warning: {ref_file} not found, skipping board {board_idx} (reference)")
#
#         board_starts.append(board_idx * channels_per_board)
#
#     fig, ax = plt.subplots(figsize=(16, 6))
#
#     for x in ref_flagged_x:
#         ax.axvspan(x - 0.5, x + 0.5, color='lightgrey', alpha=0.6, zorder=0)
#     for x in gen_flagged_x:
#         ax.axvline(x, color='red', linestyle='--', linewidth=1, alpha=0.8, zorder=2)
#
#     ax.plot(ref_x, ref_y, 'o', label=f"reference ({reference_dir.name})",
#             color='tab:gray', alpha=0.7, markersize=4)
#     ax.plot(gen_x, gen_y, 'x', label=f"generated ({generated_dir.name})",
#             color='tab:red', markersize=5)
#
#     y_top = 50
#     for board_idx, x_start in enumerate(board_starts):
#         ax.axvline(x_start - 0.5, color='black', linewidth=0.5, linestyle='--', alpha=0.4)
#         ax.text(x_start, y_top * 0.99, f"Board {board_idx}", rotation=90,
#                 va='top', ha='left', fontsize=8, color='black')
#
#     ax.set_xticks([])
#     ax.set_xlabel("Board / Channel")
#     ax.set_ylabel("1phe value")
#     ax.set_ylim(0, y_top)
#     ax.set_title("1phe comparison: generated vs reference")
#
#     handles, labels = ax.get_legend_handles_labels()
#     if ref_flagged_x:
#         handles.append(mpatches.Patch(color='lightgrey', alpha=0.6))
#         labels.append('flagged in reference')
#     if gen_flagged_x:
#         handles.append(Line2D([0], [0], color='red', linestyle='--', linewidth=1))
#         labels.append('flagged in generated')
#     ax.legend(handles, labels)
#     fig.tight_layout()
#
#     if output_path is not None:
#         output_path = Path(output_path)
#         output_path.parent.mkdir(parents=True, exist_ok=True)
#         fig.savefig(output_path, dpi=150)
#         print(f"Saved comparison plot to {output_path}")
#
#     if show:
#         plt.show()
#     plt.close(fig)


def process_pedestal_files(
    input_dir: str | Path,
    output_base: str | Path,
    stripIndices: Sequence[int],
    telescope: str,
    input_pattern: str = "ped_onephe_sk{}.txt",
    output_pattern: str = "pedestal_{}",
    num_boards: int = 16,
) -> list[Path]:
    """
    Process multiple pedestal files from input directory and write to output directory.
    Output files are placed in output_base/telescope/, since all runs belonging to the
    same telescope (NERO, BLU, ROSSO) share the same 1phe calibration values.

    Args:
        input_dir:     Directory containing input files
        output_base:   Base output directory (telescope sub-folder is created inside)
        stripIndices:  Only its length is used, as the number of channels per board (same for all boards)
        telescope:     Telescope name used as sub-folder, e.g. "NERO"
        input_pattern: Pattern for input filenames (use {} for board number, 1-indexed)
        output_pattern:Pattern for output filenames (use {} for board number, 0-indexed)
        num_boards:    Number of boards to process (default 16)

    Returns:
        List of the pedestal files actually written.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_base) / telescope
    output_dir.mkdir(parents=True, exist_ok=True)
    written_files: list[Path] = []

    for board_idx in range(num_boards):
        # Input files are 1-indexed (sk1, sk2, ..., sk16)
        input_file = input_dir / input_pattern.format(board_idx + 1)

        # Output files are 0-indexed (pedestal_0, pedestal_1, ..., pedestal_15)
        output_file = output_dir / output_pattern.format(board_idx)

        if input_file.exists():
            print(f"Processing {input_file.name} -> {output_file}")
            ped_data = read_pedestal_file(input_file, stripIndices, board_idx=board_idx)
            write_pedestal_file(output_file, ped_data, board_idx=board_idx)
            written_files.append(output_file)
        else:
            print(f"Warning: {input_file} not found, skipping board {board_idx}")

    return written_files


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Read pedestal/1phe files, keeping channels in raw file order (no strip remapping)"
    )
    parser.add_argument(
        "--input-dir",
        default=str(Path(__file__).parent / "NERO_PED_BANCHMARK"),
        help="Directory containing ped_onephe_sk*.txt files (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).parent),
        help="Base output directory; a sub-folder named after --telescope is created inside (default: %(default)s)",
    )
    parser.add_argument(
        "--telescope",
        choices=["NERO", "BLU", "ROSSO"],
        default=None,
        help="Telescope name; all runs of this telescope share the same 1phe values. "
             "Defaults to the prefix of --input-dir's name (e.g. 'NERO_PED_BANCHMARK' -> 'NERO').",
    )
    parser.add_argument(
        "--num-channels", type=int, default=32,
        help="Number of channels per board (default: %(default)s)",
    )
    parser.add_argument(
        "--num-boards", type=int, default=16,
        help="Number of boards/files to process (default: %(default)s)",
    )
    parser.add_argument(
        "--stamp-file",
        default=None,
        help="If set, write a stamp file listing the pedestal files produced "
             "(used as the Snakemake completion marker for this rule).",
    )
    args = parser.parse_args()

    strip_indices = list(range(args.num_channels))
    print(f"Channels per board: {len(strip_indices)}")

    telescope = args.telescope
    if telescope is None:
        telescope = Path(args.input_dir).name.split("_")[0]
        if telescope not in ("NERO", "BLU", "ROSSO"):
            parser.error(
                f"Could not infer telescope from --input-dir ('{args.input_dir}'); "
                "please pass --telescope explicitly."
            )
    print(f"Telescope: {telescope}")

    written_files = process_pedestal_files(
        input_dir=args.input_dir,
        output_base=args.output_dir,
        stripIndices=strip_indices,
        telescope=telescope,
        num_boards=args.num_boards,
    )
    print(f"Done. Output written to: {Path(args.output_dir) / telescope}")

    if args.stamp_file:
        stamp_path = Path(args.stamp_file)
        stamp_path.parent.mkdir(parents=True, exist_ok=True)
        with stamp_path.open("w", encoding="utf-8") as handle:
            for written_file in written_files:
                handle.write(f"{written_file}\n")
        print(f"Wrote stamp file: {stamp_path}")
