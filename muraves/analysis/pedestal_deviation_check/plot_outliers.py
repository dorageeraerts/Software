from collections import defaultdict
import matplotlib.pyplot as plt


# File containing all measurements:
# run channel value
ORIGINAL_FILE = "NERO_HISTchi2.txt"

# Output produced by the outlier-finding script:
# Run Channel Value Mean Sigma Zscore
OUTLIER_FILE = "NERO_HISTchi2_outliers.txt"

SUMMARY_FILE = "NERO_HISTchi2_outliers_SUMMARY.txt"
PLOT_FILE = "NERO_HISTchi2_outliers_PLOT.png"

CHANNEL_THRESHOLD = 10


def read_all_runs(filename):
    """Read all run numbers from the original measurement file."""
    runs = set()

    with open(filename, "r") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            columns = line.split()

            if len(columns) < 3:
                print(
                    f"Warning: ignoring malformed line "
                    f"{line_number} in {filename}"
                )
                continue

            try:
                run = int(columns[0])
            except ValueError:
                print(
                    f"Warning: invalid run number on line "
                    f"{line_number} in {filename}"
                )
                continue

            runs.add(run)

    return runs


def read_outlier_channels(filename):
    """
    Read the outlier file and store unique affected channels for each run.

    Using a set prevents the same channel from being counted twice
    in one run.
    """
    affected_channels = defaultdict(set)

    with open(filename, "r") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            columns = line.split()

            if len(columns) < 2:
                print(
                    f"Warning: ignoring malformed line "
                    f"{line_number} in {filename}"
                )
                continue

            try:
                run = int(columns[0])
                channel = int(columns[1])
            except ValueError:
                print(
                    f"Warning: invalid run or channel on line "
                    f"{line_number} in {filename}"
                )
                continue

            affected_channels[run].add(channel)

    return affected_channels


def main():
    all_runs = read_all_runs(ORIGINAL_FILE)
    affected_channels = read_outlier_channels(OUTLIER_FILE)

    if not all_runs:
        raise RuntimeError(f"No valid runs found in {ORIGINAL_FILE}")

    # Include runs appearing only in the outlier file as a safety measure.
    all_runs.update(affected_channels.keys())

    sorted_runs = sorted(all_runs)

    # Runs without outliers receive a count of zero.
    outlier_counts = [
        len(affected_channels.get(run, set()))
        for run in sorted_runs
    ]

    total_runs = len(sorted_runs)

    runs_over_threshold = [
        run
        for run, count in zip(sorted_runs, outlier_counts)
        if count > CHANNEL_THRESHOLD
    ]

    # Write run-by-run output.
    with open(SUMMARY_FILE, "w") as file:
        file.write("# Run Number_of_affected_channels\n")

        for run, count in zip(sorted_runs, outlier_counts):
            file.write(f"{run} {count}\n")

    # Create the plot.
    plt.figure(figsize=(12, 6))
    plt.plot(
        sorted_runs,
        outlier_counts,
        marker=".",
        linewidth=1,
    )

    # Show the threshold of 10 channels.
    plt.axhline(
        y=CHANNEL_THRESHOLD,
        linestyle="--",
        label=f"{CHANNEL_THRESHOLD} channels",
    )

    plt.xlabel("Run number")
    plt.ylabel("Number of affected channels")
    plt.title("Number of outlier channels per run")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig(PLOT_FILE, dpi=200)
    plt.show()

    print(f"Total number of runs: {total_runs}")
    print(
        f"Runs with more than {CHANNEL_THRESHOLD} affected channels: "
        f"{len(runs_over_threshold)}"
    )

    if runs_over_threshold:
        print("Run numbers with more than 10 affected channels:")
        print(" ".join(str(run) for run in runs_over_threshold))

    print(f"\nRun-by-run counts written to: {SUMMARY_FILE}")
    print(f"Plot saved to: {PLOT_FILE}")


if __name__ == "__main__":
    main()
