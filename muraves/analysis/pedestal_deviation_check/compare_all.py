"""
Created on Wed May  6 16:47:03 2026
@author: Gabor
"""

import os
import numpy as np

# -------------------------
# USER SETTINGS
# -------------------------
hist_file = "NERO_standard_hist_bin3.txt"

run_numbers = range(4593, 12199)

dataname = "PIEDISTALLI"
dir_ADC = "/media/gabor/Expansion/MURAVES/PARSERED/NERO"

n_boards = 16
n_events_perrun = 50000

bins = np.arange(1400, 1801, 3)

output_filename = "NERO_chi2_results_full.txt"

# -------------------------
# HELPER PRINT
# -------------------------
def log(message):
    print(message, flush=True)

# -------------------------
# OPEN OUTPUT FILE
# -------------------------
outfile = open(output_filename, "w")

# -------------------------
# LOAD SAVED HISTOGRAMS
# -------------------------
log(f"[READ ] Loading histogram file: {hist_file}")

saved = np.loadtxt(hist_file, comments="#")

log("[DONE ] Histogram file loaded")

all_centers = saved[:, 0]

# rows=bins, cols=channels
saved_counts = saved[:, 1:]

n_channels = saved_counts.shape[1]

epsilon = 1e-12

# -------------------------
# WRITE HEADER
# -------------------------
outfile.write("# Run Channel Chi2NDF\n")

# -------------------------
# LOOP OVER RUNS
# -------------------------
total_runs = len(run_numbers)

for i, run in enumerate(run_numbers, start=1):

    ADCname = f"{dataname}_run{run}.txt"
    filename = os.path.join(dir_ADC, ADCname)

    log(f"[RUN {i:02d}/{total_runs}] Processing run {run}")

    if not os.path.isfile(filename):
        log(f"[MISS ] Run {run}: file not found")
        continue

    # -------------------------
    # LOAD RUN
    # -------------------------
    try:
        data = np.genfromtxt(
            filename,
            filling_values=0,
            invalid_raise=False
        )

        # Skip empty files
        if data.size == 0:
            log(f"[SKIP ] Run {run}: empty file")
            continue

        # Ensure 2D shape
        data = np.atleast_2d(data)

        data = np.nan_to_num(data)

    except Exception as e:
        log(f"[ERROR] Failed reading run {run}: {e}")
        continue

    # -------------------------
    # BUILD ADC ARRAY
    # -------------------------
    data_ADC_run = np.zeros(
        (n_events_perrun, n_channels),
        dtype=np.int16
    )

    for sk in range(n_boards):

        data_ADC_run[:, sk*32:sk*32+32] = data[
            :,
            3+sk*39:3+sk*39+32
        ].astype(np.int16)

    # -------------------------
    # CHANNEL LOOP
    # -------------------------
    for ch in range(n_channels):

        saved_ch = saved_counts[:, ch]

        counts_run, edges = np.histogram(
            data_ADC_run[:, ch],
            bins=bins
        )

        peak = np.max(counts_run)

        if peak == 0:
            continue

        counts_run = counts_run / peak

        # Peak from saved histogram
        peak_index = np.argmax(saved_ch)
        peak_adc = all_centers[peak_index]

        cut_min = peak_adc
        cut_max = peak_adc + 40

        mask = (
            (all_centers >= cut_min) &
            (all_centers <= cut_max)
        )

        y_saved = saved_ch[mask]
        y_run = counts_run[mask]

        chi2 = np.sum(
            (y_run - y_saved)**2 /
            (y_saved + epsilon)
        )

        ndf = len(y_saved) - 1

        chi2_ndf = chi2 / ndf

        # Write result
        outfile.write(
            f"{run} "
            f"{ch} "
            f"{chi2_ndf:.6f}\n"
        )

# -------------------------
# CLOSE FILE
# -------------------------
outfile.close()

log(f"[DONE ] Results saved to {output_filename}")
