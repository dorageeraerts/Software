# Good run selection (`select_runs.py`)

This folder builds the list of **good runs** for each hodoscope, by applying a sequence of
quality criteria to the outputs of the parsing, of the pedestal check and of the monitoring
scan.

The main entrypoint is:

- `select_runs.py`

For every hodoscope it writes a plain list of run numbers, a pie chart of the cut flow and two
control plots.

---

## 1) How to run

The script needs `pandas`, `numpy` and `matplotlib`, which are not in the base environment on
the working nodes. Run it inside the container:

```bash
cd muraves/analysis
singularity exec -B /user/<you> -B /ada_mnt -B /pnfs /user/<you>/muraves-sing.sif \
    python select_runs.py --hodoscope ALL --details
```

The inputs live on the storage element, so `-B /pnfs` has to be in the bind list.

Common invocations:

```bash
python select_runs.py                              # all hodoscopes, ADC files
python select_runs.py --hodoscope BLU --details    # one hodoscope, plus the per-run csv
python select_runs.py --type BOTH                  # require ADC and PIEDISTALLI to pass
python select_runs.py --no-plot --output-dir /tmp  # numbers only, written elsewhere
```

The script prints the full cut flow as it goes, so the run log is the first place to look when
a number surprises you.

---

## 2) Inputs

Everything is read from the production data on the storage element, under

```
/pnfs/iihe/muraves/muraves_DATA          # DATA_PATH, --data-path to read a copy elsewhere
```

| File, relative to that root | Provides | Used by |
|---|---|---|
| `PARSED/<COLOR>/<version_parsing>/df_overview_<COLOR>_<version_parsing>.pkl` | one row per run and file type, from `run_manager.read_summary_file` | criteria 1, 3, 4 |
| `PEDESTAL/<COLOR>/<version_pedestal>/<COLOR>/pedestal_deviation_check/<COLOR>_HISTchi2_outliers_SUMMARY.txt` | number of pedestal outlier channels per run | criterion 5 |
| `RAW_GZ/<COLOR>/data_scan_<color>.pkl` | `trigger_rate`, `accidental_rate` and `wp` per run | criteria 6, 7 and the plots |

`<version_parsing>` and `<version_pedestal>` are the processing versions of the two production
steps. They **advance independently**, so each has its own option, `--version-parsing`
(`v0` today) and `--version-pedestal` (`v1` today). Note that the parsing version appears twice
in its path, once as the directory and once inside the file name; the script uses the same
string for both.

The `<COLOR>` segment is repeated in the pedestal path, and the scan file name is the only one
in lower case. Both are quirks of how the producers write their output, not typos here.

The three sources are **independent**. In particular the monitoring scan contains runs that
never appear in the parsing overview, because slow-control information exists for runs whose
raw files were never written. Those runs are not part of the selection.

Missing inputs degrade gracefully: if the pedestal summary or the monitoring scan is absent,
the corresponding criteria are skipped with a warning and the rest of the selection still runs.

---

## 3) The criteria

They are applied **in this order**, and a run must pass all of them.

### 1. The raw run must exist

`missing_run == False`. `read_summary_file` sets this when the list of sub-runs is empty, that
is when the raw files were not found.

### 2. All sub-runs present *(currently skipped)*

The intended cut is "exactly `N_SUBRUNS_REQUIRED` sub-runs" (4 today). The information is not
in the pkl files: `read_summary_file` only stores `missing_run`, not the sub-run count. The cut
is implemented and activates automatically as soon as the dataframe carries one of the column
names listed in `N_SUBRUNS_COLUMNS`; until then it is skipped with a warning.

To enable it: add the sub-run count in `RunValidation` / `read_summary_file` and regenerate the
parsing overview pkl files.

### 3. Parsing succeeded

`is_run_ok == True`, i.e. the summary file reports `"status": "ok"`.

### 4. Missing or shifted bits

Two conditions:

- **4a** the run has no mismatches, or all of its mismatches are recoverable;
- **4b** `mismatch_prob < MISMATCH_PROB_MAX` (0.05).

`mismatch_prob` is the probability that a line contains at least one corrupted value:

```
mismatch_prob = mismatch_counter * N_BITS_PER_LINE / (n_events * N_LINES_PER_EVENT)
```

with `n_events` 40000 for ADC and 50000 for PIEDISTALLI, 16 lines per event and 625 bits per
line. Same definition as `dev/good_runs/good_runs_stat.py`. A `mismatch_counter` of NaN means
"no mismatch at all" and is read as zero, not as unknown.

### 5. Pedestal quality

Number of affected channels `<= MAX_AFFECTED_CHANNELS` (10), from the pedestal outlier summary.
Runs absent from that summary are **kept** by default and counted in the log; pass
`--require-pedestal-info` to reject them instead.

### 6. Accidental coincidences

`accidental_rate <= MAX_ACCIDENTAL_RATE` (0), i.e. runs with any accidental coincidence are
rejected.

This cut is placed **before** the trigger rate cut on purpose. It removes the pathological
runs (see section 6 below), so the trigger rate band is computed on a sample that is already
clean.

### 7. Trigger rate stability *(applied last)*

```
median - N * sigma  <=  trigger_rate  <=  median + N * sigma
```

with `N = TRIGGER_RATE_N_SIGMA` (3) and `sigma` the **robust** spread `1.4826 * MAD`, both
computed on the runs that reach this cut.

Two choices deserve an explanation.

**Why not the standard deviation.** It is dominated by the very runs the cut has to reject. On
BLU, about one run in five reports a trigger rate of O(5e5) Hz. The standard deviation then
reaches 1.6e5 Hz against a median of 13 Hz, so the band covers the whole sample and its lower
edge is negative: a run with an abnormally **low** rate could never be rejected. The median
absolute deviation is built from medians, so a fifth of the sample can be arbitrarily wrong
without moving it. The factor 1.4826 makes it equal to the standard deviation for a gaussian
sample, which keeps the "median ± N sigma" reading of the cut.

**Why N = 3 and not 1.** The rate distribution is a narrow core with long tails, not a
gaussian. At N = 1 the cut throws away a third to a half of perfectly good runs. At N = 3 it
reproduces the usual Tukey outlier fences (`Q1 - 1.5 IQR`, `Q3 + 1.5 IQR`) almost exactly.

Runs whose rates are missing from the monitoring scan are **kept** by default for criteria 6
and 7; pass `--require-monitoring-info` to reject them. See the caveat in section 6.

---

## 4) Outputs

All written to `--output-dir` (this folder by default), one set per hodoscope and file type.

| File | Content |
|---|---|
| `<COLOR>_type<TYPE>_SELECTED_RUNS.txt` | the good runs, one run number per line |
| `<COLOR>_type<TYPE>_SELECTED_RUNS.csv` | per-run values used by the selection (`--details` only) |
| `<COLOR>_type<TYPE>_SELECTION_CUTFLOW.pdf` | pie chart of the cut flow |
| `<COLOR>_type<TYPE>_ACCIDENTAL_RATE_SELECTED.pdf` | criterion 6 control plot |
| `<COLOR>_type<TYPE>_TRIGGER_RATE_SELECTED.pdf` | criterion 7 control plot |

`--no-plot` skips the three pdf files.

### The cut flow pie chart

One slice per cut, sized by the runs **that cut removes on top of the previous ones**, plus one
green slice for the selected runs. All slices are fractions of the initial number of runs, so
they add up to 100%. Cuts that removed nothing get no slice and are listed under the chart.

### The two control plots

Both show trigger rate versus run number, with a lower panel giving the **working point
temperature** (`wp`, column 5 of the slow-control file, in °C) of the same runs. Only runs that
appear in the rate panel are drawn in the working point panel.

- **`ACCIDENTAL_RATE_SELECTED`** draws the runs reaching criterion 6, that is the survivors of
  criteria 1 to 5, split into the ones it keeps (blue) and the ones it rejects (orange). No
  band is drawn: the band is only defined at criterion 7, after this cut.
- **`TRIGGER_RATE_SELECTED`** draws the runs reaching criterion 7, split into the selected runs
  (blue) and the ones the band rejects (orange), with the median, the band and its edges.

Neither plot shows runs rejected by criteria 1 to 5, nor runs that are in the monitoring scan
but not in the parsing overview. That is deliberate: every point drawn is a run the selection
actually considered, and the two plots chain, the blue points of the first being exactly the
points of the second.

When the rates span more than a decade above the median, the rate panel is split in two: the
full range on a symlog scale, and a linear zoom on the core of the distribution. This happens
for BLU, whose O(5e5) Hz runs would otherwise squash everything else onto a single line.

---

## 5) Command line options

| Option | Default | Effect |
|---|---|---|
| `--hodoscope {BLU,NERO,ROSSO,ALL}` | `ALL` | which hodoscope to process |
| `--type {ADC,PIEDISTALLI,BOTH}` | `ADC` | file type the criteria apply to; `BOTH` requires the run to pass for both types |
| `--data-path` | `/pnfs/iihe/muraves/muraves_DATA` | root of the production data |
| `--version-parsing` | `v0` | processing version of the parsing step, used for its directory and its file name |
| `--version-pedestal` | `v1` | processing version of the pedestal step |
| `--mismatch-prob-max` | `0.05` | criterion 4b threshold |
| `--max-affected-channels` | `10` | criterion 5 threshold |
| `--n-subruns` | `4` | criterion 2 target, only used if the column exists |
| `--max-accidental-rate` | `0.0` | criterion 6 threshold |
| `--trigger-rate-n-sigma` | `3.0` | criterion 7 band half-width, in robust sigmas |
| `--require-pedestal-info` | off | reject runs absent from the pedestal summary |
| `--require-monitoring-info` | off | reject runs absent from the monitoring scan |
| `--details` | off | also write the per-run csv |
| `--no-plot` | off | do not write the pdf files |
| `--output-dir` | this folder | where to write everything |

The defaults live at the top of the script (`DATA_PATH`, `VERSION_PARSING`, `VERSION_PEDESTAL`,
`MISMATCH_PROB_MAX`, `MAX_AFFECTED_CHANNELS`, `MAX_ACCIDENTAL_RATE`, `TRIGGER_RATE_N_SIGMA`,
...), so a permanent change belongs there rather than in a habit of typing flags. When a new
production version comes out, bumping `VERSION_PARSING` or `VERSION_PEDESTAL` is the whole
change.

---

## 6) Current results and things to know

Cut flow for `--type ADC`, parsing version v0 and pedestal version v1:

| | BLU | NERO | ROSSO |
|---|---|---|---|
| runs considered | 10100 | 9700 | 14899 |
| 1. run not available | 2228 | 537 | 2210 |
| 3. parsing failed | 313 | 890 | 1294 |
| 4a. unrecoverable mismatches | 0 | 0 | 0 |
| 4b. mismatch probability | 0 | 771 | 263 |
| 5. bad pedestal | 141 | 169 | 148 |
| 6. accidental rate > 0 | 2676 | 89 | 201 |
| 7. trigger rate outside 3 sigma | 70 | 296 | 85 |
| **selected** | **4672 (46.3%)** | **6948 (71.6%)** | **10698 (71.8%)** |

Trigger rate band actually used:

| | median [Hz] | robust sigma [Hz] | band [Hz] |
|---|---|---|---|
| BLU | 12.683 | 1.236 | 8.97 – 16.39 |
| NERO | 13.100 | 1.186 | 9.54 – 16.66 |
| ROSSO | 14.467 | 0.989 | 11.50 – 17.43 |

Four points worth keeping in mind when reading these numbers.

**BLU loses a quarter of its runs to accidental coincidences.** The same runs are the ones
reporting O(5e5) Hz trigger rates. Criterion 6 removes them, which is why criterion 7 then has
so little left to do on BLU.

**ROSSO selects 1685 runs that were never measured.** Those runs, between run 261 and run 3491,
have no trigger rate in the monitoring scan. They are kept by default, so they pass criteria 6
and 7 by convention rather than by measurement, and they are 16% of the ROSSO list. Use
`--require-monitoring-info` to reject them instead, which brings ROSSO to 9013 runs (60.5%).

**The working point changes over time.** On BLU it steps from 25 °C to 15 °C around run 4300,
to 10 °C at run 6100, and back to 25 °C at run 8100, and the trigger rate baseline follows.
Part of the spread the band measures is therefore a configuration change, not detector
instability. Computing the band per working point instead of globally is the obvious refinement
if that becomes a problem.

**Criterion 2 is not active.** The sub-run count is missing from the inputs, see above.

---

## 7) Using the script as a library

`main()` is a thin wrapper; the pieces are importable and free of side effects:

```python
import select_runs as sr

df = sr.load_parsing_overview("BLU")            # parsing overview, version v0
ped = sr.load_pedestal_outliers("BLU")          # pedestal summary, version v1, or None
mon = sr.load_monitoring_scan("BLU")            # monitoring rates, or None

good_runs, details, cutflow, trigger_info = sr.select_runs(
    df, ped, mon, run_type="ADC", trigger_rate_n_sigma=3.0, label="[BLU]")
```

- `good_runs` sorted array of run numbers;
- `details` per-run dataframe of everything the selection looked at;
- `cutflow` dataframe with one row per cut (`runs_left`, `removed`, percentages), the input of
  `plot_cutflow_pie`;
- `trigger_info` dict with the band (`median`, `sigma`, `low`, `high`) and, for each population
  the plots need, the run numbers, trigger rates and working points; the input of
  `plot_trigger_rate`.

Each loader takes the version and the root it should read from, so a comparison between two
productions is one extra argument:

```python
old = sr.load_parsing_overview("BLU", version="v0")
new = sr.load_parsing_overview("BLU", version="v1", data_path="/some/other/mirror")
```

`robust_sigma(values)` is available on its own if you need the same spread elsewhere.
