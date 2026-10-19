# CSD analysis

## 1) Delta peak / trough → CSD-ready LFP (MATLAB)

**Aim:** test whether delta extrema look like laminar brain generators (sink–source) vs common-mode artifact.

| File | Role |
|------|------|
| `GENERATE_DELTA_PEAK_TROUGH_CSD_LFP.m` | One session: LPAG-center detection, Partners-session events, raw ±500 ms cuts |
| `Estimate_delta_peak_trough_csd.m` | Fixed params, loop animals, write `.npy` |
| `write_csd_npy.m` | MATLAB → NumPy `.npy` writer |

### Rules baked into GENERATE

1. Detection channel = **center of LPAG** (same geometry rule as Fig 2 PSTHs). **No LPAG → skip session/shank.**
2. **NPX1:** one folder with the full channel list. **NPX2:** one folder **per probe/shank** that has LPAG (reduced `n_chan`).
3. **No re-referencing.**
4. Window **t1 = −0.50 s, t2 = +0.50 s** (±500 ms; subsample later in Python if needed).
5. **No baseline correction** in the export (Python can enable later; `normalize` hook reserved).
6. Events restricted to **`play_session`**: rows with `Type == 'Partners session'` in the behavior txt (same as Fig 6 / phase-coupling GENERATEs).

### Output tree

`Data/Analysis results/delta peak trough CSD/`

```
peak/<animal>/                 # NPX1
peak/<animal>/probe2/          # NPX2 shank
trough/<animal>/
trough/<animal>/probe3/
```

Each leaf folder:

- `lfp.npy` `(n_events, n_chan, n_time)` raw snippets  
- `x_coords.npy`, `y_coords.npy`  
- `T0.npy` sample index of the extremum inside the snippet  
- `fs.npy`, `meta.mat`

Run on a machine with the experimentFS `Data` share:

```matlab
Estimate_delta_peak_trough_csd
```

---

## 2) Visualize evoked LFP (MATLAB, before Python CSD)

| File | Role |
|------|------|
| `Plot_delta_peak_trough_evoked.m` | Load saved npy; one figure/session with peak \| trough subplots |
| `read_csd_npy.m` | Reader for files written by `write_csd_npy` |

Edit `root_folder` at the top (defaults to your `.../delta peak trough CSD` path), then run.

## 3) 1D CSD in Python

| File | Role |
|------|------|
| `csd_1d_pipeline.py` | Load → CSD; `run_delta_csd_from_root(folder)` is the one-call entry |
| `CSD_1D_general.ipynb` | Set `DELTA_ROOT` and run |

```python
from csd_1d_pipeline import run_delta_csd_from_root
run_delta_csd_from_root(
    r"\\experimentfs.home.bccn-berlin.de\...\Data\Analysis results\delta peak trough CSD"
)
```

Or CLI (positional path is enough):

```bash
python csd_1d_pipeline.py "\\experimentfs.home.bccn-berlin.de\...\delta peak trough CSD"
```

Defaults for delta (`DELTA_CSD_DEFAULTS`, used by `run_delta_csd_from_root`):  
`ref_chan=None`, **`baseline_correct=True`**, window **−500 … +500 ms**.

Do not confuse with `DatasetConfig` class fields `t_pre=200` / `t_post=800` — those are only the generic Lombard/stim defaults and are **overridden** for delta runs.

### Saved outputs (per session folder)

Under `<session>/`:
- `csd_1d_array.npy` — stacked CSD `(n_y, n_time, n_x)`
- `csd_1d_x.npy`, `csd_1d_time.npy`
- `csd_1d/csd_x{X}.npy`, `y_x{X}.npy`, `V_x{X}.npy`, `time.npy`
- `csd_1d/figures/kcsd_x{X}.png` — one PNG per shank column

Also copied to `<DELTA_ROOT>/csd_figures/<session_name>/`.

---

## 4) Behavior / play-bout onset and offset with one shared baseline

### MATLAB export

Run:

```matlab
Estimate_behavior_aligned_csd
```

It calls `GENERATE_BEHAVIOR_ALIGNED_CSD_LFP` once per NPX session and writes:

- all 14 individual behaviors with implanted and partner events pooled;
- categories:
  - play: `Pounce`, `CC`, `Boxing`, `Evasion`, `Pin`, `Escape`, `CB`, `CD`;
  - aggression: `PounceI`, `Bite`;
  - exploration: `Grooming`, `Rearing`, `Sniffing`, `Scratching`;
- `self_play` (implanted animal) and `other_play` (partners);
- whole play bouts from `play_bout.m`.

Only complete −1…+2 s onset/offset pairs inside a `Partners session` are
exported. A pair is excluded when any genuine behavior annotation from either
animal overlaps its −1…0 s onset baseline. Exact duplicate annotations are
removed. Neighboring behaviors and overlapping onset/offset windows are
flagged in `events.mat`, while groups with fewer than five events are retained
but marked as low-count.

Every leaf contains raw `lfp.npy`, matching `event_ids.npy`, geometry/timing
arrays, `meta.mat`, and `events.mat` content inside `meta.mat`. It also contains
`evoked_lfp.npy`: the raw-voltage LFP PSTH after the event's onset baseline was
subtracted from both onset and offset trials. No band-pass filtering or
rereferencing is applied.

### Python paired kCSD

`behavior_aligned_kcsd.py` handles the paired onset and offset exports. It uses
raw, unfiltered LFP and deliberately keeps baseline preparation outside the
general kCSD estimator:

1. Match onset and offset rows using `event_ids.npy`.
2. Measure one event- and channel-specific baseline from −1000 to 0 ms
   relative to onset.
3. Subtract that same baseline from both the onset and offset epochs.
4. Average corrected events and pass the two evoked voltages independently
   to `compute_kcsd_1d_columns`.

This means offset is **not** normalized to its own pre-offset period.

Expected paired leaves can be anywhere under the export root as long as one
path component is exactly `onset` and replacing it with `offset` finds the
matching leaf:

```text
behavior/Pounce/onset/<session>/[probeN]/
behavior/Pounce/offset/<session>/[probeN]/
play_bout/onset/<session>/[probeN]/
play_bout/offset/<session>/[probeN]/
```

Each leaf uses the same `lfp.npy`, coordinate, `T0.npy`, and optional `fs.npy`
format described above. Both leaves should also contain identical
`event_ids.npy`; this prevents an offset from receiving another event's
baseline. Matching IDs are required by default. Row-order pairing is available
for old exports only when `stop_on_missing_event_ids=False` is set explicitly.

Run one pair:

```python
from behavior_aligned_kcsd import run_behavior_aligned_kcsd

result = run_behavior_aligned_kcsd(
    r"...\behavior\Pounce\onset\session1\probe2",
    r"...\behavior\Pounce\offset\session1\probe2",
    stop_on_missing_event_ids=True,
)
```

Or discover all individual behaviors, grouped behavior categories, self/other
play, and play-bout pairs below one root:

```python
from behavior_aligned_kcsd import run_behavior_aligned_kcsd_from_root

results = run_behavior_aligned_kcsd_from_root(
    r"...\behavior aligned CSD",
    stop_on_error=False,
)
```

Notebook users can open `CSD_1D_general.ipynb` and run all cells under
**Section C: Behavior/play-bout onset and offset (self-locating)**. That
section searches both experimentFS hostnames and does not depend on Jupyter's
current working directory.

Defaults are raw LFP, no rereference, −1000…+2000 ms response epochs, onset
baseline −1000…0 ms, and the same spatial kCSD parameters as the delta
peak/trough pipeline.
Outputs are written beneath each alignment leaf in
`csd_1d_onset_baseline/`. `paired_baseline_metadata.json` records the baseline
source/window and pairing method; `onset_baseline.npy` stores the actual
event-by-channel baseline values. Trial matrices are memory-mapped and averaged
in event chunks by default; pass `retain_epochs=True` only when trial-level
arrays must remain in the returned Python object.

During a batch run, the Python pipeline reports the current dataset number,
event count, onset/offset stage, probe column, matrix destination, figure path,
per-dataset duration, and final success/failure count. All selected time frames
for one probe column are solved in a single vectorized `KCSD1D` call; the
previous frame-by-frame implementation produced the same values but was much
slower for three-second epochs.

For behavior datasets, kCSD is still fitted independently to each staggered
electrode column. The column maps are then linearly interpolated onto their
shared, non-extrapolated depth range and averaged with equal column weight:
four columns for NPX1 and the two columns within each NPX2 probe. The primary
image is `kcsd_xmean.png`; onset and offset use shared symmetric color limits.
Per-column matrices remain saved for QC, while legacy per-column PNGs are
removed by default. Averaged outputs include `csd_mean_columns.npy`,
`csd_std_columns.npy`, `y_mean_columns.npy`, `V_mean_columns.npy`,
`n_columns_contributing.npy`, and `x_center.npy`.

### Sorted figure copies

If per-column kCSD matrices are already saved, do **not** rerun kCSD. Run the
final Section C notebook cell, which first averages the existing matrices:

```python
from behavior_aligned_kcsd import (
    average_saved_behavior_kcsd_from_root,
    organize_behavior_aligned_figures,
)
average_saved_behavior_kcsd_from_root(BEHAVIOR_ROOT)
organize_behavior_aligned_figures(BEHAVIOR_ROOT)
```

The averaging function loads `csd_x*.npy`, writes the mean/std arrays and
`kcsd_xmean.png`, and removes legacy individual-column PNGs. It never calls
`KCSD1D`. The organizer then rebuilds:

```text
sorted_figures/
  by_session/onset/
  by_session/offset/
  by_behavior/onset/
  by_behavior/offset/
  figure_index.csv
```

Session-first names use
`session__scope__behavior__probe__xpos.png`; behavior-first names use
`scope__behavior__session__probe__xpos.png`. Including scope and probe avoids
ambiguity and collisions between individual behaviors, pooled categories,
self/other play, play bouts, and multi-probe sessions. New runs copy only the
averaged `xmean` PNG; per-column PNGs are used only as a fallback for legacy
results that do not yet contain an averaged image.

### Session area overlap

Before pooling sessions, estimate the anatomical areas actually sampled by
each recording from `y_mean_columns.npy` and `Area_limits_GoodLooking.xlsx`.
The workbook is only read. For every pair of recordings the overlap matrix
stores, at row `i` and column `j`, the percentage of areas sampled by `i`
that are also sampled by `j`. That matrix is not symmetric.

A recording pair may be averaged together only when both directions are at
least 50%. The retained group is the largest set in which every pair passes
that test. If two such sets are the same size, the one with the higher mean
overlap is kept. Recordings outside that group stay out of the cross-session
area average. Results are written to `session_overlap/`.

```python
from area_csd_consistency import summarize_session_area_overlap
session_overlap = summarize_session_area_overlap(BEHAVIOR_ROOT, AREA_TABLE)
summarize_area_csd_across_sessions(
    BEHAVIOR_ROOT,
    AREA_TABLE,
    session_allowlist=session_overlap["selected_sessions"],
)
```

### Area consistency across sessions

After `csd_mean_columns.npy` exists, `area_csd_consistency.py` groups channels
from the same anatomical structure and pools those maps across sessions. It
does not estimate kCSD. The area table is `Data/Area_limits_GoodLooking.xlsx`,
using the animal key and depth rules from `GENERATE_PSTH_AREA_MAPS.m`:

- NPX1 depths must have `Probe_Area` PAG.
- NPX2 depths must match the `probeN` folder's probe number.
- A depth matches a row when `depth_start <= y < depth_end + 1`.
- The narrowest matching span is kept. A remaining tie is left unassigned.

Depth is converted to relative position inside the selected area,
`(y - depth_start) / (depth_end - depth_start)`. Values are interpolated only
inside depths the electrodes actually cover. Probes that contain the same area
in one session are averaged first, with equal probe weight, so each session is
one observation. Across sessions the outputs are the mean, SEM, fraction
positive, fraction negative, and sign consistency of the **saved kCSD sign**.
That sign is not relabeled as sink or source.

```python
from area_csd_consistency import summarize_area_csd_across_sessions
summarize_area_csd_across_sessions(BEHAVIOR_ROOT, AREA_TABLE)
```

Results are written to `area_consistency/<scope>/<behavior>/<onset|offset>/<area>/`
with `mean_csd.png`, `sign_consistency.png`, and the arrays behind them.
`summary.csv` lists every area map. `stable_features.csv` keeps bins that have
at least two sessions, high sign agreement, and a large absolute mean.
Consistency images are blank where fewer than two sessions overlap. The
`kCSD_env` kernel needs `openpyxl` (`python -m pip install openpyxl`) to read
the spreadsheet. The reader copies the workbook into memory and does not
save it. `Area_limits_GoodLooking.xlsx` is never modified.

Before that average, each probe segment is smoothed with a Gaussian whose
full width at half maximum is 250 ms and 100 µm. Samples are dropped where
half or more of the kernel falls outside the observed depths or the saved
time range. A segment whose unique depths span less than 200 µm is left out
of the session. With at least five sessions, every remaining depth and time
is compared with zero by a one-sample t-test. Positive and negative bins
form separate clusters. A cluster is kept when its block is at least 25 ms
long and covers at least 10% of the 0–1 depth axis. This is a size filter
after an uncorrected pointwise p < 0.05 threshold. The unsmoothed mean and
sign-consistency figures remain. `smoothed_mean_clusters.png` draws the kept
clusters on the smoothed mean, and `clusters.csv` lists them.

The same run copies those PNGs into two browsing folders. Filenames start with
the behavior, so one folder lists every merged onset map and the other lists
every merged offset map:

```text
area_consistency/by_behavior/onset/
area_consistency/by_behavior/offset/
area_consistency/by_behavior/figure_index.csv
```

`play__category__LPAG__mean_csd.png` is the across-session play map for LPAG.
Original area folders are not changed. If the maps already exist, rerun only
`organize_area_consistency_figures(BEHAVIOR_ROOT)`.

While each saved map is read, a copy is written to
`area_consistency/session_cache/`. Refreshing the area folder does not delete
that cache. On experimentFS, Windows often raises `PermissionError` when
`shutil.rmtree` tries to remove `area_consistency` (a file is open or
read-only). The run now overwrites the old files instead of stopping. The
failed run that hit that delete had already finished the read in memory, but
the old code wrote nothing, so that attempt cannot be resumed. Rerun the area
cell once. After `session_cache` exists, rebuild the means, clusters, and
`by_behavior` copies without reading `csd_mean_columns.npy` again:

```python
from area_csd_consistency import finalize_area_consistency_from_cache
finalize_area_consistency_from_cache(BEHAVIOR_ROOT)
```
