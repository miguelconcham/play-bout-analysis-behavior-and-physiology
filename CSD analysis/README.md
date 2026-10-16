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

### Sorted figure copies

After kCSD finishes, run the final Section C notebook cell or call:

```python
from behavior_aligned_kcsd import organize_behavior_aligned_figures
organize_behavior_aligned_figures(BEHAVIOR_ROOT)
```

This leaves source PNGs untouched and rebuilds:

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
self/other play, play bouts, and multi-probe sessions.
