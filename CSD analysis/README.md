# CSD analysis

## 1) Delta peak / trough → CSD-ready LFP (MATLAB)

**Aim:** test whether delta extrema look like laminar brain generators (sink–source) vs common-mode artifact.

| File | Role |
|------|------|
| `GENERATE_DELTA_PEAK_TROUGH_CSD_LFP.m` | One session: LPAG-center detection, Partners-session events, raw ±250 ms cuts |
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

Defaults: `ref_chan=None`, **`baseline_correct=True`**, window **−500 … +500 ms** (`t_pre=t_post=500`).  
Optional later: `normalize='zscore_channel'`.
