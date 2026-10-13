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

## 2) 1D CSD in Python

| File | Role |
|------|------|
| `csd_1d_pipeline.py` | Load → (optional) ref/baseline → 1D kCSD / FD CSD |
| `CSD_1D_general.ipynb` | Batch notebook |

### Delta defaults (NPX1 + NPX2 shank folders)

```python
from csd_1d_pipeline import DELTA_CSD_DEFAULTS, discover_delta_csd_datasets, run_many

root = r"\\experimentfs.bccn-berlin.pri\...\Data\Analysis results\delta peak trough CSD"
datasets = discover_delta_csd_datasets(root)
run_many(datasets, defaults=DELTA_CSD_DEFAULTS)
```

Or CLI:

```bash
python csd_1d_pipeline.py --delta-root "/path/to/delta peak trough CSD"
```

`DELTA_CSD_DEFAULTS` sets `ref_chan=None`, `baseline_correct=False`, `t_pre=t_post=250` ms.  
Set `normalize='zscore_channel'` later if you want scaling without baseline subtraction.

Both geometries load the same four files; NPX2 simply has fewer channels (one shank). The loader reports `single_shank` vs `multi_column` from unique `x`.

### Lombard / stim CSD notebooks

`EXAMPLE_DATASETS` still mirrors the older session-specific stim notebooks (those used referencing / baseline).
