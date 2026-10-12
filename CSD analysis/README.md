# CSD analysis (1D)

Generalized pipeline from the session notebooks (`CSD 1D 230516`, `230606`, `230608`, `250320`, …).

## Files

| File | Role |
|------|------|
| `csd_1d_pipeline.py` | Load → preprocess → 1D kCSD / FD CSD → save |
| `CSD_1D_general.ipynb` | Edit a `DATASETS` list and run many sessions |

## Per-folder inputs

```
lfp.npy          (n_stim, n_chan, n_time)
x_coords.npy     (n_chan,)
y_coords.npy     (n_chan,)
T0.npy           scalar onset sample
```

## Run on new datasets

1. Open `CSD_1D_general.ipynb`
2. Fill `DATASETS` with each folder + any overrides (`chan_start` / `chan_end`, `channels_remove`, `ref_chan`, …)
3. Run all cells

Or from a terminal (on a machine that can see the data):

```bash
pip install numpy scipy matplotlib kcsd
python "csd_1d_pipeline.py" --folder "/path/to/ds1" --folder "/path/to/ds2" --ref deepest
```

`EXAMPLE_DATASETS` inside `csd_1d_pipeline.py` mirrors the settings used in the uploaded Lombard session notebooks.

## Fixes vs session notebooks

- Bad-channel removal, referencing, and channel-range slicing stay index-aligned (older notebooks subtracted a reference on the full probe then sliced `x`/`y` only).
- Channels are sorted by depth before 1D kCSD.
- Finite-difference CSD is available even if `kcsd` is not installed (`compute_kcsd=False`).
