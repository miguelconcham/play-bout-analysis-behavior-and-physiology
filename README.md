# Play Bout Analysis — Figure 1

MATLAB and Python code for Figure 1 analyses: play-bout statistics, motion–USV (LDA/UMAP) space, and HMM-based play-state detection and plotting.

This repository contains **code only**. Experimental data are archived separately (see [Data availability](#data-availability)).

## Repository structure

```
Figure 1/
├── Figure 1 Play Bouts.m          # Play-bout duration / interval distributions (Fig 1G)
├── Figure 1 LDA.m                 # Motion–USV UMAP + LDA plots (Fig 1C, S1B–D)
├── Figure 1 HMM.m                 # HMM example trace + summary panels (Fig 1B, 2E, S2)
├── BACKUP_DEPENDENCIES.txt        # Full list of inputs each script loads
├── outputs/                       # Saved figures (created locally; not in git)
├── LDA analysis/
│   ├── Obtain_Feature_USV_space.m
│   ├── get_behavior_properties.m
│   ├── normalizedRadialPosition.m
│   └── Convolutional classification.ipynb
└── HMM modeling/
    ├── CreateHMMFIiles.m          # Batch: inputs → HMM fit → analysis
    ├── Create_HMM_inputs.m
    ├── Load_HMM_outputs_and_analyze.m
    ├── python hmm estimate.txt
    └── README.md
```

## Requirements

### MATLAB
- Statistics and Machine Learning Toolbox (`pca`, `fitglm`, `mdscale`, …)
- Signal Processing Toolbox (`spectrogram`, `smoothdata`)
- Image Processing Toolbox (`imfilter`, `imsharpen`, …)
- Third-party on path: `readNPY` / `writeNPY` ([npy-matlab](https://github.com/kwikteam/npy-matlab)), `run_umap`, `play_bout`, `fit_exp_heaviside3`

### Python
- `numpy`, `ssm` (HMM fitting), `torch` (CNN notebook)

## Setup after cloning

1. Clone this repository.
2. Download the **Figure 1 data bundle** (see [Data availability](#data-availability)).
3. Extract so the layout matches:

```
Codes repository/
├── Figure 1/          ← this repo
└── Data/              ← downloaded archive
    ├── Behavior backups/
    ├── Analysis results/
    ├── HMM data/
    ├── Synch data/
    └── CallDetectionBackup/
```

4. Scripts expect data at:
   `...\Codes repository\Data\`
   (already configured in the Figure 1 `.m` files).

5. Add external MATLAB functions to your path (`play_bout`, `readNPY`, `run_umap`, etc.).

## Running analyses

| Goal | Script |
|------|--------|
| Play-bout plots & fits | `Figure 1/Figure 1 Play Bouts.m` |
| LDA / UMAP figures | `Figure 1/Figure 1 LDA.m` |
| HMM figures | `Figure 1/Figure 1 HMM.m` |
| Regenerate HMM inputs (needs full raw data) | `Figure 1/HMM modeling/CreateHMMFIiles.m` |
| Fit HMMs in Python | `Figure 1/HMM modeling/python hmm estimate.txt` |

Pipeline order for HMM: see `Figure 1/HMM modeling/README.md`.

---

## Data availability

**Data are not stored in this GitHub repository** (too large for git; see `.gitignore`).

### What the bundled `Data/` folder contains (minimum for plotting)

Copied subset used by Figure 1 scripts (~141 files):

- `Behavior backups/` — manual behavior annotations (`.txt`)
- `Analysis results/locomotive behaviors 2 partners/` — `all_behavior.mat`, etc.
- `Analysis results/Behavior classification/` — trained LDA models
- `Analysis results/HMM 2 and 3 states 2 partners/` — HMM analysis structs
- `HMM data/HMM raw data/` — HMM states, model comparison, `PropAndTime.mat`
- `Synch data/B1D1 1013 Dual/` — example session sync model (Fig 1B)
- `CallDetectionBackup/` — example `.wav` + `_Stats.xlsx` (Fig 1B)

Full file list: see `Data/COPY_MANIFEST.txt` in the archived bundle.

### Where to download data

> **Update these links before leaving the lab:**

| Archive | Link | Notes |
|---------|------|--------|
| Zenodo (recommended) | `[DOI or URL]` | Permanent DOI; good for publications |
| OSF | `[URL]` | Alternative open archive |
| Lab / institutional storage | `[URL]` | If Zenodo/OSF not used |

**Suggested Zenodo record title:**  
`Play Bout Analysis — Figure 1 data bundle (behavior, HMM, LDA)`

Upload a single zip: `Figure1_Data_bundle.zip` containing the `Data/` folder.

### Full raw data (optional, for re-running pipelines from scratch)

Not included in the minimal bundle. If archived separately:

- `Traking backups/` — `traking_structure` `.mat` per session  
- `Synch data/` — all animals  
- `CallDetectionBackup/` — all sessions  
- `OLD DATA FORMAT/` or equivalent session folders  

Document the link in the same Zenodo record (second zip) or in lab handover notes.


## Citation

If you use this code, cite the associated publication and the data DOI once registered.

## Contact

> **Update before archiving:** `[your email]`
