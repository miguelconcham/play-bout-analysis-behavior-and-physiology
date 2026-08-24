# Play Bout Analysis — Behavior and Physiology

MATLAB and Python code for play-bout analyses: behavior statistics, motion–USV (LDA/UMAP) space, HMM play-state detection, and LFP physiology (PSTH power by frequency and behavior, coupling, phase).

This repository contains **code only**. Experimental data are archived separately (see [Data availability](#data-availability)).

**GitHub:** [play-bout-analysis-behavior-and-physiology](https://github.com/miguelconcham/play-bout-analysis-behavior-and-physiology)

## Repository structure

```
play bout analysis behavior and physiology/
├── Figure 1/
│   ├── Figure 1 Play Bouts.m          # Play-bout duration / interval distributions (Fig 1G)
│   ├── Figure 1 LDA.m                 # Motion–USV UMAP + LDA plots (Fig 1C, S1B–D)
│   ├── Figure 1 HMM.m                 # HMM example trace + summary panels (Fig 1B, 2E, S2)
│   ├── BACKUP_DEPENDENCIES.txt        # Full list of inputs each script loads
│   ├── outputs/                       # Saved figures (created locally; not in git)
│   ├── LDA analysis/
│   └── HMM modeling/
│       └── README.md                  # HMM pipeline order
│
├── Figure 2/
│   ├── Figure 2.m                     # Main Figure 2 plotting script
│   ├── Play partner.m
│   ├── Reciprocity.m
│   ├── Estimate percentage of R2.m
│   ├── Supporting codes/              # Original scripts (legacy; kept unchanged)
│   └── Supporting codes v2/         # Renamed, simplified scripts (recommended)
│
├── Data/
│   └── Analysis results/
│       └── psth power by frequency and behavior/   # Figure 2 PSTH / results .mat (local copy)
└── README.md
```

## Figure 1

### Requirements

**MATLAB**
- Statistics and Machine Learning Toolbox (`pca`, `fitglm`, `mdscale`, …)
- Signal Processing Toolbox (`spectrogram`, `smoothdata`)
- Image Processing Toolbox (`imfilter`, `imsharpen`, …)
- Third-party on path: `readNPY` / `writeNPY` ([npy-matlab](https://github.com/kwikteam/npy-matlab)), `run_umap`, `play_bout`, `fit_exp_heaviside3`

**Python**
- `numpy`, `ssm` (HMM fitting), `torch` (CNN notebook)

### Running Figure 1

| Goal | Script |
|------|--------|
| Play-bout plots & fits | `Figure 1/Figure 1 Play Bouts.m` |
| LDA / UMAP figures | `Figure 1/Figure 1 LDA.m` |
| HMM figures | `Figure 1/Figure 1 HMM.m` |
| Regenerate HMM inputs (needs full raw data) | `Figure 1/HMM modeling/CreateHMMFIiles.m` |
| Fit HMMs in Python | `Figure 1/HMM modeling/python hmm estimate.txt` |

Pipeline order for HMM: see `Figure 1/HMM modeling/README.md`.

---

## Figure 2

Neural analyses around play bouts: band-limited LFP power PSTHs, area-resolved maps, behavior-resolved PSTHs, cross-frequency coupling, phase at exploratory-bout onset, speed–frequency relations, and call-locked PSTHs.

### Supporting codes v2 (recommended)

Use **`Figure 2/Supporting codes v2/`** for new work. Scripts follow a consistent naming convention:

| Prefix | Role | Example |
|--------|------|---------|
| `GENERATE_*` | Build per-session structs from raw/preprocessed NPX data | `GENERATE_PSTH_AREA_MAPS.m` |
| `Estimate_*` | Batch drivers: loop animals, call `GENERATE_*`, merge results | `Estimate_psth_area_maps.m` |
| `Analyze_*` | Load saved `.mat` files, statistics and figures | `Analyze_psth_area_maps.m` |
| `Compute_*` | Shared utilities (e.g. PLI/wPLI) | `Compute_pli_wpli.m` |

Only **`GENERATE_*`** filenames use all capitals; other prefixes use mixed case.

**Typical pipeline**

```
Estimate_*.m  →  GENERATE_*  →  save .mat  →  Analyze_*.m
```

| Analysis | Estimate (generate) | Analyze (plot) |
|----------|---------------------|----------------|
| Area-resolved play-bout PSTH maps | `Estimate_psth_area_maps.m` | `Analyze_psth_area_maps.m` |
| PSTH per behavior (self / partner) | `Estimate_psth_all_behaviors.m` | `Analyze_psth_all_behaviors.m` |
| Full-spectrum PSTH | `Estimate_psth_full_spectrogram.m` | `Analyze_psth_full_spectrogram.m` |
| Single-channel band power PSTH | `Estimate_psth_band_power.m` | — |
| PSTH + kinematic/call regressors | `Estimate_psth_all_regressors.m` | — |
| Call-locked PSTH | `Estimate_psth_calls_all_animals.m` | (includes plots in same script) |
| Cross-frequency coupling (mid-PAG) | `Estimate_freq_coupling.m` | `Analyze_freq_coupling.m` |
| Cross-frequency coupling maps | `Estimate_freq_coupling_maps.m` | `Analyze_freq_coupling_maps.m` |
| Phase PSTH (band/behavior via GENERATE) | `Estimate_phase_psth.m` | `Analyze_phase_psth.m` |
| Speed vs band-limited power | `Estimate_speed_freq_relation.m` | `Analyze_speed_freq_relation.m` |
| Coherence / PLI maps | `Estimate_coherence_maps.m` | — |

Main panel script: **`Figure 2/Figure 2.m`** (loads precomputed PSTH structs and assembles publication figures).

### Supporting codes (legacy)

**`Figure 2/Supporting codes/`** holds the original scripts (mixed naming, typos such as `Estiamte_*` / `Abalyze_*`). These are **not modified**; use v2 for cleaner names and simplified code.

### Figure 2 requirements

**MATLAB**
- Signal Processing Toolbox (`spectrogram`, `filtfilt`, `hilbert`, `designfilt`)
- Statistics and Machine Learning Toolbox (`fitlme`, `fitlm`, `ttest`, …)
- Third-party on path: `play_bout`, `wrap_probe_values`, `generateDistinctColors`, `align_by_area`, circular statistics helpers as used in the original pipeline

**Data** (not in git): NPX raw/preprocessed LFP, behavior backups, synch models, channel maps, area limits table, and precomputed results under `Analysis results/psth power by frequency and behavior/`.

---

## Setup after cloning

1. Clone this repository.
2. Download the data bundle(s) (see [Data availability](#data-availability)).
3. Extract so the layout matches:

```
play bout analysis behavior and physiology/
├── Figure 1/
├── Figure 2/
└── Data/                    ← downloaded archive
    ├── Behavior backups/
    ├── Analysis results/
    ├── HMM data/
    ├── Synch data/
    ├── CallDetectionBackup/
    └── …
```

4. **Figure 1** scripts are configured for:
   `...\play bout analysis behavior and physiology\Data\`

5. **Figure 2** v2 scripts currently point to lab network paths under `DataSets\` on the BCCN share. To run locally, update paths in the relevant `Estimate_*` / `Analyze_*` scripts or add a shared path config.

6. Add external MATLAB functions to your path (`play_bout`, `readNPY`, `run_umap`, `wrap_probe_values`, etc.).

7. Add **`Figure 2/Supporting codes v2`** to the MATLAB path before running Figure 2 estimate/analyze scripts.

---

## Data availability

**Data are not stored in this GitHub repository** (too large for git; see `.gitignore`).

### Figure 1 data bundle (minimum for plotting)

Copied subset used by Figure 1 scripts (~141 files):

- `Behavior backups/` — manual behavior annotations (`.txt`)
- `Analysis results/locomotive behaviors 2 partners/` — `all_behavior.mat`, etc.
- `Analysis results/Behavior classification/` — trained LDA models
- `Analysis results/HMM 2 and 3 states 2 partners/` — HMM analysis structs
- `HMM data/HMM raw data/` — HMM states, model comparison, `PropAndTime.mat`
- `Synch data/B1D1 1013 Dual/` — example session sync model (Fig 1B)
- `CallDetectionBackup/` — example `.wav` + `_Stats.xlsx` (Fig 1B)

Full file list: see `Data/COPY_MANIFEST.txt` in the archived bundle.

### Figure 2 data (local copy)

Precomputed PSTH, mixed-model results, and CV structs for Figure 2 live in:

`Data/Analysis results/psth power by frequency and behavior/`

File combinations (frequency band × behavior / calls / CV) are documented in:

`Figure 2/Figure 2 Psth animal names and result combinations.txt`

See `Data/Analysis results/psth power by frequency and behavior/COPY_MANIFEST.txt` for the full file list.

### Where to download data

> **Update these links before leaving the lab:**

| Archive | Link | Notes |
|---------|------|--------|
| Zenodo (recommended) | `[DOI or URL]` | Permanent DOI; good for publications |
| OSF | `[URL]` | Alternative open archive |
| Lab / institutional storage | `[URL]` | If Zenodo/OSF not used |

**Suggested Zenodo record titles:**
- `Play Bout Analysis — Behavior and Physiology (Figure 1 data: behavior, HMM, LDA)`
- `Play Bout Analysis — Behavior and Physiology (Figure 2 data: PSTH power by frequency and behavior)`

### Full raw data (optional, for re-running pipelines from scratch)

Not included in minimal bundles. If archived separately:

- `NPX data/NPX raw data/` — session folders with LFP and channel maps
- `Traking backups/` — `traking_structure` `.mat` per session
- `Synch data/` — all animals
- `CallDetectionBackup/` — all sessions

Document links in the same Zenodo record (additional zip) or in lab handover notes.

---

## Citation

If you use this code, cite the associated publication and the data DOI once registered.

## Contact

> **Update before archiving:** `[your email]`
