# Analysis catalog

Which script to run for each analysis in this repository. Paper-panel letters live in [`README.md`](README.md). This file lists **every analysis that exists**, including supporting and complementary work that is not a numbered main or supplementary figure.

Run `add_repo_paths.m` once per MATLAB session before Estimate / Analyze / GENERATE scripts.

**Pipeline convention**

| Prefix | What it is | You usually run |
|--------|------------|-----------------|
| `GENERATE_*` | Per-session function (one animal / one NPX folder) | Called by Estimate; do not start here |
| `Estimate_*` | Loop sessions, call GENERATE, save `.mat` | **First** (rebuild data) |
| `Analyze_*` or `Figure *.m` | Load `.mat`, statistics and figures | **Second** (plot / test) |

Figure 2 supporting code: use **`Figure 2/Supporting codes v2/`**. The older `Figure 2/Supporting codes/` copies are legacy.

---

## All analyses (run this)

| Analysis | Paper? | Rebuild (Estimate) | Plot / test (run this) |
|----------|--------|--------------------|------------------------|
| Play-bout intervals, duration, speed–call CCG | Fig 1g–h | — | `Figure 1/Figure 1 Play Bouts.m` |
| Motion–USV features for LDA / UMAP | Fig 1c–d, Supp 1 | `Figure 1/LDA analysis/Obtain_Feature_USV_space.m` | `Figure 1/Figure 1 LDA.m` |
| CNN on locomotion features | extra | `Obtain_Feature_USV_space.m` (features) | `Figure 1/LDA analysis/Convolutional classification.ipynb` |
| HMM inputs from speed + calls | Fig 1, Supp 2 | `Figure 1/HMM modeling/CreateHMMFIiles.m` (calls `Create_HMM_inputs.m`) | Python: `Figure 1/HMM modeling/python hmm estimate.txt` |
| HMM states, play score, engaged vs unengaged | Fig 1e–f, Supp 2 | after Python HMM | `Figure 1/HMM modeling/Load_HMM_outputs_and_analyze.m` then `Figure 1/Figure 1 HMM.m` |
| Play-bout LFP band-power PSTH (delta / theta / gamma) | Fig 2c–f, Supp 4a,d | `Figure 2/Supporting codes v2/Estimate_psth_band_power.m` | `Figure 2/Figure 2.m` |
| Unique \(R^2\) / GLM of band power | Fig 2g, Supp 4c,f | `Figure 2/Supporting codes v2/Estimate_psth_all_regressors.m` | `Figure 2/Figure 2.m` (`%% 7`) |
| Reciprocal vs non-reciprocal play vs delta | Fig 2h, Supp 3a–f | `Estimate_psth_all_regressors.m` | `Figure 2/Reciprocity.m` |
| First vs second play partner vs delta | Supp 3g–i | play-bout PSTH (`Estimate_psth_band_power.m`) | `Figure 2/Play partner.m` |
| Delta PSTH per annotated behavior | extra (feeds Supp 3) | `Figure 2/Supporting codes v2/Estimate_psth_all_behaviors.m` | `Figure 2/Supporting codes v2/Analyze_psth_all_behaviors.m` |
| Call-locked LFP power PSTH | Supp 4e | `Figure 2/Supporting codes v2/Estimate_psth_calls_all_animals.m` | `Figure 2/Figure 2.m` (`%% 13`–`%% 17`); plots also inside Estimate |
| Speed vs LFP band power | Supp 4b; extra | `Figure 2/Supporting codes v2/Estimate_speed_freq_relation.m` | `Figure 2/Supporting codes v2/Analyze_speed_freq_relation.m`; also `Figure 2.m` `%% 18`–`%% 21` |
| Cross-frequency coupling (mid-PAG) | extra / CFC panel | `Figure 2/Supporting codes v2/Estimate_freq_coupling.m` | `Figure 2/Supporting codes v2/Analyze_freq_coupling.m`; `Figure 2/CrossFreq coupling.m` |
| CFC maps along the probe | extra | `Figure 2/Supporting codes v2/Estimate_freq_coupling_maps.m` | `Figure 2/Supporting codes v2/Analyze_freq_coupling_maps.m` |
| Inter-area PLI / wPLI / coherence maps | extra | `Figure 2/Supporting codes v2/Estimate_coherence_maps.m` | same script (Estimate also plots) |
| LFP phase reset at bout onset | extra / phase-reset panels | `Figure 2/Supporting codes v2/Estimate_phase_psth.m` | `Figure 2/Supporting codes v2/Analyze_phase_psth.m` |
| Phase reset (supplementary layout) | extra | `Estimate_phase_psth.m` | `Figure 2/Supporting codes v2/Analyze_phase_psth_SuppFig.m` |
| Play-bout PSTH maps by channel / area | extra | `Figure 2/Supporting codes v2/Estimate_psth_area_maps.m` | `Figure 2/Supporting codes v2/Analyze_psth_area_maps.m` |
| Full-spectrum spectrogram at play onset | extra | `Figure 2/Supporting codes v2/Estimate_psth_full_spectrogram.m` | `Figure 2/Supporting codes v2/Analyze_psth_full_spectrogram.m` |
| Neuron–LFP phase locking (PPC, preferred phase, PSTH) | Fig 3a–e | `Figure 3/Supporting codes/Estiamte_phase_copuling.m` | `Figure 3/Figure 3a-e.m` |
| Play vs non-play delta peak/trough PSTH | extra | after phase-coupling Estimate | `Figure 3/Supporting codes/LOAD_AND_PLOT_psth_deltapeak_play_noplay.m` |
| Activation index and bout-aligned rates | Fig 3f–h, Supp 5a–d | `Figure 3/Supporting codes/Estiamte_activation_index_and_other_variables.m` | `Figure 3/Figure 3g-h Activation Index.m`; example: `Figure 3/SuppFig4 ActivationIndex Example.m` |
| Behavior-resolved rate modulation (A1–A7) | extra / extended | `Estiamte_activation_index_and_other_variables.m` | `Figure 3/Fig A1-A7 behavior modulation.m` |
| Activation index vs delta power | extra (former supp) | activation-index Estimate | `Figure 3/Former Supplementary Figure Correlation between power and activation index.m` |
| Pairwise spike cross-correlograms | Fig 3i, Supp 5f–g | `Figure 3/Supporting codes/Estimate cross correlogram all aniamls.m` | `Figure 3/Supporting codes/ANALYZE cross correlogram all aniamls.m`; `Figure 3/Figure 3i cross correlogram.m` |
| Coincidence / surprise around lag 0 | Fig 3k, Supp 5h | `Figure 3/Supporting codes/Estimate suprise statistics  all aniamls.m` | `Figure 3/Figure 3k Coincident events.m`; `Figure 3/Supporting codes/ANALYSE suprise statistics  all aniamls v2.m` |
| LFP phase of coincident spikes | Fig 3j | `Figure 3/Supporting codes/Estiamte_phase_coincidence.m` | `Figure 3/Supporting codes/ANALYZE_phase_coincidence.m`; `Figure 3/Fig 3j.m` |
| Coincidence-phase extra plots | extra (former supp) | `Estiamte_phase_coincidence.m` | `Figure 3/Former Supplementary Figure Phase of coincident events.m` |
| Rate + coincidence phase (shift-o-gram) | extra | `Figure 3/Supporting codes/Estimate rate coincidence phase.m` | same script (Estimate also plots) |
| Dual-animal delta mutual information | Fig 6b–e,g–h,i | `Figure 6/Supplementary codes/Estiamte animal aynch all animals.m` | `Figure 6/Figure 6 Basic MI.m`; `Figure 6/Fig 6i.m` |
| Dual-animal LFP cross-correlation | Fig 6f | `Figure 6/Supplementary codes/Estiamte animal aynch all animals_cross_correlograms.m` | `Figure 6/Fig 6f.m` |
| Acute NPX: LFP × breathing phase locking | extra | `Complementary figures/Estimate_acute_phase_locking.m` | `Complementary figures/Analyze_acute_phase_locking.m` |
| Acute NPX: call and green-window spike PSTHs | extra | `Complementary figures/Estimate_stim_call_responses.m` | outputs in `stim_call_responses_*.mat` (no Analyze yet) |
| Firing-run parsing, per-run PPC, run ACGs | extra | `Complementary figures/Estimate_spike_train_parsing.m` | `Complementary figures/Analyze_spike_train_parsing.m` |
| LIF network with shared vs private oscillation | extra | `Complementary figures/LIF neuron model.ipynb` (Estimate grid cells) | same notebook (Analyze cells) |

---

## Extra analyses (not a numbered paper panel)

These are the supporting / complementary pipelines. Rebuild with Estimate, then run Analyze.

### LFP power and behavior (Figure 2 supporting)

| What | Estimate | Analyze |
|------|----------|---------|
| Delta PSTH for every annotated behavior | `Figure 2/Supporting codes v2/Estimate_psth_all_behaviors.m` | `Analyze_psth_all_behaviors.m` |
| Channel-wise play-bout maps (LPAG vs VLPAG, etc.) | `Estimate_psth_area_maps.m` | `Analyze_psth_area_maps.m` |
| Full spectrogram at play onset | `Estimate_psth_full_spectrogram.m` | `Analyze_psth_full_spectrogram.m` |
| Speed vs band-limited power | `Estimate_speed_freq_relation.m` | `Analyze_speed_freq_relation.m` |
| Delta–gamma (or other) phase–amplitude coupling | `Estimate_freq_coupling.m` | `Analyze_freq_coupling.m` or `Figure 2/CrossFreq coupling.m` |
| CFC as depth maps | `Estimate_freq_coupling_maps.m` | `Analyze_freq_coupling_maps.m` |
| PLI / wPLI coherence maps | `Estimate_coherence_maps.m` | same file |
| Hilbert phase at exploratory / play onset | `Estimate_phase_psth.m` | `Analyze_phase_psth.m` or `Analyze_phase_psth_SuppFig.m` |

GENERATE functions live next to those Estimates (`GENERATE_PSTH_PLAY_BOUT.m`, `GENERATE_PSTH_ALL_BEHAVIORS.m`, `GENERATE_PSTH_CALLS.m`, `GENERATE_PSTH_CALL_LOC_REGRESSORS.m`, `GENERATE_SPEED_FREQ_RELATION.m`, `GENERATE_FREQ_COUPLING_STRUCT.m`, `GENERATE_FREQ_COUPLING_MAPS_STRUCT.m`, `GENERATE_COHERNECE_MAPS_STRUCT.m`, `GENERATE_PHASE_EXPLORATORY_ONSET.m`, `GENERATE_PSTH_AREA_MAPS.m`, `GENERATE_PSTH_SPECTROGRAM.m`).

### Spikes and locking (Figure 3 supporting)

| What | Estimate | Analyze |
|------|----------|---------|
| Play vs pre-play delta-aligned PSTHs | `Estiamte_phase_copuling.m` | `LOAD_AND_PLOT_psth_deltapeak_play_noplay.m` |
| A1–A7 behavior modulation of rates | `Estiamte_activation_index_and_other_variables.m` | `Figure 3/Fig A1-A7 behavior modulation.m` |
| Activation index vs LFP power | same Estimate | `Former Supplementary Figure Correlation between power and activation index.m` |
| Population CCGs by area and lock type | `Estimate cross correlogram all aniamls.m` | `ANALYZE cross correlogram all aniamls.m` |
| Surprise / coincidence time course | `Estimate suprise  all aniamls.m` or `Estimate suprise statistics  all aniamls.m` | `ANALYSE suprise statistics  all aniamls v2.m` |
| Phase of coincident spikes (extended) | `Estiamte_phase_coincidence.m` | `ANALYZE_phase_coincidence.m`; former plots: `Former Supplementary Figure Phase of coincident events.m` |
| Rate-coincidence phase | `Estimate rate coincidence phase.m` | same file |

BC / dithering variants (same question, different shuffle): `Estimate suprise statistics  all aniamls BC.m`, `GENERATE_SURPIRSE_DYNAMICS_STATS_DITTERING.m`.

### Complementary figures (acute, parsing, LIF)

| What | Estimate | Analyze |
|------|----------|---------|
| Acute green-window LFP + breathing lock, 2×2 PSTHs, region split (SC, PAG, Raphe) | `Complementary figures/Estimate_acute_phase_locking.m` | `Complementary figures/Analyze_acute_phase_locking.m` |
| Acute call PSTHs (drop onsets within 100 ms after StimEnd) and green-window PSTHs; session-rate z-score with StimEnd NaNs | `Complementary figures/Estimate_stim_call_responses.m` | no Analyze yet; inspect `Data/Analysis results/acute phase locking/stim_call_responses_*.mat` |
| Firing runs, per-run PPC, edge-corrected autocorrelograms (locked low/high PPC vs unlocked) | `Complementary figures/Estimate_spike_train_parsing.m` | `Complementary figures/Analyze_spike_train_parsing.m` |
| LIF with shared vs independent oscillatory drive | notebook Estimate grid | `Complementary figures/LIF neuron model.ipynb` |

Older one-file versions (do not start here unless you need the original session walkthrough): `Complementary figures/Phase_locking in acute experiments.m`, `Complementary figures/spike_train_parsing_phase_lock.m`.

---

## Figure 1 supporting (not a panel by itself)

| What | Run |
|------|-----|
| Feature table for LDA | `Figure 1/LDA analysis/Obtain_Feature_USV_space.m` (uses `get_behavior_properties.m`) |
| Conv-net behavior classifier | `Figure 1/LDA analysis/Convolutional classification.ipynb` |
| HMM input files | `Figure 1/HMM modeling/CreateHMMFIiles.m` |
| HMM fit | `Figure 1/HMM modeling/python hmm estimate.txt` (see `Figure 1/HMM modeling/README.md`) |
| HMM post-processing | `Figure 1/HMM modeling/Load_HMM_outputs_and_analyze.m` |

---

## Dual-animal (Figure 6 supporting)

| What | Estimate | Plot |
|------|----------|------|
| Per-bout and time-resolved MI | `Figure 6/Supplementary codes/Estiamte animal aynch all animals.m` | `Figure 6/Figure 6 Basic MI.m`, `Fig 6i.m` |
| Two-animal LFP CCG | `Estiamte animal aynch all animals_cross_correlograms.m` | `Figure 6/Fig 6f.m` |

GENERATE: `GENERATE_STRUCUTRE_animal_synch.m`, `GENERATE_STRUCUTRE_animal_synch_cross_correlogram.m`.
