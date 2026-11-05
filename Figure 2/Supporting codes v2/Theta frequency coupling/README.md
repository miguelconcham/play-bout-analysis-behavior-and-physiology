# Theta frequency coupling

Behavior-band theta/gamma power and cross-frequency coupling for Figure 2.

- `Estimate_psth_all_behaviors_theta_gamma.m` and
  `GENERATE_PSTH_ALL_BEHAVIORS_BANDS.m`: estimate behavior-aligned theta
  power (plus the optional delta and gamma bands).
- `Analyze_psth_all_behaviors_theta_gamma.m`: grand averages and condition
  differences for those PSTHs.
- `Analyze_psth_preonset_spectrum.m`: test the 2–3 Hz oscillation of theta
  power before behavior onset, including play/non-play summaries,
  animal-level checks, and trial-wise circular-shift nulls.
- `Estimate_frequency_resolved_theta_pac.m`: primary frequency-resolved
  slow-phase/theta-power analysis, independent block-shift nulls, harmonic
  controls, play/non-play interaction, and hierarchical inference.
- `Estimate_freq_coupling.m`, `GENERATE_FREQ_COUPLING_STRUCT.m`, and
  `Analyze_freq_coupling.m`: legacy mid-PAG phase-amplitude analysis.
- `Estimate_freq_coupling_maps.m`,
  `GENERATE_FREQ_COUPLING_MAPS_STRUCT.m`, and
  `Analyze_freq_coupling_maps.m`: legacy probe-depth maps.
- `CrossFreq coupling.m`: legacy combined analysis/figure script.

The theta-amplitude oscillation pipeline is in the sibling folder
`Theta envelope` (estimate, generate, analyze, and the synthetic test).
Its surrogate, group, and plotting helpers are in
`Theta envelope statistics`.

Run the repository-level `add_repo_paths.m` once per MATLAB session before
using the Estimate/GENERATE scripts. That call adds every Supporting codes
v2 theme folder, including `Helper functions`.
