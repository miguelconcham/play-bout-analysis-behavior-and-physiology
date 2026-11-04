# Theta-power coupling

This folder contains the Figure 2 theta-power oscillation and
cross-frequency coupling analyses.

- `Estimate_psth_all_behaviors_theta_gamma.m` and
  `GENERATE_PSTH_ALL_BEHAVIORS_BANDS.m`: estimate behavior-aligned theta
  power (plus the optional delta and gamma bands).
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

Run the repository-level `add_repo_paths.m` once per MATLAB session before
using the Estimate/GENERATE scripts.
