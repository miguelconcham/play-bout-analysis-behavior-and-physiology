# Theta-power coupling

This folder contains the Figure 2 cross-frequency coupling analyses.

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
using the legacy Estimate/GENERATE scripts.
