function summary_table = delta_episode_summary_table(recordings)
%DELTA_EPISODE_SUMMARY_TABLE  One row per mid-LPAG channel.

n_recordings = numel(recordings);
animal_name = cell(n_recordings, 1);
interval_mode = cell(n_recordings, 1);
probe_number = nan(n_recordings, 1);
lpag_channel = nan(n_recordings, 1);
lpag_x_um = nan(n_recordings, 1);
lpag_y_um = nan(n_recordings, 1);
depth_start = nan(n_recordings, 1);
depth_end = nan(n_recordings, 1);
n_intervals = nan(n_recordings, 1);
analyzed_seconds = nan(n_recordings, 1);
median_amplitude_threshold = nan(n_recordings, 1);
n_threshold_runs = nan(n_recordings, 1);
n_episodes = nan(n_recordings, 1);
n_rejected_runs = nan(n_recordings, 1);
episode_time_fraction = nan(n_recordings, 1);
median_duration_s = nan(n_recordings, 1);
iqr_duration_s = nan(n_recordings, 1);
median_n_cycles = nan(n_recordings, 1);
median_cycle_period_s = nan(n_recordings, 1);
median_cycle_period_iqr_s = nan(n_recordings, 1);
median_cycle_period_cv = nan(n_recordings, 1);
median_fraction_cycles_in_delta_band = nan(n_recordings, 1);
pooled_n_cycles = nan(n_recordings, 1);
pooled_median_cycle_period_s = nan(n_recordings, 1);
pooled_fraction_cycles_in_delta_band = nan(n_recordings, 1);
median_n_trough_cycles = nan(n_recordings, 1);
median_trough_period_s = nan(n_recordings, 1);
median_trough_period_cv = nan(n_recordings, 1);
median_cycle_amplitude = nan(n_recordings, 1);
median_cycle_amplitude_cv = nan(n_recordings, 1);
median_rdsym = nan(n_recordings, 1);
median_ptsym = nan(n_recordings, 1);
median_monotonicity = nan(n_recordings, 1);
median_amplitude_consistency = nan(n_recordings, 1);
median_period_consistency = nan(n_recordings, 1);
n_episodes_with_sustained_sequence = nan(n_recordings, 1);
pooled_n_trough_cycles = nan(n_recordings, 1);
pooled_n_sustained_cycles = nan(n_recordings, 1);
pooled_median_trough_period_s = nan(n_recordings, 1);
pooled_median_rdsym = nan(n_recordings, 1);
pooled_median_ptsym = nan(n_recordings, 1);

for recording_n = 1:n_recordings
    summary = summarize_delta_episode_recording(recordings(recording_n));
    animal_name{recording_n} = summary.animal_name;
    interval_mode{recording_n} = summary.interval_mode;
    probe_number(recording_n) = summary.probe_number;
    lpag_channel(recording_n) = summary.lpag_channel;
    lpag_x_um(recording_n) = summary.lpag_x_um;
    lpag_y_um(recording_n) = summary.lpag_y_um;
    depth_start(recording_n) = summary.depth_start;
    depth_end(recording_n) = summary.depth_end;
    n_intervals(recording_n) = summary.n_intervals;
    analyzed_seconds(recording_n) = summary.analyzed_seconds;
    median_amplitude_threshold(recording_n) = summary.median_amplitude_threshold;
    n_threshold_runs(recording_n) = summary.n_threshold_runs;
    n_episodes(recording_n) = summary.n_episodes;
    n_rejected_runs(recording_n) = summary.n_rejected_runs;
    episode_time_fraction(recording_n) = summary.episode_time_fraction;
    median_duration_s(recording_n) = summary.median_duration_s;
    iqr_duration_s(recording_n) = summary.iqr_duration_s;
    median_n_cycles(recording_n) = summary.median_n_cycles;
    median_cycle_period_s(recording_n) = summary.median_cycle_period_s;
    median_cycle_period_iqr_s(recording_n) = summary.median_cycle_period_iqr_s;
    median_cycle_period_cv(recording_n) = summary.median_cycle_period_cv;
    median_fraction_cycles_in_delta_band(recording_n) = ...
        summary.median_fraction_cycles_in_delta_band;
    pooled_n_cycles(recording_n) = summary.pooled_n_cycles;
    pooled_median_cycle_period_s(recording_n) = summary.pooled_median_cycle_period_s;
    pooled_fraction_cycles_in_delta_band(recording_n) = ...
        summary.pooled_fraction_cycles_in_delta_band;
    median_n_trough_cycles(recording_n) = summary.median_n_trough_cycles;
    median_trough_period_s(recording_n) = summary.median_trough_period_s;
    median_trough_period_cv(recording_n) = summary.median_trough_period_cv;
    median_cycle_amplitude(recording_n) = summary.median_cycle_amplitude;
    median_cycle_amplitude_cv(recording_n) = summary.median_cycle_amplitude_cv;
    median_rdsym(recording_n) = summary.median_rdsym;
    median_ptsym(recording_n) = summary.median_ptsym;
    median_monotonicity(recording_n) = summary.median_monotonicity;
    median_amplitude_consistency(recording_n) = summary.median_amplitude_consistency;
    median_period_consistency(recording_n) = summary.median_period_consistency;
    n_episodes_with_sustained_sequence(recording_n) = ...
        summary.n_episodes_with_sustained_sequence;
    pooled_n_trough_cycles(recording_n) = summary.pooled_n_trough_cycles;
    pooled_n_sustained_cycles(recording_n) = summary.pooled_n_sustained_cycles;
    pooled_median_trough_period_s(recording_n) = summary.pooled_median_trough_period_s;
    pooled_median_rdsym(recording_n) = summary.pooled_median_rdsym;
    pooled_median_ptsym(recording_n) = summary.pooled_median_ptsym;
end

summary_table = table(animal_name, probe_number, lpag_channel, lpag_x_um, ...
    lpag_y_um, depth_start, depth_end, interval_mode, n_intervals, ...
    analyzed_seconds, median_amplitude_threshold, n_threshold_runs, ...
    n_episodes, n_rejected_runs, episode_time_fraction, median_duration_s, ...
    iqr_duration_s, median_n_cycles, median_cycle_period_s, ...
    median_cycle_period_iqr_s, median_cycle_period_cv, ...
    median_fraction_cycles_in_delta_band, pooled_n_cycles, ...
    pooled_median_cycle_period_s, pooled_fraction_cycles_in_delta_band, ...
    median_n_trough_cycles, median_trough_period_s, median_trough_period_cv, ...
    median_cycle_amplitude, median_cycle_amplitude_cv, median_rdsym, ...
    median_ptsym, median_monotonicity, median_amplitude_consistency, ...
    median_period_consistency, n_episodes_with_sustained_sequence, ...
    pooled_n_trough_cycles, pooled_n_sustained_cycles, ...
    pooled_median_trough_period_s, pooled_median_rdsym, pooled_median_ptsym);
end
