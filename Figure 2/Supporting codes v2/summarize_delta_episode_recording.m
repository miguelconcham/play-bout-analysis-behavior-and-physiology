function summary = summarize_delta_episode_recording(recording)
%SUMMARIZE_DELTA_EPISODE_RECORDING  Scalar description of one mid-LPAG channel.
%
% Durations, cycle counts, median periods, and period variability are
% summarized across episodes by their median. Variability of an episode
% is the interquartile range of its peak-to-peak periods. Pooled cycle
% statistics use every peak-to-peak interval from the kept episodes.
%
% Trough-to-trough fields are the Cole & Voytek cycle measures. They are
% left NaN when an older episode struct does not carry them. Period CV,
% amplitude CV, and symmetry are continuous summaries, not cutoffs.

episodes = recording.episodes;
if isempty(recording.intervals)
    analyzed_seconds = 0;
    n_runs = 0;
    median_threshold = NaN;
else
    analyzed_seconds = sum([recording.intervals.duration_s]);
    n_runs = sum([recording.intervals.n_threshold_runs]);
    median_threshold = median([recording.intervals.amplitude_threshold]);
end
n_episodes = numel(episodes);

summary = struct();
summary.animal_name = recording.animal_name;
summary.probe_number = recording.lpag_probe_number;
summary.lpag_channel = recording.lpag_channel;
summary.lpag_x_um = recording.lpag_x_um;
summary.lpag_y_um = recording.lpag_y_um;
summary.depth_start = recording.lpag_depth_range(1);
summary.depth_end = recording.lpag_depth_range(2);
summary.interval_mode = recording.interval_mode;
summary.n_intervals = numel(recording.intervals);
summary.analyzed_seconds = analyzed_seconds;
summary.median_amplitude_threshold = median_threshold;
summary.n_threshold_runs = n_runs;
summary.n_episodes = n_episodes;
summary.n_rejected_runs = n_runs - n_episodes;

summary.median_duration_s = NaN;
summary.iqr_duration_s = NaN;
summary.median_n_cycles = NaN;
summary.median_cycle_period_s = NaN;
summary.median_cycle_period_iqr_s = NaN;
summary.median_cycle_period_cv = NaN;
summary.median_fraction_cycles_in_delta_band = NaN;
summary.pooled_n_cycles = 0;
summary.pooled_median_cycle_period_s = NaN;
summary.pooled_fraction_cycles_in_delta_band = NaN;
summary.episode_time_fraction = episode_fraction(0, analyzed_seconds);
summary.median_n_trough_cycles = NaN;
summary.median_trough_period_s = NaN;
summary.median_trough_period_cv = NaN;
summary.median_cycle_amplitude = NaN;
summary.median_cycle_amplitude_cv = NaN;
summary.median_rdsym = NaN;
summary.median_ptsym = NaN;
summary.median_monotonicity = NaN;
summary.median_amplitude_consistency = NaN;
summary.median_period_consistency = NaN;
summary.n_episodes_with_sustained_sequence = NaN;
summary.pooled_n_trough_cycles = NaN;
summary.pooled_n_sustained_cycles = NaN;
summary.pooled_median_trough_period_s = NaN;
summary.pooled_median_rdsym = NaN;
summary.pooled_median_ptsym = NaN;

if n_episodes == 0
    summary = add_trough_cycle_summary(summary, recording);
    return
end

durations = [episodes.duration_s]';
n_cycles = [episodes.n_cycles]';
median_periods = [episodes.median_cycle_period_s]';
period_iqr = [episodes.cycle_period_iqr_s]';
period_cv = [episodes.cycle_period_cv]';
fraction_in_band = [episodes.fraction_cycles_in_delta_band]';

summary.median_duration_s = median(durations);
summary.median_n_cycles = median(n_cycles);
summary.median_cycle_period_s = median(median_periods);
summary.median_cycle_period_iqr_s = median_finite(period_iqr);
summary.median_cycle_period_cv = median_finite(period_cv);
summary.median_fraction_cycles_in_delta_band = median_finite(fraction_in_band);
summary.episode_time_fraction = episode_fraction(sum(durations), analyzed_seconds);
if n_episodes >= 2
    summary.iqr_duration_s = linear_iqr(durations);
end

periods = zeros(0, 1);
for episode_n = 1:n_episodes
    periods = [periods; episodes(episode_n).cycle_periods_s(:)]; %#ok<AGROW>
end
summary.pooled_n_cycles = numel(periods);
if ~isempty(periods)
    summary.pooled_median_cycle_period_s = median(periods);
    band = [1 / recording.parameters.delta_freq(2), ...
        1 / recording.parameters.delta_freq(1)];
    summary.pooled_fraction_cycles_in_delta_band = ...
        mean(periods >= band(1) & periods <= band(2));
end
summary = add_trough_cycle_summary(summary, recording);
end

function summary = add_trough_cycle_summary(summary, recording)
episodes = recording.episodes;
if ~isempty(episodes) && isfield(episodes, 'n_trough_cycles')
    summary.median_n_trough_cycles = median([episodes.n_trough_cycles]');
    summary.median_trough_period_s = median_finite([episodes.median_trough_period_s]');
    summary.median_trough_period_cv = median_finite([episodes.trough_period_cv]');
    summary.median_cycle_amplitude = median_finite([episodes.median_cycle_amplitude]');
    summary.median_cycle_amplitude_cv = median_finite([episodes.cycle_amplitude_cv]');
    summary.median_rdsym = median_finite([episodes.median_rdsym]');
    summary.median_ptsym = median_finite([episodes.median_ptsym]');
    summary.median_monotonicity = median_finite([episodes.median_monotonicity]');
    summary.median_amplitude_consistency = ...
        median_finite([episodes.median_amplitude_consistency]');
    summary.median_period_consistency = ...
        median_finite([episodes.median_period_consistency]');
    summary.n_episodes_with_sustained_sequence = ...
        nnz([episodes.has_sustained_sequence]);
end
if ~isfield(recording, 'cycles')
    return
end
cycles = recording.cycles;
summary.pooled_n_trough_cycles = numel(cycles);
if isempty(cycles)
    summary.pooled_n_sustained_cycles = 0;
    return
end
summary.pooled_n_sustained_cycles = nnz([cycles.in_sustained_sequence]);
summary.pooled_median_trough_period_s = median([cycles.period_s]');
summary.pooled_median_rdsym = median_finite([cycles.rdsym]');
summary.pooled_median_ptsym = median_finite([cycles.ptsym]');
end

function value = median_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function fraction = episode_fraction(episode_seconds, analyzed_seconds)
if analyzed_seconds > 0
    fraction = episode_seconds / analyzed_seconds;
else
    fraction = NaN;
end
end
