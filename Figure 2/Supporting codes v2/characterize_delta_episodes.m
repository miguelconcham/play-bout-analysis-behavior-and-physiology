function detected = characterize_delta_episodes(filtered, amplitude, sampling_rate, options)
%CHARACTERIZE_DELTA_EPISODES  Episodes of elevated delta amplitude.
%
% filtered and amplitude are one continuous interval (a partner-session
% block, or a whole recording). amplitude is the Hilbert envelope of the
% 1-5 Hz signal, in the same samples as filtered.
%
% A candidate run is a stretch where the envelope is at or above the
% amplitude threshold. Positive peaks of the filtered signal inside that
% run give peak-to-peak cycle durations. A run is kept when it has at
% least min_cycles peak-to-peak intervals and is at least
% min_duration_seconds long. Trough-to-trough shape is added afterwards
% by characterize_cycle_features; the fields below stay peak-to-peak.
%
% detected.episodes uses sample indices into this interval. The caller
% maps them onto recording time.

if nargin < 4 || isempty(options)
    options = struct();
end
options = delta_episode_options(options);

filtered = double(filtered(:).');
amplitude = double(amplitude(:).');
if numel(filtered) ~= numel(amplitude)
    error('characterize_delta_episodes:size', ...
        'Filtered LFP and its amplitude must have the same length.');
end
if sampling_rate <= 0
    error('characterize_delta_episodes:rate', 'Sampling rate must be positive.');
end
n_samples = numel(filtered);

if options.amplitude_smooth_seconds > 0
    smooth_samples = max(1, round(options.amplitude_smooth_seconds * sampling_rate));
    amplitude_used = movmean(amplitude, smooth_samples);
else
    amplitude_used = amplitude;
end

if isempty(options.edge_margin_seconds)
    edge_margin_seconds = options.delta_filter_order / sampling_rate;
else
    edge_margin_seconds = options.edge_margin_seconds;
end
edge_margin_samples = round(edge_margin_seconds * sampling_rate);
eligible = true(1, n_samples);
if edge_margin_samples > 0
    eligible(1:min(edge_margin_samples, n_samples)) = false;
    eligible(max(1, n_samples - edge_margin_samples + 1):n_samples) = false;
end

detected = struct();
detected.threshold = NaN;
detected.amplitude_mean = NaN;
detected.amplitude_std = NaN;
detected.amplitude_median = NaN;
detected.edge_margin_seconds = edge_margin_seconds;
detected.edge_margin_samples = edge_margin_samples;
detected.n_eligible_samples = nnz(eligible);
detected.n_threshold_runs = 0;
detected.n_episodes = 0;
detected.episodes = empty_local_episodes();

if nnz(eligible) < 2 || ~any(isfinite(amplitude_used(eligible)))
    return
end

used = amplitude_used(eligible);
used = used(isfinite(used));
detected.amplitude_mean = mean(used);
detected.amplitude_std = std(used, 0);
detected.amplitude_median = median(used);
detected.threshold = amplitude_threshold(used, detected, options);

above = eligible & isfinite(amplitude_used) & amplitude_used >= detected.threshold;
runs = find_beg_end(double(above));
runs = merge_close_runs(runs, round(options.max_gap_seconds * sampling_rate));
detected.n_threshold_runs = size(runs, 1);
if isempty(runs)
    return
end

min_peak_distance = max(1, round(options.min_peak_distance_seconds * sampling_rate));
band_min = 1 / options.delta_freq(2);
band_max = 1 / options.delta_freq(1);
kept = empty_local_episodes();

for run_n = 1:size(runs, 1)
    episode = episode_from_run(filtered, runs(run_n, :), sampling_rate, ...
        min_peak_distance, band_min, band_max);
    if episode.duration_s >= options.min_duration_seconds && ...
            episode.n_cycles >= options.min_cycles
        if isempty(kept)
            kept = episode;
        else
            kept(end + 1) = episode; %#ok<AGROW>
        end
    end
end
detected.episodes = kept;
detected.n_episodes = numel(kept);
end

function threshold = amplitude_threshold(used, detected, options)
switch options.amplitude_threshold_mode
    case 'std'
        threshold = detected.amplitude_mean + ...
            options.amplitude_threshold * detected.amplitude_std;
    case 'median'
        threshold = options.amplitude_threshold * detected.amplitude_median;
    case 'absolute'
        threshold = options.amplitude_threshold;
    otherwise
        error('characterize_delta_episodes:threshold', ...
            'Unknown amplitude_threshold_mode.');
end
end

function episode = episode_from_run(filtered, run, sampling_rate, ...
        min_peak_distance, band_min, band_max)
start_sample = run(1);
end_sample = run(2);
segment = filtered(start_sample:end_sample);
peak_samples = positive_peaks(segment, min_peak_distance);
peak_samples = peak_samples + start_sample - 1;
cycle_periods = diff(peak_samples(:)) / sampling_rate;
n_cycles = numel(cycle_periods);

episode = struct();
episode.start_sample = start_sample;
episode.end_sample = end_sample;
episode.duration_s = (end_sample - start_sample + 1) / sampling_rate;
episode.n_peaks = numel(peak_samples);
episode.n_cycles = n_cycles;
episode.peak_samples = peak_samples(:);
episode.cycle_periods_s = cycle_periods(:);
if n_cycles == 0
    episode.median_cycle_period_s = NaN;
    episode.cycle_period_iqr_s = NaN;
    episode.cycle_period_std_s = NaN;
    episode.cycle_period_cv = NaN;
    episode.fraction_cycles_in_delta_band = NaN;
    return
end
episode.median_cycle_period_s = median(cycle_periods);
if n_cycles >= 2
    episode.cycle_period_iqr_s = linear_iqr(cycle_periods);
    episode.cycle_period_std_s = std(cycle_periods, 0);
else
    episode.cycle_period_iqr_s = NaN;
    episode.cycle_period_std_s = NaN;
end
period_mean = mean(cycle_periods);
if period_mean > 0 && n_cycles >= 2
    episode.cycle_period_cv = episode.cycle_period_std_s / period_mean;
else
    episode.cycle_period_cv = NaN;
end
in_band = cycle_periods >= band_min & cycle_periods <= band_max;
episode.fraction_cycles_in_delta_band = mean(in_band);
end

function peak_samples = positive_peaks(segment, min_peak_distance)
% Local maxima above zero. A shorter neighbor within min_peak_distance
% loses to the taller one, so one cycle is not counted twice.
segment = segment(:);
peak_samples = zeros(0, 1);
n_samples = numel(segment);
if n_samples < 3
    return
end
is_peak = segment(2:n_samples - 1) > segment(1:n_samples - 2) & ...
    segment(2:n_samples - 1) >= segment(3:n_samples) & ...
    segment(2:n_samples - 1) > 0;
candidates = find(is_peak) + 1;
if isempty(candidates)
    return
end
[~, order] = sort(segment(candidates), 'descend');
candidates = candidates(order);
blocked = false(n_samples, 1);
kept = zeros(numel(candidates), 1);
n_kept = 0;
for index = 1:numel(candidates)
    sample = candidates(index);
    if blocked(sample)
        continue
    end
    n_kept = n_kept + 1;
    kept(n_kept) = sample;
    left = max(1, sample - (min_peak_distance - 1));
    right = min(n_samples, sample + (min_peak_distance - 1));
    blocked(left:right) = true;
end
peak_samples = sort(kept(1:n_kept));
end

function runs = merge_close_runs(runs, max_gap_samples)
if isempty(runs) || max_gap_samples < 1
    return
end
merged = runs(1, :);
for run_n = 2:size(runs, 1)
    gap = runs(run_n, 1) - merged(end, 2) - 1;
    if gap <= max_gap_samples
        merged(end, 2) = runs(run_n, 2);
    else
        merged = [merged; runs(run_n, :)]; %#ok<AGROW>
    end
end
runs = merged;
end

function episodes = empty_local_episodes()
episodes = struct( ...
    'start_sample', {}, ...
    'end_sample', {}, ...
    'duration_s', {}, ...
    'n_peaks', {}, ...
    'n_cycles', {}, ...
    'peak_samples', {}, ...
    'cycle_periods_s', {}, ...
    'median_cycle_period_s', {}, ...
    'cycle_period_iqr_s', {}, ...
    'cycle_period_std_s', {}, ...
    'cycle_period_cv', {}, ...
    'fraction_cycles_in_delta_band', {});
end
