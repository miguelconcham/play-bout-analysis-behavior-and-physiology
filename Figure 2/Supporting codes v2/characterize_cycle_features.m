function [episodes, cycles] = characterize_cycle_features(broad, narrow, ...
        sampling_rate, episodes, options)
%CHARACTERIZE_CYCLE_FEATURES  Trough cycles inside detected delta episodes.
%
% Cole & Voytek (2019) J Neurophysiol 122:849-861. The narrowband signal
% supplies rise and decay zero-crossings. Peaks and troughs are the
% maximum and minimum of the low-pass (broad) signal between those
% crossings. Voltages, durations, and monotonicity are read on the broad
% signal, not on the 1-5 Hz trace.
%
% A cycle starts and ends at consecutive troughs and is kept only when
% exactly one peak lies between them. Both troughs must fall inside a
% candidate episode. A trough pair cut by the episode edge is incomplete
% and is omitted from the cycle statistics. Samples inside the filter
% edge margin are not used as extrema.
%
% Cycle period is the time between those troughs. Frequency is its
% inverse. Amplitude is the mean of the rise voltage (trough to peak)
% and the decay voltage (peak to next trough). Rise-decay symmetry is
% the rise duration divided by the period. Peak-trough symmetry is the
% peak duration divided by the peak duration plus the previous trough
% duration. The peak duration runs from the rise midpoint to the next
% decay midpoint. The trough duration runs from the previous decay
% midpoint to the rise midpoint. A midpoint is the median sample at
% which the broad signal crosses the voltage halfway between the two
% extrema (the paper's flank midpoint; bycycle find_zerox).
%
% Amplitude consistency is the smallest adjacent-flank ratio
% min(|voltage_a|, |voltage_b|) / max(|voltage_a|, |voltage_b|) among the
% current rise and decay, the previous decay versus the current rise,
% and the current decay versus the next rise. Period consistency is the
% smaller adjacent-period ratio against the previous and next cycles.
% The first and last complete cycles of the interval have no neighbor
% pair, so those consistencies are undefined and those cycles cannot
% meet the sequence test. Neighbors may sit just outside the episode;
% only cycles inside the episode are saved.
%
% Monotonicity is the fraction of sample-to-sample steps that rise on
% the rise flank and fall on the decay flank, pooled across both flanks.
% A flat step is not counted as correct. Bycycle instead averages the
% two flank fractions without weighting by length; the two agree when
% the flanks have the same number of steps.
%
% A cycle meets the published sequence test when its amplitude is at or
% above the chosen percentile of all complete cycles in the interval
% (no percentile on the simulation preset), and its amplitude
% consistency, period consistency, and monotonicity are at least the
% configured thresholds. At least min_consecutive_cycles such cycles
% must be consecutive inside one episode. Period CV, amplitude CV, and
% either symmetry stay continuous: they do not gate the sequence.
%
% The motor-cortex thresholds (0.3, 0.5, 0.6, 20th percentile) and the
% simulation thresholds (0.6, 0.6, 0.9, no percentile) were chosen by
% Cole and Voytek by visual inspection. The beta values are not a
% delta-specific recommendation. Bycycle's later defaults (0.5, 0.5,
% 0.8) are not used. Comparison is >= so a value sitting on a threshold
% passes. Equal amplitudes all pass the percentile rule; a strict
% "above" test would reject a constant-amplitude oscillation.

if nargin < 5 || isempty(options)
    options = struct();
end
options = delta_episode_options(options);

broad = double(broad(:));
narrow = double(narrow(:));
if numel(broad) ~= numel(narrow)
    error('characterize_cycle_features:size', ...
        'Broad and narrowband signals must have the same length.');
end
if sampling_rate <= 0
    error('characterize_cycle_features:rate', 'Sampling rate must be positive.');
end
if isempty(episodes)
    cycles = empty_cycles();
    return
end
episodes = episodes(:);

[peaks, troughs] = locate_extrema(broad, narrow, edge_limits(numel(broad), ...
    sampling_rate, options));
cycles = cycles_from_extrema(broad, peaks, troughs, sampling_rate);
cycles = add_consistency(cycles);
cutoff = amplitude_cutoff(cycles, options.cycle_amplitude_percentile);
cycles = add_criteria(cycles, cutoff, options);
cycles = assign_episodes(cycles, episodes);
cycles = mark_sustained_sequences(cycles, options.min_consecutive_cycles);
episodes = summarize_episodes(episodes, cycles, cutoff);
cycles = cycles_inside_episodes(cycles);
end

function limits = edge_limits(n_samples, sampling_rate, options)
if isempty(options.edge_margin_seconds)
    edge_margin_seconds = options.delta_filter_order / sampling_rate;
else
    edge_margin_seconds = options.edge_margin_seconds;
end
edge_samples = round(edge_margin_seconds * sampling_rate);
limits = [min(n_samples, edge_samples + 1), max(1, n_samples - edge_samples)];
end

function [peaks, troughs] = locate_extrema(broad, narrow, limits)
rise_xs = flank_crossings(narrow, 0, 'rise');
decay_xs = flank_crossings(narrow, 0, 'decay');
peaks = extrema_between(broad, rise_xs, decay_xs, 'peak');
troughs = extrema_between(broad, decay_xs, rise_xs, 'trough');
peaks = peaks(peaks >= limits(1) & peaks <= limits(2));
troughs = troughs(troughs >= limits(1) & troughs <= limits(2));
end

function indexes = flank_crossings(signal, midpoint, flank)
if strcmp(flank, 'rise')
    before_cross = signal <= midpoint;
else
    before_cross = signal > midpoint;
end
indexes = find(before_cross(1:end - 1) & ~before_cross(2:end));
indexes = indexes(:);
end

function extrema = extrema_between(broad, start_xs, end_xs, kind)
extrema = zeros(0, 1);
end_cursor = 1;
n_ends = numel(end_xs);
for start_n = 1:numel(start_xs)
    while end_cursor <= n_ends && end_xs(end_cursor) <= start_xs(start_n)
        end_cursor = end_cursor + 1;
    end
    if end_cursor > n_ends
        break
    end
    segment = start_xs(start_n):(end_xs(end_cursor) - 1);
    if isempty(segment)
        continue
    end
    if strcmp(kind, 'peak')
        [~, relative] = max(broad(segment));
    else
        [~, relative] = min(broad(segment));
    end
    extrema(end + 1, 1) = segment(relative); %#ok<AGROW>
end
extrema = unique(extrema);
end

function cycles = cycles_from_extrema(broad, peaks, troughs, sampling_rate)
cycles = empty_cycles();
if numel(troughs) < 2 || isempty(peaks)
    return
end
peaks = sort(peaks(:));
troughs = sort(troughs(:));
found = cell(numel(troughs) - 1, 1);
n_found = 0;
for trough_n = 1:numel(troughs) - 1
    trough_start = troughs(trough_n);
    trough_end = troughs(trough_n + 1);
    between = peaks(peaks > trough_start & peaks < trough_end);
    if numel(between) ~= 1
        continue
    end
    peak = between(1);
    if ~all(isfinite(broad([trough_start, peak, trough_end])))
        continue
    end
    previous_peak = peaks(peaks < trough_start);
    if isempty(previous_peak)
        previous_peak = NaN;
    else
        previous_peak = previous_peak(end);
    end
    n_found = n_found + 1;
    found{n_found} = one_cycle(broad, trough_start, peak, trough_end, ...
        previous_peak, sampling_rate);
end
if n_found == 0
    return
end
cycles = vertcat(found{1:n_found});
end

function cycle = one_cycle(broad, trough_start, peak, trough_end, ...
        previous_peak, sampling_rate)
cycle = blank_cycle();
cycle.start_sample = trough_start;
cycle.end_sample = trough_end;
cycle.peak_sample = peak;
cycle.period_s = (trough_end - trough_start) / sampling_rate;
if cycle.period_s > 0
    cycle.frequency_hz = 1 / cycle.period_s;
end
cycle.rise_duration_s = (peak - trough_start) / sampling_rate;
cycle.decay_duration_s = (trough_end - peak) / sampling_rate;
cycle.rise_amplitude = broad(peak) - broad(trough_start);
cycle.decay_amplitude = broad(peak) - broad(trough_end);
cycle.amplitude = (cycle.rise_amplitude + cycle.decay_amplitude) / 2;
cycle.rdsym = cycle.rise_duration_s / cycle.period_s;
cycle.monotonicity = flank_monotonicity(broad, trough_start, peak, trough_end);

cycle.rise_midpoint_sample = flank_midpoint(broad, trough_start, peak, 'rise');
cycle.decay_midpoint_sample = flank_midpoint(broad, peak, trough_end, 'decay');
if isfinite(previous_peak)
    cycle.previous_decay_midpoint_sample = flank_midpoint( ...
        broad, previous_peak, trough_start, 'decay');
end
peak_samples = cycle.decay_midpoint_sample - cycle.rise_midpoint_sample;
trough_samples = cycle.rise_midpoint_sample - cycle.previous_decay_midpoint_sample;
if isfinite(peak_samples) && isfinite(trough_samples) && ...
        (peak_samples + trough_samples) ~= 0
    cycle.peak_duration_s = peak_samples / sampling_rate;
    cycle.trough_duration_s = trough_samples / sampling_rate;
    cycle.ptsym = peak_samples / (peak_samples + trough_samples);
end
end

function sample = flank_midpoint(signal, start_sample, end_sample, flank)
segment = signal(start_sample:end_sample);
if numel(segment) < 2 || all(segment == 0) || flank_is_reversed(segment, flank)
    sample = start_sample + floor(numel(segment) / 2);
    return
end
midpoint = (segment(1) + segment(end)) / 2;
crossings = flank_crossings(segment, midpoint, flank);
if isempty(crossings)
    sample = start_sample + floor(numel(segment) / 2);
    return
end
sample = start_sample + floor(median(crossings)) - 1;
end

function reversed = flank_is_reversed(segment, flank)
if strcmp(flank, 'rise')
    reversed = segment(1) > segment(end);
else
    reversed = segment(1) < segment(end);
end
end

function fraction = flank_monotonicity(signal, trough_start, peak, trough_end)
rise_steps = diff(signal(trough_start:peak));
decay_steps = diff(signal(peak:trough_end));
n_steps = numel(rise_steps) + numel(decay_steps);
if n_steps == 0
    fraction = NaN;
    return
end
n_correct = nnz(rise_steps > 0) + nnz(decay_steps < 0);
fraction = n_correct / n_steps;
end

function cycles = add_consistency(cycles)
n_cycles = numel(cycles);
if n_cycles == 0
    return
end
rises = [cycles.rise_amplitude]';
decays = [cycles.decay_amplitude]';
periods = [cycles.period_s]';
for cycle_n = 1:n_cycles
    if cycle_n == 1 || cycle_n == n_cycles
        cycles(cycle_n).amplitude_consistency = NaN;
        cycles(cycle_n).period_consistency = NaN;
        continue
    end
    current_ratio = positive_ratio(rises(cycle_n), decays(cycle_n));
    previous_ratio = positive_ratio(rises(cycle_n), decays(cycle_n - 1));
    next_ratio = positive_ratio(rises(cycle_n + 1), decays(cycle_n));
    cycles(cycle_n).amplitude_consistency = min_finite([current_ratio, ...
        previous_ratio, next_ratio]);
    previous_period = positive_ratio(periods(cycle_n), periods(cycle_n - 1));
    next_period = positive_ratio(periods(cycle_n + 1), periods(cycle_n));
    cycles(cycle_n).period_consistency = min(previous_period, next_period);
end
end

function value = min_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = min(x);
end
end

function ratio = positive_ratio(left, right)
% Absolute voltage changes. "Twice as big" does not depend on sign.
if ~all(isfinite([left, right]))
    ratio = NaN;
    return
end
left = abs(left);
right = abs(right);
largest = max(left, right);
if largest == 0
    ratio = NaN;
    return
end
ratio = min(left, right) / largest;
end

function cutoff = amplitude_cutoff(cycles, percentile)
cutoff = NaN;
if isempty(cycles) || ~isfinite(percentile)
    return
end
amplitudes = [cycles.amplitude]';
amplitudes = amplitudes(isfinite(amplitudes));
if isempty(amplitudes)
    return
end
cutoff = linear_percentile(amplitudes, percentile);
end

function cycles = add_criteria(cycles, cutoff, options)
for cycle_n = 1:numel(cycles)
    if ~isfinite(cutoff)
        meets_amplitude = isfinite(cycles(cycle_n).amplitude);
    else
        meets_amplitude = cycles(cycle_n).amplitude >= cutoff;
    end
    meets_amplitude_consistency = ...
        cycles(cycle_n).amplitude_consistency >= options.amp_consistency_threshold;
    meets_period_consistency = ...
        cycles(cycle_n).period_consistency >= options.period_consistency_threshold;
    meets_monotonicity = ...
        cycles(cycle_n).monotonicity >= options.monotonicity_threshold;
    cycles(cycle_n).meets_amplitude = meets_amplitude;
    cycles(cycle_n).meets_amplitude_consistency = meets_amplitude_consistency;
    cycles(cycle_n).meets_period_consistency = meets_period_consistency;
    cycles(cycle_n).meets_monotonicity = meets_monotonicity;
    cycles(cycle_n).meets_cycle_criteria = meets_amplitude && ...
        meets_amplitude_consistency && meets_period_consistency && ...
        meets_monotonicity;
end
end

function cycles = assign_episodes(cycles, episodes)
for cycle_n = 1:numel(cycles)
    for episode_n = 1:numel(episodes)
        if cycles(cycle_n).start_sample >= episodes(episode_n).start_sample && ...
                cycles(cycle_n).end_sample <= episodes(episode_n).end_sample
            cycles(cycle_n).episode_index = episode_n;
            break
        end
    end
end
end

function cycles = mark_sustained_sequences(cycles, min_consecutive)
n_cycles = numel(cycles);
if n_cycles == 0
    return
end
inside = [cycles.episode_index]' > 0;
meets = [cycles.meets_cycle_criteria]';
qualifying = enforce_min_run(inside & meets, min_consecutive);
for cycle_n = 1:n_cycles
    cycles(cycle_n).in_sustained_sequence = qualifying(cycle_n);
end
end

function mask = enforce_min_run(mask, min_consecutive)
mask = mask(:).';
if isempty(mask)
    return
end
changes = diff([false, mask, false]);
starts = find(changes == 1);
stops = find(changes == -1);
for run_n = 1:numel(starts)
    if (stops(run_n) - starts(run_n)) < min_consecutive
        mask(starts(run_n):stops(run_n) - 1) = false;
    end
end
mask = mask(:);
end

function episodes = summarize_episodes(episodes, cycles, cutoff)
episodes = blank_cycle_summaries(episodes, cutoff);
if isempty(cycles)
    return
end
for episode_n = 1:numel(episodes)
    selected = cycles([cycles.episode_index]' == episode_n);
    if isempty(selected)
        continue
    end
    periods = [selected.period_s]';
    amplitudes = [selected.amplitude]';
    episodes(episode_n).n_trough_cycles = numel(selected);
    episodes(episode_n).trough_periods_s = periods;
    episodes(episode_n).median_trough_period_s = median(periods);
    episodes(episode_n).trough_period_iqr_s = variability_iqr(periods);
    episodes(episode_n).trough_period_cv = coefficient_of_variation(periods);
    finite_frequency = [selected.frequency_hz]';
    episodes(episode_n).median_trough_frequency_hz = finite_median(finite_frequency);
    episodes(episode_n).median_cycle_amplitude = finite_median(amplitudes);
    episodes(episode_n).cycle_amplitude_cv = coefficient_of_variation(amplitudes);
    episodes(episode_n).cycle_amplitude_iqr = variability_iqr(amplitudes);
    episodes(episode_n).median_rdsym = finite_median([selected.rdsym]');
    episodes(episode_n).median_ptsym = finite_median([selected.ptsym]');
    episodes(episode_n).median_monotonicity = finite_median([selected.monotonicity]');
    episodes(episode_n).median_amplitude_consistency = ...
        finite_median([selected.amplitude_consistency]');
    episodes(episode_n).median_period_consistency = ...
        finite_median([selected.period_consistency]');
    episodes(episode_n).n_cycles_meeting_criteria = nnz([selected.meets_cycle_criteria]);
    episodes(episode_n).n_sustained_cycles = nnz([selected.in_sustained_sequence]);
    episodes(episode_n).n_sustained_runs = count_true_runs([selected.in_sustained_sequence]');
    episodes(episode_n).has_sustained_sequence = ...
        episodes(episode_n).n_sustained_runs > 0;
end
end

function episodes = blank_cycle_summaries(episodes, cutoff)
for episode_n = 1:numel(episodes)
    episodes(episode_n).n_trough_cycles = 0;
    episodes(episode_n).trough_periods_s = zeros(0, 1);
    episodes(episode_n).median_trough_period_s = NaN;
    episodes(episode_n).trough_period_iqr_s = NaN;
    episodes(episode_n).trough_period_cv = NaN;
    episodes(episode_n).median_trough_frequency_hz = NaN;
    episodes(episode_n).median_cycle_amplitude = NaN;
    episodes(episode_n).cycle_amplitude_cv = NaN;
    episodes(episode_n).cycle_amplitude_iqr = NaN;
    episodes(episode_n).median_rdsym = NaN;
    episodes(episode_n).median_ptsym = NaN;
    episodes(episode_n).median_monotonicity = NaN;
    episodes(episode_n).median_amplitude_consistency = NaN;
    episodes(episode_n).median_period_consistency = NaN;
    episodes(episode_n).n_cycles_meeting_criteria = 0;
    episodes(episode_n).n_sustained_cycles = 0;
    episodes(episode_n).n_sustained_runs = 0;
    episodes(episode_n).has_sustained_sequence = false;
    episodes(episode_n).amplitude_percentile_cutoff = cutoff;
end
end

function value = variability_iqr(x)
x = x(isfinite(x));
if numel(x) < 2
    value = NaN;
else
    value = linear_iqr(x);
end
end

function value = coefficient_of_variation(x)
x = x(isfinite(x));
if numel(x) < 2 || mean(x) == 0
    value = NaN;
    return
end
value = std(x, 0) / mean(x);
end

function value = finite_median(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function n_runs = count_true_runs(mask)
mask = logical(mask(:).');
if isempty(mask) || ~any(mask)
    n_runs = 0;
    return
end
changes = diff([false, mask, false]);
n_runs = nnz(changes == 1);
end

function cycles = cycles_inside_episodes(cycles)
if isempty(cycles)
    cycles = empty_cycles();
    return
end
cycles = cycles([cycles.episode_index]' > 0);
if isempty(cycles)
    cycles = empty_cycles();
end
end

function cycles = empty_cycles()
cycles = blank_cycle();
cycles(1) = [];
end

function cycle = blank_cycle()
cycle = struct( ...
    'interval_index', NaN, ...
    'episode_index', 0, ...
    'start_sample', NaN, ...
    'end_sample', NaN, ...
    'peak_sample', NaN, ...
    'rise_midpoint_sample', NaN, ...
    'decay_midpoint_sample', NaN, ...
    'previous_decay_midpoint_sample', NaN, ...
    'period_s', NaN, ...
    'frequency_hz', NaN, ...
    'rise_duration_s', NaN, ...
    'decay_duration_s', NaN, ...
    'rise_amplitude', NaN, ...
    'decay_amplitude', NaN, ...
    'amplitude', NaN, ...
    'rdsym', NaN, ...
    'peak_duration_s', NaN, ...
    'trough_duration_s', NaN, ...
    'ptsym', NaN, ...
    'monotonicity', NaN, ...
    'amplitude_consistency', NaN, ...
    'period_consistency', NaN, ...
    'meets_amplitude', false, ...
    'meets_amplitude_consistency', false, ...
    'meets_period_consistency', false, ...
    'meets_monotonicity', false, ...
    'meets_cycle_criteria', false, ...
    'in_sustained_sequence', false);
end
