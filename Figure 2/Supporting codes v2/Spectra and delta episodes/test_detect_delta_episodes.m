function test_detect_delta_episodes()
%TEST_DETECT_DELTA_EPISODES  Checks that do not need the experiment files.
%
% Channel choice, amplitude runs, peak-to-peak periods, and the session
% summary. Run this file by name after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(fileparts(this_dir)), 'Custom functions'));
end

options = delta_episode_options(struct( ...
    'max_gap_seconds', 0, ...
    'min_duration_seconds', 0, ...
    'amplitude_threshold', 1));
assert(options.max_gap_seconds == 0, 'A gap of 0 s was replaced by the default.');
assert(options.min_duration_seconds == 0, 'A duration minimum of 0 was replaced.');
assert(options.min_cycles == 3, 'The default cycle minimum was lost.');
assert(isequal(options.delta_freq, [1 5]), 'The delta band is not 1-5 Hz.');
assert(options.delta_filter_order == 2000, 'The delta FIR order is not 2000.');
assert(isequal(options.histogram_ms, [50 1000]), ...
    'The cycle-period axis is not 50-1000 ms.');
assert(strcmp(options.cycle_threshold_preset, 'motor_beta'), ...
    'The default cycle preset is not the published motor-cortex set.');
assert(options.amp_consistency_threshold == 0.3, ...
    'Motor-cortex amplitude consistency is not 0.3.');
assert(options.period_consistency_threshold == 0.5, ...
    'Motor-cortex period consistency is not 0.5.');
assert(options.monotonicity_threshold == 0.6, ...
    'Motor-cortex monotonicity is not 0.6.');
assert(options.cycle_amplitude_percentile == 20, ...
    'The motor-cortex amplitude gate is not the 20th percentile.');
assert(options.min_consecutive_cycles == 3, ...
    'The published minimum is not 3 consecutive cycles.');
assert(options.cycle_lowpass_hz == 20, ...
    'The cycle low-pass is not 4 times the 5 Hz upper edge.');
assert(options.cycle_lowpass_order == round(0.1 * options.sampling_rate), ...
    'The cycle low-pass is not the published 100 ms duration.');

simulation = delta_episode_options(struct('cycle_threshold_preset', 'simulation'));
assert(simulation.amp_consistency_threshold == 0.6 && ...
    simulation.period_consistency_threshold == 0.6 && ...
    simulation.monotonicity_threshold == 0.9, ...
    'The simulation preset is not 0.6, 0.6, 0.9.');
assert(isnan(simulation.cycle_amplitude_percentile), ...
    'The simulation preset should not apply an amplitude percentile.');
overridden = delta_episode_options(struct( ...
    'cycle_threshold_preset', 'simulation', ...
    'amp_consistency_threshold', 0.25));
assert(overridden.amp_consistency_threshold == 0.25 && ...
    overridden.period_consistency_threshold == 0.6, ...
    'An explicit consistency threshold was replaced by the preset.');
wide = delta_episode_options(struct('delta_freq', [2 8], 'sampling_rate', 1000));
assert(wide.cycle_lowpass_hz == 32 && wide.cycle_lowpass_order == 100, ...
    'The low-pass did not follow 4x the upper cutoff and 100 ms.');

failed = false;
try
    delta_episode_options(struct('delta_freq', [5 1]));
catch
    failed = true;
end
assert(failed, 'A reversed delta band was accepted.');

test_channel_selection();
test_regular_episode();
test_gap_and_rejection();
test_summary();
test_cycle_shape();
test_cycle_boundaries_and_consistency();
if exist('designfilt', 'file') == 2 || exist('designfilt', 'builtin') == 5
    test_filter_on_a_sine([]);
elseif exist('fir1', 'file') == 2
    sampling_rate = 2500;
    fir_coefficients = fir1(2000, [1 5] / (sampling_rate / 2), 'bandpass', hamming(2001));
    test_filter_on_a_sine(fir_coefficients);
else
    disp('No FIR designer is available here; the filter integration check was skipped.');
end
disp('test_detect_delta_episodes passed');
end

function test_channel_selection()
xcoords = [16; 16; 16];
ycoords = [100; 110; 200];
chan_map = [4; 5; 6];
selected = select_mid_lpag_channels(xcoords, ycoords, chan_map, [1 50 250], 1);
assert(selected.channel == 5 && selected.lfp_row == 5, ...
    'The mid-LPAG site was not the channel nearest the LPAG centroid.');
assert(selected.y_um == 110, 'The selected depth is not the centroid depth.');

xcoords = [0; 100; 50];
ycoords = [0; 0; 11];
chan_map = [1; 2; 3];
selected = select_mid_lpag_channels(xcoords, ycoords, chan_map, [1 -10 10], 1);
assert(selected.channel ~= 3, ...
    'A site outside the LPAG depth was closer and was selected.');

xcoords = [8; 40; 8; 258];
ycoords = [100; 100; 400; 100];
chan_map = [10; 20; 30; 40];
selected = select_mid_lpag_channels(xcoords, ycoords, chan_map, [1 50 450], 2);
assert(selected.channel == 10 && selected.lfp_row == 10, ...
    'The NPX2 mid-LPAG site was not the nearest site on that shank.');
assert(selected.probe_number == 1, 'The LPAG probe number was not kept.');

chan_map = [(1:191).'; (193:384).'];
xcoords = ones(numel(chan_map), 1);
ycoords = (1:numel(chan_map)).';
selected = select_mid_lpag_channels(xcoords, ycoords, chan_map, [1 100 120], 1);
assert(selected.channel == 110 && selected.y_um == 110, ...
    'The NPX1 map that skips channel 192 did not keep the mid-LPAG site.');
end

function test_regular_episode()
sampling_rate = 2500;
n_samples = 8 * sampling_rate;
time_s = (0:n_samples - 1) / sampling_rate;
filtered = sin(2 * pi * 3 * time_s);
amplitude = ones(1, n_samples);
options = detection_options(struct( ...
    'amplitude_threshold_mode', 'absolute', ...
    'amplitude_threshold', 0.5));
detected = characterize_delta_episodes(filtered, amplitude, sampling_rate, options);

assert(detected.n_episodes == 1, 'A continuous 3 Hz sine was not one episode.');
episode = detected.episodes(1);
assert(episode.n_cycles >= 20 && episode.n_cycles <= 26, ...
    'The 8 s sine did not yield about 24 peak-to-peak cycles.');
assert(abs(episode.median_cycle_period_s - (1 / 3)) < 0.002, ...
    'The median peak-to-peak period was not 1/3 s.');
assert(episode.fraction_cycles_in_delta_band == 1, ...
    'A regular 3 Hz run had cycles outside 1-5 Hz.');
assert(episode.cycle_period_iqr_s < 0.002, ...
    'A regular sine had a large cycle-period IQR.');
assert(abs(episode.duration_s - 8) < 1e-9, ...
    'The episode duration was not the full interval.');

std_options = detection_options(struct('amplitude_threshold_mode', 'std'));
std_detected = characterize_delta_episodes(filtered, amplitude, ...
    sampling_rate, std_options);
assert(abs(std_detected.threshold - 1) < 1e-9, ...
    'A flat envelope did not threshold at its mean.');
assert(std_detected.n_episodes == 1, ...
    'A flat envelope above its mean was not one episode.');
end

function test_gap_and_rejection()
sampling_rate = 2500;
n_samples = 10 * sampling_rate;
time_s = (0:n_samples - 1) / sampling_rate;
filtered = sin(2 * pi * 3 * time_s);
gapped = filtered;
gapped((4 * sampling_rate + 1):(6 * sampling_rate)) = 0;
amplitude = 0.1 * ones(1, n_samples);
amplitude((2 * sampling_rate + 1):(4 * sampling_rate)) = 1;
amplitude((6 * sampling_rate + 1):(8 * sampling_rate)) = 1;

separate = characterize_delta_episodes(gapped, amplitude, sampling_rate, ...
    detection_options(struct('max_gap_seconds', 0, ...
    'amplitude_threshold_mode', 'absolute', 'amplitude_threshold', 0.5)));
assert(separate.n_episodes == 2, 'Two amplitude bursts were merged with no gap fill.');
assert(all(arrayfun(@(episode) max(episode.cycle_periods_s) < 0.5, separate.episodes)), ...
    'A separated burst included the silent gap as a cycle.');

merged = characterize_delta_episodes(gapped, amplitude, sampling_rate, ...
    detection_options(struct('max_gap_seconds', 3, ...
    'amplitude_threshold_mode', 'absolute', 'amplitude_threshold', 0.5)));
assert(merged.n_episodes == 1, 'A 2 s gap was not joined when the limit was 3 s.');
assert(max(merged.episodes(1).cycle_periods_s) > 1, ...
    'The joined episode did not record the long cycle across the silent gap.');

short_amplitude = 0.1 * ones(1, n_samples);
short_amplitude((5 * sampling_rate + 1):(5 * sampling_rate + round(0.5 * sampling_rate))) = 1;
short = characterize_delta_episodes(filtered, short_amplitude, sampling_rate, ...
    detection_options(struct('amplitude_threshold_mode', 'absolute', ...
    'amplitude_threshold', 0.5, 'min_cycles', 3)));
assert(short.n_threshold_runs == 1, 'The 0.5 s burst was not a threshold run.');
assert(short.n_episodes == 0, 'A burst shorter than 3 cycles was kept.');

edge_amplitude = 0.1 * ones(1, n_samples);
edge_amplitude(1:round(0.3 * sampling_rate)) = 1;
edge = characterize_delta_episodes(filtered, edge_amplitude, sampling_rate, ...
    detection_options(struct('amplitude_threshold_mode', 'absolute', ...
    'amplitude_threshold', 0.5, 'edge_margin_seconds', 1, 'min_cycles', 1)));
assert(edge.n_threshold_runs == 0, ...
    'A burst inside the filter-edge margin was still searched.');
end

function test_summary()
recording = struct();
recording.animal_name = 'Test Animal';
recording.lpag_probe_number = 1;
recording.lpag_channel = 10;
recording.lpag_x_um = 8;
recording.lpag_y_um = 100;
recording.lpag_depth_range = [50 250];
recording.interval_mode = 'partners_session';
recording.parameters = struct('delta_freq', [1 5]);
recording.intervals = struct( ...
    'duration_s', 10, ...
    'n_threshold_runs', 3, ...
    'amplitude_threshold', 2);
recording.episodes = [ ...
    fake_episode(2, 4, 0.3, 0.02, 0.1, 1, [0.3; 0.3; 0.3; 0.3]), ...
    fake_episode(4, 6, 0.5, 0.04, 0.2, 2 / 3, [0.5; 0.5; 1.5])];

summary = summarize_delta_episode_recording(recording);
assert(summary.n_episodes == 2, 'The episode count was wrong.');
assert(summary.n_rejected_runs == 1, 'Rejected threshold runs were not counted.');
assert(abs(summary.episode_time_fraction - 0.6) < 1e-12, ...
    'Episode time fraction was not 0.6.');
assert(abs(summary.median_duration_s - 3) < 1e-12, 'Median duration was not 3 s.');
assert(abs(summary.median_n_cycles - 5) < 1e-12, 'Median cycle count was not 5.');
assert(abs(summary.median_cycle_period_s - 0.4) < 1e-12, ...
    'Median of median cycle periods was not 0.4 s.');
assert(abs(summary.median_cycle_period_iqr_s - 0.03) < 1e-12, ...
    'Median cycle-period IQR was not 0.03 s.');
assert(summary.pooled_n_cycles == 7, 'Pooled cycle count was not 7.');
assert(abs(summary.pooled_median_cycle_period_s - 0.3) < 1e-12, ...
    'Pooled median period was not 0.3 s.');
assert(abs(summary.pooled_fraction_cycles_in_delta_band - (6 / 7)) < 1e-12, ...
    'The 1.5 s cycle was counted inside the 1-5 Hz band.');
assert(abs(summary.iqr_duration_s - 1) < 1e-12, ...
    'The IQR of episode durations was not 1 s.');
assert(abs(linear_iqr([1 2 3 4]) - 1.5) < 1e-12, ...
    'Linear IQR of 1:4 was not 1.5.');

if exist('table', 'file') == 2 || exist('table', 'builtin') == 5
    summary_table = delta_episode_summary_table(recording);
    assert(height(summary_table) == 1 && summary_table.n_episodes == 2, ...
        'The session table did not keep the episode count.');
end
end

function test_filter_on_a_sine(fir_coefficients)
sampling_rate = 2500;
time_s = (0:20 * sampling_rate - 1) / sampling_rate;
lfp = sin(2 * pi * 3 * time_s);
lfp(1:40) = 20;
options = detection_options(struct( ...
    'edge_margin_seconds', 2, ...
    'amplitude_threshold_mode', 'absolute', ...
    'amplitude_threshold', 0.4, ...
    'min_cycles', 3));
if ~isempty(fir_coefficients)
    options.filter_coefficients = fir_coefficients;
end
detected = detect_delta_episodes(lfp, sampling_rate, options);
assert(detected.n_episodes >= 1, 'The filtered 3 Hz sine produced no episode.');
assert(abs(detected.episodes(1).median_cycle_period_s - (1 / 3)) < 0.01, ...
    'The filtered sine median period was not near 1/3 s.');
assert(max(abs(detected.lfp(1:40))) > 10, ...
    'The minimally processed trace no longer shows the transient.');
assert(detected.episodes(1).n_trough_cycles >= 10, ...
    'The filtered sine did not yield trough-to-trough cycles.');
assert(abs(detected.episodes(1).median_trough_period_s - (1 / 3)) < 0.02, ...
    'The trough-to-trough period was not near 1/3 s.');
end

function test_cycle_shape()
sampling_rate = 250;
n_nominal_cycles = 7;
[broad, n_samples] = asymmetric_train(n_nominal_cycles, sampling_rate);
episodes = struct('start_sample', 1, 'end_sample', n_samples);
[episodes, cycles] = characterize_cycle_features(broad, broad, ...
    sampling_rate, episodes, cycle_options(struct()));

assert(numel(cycles) == 5, ...
    'A 7-cycle train did not leave 5 cycles with a zero-crossing on both sides.');
assert(abs(episodes.median_trough_period_s - 0.4) < 1e-9, ...
    'The trough-to-trough period was not 0.4 s.');
assert(abs(episodes.median_rdsym - 0.2) < 1e-9, ...
    'Rise-decay symmetry was not the 20-sample rise over a 100-sample period.');
assert(abs(episodes.median_cycle_amplitude - 4) < 1e-9, ...
    'Cycle amplitude was not the mean of the two flank voltages.');
assert(abs(episodes.median_ptsym - 0.49) < 1e-9, ...
    'Peak-trough symmetry was not peak/(peak + previous trough).');
assert(episodes.median_monotonicity == 1, ...
    'A straight rise and decay was not fully monotonic.');
assert(episodes.trough_period_cv == 0, ...
    'Identical periods did not have CV 0.');
assert(episodes.n_sustained_cycles == 3 && episodes.has_sustained_sequence, ...
    'Three interior cycles were not a sustained sequence.');
assert(~cycles(1).in_sustained_sequence && ~cycles(end).in_sustained_sequence, ...
    'The first or last cycle was marked sustained without a neighbor.');
assert(abs(cycles(3).frequency_hz - 2.5) < 1e-9, ...
    'Cycle frequency was not 1/period.');

shifted = broad;
for cycle_n = 0:(n_nominal_cycles - 1)
    shifted(1 + 100 * cycle_n + 40) = 3;
end
[shifted_episodes, shifted_cycles] = characterize_cycle_features( ...
    shifted, broad, sampling_rate, episodes, cycle_options(struct()));
assert(all(abs([shifted_cycles.rdsym] - 0.4) < 1e-9), ...
    'The peak was not taken from the low-pass signal.');
assert(abs(shifted_episodes.median_cycle_amplitude - 5) < 1e-9, ...
    'The shifted peak did not set the flank amplitude.');
end

function test_cycle_boundaries_and_consistency()
sampling_rate = 250;
n_nominal_cycles = 7;
[broad, n_samples] = asymmetric_train(n_nominal_cycles, sampling_rate);
partial = struct('start_sample', 150, 'end_sample', 350);
[partial_episodes, partial_cycles] = characterize_cycle_features( ...
    broad, broad, sampling_rate, partial, cycle_options(struct()));
assert(partial_episodes.n_trough_cycles == 1 && numel(partial_cycles) == 1, ...
    'A cycle cut by the episode boundary was kept.');
assert(~partial_episodes.has_sustained_sequence, ...
    'One qualifying cycle was called a sustained sequence.');

damaged = broad;
% One decay step reversed on a cycle that would otherwise be interior.
damaged(331) = damaged(330) + 0.2;
damaged_options = cycle_options(struct('monotonicity_threshold', 0.995));
[damaged_episodes, damaged_cycles] = characterize_cycle_features( ...
    damaged, broad, sampling_rate, ...
    struct('start_sample', 1, 'end_sample', n_samples), damaged_options);
damaged_cycle = damaged_cycles([damaged_cycles.peak_sample] > 320 & ...
    [damaged_cycles.peak_sample] < 360);
assert(numel(damaged_cycle) == 1, 'The damaged cycle was not found.');
assert(abs(damaged_cycle.monotonicity - 0.99) < 1e-9, ...
    'Monotonicity was not the fraction of steps across both flanks.');
assert(~damaged_cycle.meets_monotonicity, ...
    'A 0.99 monotonicity passed a 0.995 threshold.');
assert(~damaged_episodes.has_sustained_sequence, ...
    'Breaking one interior cycle still left three consecutive cycles.');

long = asymmetric_train_with_gap(sampling_rate);
[long_episodes, long_cycles] = characterize_cycle_features( ...
    long, long, sampling_rate, ...
    struct('start_sample', 1, 'end_sample', numel(long)), cycle_options(struct()));
long_cycle = long_cycles([long_cycles.period_s] > 0.6);
assert(~isempty(long_cycle), 'The long cycle was not measured.');
assert(long_cycle(1).period_consistency < 0.5, ...
    'A neighbor twice as long still met period consistency 0.5.');
assert(long_episodes.trough_period_cv > 0, ...
    'Unequal periods did not produce a period CV.');
assert(~long_episodes.has_sustained_sequence, ...
    'Failed period consistency was still labeled a sustained sequence.');
end

function [broad, n_samples] = asymmetric_train(n_nominal_cycles, sampling_rate)
% Piecewise-linear 2.5 Hz cycles at sampling_rate 250: 20 samples up,
% 80 samples down, from -2 to +2. sampling_rate is accepted so the
% caller and the signal stay tied together; the geometry is in samples.
if sampling_rate ~= 250
    error('This fixture is defined at 250 Hz.');
end
n_samples = 1 + 100 * n_nominal_cycles;
broad = zeros(1, n_samples);
for cycle_n = 0:(n_nominal_cycles - 1)
    trough = 1 + 100 * cycle_n;
    peak = trough + 20;
    next_trough = trough + 100;
    broad(trough:peak) = linspace(-2, 2, peak - trough + 1);
    broad(peak:next_trough) = linspace(2, -2, next_trough - peak + 1);
end
end

function broad = asymmetric_train_with_gap(sampling_rate)
if sampling_rate ~= 250
    error('This fixture is defined at 250 Hz.');
end
% Five 100-sample cycles, then one 250-sample cycle, then two more.
lengths = [100 100 100 100 100 250 100 100];
n_samples = 1 + sum(lengths);
broad = zeros(1, n_samples);
cursor = 1;
for length_n = 1:numel(lengths)
    span = lengths(length_n);
    peak = cursor + round(0.2 * span);
    next_trough = cursor + span;
    broad(cursor:peak) = linspace(-2, 2, peak - cursor + 1);
    broad(peak:next_trough) = linspace(2, -2, next_trough - peak + 1);
    cursor = next_trough;
end
end

function options = cycle_options(overrides)
defaults = struct('edge_margin_seconds', 0);
names = fieldnames(overrides);
for name_n = 1:numel(names)
    defaults.(names{name_n}) = overrides.(names{name_n});
end
options = delta_episode_options(defaults);
end

function options = detection_options(overrides)
defaults = struct( ...
    'edge_margin_seconds', 0, ...
    'min_cycles', 3, ...
    'min_duration_seconds', 0, ...
    'max_gap_seconds', 0, ...
    'min_peak_distance_seconds', 0.05, ...
    'amplitude_smooth_seconds', 0, ...
    'delta_freq', [1 5]);
names = fieldnames(overrides);
for name_n = 1:numel(names)
    defaults.(names{name_n}) = overrides.(names{name_n});
end
options = delta_episode_options(defaults);
end

function episode = fake_episode(duration_s, n_cycles, median_period, period_iqr, ...
        period_cv, fraction_in_band, periods_s)
episode = struct();
episode.duration_s = duration_s;
episode.n_cycles = n_cycles;
episode.median_cycle_period_s = median_period;
episode.cycle_period_iqr_s = period_iqr;
episode.cycle_period_cv = period_cv;
episode.fraction_cycles_in_delta_band = fraction_in_band;
episode.cycle_periods_s = periods_s;
end
