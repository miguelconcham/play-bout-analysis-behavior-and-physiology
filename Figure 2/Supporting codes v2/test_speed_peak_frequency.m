function test_speed_peak_frequency()
%TEST_SPEED_PEAK_FREQUENCY  Checks that do not need the experiment files.
%
% Speed bins, band-peak selection, within-session
% correlations, and the session table. Run this file by name after
% add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

options = speed_peak_options();
assert(isequal(options.speed_bin_fraction_limits, [0.05 0.10]), ...
    'Speed bins are not 5-10% of each session.');
assert(isequal(options.delta_band, [1 5]), 'The delta band is not 1-5 Hz.');
assert(isequal(options.theta_band, [6 12]), 'The theta band is not 6-12 Hz.');
assert(options.window_seconds == 10, 'The spectral window is not 10 s.');
assert(options.window_samples == 25000, 'A 10 s window is not 25000 samples.');
assert(abs(options.frequency_resolution - 0.1) < 1e-12, ...
    'The frequency grid is not 0.1 Hz.');
assert(abs(options.frequency(1) - 0.1) < 1e-12 && ...
    abs(options.frequency(end) - 40) < 1e-9, ...
    'The frequency grid does not run from 0.1 to 40 Hz.');
assert(strcmp(options.fooof_aperiodic_mode, 'fixed'), ...
    'FOOOF is not in fixed aperiodic mode.');
assert(isequal(options.fooof_peak_width_limits, [0.5 12]), ...
    'FOOOF peak-width limits are not [0.5 12].');
assert(options.min_speed_samples_per_window == 240, ...
    '80% of a 10 s, 30 Hz window is not 240 speed samples.');
assert(isfield(options, 'speed_definition') && isfield(options, 'session_definition'), ...
    'The saved parameters do not describe speed and sessions.');

changed = speed_peak_options(struct('window_seconds', 4, 'frequency_max', 20));
assert(changed.window_samples == 10000, 'A 4 s window was not honored.');
assert(abs(changed.frequency_resolution - 0.25) < 1e-12, ...
    'A 4 s window does not give 0.25 Hz bins.');
assert(abs(changed.frequency(end) - 20) < 1e-9, ...
    'frequency_max was not applied to the derived grid.');

bins = speed_quantile_bins((1:100).', [0.05 0.10]);
assert(bins.n_bins == 10, '100 samples were not split into 10 bins.');
assert(all(abs(bins.fraction - 0.10) < 1e-12), ...
    'Each bin of 100 samples is not 10% of the data.');
assert(bins.within_requested_fractions, ...
    'Exact 10% bins were flagged as outside 5-10%.');
assert(numel(unique(bins.bin_id)) == 10, 'Not every bin was used.');
assert(all(diff(bins.speed_median) > 0), 'Bin medians do not increase with speed.');
assert(sum(bins.count) == 100, 'Bin counts do not cover every sample.');

scrambled = (1:100).';
scrambled = scrambled(randperm(100));
scrambled_bins = speed_quantile_bins(scrambled, [0.05 0.10]);
assert(isequal(sort(scrambled_bins.bin_id), sort(bins.bin_id)), ...
    'Equal-count bins changed after a permutation.');

short_bins = speed_quantile_bins((1:25).', [0.05 0.10]);
assert(all(short_bins.fraction <= 0.10 + 1e-12), ...
    'A 25-sample session put more than 10% in a bin.');
assert(~short_bins.within_lower_fraction, ...
    'The unavoidable bin below 5% was reported as inside the range.');

tiny_bins = speed_quantile_bins((1:8).', [0.05 0.10]);
assert(tiny_bins.n_bins == 8, 'Eight samples were not kept one per bin.');
assert(~tiny_bins.within_upper_fraction, ...
    'Bins above 10% were reported as inside the upper limit.');
assert(all(isfinite(tiny_bins.bin_id)), 'A finite sample was left unassigned.');

with_nan = [1 2 NaN 4 5 Inf 7 8 9 10].';
nan_bins = speed_quantile_bins(with_nan, [0.05 0.10]);
assert(nan_bins.n_samples == 8, 'Non-finite speeds were counted.');
assert(all(isnan(nan_bins.bin_id(~isfinite(with_nan)))), ...
    'A non-finite speed received a bin.');
assert(sum(nan_bins.count) == 8, 'Finite samples were dropped.');

tied = zeros(20, 1);
tied_bins = speed_quantile_bins(tied, [0.05 0.10]);
assert(tied_bins.n_bins == 10, 'Identical speeds were not split into 10% bins.');
assert(all(tied_bins.speed_median == 0), 'A tied bin changed the speed value.');

peaks = [1.5 0.2 1; 3 0.9 1.2; 7 0.4 2; 8 1.1 2; 4 0.9 1];
delta_peak = select_band_peak(peaks, [1 5]);
theta_peak = select_band_peak(peaks, [6 12]);
assert(delta_peak.present && delta_peak.center_frequency == 3, ...
    'The stronger delta peak was not selected.');
assert(delta_peak.n_peaks_in_band == 3, 'The delta peak count is wrong.');
assert(theta_peak.center_frequency == 8 && theta_peak.power == 1.1, ...
    'The stronger theta peak was not selected.');
none = select_band_peak(peaks, [20 30]);
assert(~none.present && isnan(none.center_frequency), ...
    'An empty band returned a peak.');
empty_peak = select_band_peak(zeros(0, 3), [1 5]);
assert(~empty_peak.present, 'An empty FOOOF peak matrix returned a peak.');
tied_power = select_band_peak([4 0.5 1; 2 0.5 1], [1 5]);
assert(tied_power.center_frequency == 2, ...
    'Equal peak powers did not keep the lower center frequency.');

time = (0:9).' / 30;
position = [(0:9).' * 2, zeros(10, 1)];
[speed, speed_time] = focal_animal_speed(position, time);
assert(all(abs(speed - 60) < 1e-9), ...
    'Two pixels per frame at 30 Hz is not 60 pixels/s.');
assert(numel(speed_time) == 9 && abs(speed_time(1) - time(2)) < 1e-12, ...
    'Speed was not timed at the end of each step.');

speed_axis = (1:10).';
delta_cf = 2 + 0.1 * speed_axis;
theta_cf = 8 - 0.05 * speed_axis;
stats = speed_peak_correlation(speed_axis, delta_cf, theta_cf, 3);
assert(abs(stats.pearson_delta - 1) < 1e-12, 'A perfect delta-speed line is not r = 1.');
assert(abs(stats.pearson_theta + 1) < 1e-12, 'A perfect negative theta-speed line is not r = -1.');
assert(abs(stats.pearson_delta_theta + 1) < 1e-12, ...
    'Opposite delta and theta lines are not r = -1.');
assert(stats.n_delta == 10 && stats.n_theta == 10 && stats.n_delta_theta == 10, ...
    'Finite bins were dropped from the correlation.');

missing = delta_cf;
missing(1:2) = NaN;
missing_stats = speed_peak_correlation(speed_axis, missing, theta_cf, 3);
assert(missing_stats.n_delta == 8, 'NaN peak frequencies were correlated.');
short_stats = speed_peak_correlation(speed_axis(1:2), delta_cf(1:2), theta_cf(1:2), 3);
assert(isnan(short_stats.pearson_delta), 'Two bins produced a correlation.');
flat_stats = speed_peak_correlation(ones(6, 1), delta_cf(1:6), theta_cf(1:6), 3);
assert(isnan(flat_stats.pearson_delta) && isnan(flat_stats.spearman_theta), ...
    'A constant speed produced a correlation.');

recording = fake_recording();
[peak_table, correlation_table] = speed_peak_session_table(recording);
assert(height(peak_table) == 2, 'The session table does not have one row per speed bin.');
assert(peak_table.SpectrumSpeed(2) == 20, 'Spectrum speed was not copied.');
assert(peak_table.WelchThetaCF(1) == 7, 'The Welch theta peak was not copied.');
assert(peak_table.FourierDeltaCF(2) == 3, 'The Fourier delta peak was not copied.');
assert(height(correlation_table) == 2, ...
    'Welch and Fourier were not kept on separate correlation rows.');
assert(correlation_table.PearsonTheta(strcmp(correlation_table.Method, 'Welch')) == 0.5, ...
    'The stored within-session correlation was not copied.');
[empty_peaks, empty_correlations] = speed_peak_session_table(struct([]));
assert(isempty(empty_peaks) && isempty(empty_correlations) && ...
    ismember('SpectrumSpeed', empty_peaks.Properties.VariableNames), ...
    'An empty recording did not return the peak-table columns.');

disp('test_speed_peak_frequency passed');
end

function recording = fake_recording()
peak = select_band_peak(zeros(0, 3), [1 5]);
delta_fast = peak;
delta_fast.present = true;
delta_fast.center_frequency = 3;
delta_fast.power = 0.4;
delta_fast.bandwidth = 1;
theta_slow = peak;
theta_slow.present = true;
theta_slow.center_frequency = 7;
theta_slow.power = 0.8;
theta_slow.bandwidth = 2;
theta_fast = theta_slow;
theta_fast.center_frequency = 8;

speed_bin = struct( ...
    'n_speed_samples', 10, ...
    'fraction_of_speed_samples', 0.10, ...
    'speed_min', 0, ...
    'speed_max', 1, ...
    'speed_median', 0.5, ...
    'speed_mean', 0.5, ...
    'n_windows', 1, ...
    'window_index', 1, ...
    'spectrum_speed', 10);
speed_bin(2) = speed_bin(1);
speed_bin(2).spectrum_speed = 20;
speed_bin(2).window_index = 2;

method = struct();
method.delta_peak = [delta_fast, delta_fast];
method.theta_peak = [theta_slow, theta_fast];
method.r_squared = [0.9 0.8];

stats = speed_peak_correlation([10; 20], [3; 3], [7; 8], 3);
stats.pearson_theta = 0.5;

session = struct( ...
    'partner_index', 2, ...
    'session_index', 1, ...
    'session_source', 'partners session', ...
    'bins_within_requested_fractions', true, ...
    'speed_bins', speed_bin, ...
    'welch', method, ...
    'fourier', method, ...
    'correlation_welch', stats, ...
    'correlation_fourier', stats);
recording = struct( ...
    'animal_name', 'B1D1 1013 Dual', ...
    'animal_area_table_name', 'Batch1Dual1', ...
    'lpag_probe_number', 1, ...
    'lpag_channel', 100, ...
    'sessions', session);
end
