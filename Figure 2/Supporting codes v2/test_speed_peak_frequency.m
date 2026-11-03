function test_speed_peak_frequency()
%TEST_SPEED_PEAK_FREQUENCY  Checks that do not need the experiment files.
%
% Speed bins, band-peak selection, within-session
% correlations, peak-frequency regressions, and the session table.
% Run this file by name after add_repo_paths.

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

x = (1:4).';
y = [2; 2; 4; 4];
line_fit = speed_peak_ols(x, y, 3);
assert(line_fit.ok && abs(line_fit.slope - 0.8) < 1e-12 && ...
    abs(line_fit.intercept - 1) < 1e-12 && ...
    abs(line_fit.slope_se - sqrt(0.08)) < 1e-12, ...
    'The closed-form slope is not 0.8 with intercept 1.');
perfect = speed_peak_ols((1:6).', 2 + 3 * (1:6).', 3);
assert(abs(perfect.slope - 3) < 1e-10 && abs(perfect.intercept - 2) < 1e-10 && ...
    abs(perfect.r_squared - 1) < 1e-12 && perfect.slope_p == 0, ...
    'A perfect line was not recovered.');
flat_peak = speed_peak_ols((1:5).', 4 * ones(5, 1), 3);
assert(flat_peak.ok && flat_peak.slope == 0 && flat_peak.slope_p == 1, ...
    'A constant peak frequency did not give a zero slope.');
too_short = speed_peak_ols((1:2).', (1:2).', 3);
assert(~too_short.ok && isnan(too_short.slope), ...
    'Two rows produced a regression.');
assert(~speed_peak_ols(ones(5, 1), (1:5).', 3).ok, ...
    'A constant predictor produced a regression.');

speed_for_models = [0; (1:6).'; 8];
log_peak = [2; exp(0.2 + 0.3 * (1:6).'); -1];
peak_models = speed_peak_models(speed_for_models, log_peak, 3);
log_frequency = peak_models(strcmp({peak_models.name}, 'log_frequency'));
linear_model = peak_models(strcmp({peak_models.name}, 'linear'));
assert(log_frequency.n == 6, 'Non-positive speed or peak entered the regression.');
assert(abs(log_frequency.slope - 0.3) < 1e-8 && ...
    abs(log_frequency.intercept - 0.2) < 1e-8, ...
    'log(PeakHz) ~ Speed did not recover the log-frequency line.');
assert(log_frequency.aic_hz < linear_model.aic_hz, ...
    'The log-frequency model lost to a linear model on log-linear data.');
jacobian_speed = (1:8).';
jacobian_peak = exp(0.1 + 0.2 * jacobian_speed) + 0.15 * sin(jacobian_speed);
jacobian_models = speed_peak_models(jacobian_speed, jacobian_peak, 3);
jacobian_log = jacobian_models(strcmp({jacobian_models.name}, 'log_frequency'));
assert(isfinite(jacobian_log.aic) && abs(jacobian_log.aic_hz - ...
    (jacobian_log.aic + 2 * sum(log(jacobian_peak)))) < 1e-6, ...
    'The hertz AIC is missing the log Jacobian.');

log_speed_axis = [0; exp((1:6).'); 4];
log_speed_peak = [NaN; 4 + 1.5 * log(log_speed_axis(2:7)); 0];
log_speed_models = speed_peak_models(log_speed_axis, log_speed_peak, 3);
log_speed = log_speed_models(strcmp({log_speed_models.name}, 'log_speed'));
assert(log_speed.n == 6 && abs(log_speed.slope - 1.5) < 1e-8 && ...
    abs(log_speed.intercept - 4) < 1e-8, ...
    'PeakHz ~ log(Speed) did not recover the log-speed line.');
assert(log_speed.aic_hz < log_speed_models(strcmp({log_speed_models.name}, 'linear')).aic_hz, ...
    'The log-speed model lost to a linear model on logarithmic data.');

animal = {'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'A'; 'B'; 'B'; 'A'};
session_key = {'A|s1'; 'A|s1'; 'A|s1'; 'A|s1'; 'A|s1'; 'A|s2'; 'A|s2'; 'A|s1'; 'A|s1'; 'B|s1'; 'B|s1'; 'A|s1'};
bin_id = [1; 1; 2; 2; 3; 1; 1; 4; 5; 1; 2; 6];
animal_speed = [1; 1; 2; 2; 3; 6; 6; 0; 4; 1; 2; NaN];
animal_peak = [8; 8; 7; 7; 6; 5; 5; 9; NaN; 4; 3; 3];
[speed_z, speed_reference] = within_animal_speed_z( ...
    animal, session_key, bin_id, animal_speed);
a_reference = speed_reference(strcmp({speed_reference.animal}, 'A'));
a_unique = [0 1 2 3 4 6];
assert(a_reference.n_unique == numel(a_unique), ...
    'Duplicate channels changed the number of speeds in the z-score.');
assert(abs(a_reference.mean - mean(a_unique)) < 1e-12 && ...
    abs(a_reference.sd - std(a_unique, 0)) < 1e-12, ...
    'The within-animal SD was not the SD of the unique session bins.');
assert(abs(speed_z(1) - (1 - mean(a_unique)) / std(a_unique, 0)) < 1e-12 && ...
    speed_z(1) == speed_z(2), ...
    'The two channels of one bin did not share a z-score.');
assert(isfinite(speed_z(9)), ...
    'A bin with a missing peak was left out of the speed scale.');
assert(isnan(speed_z(12)), 'A missing speed received a z-score.');
b_unique = [1 2];
b_z = (1 - mean(b_unique)) / std(b_unique, 0);
assert(abs(speed_z(10) - b_z) < 1e-12, ...
    'Animal B was z-scored on animal A''s speeds.');

obs = speed_peak_pooled_observations(animal, session_key, bin_id, ...
    animal_speed, animal_peak);
assert(obs.n == 9, 'The pooled rows kept a non-positive speed or a missing peak.');
assert(abs(obs.speed_z(1) - speed_z(1)) < 1e-12, ...
    'A missing peak changed the z-score of the bins that remain.');
assert(~any(obs.speed == 0) && all(obs.peak > 0), ...
    'The shared comparison rows include a non-positive speed or peak.');
linear_rows = speed_peak_ols(obs.speed_z, obs.peak, 3);
log_frequency_rows = speed_peak_ols(obs.speed_z, obs.log_peak, 3);
log_speed_rows = speed_peak_ols(obs.log_speed, obs.peak, 3);
assert(linear_rows.n == obs.n && log_frequency_rows.n == obs.n && ...
    log_speed_rows.n == obs.n, ...
    'The three pooled functional forms did not use the same rows.');

two_sessions = speed_peak_pooled_formula( ...
    {'A'; 'A'; 'B'}, {'A|s1'; 'A|s2'; 'B|s1'});
assert(strcmp(two_sessions.random_terms, '(1|Animal) + (1|Session)'), ...
    'Repeated sessions within animal did not keep both intercepts.');
one_session_each = speed_peak_pooled_formula({'A'; 'B'}, {'A|s1'; 'B|s1'});
assert(strcmp(one_session_each.random_terms, '(1|Animal)'), ...
    'One session per animal kept a separate session intercept.');
one_animal = speed_peak_pooled_formula({'A'; 'A'}, {'A|s1'; 'A|s2'});
assert(strcmp(one_animal.random_terms, '(1|Session)'), ...
    'A single animal did not use its sessions as the repeated measure.');

spec = speed_peak_model_specification();
assert(any(strcmp(spec.per_session_models, 'PeakHz ~ Speed')) && ...
    any(strcmp(spec.per_session_models, 'log(PeakHz) ~ Speed')) && ...
    any(strcmp(spec.per_session_models, 'PeakHz ~ log(Speed)')), ...
    'The saved specification is missing a peak-frequency model.');
assert(~isempty(strfind(spec.why_zscore_linear_speed, 'within-animal z-score')) && ...
    ~isempty(strfind(spec.log_speed_normalization, 'animal intercept')), ...
    'The saved specification does not record the speed-scale decision.');

if exist('cell2table', 'file')
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
end

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
