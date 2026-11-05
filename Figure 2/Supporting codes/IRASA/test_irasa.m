function test_irasa()
%TEST_IRASA  IRASA checks that do not need the experiment files.
%
% Defaults from Wen & Liu (2016) section 2.4, recovery of a power-law
% spectrum, a 10 Hz peak left in the residual, and a 16 Hz harmonic that
% stays in the residual. Run this file by name after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

options = irasa_options();
assert(numel(options.h_values) == 17, 'The published h set does not have 17 values.');
assert(abs(options.h_values(1) - 1.1) < 1e-12 && ...
    abs(options.h_values(end) - 1.9) < 1e-12, ...
    'h does not run from 1.1 to 1.9.');
assert(max(abs(diff(options.h_values) - 0.05)) < 1e-12, ...
    'The h increment is not 0.05.');
assert(options.n_segments == 10, 'The default is not 10 windows.');
assert(abs(options.segment_fraction - 0.9) < 1e-12, ...
    'The default window is not 90% of the epoch.');
assert(options.detrend && options.antialias, ...
    'Detrend and anti-alias are not on by default.');
assert(isempty(options.frequency_range), ...
    'The frequency crop is not left to the sampling rate.');

overridden = irasa_options(struct('n_segments', 4, 'h_values', [1.2 1.5]));
assert(overridden.n_segments == 4 && numel(overridden.h_values) == 2, ...
    'An explicit h set or window count was replaced.');
assert(abs(overridden.segment_fraction - 0.9) < 1e-12, ...
    'An unset window fraction did not keep 0.9.');

failed = false;
try
    irasa_options(struct('h_values', [0.5 1.2]));
    failed = true;
catch
end
assert(~failed, 'An h value below 1 was accepted.');

rng(7);
sampling_rate = 200;
n_samples = 4000;
time = (0:n_samples - 1).' / sampling_rate;
fractal = irasa_simulate_fractal(n_samples, sampling_rate, 1.5);
pure = irasa_separate(fractal, sampling_rate);
assert_separation_identity(pure);
assert(pure.n_segments == 10, 'The separation did not use 10 windows.');
assert(pure.segment_length == round(0.9 * n_samples), ...
    'The window length is not 90% of the epoch.');
assert(pure.segment_starts(1) == 1, 'The first window does not start at sample 1.');
assert(pure.segment_starts(end) == n_samples - pure.segment_length + 1, ...
    'The last window does not end on the last sample.');
assert(max(abs(diff(pure.segment_starts) - median(diff(pure.segment_starts)))) <= 1, ...
    'The windows are not evenly spaced.');
expected_nfft = 2 * 2 ^ ceil(log2(pure.segment_length));
assert(pure.nfft == expected_nfft, 'NFFT is not twice the next power of two.');
assert(all(diff(pure.frequency) > 0), 'Frequencies are not increasing.');
assert(pure.frequency(1) > 0, 'DC was kept in the default spectrum.');
assert(pure.frequency(end) <= sampling_rate / (2 * max(pure.h_values)) + 1e-9, ...
    'The spectrum extends past the largest-h bandwidth.');
row_input = irasa_separate(fractal.', sampling_rate);
assert(max(abs(row_input.fractal - pure.fractal)) < 1e-10, ...
    'A row vector did not match the column vector.');

fit = irasa_power_law_fit(pure, [2 40]);
assert(abs(fit.exponent - 1.5) < 0.25, ...
    'A pure power law was not recovered near exponent 1.5.');
in_fit = pure.frequency >= 2 & pure.frequency <= 40;
fractal_log_rms = sqrt(mean((log10(pure.fractal(in_fit)) - ...
    log10(fit.fitted_power)) .^ 2));
mixed_log_rms = sqrt(mean((log10(pure.mixed(in_fit)) - ...
    log10(fit.fitted_power)) .^ 2));
assert(fractal_log_rms < 0.2, ...
    'The fractal spectrum is not close to its power-law line.');
assert(fractal_log_rms < mixed_log_rms / 2, ...
    'The fractal estimate is not smoother than the mixed periodogram.');
scale_ratio = median(pure.fractal(in_fit) ./ pure.mixed(in_fit));
assert(scale_ratio > 0.5 && scale_ratio < 1.5, ...
    'The fractal spectrum is not on the same scale as the mixed spectrum.');

oscillation = 0.8 * sin(2 * pi * 10 * time);
mixed_signal = fractal + oscillation;
separated = irasa_separate(mixed_signal, sampling_rate);
assert_separation_identity(separated);
assert(local_peak_ratio(separated, 10) > 4, ...
    'The 10 Hz rhythm did not remain a peak after subtraction.');
assert(fractal_peak_ratio(separated, 10) < 1.5, ...
    'The fractal estimate still has a 10 Hz peak.');
mixed_ratio = band_ratio(separated.frequency, separated.mixed, 10);
fractal_ratio = band_ratio(separated.frequency, separated.fractal, 10);
assert(mixed_ratio > 3 * fractal_ratio, ...
    'The fractal estimate did not reject the 10 Hz peak relative to the mixed spectrum.');

harmonic = sin(2 * pi * 8 * time) + 0.55 * sin(2 * pi * 16 * time);
harmonic_signal = 0.35 * fractal + harmonic;
harmonic_spec = irasa_separate(harmonic_signal, sampling_rate);
assert(local_peak_ratio(harmonic_spec, 8) > 4, ...
    'The 8 Hz peak did not remain in the residual.');
assert(local_peak_ratio(harmonic_spec, 16) > 4, ...
    'The 16 Hz harmonic did not remain in the residual.');

fast = sin(2 * pi * 90 * time);
aliased = irasa_separate(0.5 * fractal + fast, sampling_rate, ...
    struct('frequency_range', [1 95]));
low = aliased.frequency >= 1 & aliased.frequency <= 30;
high = nearest_bin(aliased.frequency, 90);
assert(abs(aliased.oscillatory(high)) > 5 * max(abs(aliased.oscillatory(low))), ...
    'The anti-alias filter left the 90 Hz peak inside 1-30 Hz.');

two_channels = irasa_separate([fractal, fractal], sampling_rate);
assert(max(abs(two_channels.fractal(:, 1) - two_channels.fractal(:, 2))) < 1e-10, ...
    'Identical channels did not produce identical fractal spectra.');

cropped = irasa_separate(fractal, sampling_rate, ...
    struct('frequency_range', [2 40]));
session_summary = irasa_session_summary([cropped cropped], [2 20]);
assert(session_summary.n_windows == 2, 'The session summary dropped a window.');
assert(max(abs(session_summary.mixed - cropped.mixed)) < 1e-10, ...
    'Averaging identical windows changed the mixed spectrum.');
assert(max(abs(session_summary.oscillatory - ...
    (session_summary.mixed - session_summary.fractal))) < 1e-10, ...
    'The session residual is not mixed minus fractal.');
direct_fit = irasa_power_law_fit(cropped, [2 20]);
assert(abs(session_summary.exponent - direct_fit.exponent) < 1e-8, ...
    'The session-mean exponent does not match the window fit.');
assert(isequal(session_summary.fit_range, [2 20]), ...
    'The session fit range was not stored.');
failed = false;
try
    mismatched = cropped;
    mismatched.frequency = cropped.frequency(1:end - 1);
    irasa_session_summary([cropped mismatched], [2 20]);
    failed = true;
catch
end
assert(~failed, 'Windows with different frequency grids were averaged.');

sampling_rate_raw = 250;
n_raw = sampling_rate_raw * 10;
raw_time = (0:n_raw - 1).' / sampling_rate_raw;
raw = irasa_raw_window_spectrum(sin(2 * pi * 10 * raw_time), sampling_rate_raw);
assert(abs(raw.frequency(1) - 0.1) < 1e-12 && ...
    abs(raw.frequency(end) - 40) < 1e-9, ...
    'The raw spectrum is not the 0.1-40 Hz FOOOF grid.');
[~, welch_peak] = max(raw.welch);
[~, fourier_peak] = max(raw.fourier);
assert(abs(raw.frequency(welch_peak) - 10) < 1e-9 && ...
    abs(raw.frequency(fourier_peak) - 10) < 1e-9, ...
    'The raw Welch or Fourier spectrum missed the 10 Hz sinusoid.');
components = irasa_session_components(session_summary, ...
    raw.frequency, [raw.welch.'; raw.welch.'], [raw.fourier.'; raw.fourier.']);
assert(max(abs(components.raw.welch.mean_power - raw.welch)) < 1e-12, ...
    'The raw Welch mean is not the average of the windows.');
assert(isequal(size(components.raw.welch.window_power), [2, numel(raw.frequency)]), ...
    'The raw Welch spectrogram is not windows by frequency.');
assert(max(abs(components.aperiodic.mean_power - session_summary.fractal)) < 1e-12, ...
    'The aperiodic spectrum is not the IRASA fractal spectrum.');
assert(max(abs(components.periodic.mean_power - session_summary.oscillatory)) < 1e-12, ...
    'The periodic spectrum is not the IRASA oscillatory spectrum.');

fprintf(['IRASA recovered exponent %.3f (target 1.5), fractal log-rms %.3f, ' ...
    '10 Hz oscillatory peak ratio %.2f.\n'], ...
    fit.exponent, fractal_log_rms, local_peak_ratio(separated, 10));
disp('test_irasa passed');
end

function assert_separation_identity(spec)
difference = spec.mixed - spec.fractal - spec.oscillatory;
assert(max(abs(difference(:))) < 1e-10, ...
    'oscillatory is not mixed minus fractal.');
end

function ratio = local_peak_ratio(spec, center_hz)
index = nearest_bin(spec.frequency, center_hz);
neighbors = spec.frequency > center_hz - 2 & spec.frequency < center_hz + 2 & ...
    abs(spec.frequency - center_hz) > 0.6;
ratio = spec.oscillatory(index) / median(abs(spec.oscillatory(neighbors)));
end

function ratio = fractal_peak_ratio(spec, center_hz)
index = nearest_bin(spec.frequency, center_hz);
neighbors = spec.frequency > center_hz - 2 & spec.frequency < center_hz + 2 & ...
    abs(spec.frequency - center_hz) > 0.6;
ratio = spec.fractal(index) / median(spec.fractal(neighbors));
end

function ratio = band_ratio(frequency, power, center_hz)
index = nearest_bin(frequency, center_hz);
neighbors = frequency > center_hz - 2 & frequency < center_hz + 2 & ...
    abs(frequency - center_hz) > 0.6;
ratio = power(index) / median(power(neighbors));
end

function index = nearest_bin(frequency, center_hz)
[~, index] = min(abs(frequency - center_hz));
end
