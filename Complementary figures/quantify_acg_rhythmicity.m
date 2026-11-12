function R = quantify_acg_rhythmicity(acg, lags_s, varargin)
%QUANTIFY_ACG_RHYTHMICITY  Oscillatory structure of a positive-lag ACG.
%
% acg is one autocorrelogram. lags_s is the matching positive lags in
% seconds. The spectral score adapts the Muresan et al. (2008) oscillation
% score. It is not an exact reproduction.
%
% Name-value parameters:
%   TroughWindows      Nx2 trough windows, seconds
%   PeakWindows        Nx2 peak windows, seconds
%   BaselineWindow     tail interval for the baseline, seconds
%   FrequencyBand      spectral band, Hz
%   FastSmoothSD       Gaussian smoothing SD, seconds
%   CalculateSpectrum  true to compute the frequency-domain score

p = inputParser;
addParameter(p, 'TroughWindows', [0.05 0.15; 0.35 0.45]);
addParameter(p, 'PeakWindows', [0.25 0.35; 0.55 0.65]);
addParameter(p, 'BaselineWindow', [0.8 1.0]);
addParameter(p, 'FrequencyBand', [2 5]);
addParameter(p, 'FastSmoothSD', 0.002);
addParameter(p, 'CalculateSpectrum', true);
parse(p, varargin{:});
opts = p.Results;

acg = acg(:);
lags_s = lags_s(:);
if numel(acg) ~= numel(lags_s)
    error('quantify_acg_rhythmicity:axis', ...
        'ACG and lags must have the same number of elements.');
end
if numel(lags_s) < 10
    error('quantify_acg_rhythmicity:axis', 'Insufficient number of ACG bins.');
end
if ~all(isfinite(lags_s))
    error('quantify_acg_rhythmicity:axis', 'Lags must contain finite values.');
end
if ~all(diff(lags_s) > 0)
    error('quantify_acg_rhythmicity:axis', 'Lags must be strictly increasing.');
end
if ~(lags_s(1) > 0)
    error('quantify_acg_rhythmicity:axis', 'Lags must be positive.');
end
if size(opts.TroughWindows, 2) ~= 2 || size(opts.PeakWindows, 2) ~= 2
    error('quantify_acg_rhythmicity:window', 'Windows must have two columns.');
end
if size(opts.TroughWindows, 1) ~= size(opts.PeakWindows, 1) ...
        || size(opts.TroughWindows, 1) < 2
    error('quantify_acg_rhythmicity:window', ...
        'At least two matching peak and trough windows are required.');
end
if ~all(opts.TroughWindows(:, 1) < opts.TroughWindows(:, 2)) ...
        || ~all(opts.PeakWindows(:, 1) < opts.PeakWindows(:, 2))
    error('quantify_acg_rhythmicity:window', ...
        'Window start must be smaller than window end.');
end
if ~(opts.FrequencyBand(1) > 0 && opts.FrequencyBand(2) > opts.FrequencyBand(1))
    error('quantify_acg_rhythmicity:band', 'Invalid frequency band.');
end

baseline_idx = lags_s >= opts.BaselineWindow(1) ...
    & lags_s <= opts.BaselineWindow(2) & isfinite(acg);
if any(baseline_idx)
    R.Baseline = median(acg(baseline_idx));
else
    R.Baseline = NaN;
end

n_cycles = size(opts.TroughWindows, 1);
R.PeakValues = nan(n_cycles, 1);
R.TroughValues = nan(n_cycles, 1);
R.PeakLags = nan(n_cycles, 1);
R.TroughLags = nan(n_cycles, 1);
for k = 1:n_cycles
    trough_idx = lags_s > opts.TroughWindows(k, 1) ...
        & lags_s < opts.TroughWindows(k, 2) & isfinite(acg);
    peak_idx = lags_s > opts.PeakWindows(k, 1) ...
        & lags_s < opts.PeakWindows(k, 2) & isfinite(acg);
    if any(trough_idx)
        trough_values = acg(trough_idx);
        trough_lags = lags_s(trough_idx);
        [R.TroughValues(k), i_min] = min(trough_values);
        R.TroughLags(k) = trough_lags(i_min);
    end
    if any(peak_idx)
        peak_values = acg(peak_idx);
        peak_lags = lags_s(peak_idx);
        [R.PeakValues(k), i_max] = max(peak_values);
        R.PeakLags(k) = peak_lags(i_max);
    end
end

R.ModulationPerCycle = R.PeakValues - R.TroughValues;
R.FirstModulation = R.ModulationPerCycle(1);
R.SecondModulation = R.ModulationPerCycle(2);
if abs(R.FirstModulation) > eps
    R.Persistence = R.SecondModulation / R.FirstModulation;
else
    R.Persistence = NaN;
end

if isfinite(R.Baseline) && abs(R.Baseline) > eps
    R.BaselineNormalizedFirstModulation = R.FirstModulation / abs(R.Baseline);
    R.BaselineNormalizedSecondModulation = R.SecondModulation / abs(R.Baseline);
else
    R.BaselineNormalizedFirstModulation = NaN;
    R.BaselineNormalizedSecondModulation = NaN;
end

R.PeakInterval = R.PeakLags(2) - R.PeakLags(1);
R.TroughInterval = R.TroughLags(2) - R.TroughLags(1);
if R.PeakInterval > 0
    R.PeakBasedFrequencyHz = 1 / R.PeakInterval;
else
    R.PeakBasedFrequencyHz = NaN;
end
if R.TroughInterval > 0
    R.TroughBasedFrequencyHz = 1 / R.TroughInterval;
else
    R.TroughBasedFrequencyHz = NaN;
end

R.OscillationScore = NaN;
R.DominantFrequencyHz = NaN;
R.FrequencyResolutionHz = NaN;
R.SpectrumFrequencyHz = [];
R.SpectrumMagnitude = [];
R.PeaklessACG = [];
R.MirroredLags = [];
if ~opts.CalculateSpectrum
    return
end

dt = median(diff(lags_s));
uniform = max(abs(diff(lags_s) - dt)) < max(1e-10, dt * 1e-4);
valid_idx = find(isfinite(acg));
if numel(valid_idx) < 10
    return
end
first_idx = valid_idx(1);
last_idx = valid_idx(end);
y = acg(first_idx:last_idx);
t = lags_s(first_idx:last_idx);
if ~(uniform && all(isfinite(y)))
    return
end

fs = 1 / dt;
if ~(opts.FrequencyBand(2) < fs / 2)
    error('quantify_acg_rhythmicity:band', ...
        'Frequency band must be below the Nyquist frequency.');
end
y_fast = gaussian_smooth(y, opts.FastSmoothSD / dt);
y_mirror = [flipud(y_fast(2:end)); y_fast];
t_mirror = [-flipud(t(2:end)); t];
cut_lag = R.TroughLags(1);
cut_value = interp1(t, y_fast, cut_lag, 'linear');
y_peakless = y_mirror;
if isfinite(cut_lag) && isfinite(cut_value)
    y_peakless(abs(t_mirror) < cut_lag) = cut_value;
end

n_mirror = numel(y_peakless);
n = (0:n_mirror - 1)';
if n_mirror > 1
    window = 0.42 - 0.5 * cos(2 * pi * n / (n_mirror - 1)) ...
        + 0.08 * cos(4 * pi * n / (n_mirror - 1));
else
    window = 1;
end
spectrum = fft(y_peakless .* window);
n_freq = floor(n_mirror / 2) + 1;
magnitude = abs(spectrum(1:n_freq));
frequencies = (0:n_freq - 1)' * fs / n_mirror;
band_idx = frequencies >= opts.FrequencyBand(1) ...
    & frequencies <= opts.FrequencyBand(2);
if any(band_idx)
    band_magnitude = magnitude(band_idx);
    band_frequencies = frequencies(band_idx);
    [max_magnitude, i_max] = max(band_magnitude);
    R.DominantFrequencyHz = band_frequencies(i_max);
    R.FrequencyResolutionHz = fs / n_mirror;
    mean_magnitude = mean(magnitude);
    if mean_magnitude > 0
        R.OscillationScore = max_magnitude / mean_magnitude;
    end
end
R.SpectrumFrequencyHz = frequencies;
R.SpectrumMagnitude = magnitude;
R.PeaklessACG = y_peakless;
R.MirroredLags = t_mirror;
end


function y_smooth = gaussian_smooth(y, sigma_bins)
if sigma_bins <= 0
    y_smooth = y;
    return
end
half_width = max(1, ceil(4 * sigma_bins));
x = (-half_width:half_width)';
kernel = exp(-0.5 * (x / sigma_bins) .^ 2);
kernel = kernel / sum(kernel);
y_smooth = conv(y, kernel, 'same');
weight = conv(ones(size(y)), kernel, 'same');
y_smooth = y_smooth ./ weight;
end
