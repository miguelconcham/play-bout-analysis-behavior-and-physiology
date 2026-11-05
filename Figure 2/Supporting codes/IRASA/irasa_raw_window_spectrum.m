function raw = irasa_raw_window_spectrum(signal, sampling_rate)
%IRASA_RAW_WINDOW_SPECTRUM  Raw 10 s spectrum used before FOOOF.
%
% Same two estimators as GENERATE_PLAY_SESSION_LPAG_SPECTRA, computed
% on one non-overlapping window before any aperiodic fit:
%   welch   - one periodic Hann segment, Welch density scaling
%   fourier - rectangular-window periodogram
% Both are one-sided power spectral densities in power per Hz. DC is
% removed by a constant detrend. The returned grid is df:df:40 Hz, with
% df = sampling_rate / numel(signal). A 10 s window at 2500 Hz is
% 0.1:0.1:40 Hz, the grid passed to fooof_matlab.

if nargin < 2
    error('irasa_raw_window_spectrum:input', ...
        'signal and sampling_rate are required.');
end
if ~isscalar(sampling_rate) || ~isfinite(sampling_rate) || sampling_rate <= 0
    error('irasa_raw_window_spectrum:rate', ...
        'sampling_rate must be a positive finite scalar.');
end
signal = double(signal(:));
n_samples = numel(signal);
if n_samples < 4 || mod(n_samples, 2) ~= 0 || any(~isfinite(signal))
    error('irasa_raw_window_spectrum:signal', ...
        'signal must be a finite even-length vector of at least 4 samples.');
end

signal = signal - mean(signal);
df = sampling_rate / n_samples;
n_frequency = round(40 / df);
if n_frequency < 1 || n_frequency >= n_samples / 2
    error('irasa_raw_window_spectrum:grid', ...
        'The window cannot represent a spectrum through 40 Hz.');
end
frequency = (1:n_frequency).' * df;

welch_window = 0.5 * (1 - cos(2 * pi * (0:n_samples - 1).' / n_samples));
welch = onesided_density(signal .* welch_window, sampling_rate, ...
    sum(welch_window .^ 2), n_frequency);
fourier = onesided_density(signal, sampling_rate, n_samples, n_frequency);

raw = struct();
raw.frequency = frequency;
raw.welch = welch;
raw.fourier = fourier;
end

function power = onesided_density(windowed, sampling_rate, window_energy, n_frequency)
transformed = fft(windowed);
n_positive = numel(windowed) / 2 + 1;
power = abs(transformed(1:n_positive)) .^ 2 / (sampling_rate * window_energy);
power(2:end - 1) = 2 * power(2:end - 1);
% Sample 1 is DC. The FOOOF grid starts at the first positive bin.
power = power(2:n_frequency + 1);
end
