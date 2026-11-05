function spec = irasa_separate(signal, sampling_rate, options)
%IRASA_SEPARATE  Fractal and oscillatory spectra by IRASA.
%
% Wen H, Liu Z (2016) Separating fractal and oscillatory components in
% the power spectrum of neurophysiological signal. Brain Topography
% 29:13-26. doi:10.1007/s10548-015-0448-0
%
% The measured signal is treated as y(t) = f(t) + x(t), their equation (1).
% For each resampling factor h the geometric mean of the auto-spectra of
% the upsampled and downsampled window is equation (6). The median across
% h is the fractal spectrum, equation (8). Subtracting that median from
% the original spectrum, after averaging both across windows, is the
% oscillatory spectrum implied by equation (9).
%
% Section 2.4, in order:
%   1. Ten windows, each 90% of the epoch, evenly spaced.
%   2. Hanning-tapered FFT. NFFT is twice the smallest power of two
%      greater than the window length, and is increased if a requested
%      h would make the upsampled window longer than that.
%   3. Upsample by cubic spline. Downsample by an anti-alias low-pass
%      and then a cubic spline. Same Hanning FFT and same NFFT.
%   4. Geometric mean of each pair, then the median across h, per window.
%   5. Average the fractal spectra and the original spectra across windows.
%   6. Subtract.
%
% SIGNAL is one channel (row or column) or a time-by-channel matrix.
% SAMPLING_RATE is in Hz. OPTIONS is an irasa_options struct.
%
% The spectra are one-sided periodograms on the shared NFFT grid. Bins
% other than DC and Nyquist are doubled, and each spectrum is divided by
% the sum of squared Hann weights. The same definition is used for the
% original window and both resampled windows, so the subtraction is on
% one power scale. Negative values in the oscillatory spectrum are kept.
%
% Frequencies above sampling_rate / (2 * max(h)) are dropped unless
% options.frequency_range says otherwise. Above that limit the upsampled
% branch is being read past the original Nyquist (Discussion, "Choosing
% the resampling factors").

if nargin < 2
    error('irasa_separate:input', ...
        'signal and sampling_rate are required.');
end
if nargin < 3
    options = struct();
end
options = irasa_options(options);
signal = orient_signal(signal);
sampling_rate = validate_sampling_rate(sampling_rate);

[n_samples, n_channels] = size(signal);
h_values = options.h_values(:).';
[segment_length, segment_starts] = segment_layout( ...
    n_samples, options.n_segments, options.segment_fraction);
nfft = choose_nfft(segment_length, h_values);
n_freq = nfft / 2 + 1;
frequency = (0:n_freq - 1).' * (sampling_rate / nfft);

mixed_accum = zeros(n_freq, n_channels);
fractal_accum = zeros(n_freq, n_channels);
n_segments = numel(segment_starts);
for i_segment = 1:n_segments
    start = segment_starts(i_segment);
    segment = signal(start:start + segment_length - 1, :);
    if options.detrend
        segment = linear_detrend(segment);
    end
    mixed_accum = mixed_accum + irasa_periodogram(segment, nfft);
    fractal_accum = fractal_accum + fractal_spectrum( ...
        segment, sampling_rate, h_values, nfft, logical(options.antialias));
end
mixed = mixed_accum / n_segments;
fractal = fractal_accum / n_segments;
oscillatory = mixed - fractal;

[frequency, mixed, fractal, oscillatory] = crop_frequency( ...
    frequency, mixed, fractal, oscillatory, sampling_rate, h_values, ...
    options.frequency_range);

spec = struct();
spec.frequency = frequency;
spec.sampling_rate = sampling_rate;
spec.mixed = mixed;
spec.fractal = fractal;
spec.oscillatory = oscillatory;
spec.h_values = h_values(:);
spec.n_segments = n_segments;
spec.segment_fraction = options.segment_fraction;
spec.segment_length = segment_length;
spec.segment_starts = segment_starts(:);
spec.nfft = nfft;
spec.detrend = logical(options.detrend);
spec.antialias = logical(options.antialias);
spec.reference = ['Wen H, Liu Z (2016) Brain Topography 29:13-26. ' ...
    'Section 2.4, equations (6), (8) and (9).'];
end

function fractal = fractal_spectrum(segment, sampling_rate, h_values, nfft, antialias)
n_h = numel(h_values);
n_freq = nfft / 2 + 1;
n_channels = size(segment, 2);
pair_spectra = zeros(n_freq, n_channels, n_h);
for i_h = 1:n_h
    upsampled = resample_cubic(segment, h_values(i_h), sampling_rate, ...
        'up', antialias);
    downsampled = resample_cubic(segment, h_values(i_h), sampling_rate, ...
        'down', antialias);
    up_power = irasa_periodogram(upsampled, nfft);
    down_power = irasa_periodogram(downsampled, nfft);
    pair_spectra(:, :, i_h) = sqrt(up_power .* down_power);
end
fractal = median(pair_spectra, 3);
end

function resampled = resample_cubic(segment, factor, sampling_rate, direction, antialias)
% Cubic spline on the original sample index. Upsampling stretches the
% series by h. Downsampling compresses it by h after the anti-alias
% low-pass, whose stopband is the Nyquist rate of the compressed series.
n_samples = size(segment, 1);
if strcmp(direction, 'down')
    if antialias
        segment = lowpass_before_downsample(segment, sampling_rate, ...
            sampling_rate / (2 * factor));
    end
    n_out = max(4, round(n_samples / factor));
else
    n_out = max(4, round(n_samples * factor));
end
source_time = (0:n_samples - 1).';
resampled_time = linspace(0, n_samples - 1, n_out).';
resampled = interp1(source_time, segment, resampled_time, 'spline');
end

function filtered = lowpass_before_downsample(segment, sampling_rate, stop_hz)
% Raised-cosine low-pass. Full pass below stop_hz/1.15, full stop at
% stop_hz, so the cubic spline is not asked to represent frequencies
% the downsampled series cannot carry. The 15% transition matches the
% width used by the authors' released FFT filter.
[n_samples, ~] = size(segment);
nfft = 2 ^ nextpow2(n_samples);
frequencies = (0:nfft - 1).' * (sampling_rate / nfft);
folded = frequencies;
above_nyquist = folded > sampling_rate / 2;
folded(above_nyquist) = sampling_rate - folded(above_nyquist);
pass_hz = stop_hz / 1.15;
gain = ones(nfft, 1);
gain(folded >= stop_hz) = 0;
in_transition = folded > pass_hz & folded < stop_hz;
weight = (folded(in_transition) - pass_hz) / (stop_hz - pass_hz);
gain(in_transition) = 0.5 * (1 + cos(pi * weight));
transformed = fft(segment, nfft, 1);
filtered = real(ifft(bsxfun(@times, transformed, gain), nfft, 1));
filtered = filtered(1:n_samples, :);
end

function power = irasa_periodogram(segment, nfft)
n_samples = size(segment, 1);
weights = 0.5 * (1 - cos(2 * pi * (0:n_samples - 1).' / n_samples));
windowed = bsxfun(@times, segment, weights);
transformed = fft(windowed, nfft, 1);
n_freq = nfft / 2 + 1;
power = abs(transformed(1:n_freq, :)) .^ 2;
if n_freq > 2
    power(2:n_freq - 1, :) = 2 * power(2:n_freq - 1, :);
end
power = power / sum(weights .^ 2);
end

function [frequency, mixed, fractal, oscillatory] = crop_frequency( ...
        frequency, mixed, fractal, oscillatory, sampling_rate, h_values, ...
        frequency_range)
if isempty(frequency_range)
    fmin = frequency(2);
    fmax = sampling_rate / (2 * max(h_values));
else
    fmin = frequency_range(1);
    fmax = frequency_range(2);
end
keep = frequency >= fmin & frequency <= fmax;
if ~any(keep)
    error('irasa_separate:frequency', ...
        'No FFT bin falls inside the requested frequency range.');
end
frequency = frequency(keep);
mixed = mixed(keep, :);
fractal = fractal(keep, :);
oscillatory = oscillatory(keep, :);
end

function nfft = choose_nfft(segment_length, h_values)
% Section 2.4 step 2. "Greater than" excludes a window whose length is
% already a power of two. The extra factor of two keeps the upsampled
% window inside NFFT for every published h < 2.
next_power = 2 ^ ceil(log2(segment_length));
if next_power <= segment_length
    next_power = next_power * 2;
end
nfft = 2 * next_power;
while nfft <= max(h_values) * segment_length
    nfft = nfft * 2;
end
end

function [segment_length, segment_starts] = segment_layout( ...
        n_samples, n_segments, segment_fraction)
segment_length = round(segment_fraction * n_samples);
if segment_length < 32
    error('irasa_separate:length', ...
        'Each IRASA window needs at least 32 samples.');
end
if segment_length > n_samples
    error('irasa_separate:length', ...
        'The window is longer than the signal.');
end
segment_starts = zeros(n_segments, 1);
if n_segments == 1
    segment_starts(1) = floor((n_samples - segment_length) / 2) + 1;
    return
end
span = n_samples - segment_length;
for i_segment = 1:n_segments
    segment_starts(i_segment) = round( ...
        1 + (i_segment - 1) * span / (n_segments - 1));
end
segment_starts = max(1, min(segment_starts, n_samples - segment_length + 1));
end

function signal = orient_signal(signal)
if ~isnumeric(signal) || ~isreal(signal) || any(~isfinite(signal(:)))
    error('irasa_separate:signal', ...
        'signal must be finite and real.');
end
signal = double(signal);
if isvector(signal)
    signal = signal(:);
end
if size(signal, 1) < 32
    error('irasa_separate:signal', ...
        'signal needs at least 32 samples along the first dimension.');
end
end

function sampling_rate = validate_sampling_rate(sampling_rate)
if ~isscalar(sampling_rate) || ~isnumeric(sampling_rate) || ...
        ~isfinite(sampling_rate) || sampling_rate <= 0
    error('irasa_separate:rate', ...
        'sampling_rate must be a positive finite scalar in Hz.');
end
sampling_rate = double(sampling_rate);
end

function detrended = linear_detrend(segment)
n_samples = size(segment, 1);
time = (0:n_samples - 1).';
time = time - mean(time);
slope = (time' * segment) / (time' * time);
detrended = segment - time * slope;
detrended = detrended - mean(detrended, 1);
end
