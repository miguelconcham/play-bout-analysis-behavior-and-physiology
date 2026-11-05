function power = theta_envelope_band_power(signal, rate, band, window_seconds)
%THETA_ENVELOPE_BAND_POWER  Z-scored theta amplitude from a short-time spectrum.
%
% Each sample is the mean log10 magnitude of a Hamming spectrum inside
% `band`, using a centered window of window_seconds. Samples without a full
% window are NaN. The finite center is z-scored inside this interval, on the
% same scale as the Hilbert envelope.

signal = double(signal(:));
n_samples = numel(signal);
power = nan(n_samples, 1);
n_window = round(window_seconds * rate);
if n_window < 8 || n_samples < n_window
    return
end
if mod(n_window, 2) == 0
    n_window = n_window + 1;
end
half = floor(n_window / 2);
window = theta_envelope_hamming_window(n_window);
nfft = 2 ^ nextpow2(max(n_window, round(rate / 0.5)));
frequencies = (0:nfft - 1).' * (rate / nfft);
keep = frequencies >= band(1) & frequencies <= band(2) & frequencies <= rate / 2;
if ~any(keep)
    return
end
for sample_n = (half + 1):(n_samples - half)
    snippet = signal(sample_n - half:sample_n + half);
    if any(~isfinite(snippet))
        continue
    end
    snippet = (snippet - mean(snippet)) .* window;
    spectrum = abs(fft(snippet, nfft));
    power(sample_n) = mean(log10(spectrum(keep) + realmin));
end
finite = isfinite(power);
center = power(finite);
scale = std(center);
if numel(center) < 2 || ~isfinite(scale) || scale < 1e-6
    power(finite) = 0;
    return
end
power(finite) = (center - mean(center)) / scale;
end
