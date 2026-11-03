function result = theta_envelope_session_spectrum(envelopes, rate, grids, welch)
%THETA_ENVELOPE_SESSION_SPECTRUM  Welch spectrum of theta-envelope segments.
%
% `envelopes` is a cell of column vectors, one partner-session interval.
% Interval spectra are averaged with weights equal to the number of Welch
% segments. Band power is the mean of that spectrum inside each grid window.
% Observed and surrogate spectra use the same scaling, so the z-score does
% not depend on the Fourier normalization.

if ~iscell(envelopes)
    envelopes = {envelopes};
end
settings = welch_settings(rate, welch);
spectra = [];
frequencies = [];
weights = [];
resolution = [];
for k = 1:numel(envelopes)
    x = longest_finite(envelopes{k});
    if numel(x) < round(4 * rate)
        continue
    end
    window = min(settings.target_window, numel(x));
    noverlap = min(floor(window * welch.overlap), window - 1);
    x = x - mean(x);
    [psd, f] = welch_psd(x, window, noverlap, settings.nfft, rate);
    if isempty(frequencies)
        frequencies = f(:);
    elseif numel(f) ~= numel(frequencies) || any(abs(f(:) - frequencies) > 1e-9)
        error('theta_envelope_session_spectrum:grid', ...
            'Welch frequency grids differ across intervals.');
    end
    spectra = [spectra, psd(:)]; %#ok<AGROW>
    weights(end + 1) = segment_count(numel(x), window, noverlap); %#ok<AGROW>
    resolution(end + 1) = rate / window; %#ok<AGROW>
end

result.frequencies = zeros(0, 1);
result.psd = zeros(0, 1);
result.band_power_test = nan(1, numel(grids.test.centers));
result.band_power_fine = nan(1, numel(grids.fine.centers));
result.resolution_hz = NaN;
result.n_segments = 0;
if isempty(weights)
    return
end
psd = spectra * weights(:) / sum(weights);
keep = frequencies <= max(8, grids.fine.band_hz(2));
result.frequencies = frequencies(keep);
result.psd = psd(keep);
result.band_power_test = theta_envelope_band_mean(psd, frequencies, grids.test);
result.band_power_fine = theta_envelope_band_mean(psd, frequencies, grids.fine);
result.resolution_hz = max(resolution);
result.n_segments = sum(weights);
end

function settings = welch_settings(rate, welch)
settings.nfft = max(8, round(rate / welch.resolution_hz));
settings.target_window = max(8, round(welch.window_seconds * rate));
end

function [psd, frequencies] = welch_psd(x, window, noverlap, nfft, rate)
weights = theta_envelope_hamming_window(window);
step = window - noverlap;
starts = 1:step:(numel(x) - window + 1);
n_one = floor(nfft / 2) + 1;
accum = zeros(n_one, 1);
for s = 1:numel(starts)
    idx = starts(s):(starts(s) + window - 1);
    segment = (x(idx) - mean(x(idx))) .* weights;
    spectrum = fft(segment, nfft);
    accum = accum + abs(spectrum(1:n_one)).^2;
end
psd = accum / max(numel(starts), 1);
frequencies = (0:n_one - 1).' * rate / nfft;
end

function n = segment_count(n_samples, window, noverlap)
step = window - noverlap;
n = floor((n_samples - window) / step) + 1;
n = max(n, 1);
end

function x = longest_finite(x)
x = double(x(:));
finite = isfinite(x);
if ~any(finite)
    x = [];
    return
end
edges = diff([false; finite; false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
[~, best] = max(stops - starts);
x = x(starts(best):stops(best));
end
