function result = theta_envelope_session_coherence(envelopes, lfps, rate, grids, welch, randomize_segments)
%THETA_ENVELOPE_SESSION_COHERENCE  Magnitude-squared coherence, envelope vs LFP.
%
% Segments from every partner-session interval are pooled before the
% coherence ratio is formed. When randomize_segments is true, each LFP
% segment is circularly shifted by its own lag. One global lag would leave
% magnitude-squared coherence unchanged, so it is not the null.

if nargin < 6 || isempty(randomize_segments)
    randomize_segments = false;
end
if ~iscell(envelopes)
    envelopes = {envelopes};
end
if ~iscell(lfps)
    lfps = {lfps};
end
if numel(envelopes) ~= numel(lfps)
    error('theta_envelope_session_coherence:pair', ...
        'Envelope and LFP interval counts differ.');
end

settings = welch_settings(rate, welch);
pxx = [];
pyy = [];
pxy = [];
n_segments = 0;
resolution = [];
for k = 1:numel(envelopes)
    [x, y] = longest_finite_pair(envelopes{k}, lfps{k});
    if numel(x) < round(4 * rate)
        continue
    end
    window = min(settings.target_window, numel(x));
    noverlap = min(floor(window * welch.overlap), window - 1);
    [seg_xx, seg_yy, seg_xy, n_seg] = segment_periods( ...
        x, y, window, noverlap, settings.nfft, randomize_segments);
    if isempty(pxx)
        pxx = zeros(size(seg_xx));
        pyy = zeros(size(seg_yy));
        pxy = zeros(size(seg_xy));
    end
    pxx = pxx + seg_xx;
    pyy = pyy + seg_yy;
    pxy = pxy + seg_xy;
    n_segments = n_segments + n_seg;
    resolution(end + 1) = rate / window; %#ok<AGROW>
end

n_freq = floor(settings.nfft / 2) + 1;
frequencies = (0:n_freq - 1).' * rate / settings.nfft;
result.frequencies = zeros(0, 1);
result.coherence = zeros(0, 1);
result.band_coherence_test = nan(1, numel(grids.test.centers));
result.band_coherence_fine = nan(1, numel(grids.fine.centers));
result.resolution_hz = NaN;
result.n_segments = 0;
if n_segments == 0
    return
end
coherence = (abs(pxy).^2) ./ (pxx .* pyy);
coherence(~isfinite(coherence)) = NaN;
keep = frequencies <= max(8, grids.fine.band_hz(2));
result.frequencies = frequencies(keep);
result.coherence = coherence(keep);
result.band_coherence_test = theta_envelope_band_mean(coherence, frequencies, grids.test);
result.band_coherence_fine = theta_envelope_band_mean(coherence, frequencies, grids.fine);
result.resolution_hz = max(resolution);
result.n_segments = n_segments;
end

function [pxx, pyy, pxy, n_seg] = segment_periods(x, y, window, noverlap, nfft, randomize_segments)
step = window - noverlap;
starts = 1:step:(numel(x) - window + 1);
n_seg = numel(starts);
weights = theta_envelope_hamming_window(window);
n_one = floor(nfft / 2) + 1;
pxx = zeros(n_one, 1);
pyy = zeros(n_one, 1);
pxy = zeros(n_one, 1);
if n_seg < 1
    return
end
if randomize_segments
    lags = randi([1, window - 1], n_seg, 1);
else
    lags = zeros(n_seg, 1);
end
for s = 1:n_seg
    idx = starts(s):(starts(s) + window - 1);
    xs = x(idx);
    ys = y(idx);
    if lags(s) ~= 0
        ys = circshift(ys, lags(s));
    end
    xs = (xs - mean(xs)) .* weights;
    ys = (ys - mean(ys)) .* weights;
    spec_x = fft(xs, nfft);
    spec_y = fft(ys, nfft);
    spec_x = spec_x(1:n_one);
    spec_y = spec_y(1:n_one);
    pxx = pxx + abs(spec_x).^2;
    pyy = pyy + abs(spec_y).^2;
    pxy = pxy + spec_x .* conj(spec_y);
end
end

function settings = welch_settings(rate, welch)
settings.nfft = max(8, round(rate / welch.resolution_hz));
settings.target_window = max(8, round(welch.window_seconds * rate));
end

function [x, y] = longest_finite_pair(envelope, lfp)
envelope = double(envelope(:));
lfp = double(lfp(:));
n = min(numel(envelope), numel(lfp));
envelope = envelope(1:n);
lfp = lfp(1:n);
finite = isfinite(envelope) & isfinite(lfp);
if ~any(finite)
    x = [];
    y = [];
    return
end
edges = diff([false; finite; false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
[~, best] = max(stops - starts);
idx = starts(best):stops(best);
x = envelope(idx);
y = lfp(idx);
end
