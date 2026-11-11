function powers = channel_window_delta(signal, sample_rate, starts, ends, options)
%CHANNEL_WINDOW_DELTA  FOOOF delta and Hilbert amplitude in each window.
%
% The Hilbert envelope is the absolute value of the delta-band analytic
% signal, filtered once on the whole channel. Each window then contributes
% its mean envelope. FOOOF delta is fooof_band_power of that window's
% spectrum. Sample k is treated as time k/sample_rate, matching the PAG
% LFP time base. Windows are [start, end).

signal = double(signal(:));
starts = double(starts(:));
ends = double(ends(:));
n_windows = numel(starts);
powers.fooof_delta = nan(1, n_windows);
powers.hilbert_amp = nan(1, n_windows);
if n_windows < 1 || isempty(signal) || numel(starts) ~= numel(ends)
    return
end
band = double(options.delta_band(:)');
order = double(options.filter_order);
if numel(signal) <= 3 * order
    order = max(8, 2 * floor((numel(signal) - 1) / 6));
end
filter_design = designfilt('bandpassfir', ...
    'FilterOrder', order, ...
    'CutoffFrequency1', band(1), ...
    'CutoffFrequency2', band(2), ...
    'SampleRate', sample_rate, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');
filtered = filtfilt(filter_design.Coefficients, 1, signal);
envelope = abs(hilbert(filtered));
sample_start = max(1, ceil(starts * sample_rate - 1e-9));
sample_end = min(numel(signal), floor(ends * sample_rate - 1e-9));
for w = 1:n_windows
    if sample_end(w) < sample_start(w)
        continue
    end
    snippet = signal(sample_start(w):sample_end(w));
    powers.hilbert_amp(w) = mean(envelope(sample_start(w):sample_end(w)));
    powers.fooof_delta(w) = snippet_fooof_delta(snippet, sample_rate, options);
end
end


function delta = snippet_fooof_delta(snippet, sample_rate, options)
delta = NaN;
snippet = snippet(isfinite(snippet));
if numel(snippet) < sample_rate || all(snippet == snippet(1))
    return
end
segment = min(numel(snippet), round(4 * sample_rate));
if segment < round(sample_rate)
    return
end
window = hamming(segment);
[power, frequency] = pwelch(snippet, window, floor(segment / 2), segment, sample_rate);
delta = fooof_band_power(frequency, power, options.delta_band, options);
end
