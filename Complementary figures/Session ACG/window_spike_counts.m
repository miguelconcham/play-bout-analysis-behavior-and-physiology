function counts = window_spike_counts(spikes, t0, window_s, n_windows)
%WINDOW_SPIKE_COUNTS  Spikes in each adjoining window of length window_s.
%
% Window k is [t0 + (k-1)T, t0 + kT). A spike at the right edge belongs to
% the next window. Spikes outside the covered range are not counted.

n_windows = double(n_windows);
counts = zeros(n_windows, 1);
if n_windows < 1 || ~isfinite(t0) || ~isfinite(window_s) || window_s <= 0
    return
end
spikes = double(spikes(:));
covered_end = t0 + n_windows * window_s;
spikes = spikes(isfinite(spikes) & spikes >= t0 & spikes < covered_end);
if isempty(spikes)
    return
end
bin = floor((spikes - t0) / window_s) + 1;
bin = bin(bin >= 1 & bin <= n_windows);
if isempty(bin)
    return
end
counts = accumarray(bin, 1, [n_windows, 1]);
end
