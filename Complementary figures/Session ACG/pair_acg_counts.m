function counts = pair_acg_counts(spikes, bin_edges)
%PAIR_ACG_COUNTS  Pair counts of a spike train on fixed autocorrelogram bins.
%
% spikes are in seconds. Each unordered pair is counted once at +lag and
% once at -lag, which is the same count autocorrelogram returns. A lag of
% exactly zero is left out. bin_edges must be evenly spaced and cover both
% the positive and the negative lags.

spikes = double(spikes(:));
spikes = spikes(isfinite(spikes));
bin_edges = double(bin_edges(:))';
n_bins = numel(bin_edges) - 1;
counts = zeros(1, n_bins);
if numel(spikes) < 2 || n_bins < 1
    return
end
gaps = diff(bin_edges);
bin_w = gaps(1);
if any(abs(gaps - bin_w) > 1e-9) || ~(bin_w > 0)
    error('pair_acg_counts:edges', 'ACG bin edges must be evenly spaced and increasing.');
end
spikes = sort(spikes);
limit = bin_edges(end);
n = numel(spikes);
lags = zeros(n, 1);
n_lags = 0;
right = 1;
for i = 1:n - 1
    if right < i
        right = i;
    end
    t_i = spikes(i);
    while right < n && spikes(right + 1) - t_i <= limit
        right = right + 1;
    end
    n_new = right - i;
    if n_new < 1
        continue
    end
    if n_lags + n_new > numel(lags)
        lags(max(n_lags + n_new, 2 * numel(lags)), 1) = 0; %#ok<AGROW>
    end
    lags(n_lags + 1:n_lags + n_new) = spikes(i + 1:right) - t_i;
    n_lags = n_lags + n_new;
end
if n_lags < 1
    return
end
lags = lags(1:n_lags);
lags = lags(lags > 0);
if isempty(lags)
    return
end
counts = histcounts(lags, bin_edges) + histcounts(-lags, bin_edges);
end
