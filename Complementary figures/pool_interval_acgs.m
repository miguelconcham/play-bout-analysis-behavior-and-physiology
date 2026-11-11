function rate = pool_interval_acgs(spikes, starts, ends, bin_edges)
%POOL_INTERVAL_ACGS  Overlap-weighted autocorrelogram across intervals.
%
% Each interval contributes its pair counts and its overlap (T - |lag|).
% The pooled bin is the sum of the counts divided by the sum of the
% overlaps, so a long interval weighs more than a short one. Intervals are
% [start, end) and must not overlap. Spikes on a shared boundary belong to
% the later interval.

bin_edges = double(bin_edges(:))';
n_bins = max(numel(bin_edges) - 1, 0);
rate = nan(1, n_bins);
if n_bins < 1
    return
end
starts = double(starts(:));
ends = double(ends(:));
if numel(starts) ~= numel(ends)
    error('pool_interval_acgs:windows', 'Window starts and ends must have the same length.');
end
if isempty(starts)
    return
end
[starts, order] = sort(starts);
ends = ends(order);
spikes = double(spikes(:));
spikes = sort(spikes(isfinite(spikes)));
bin_w = bin_edges(2) - bin_edges(1);
centers = bin_edges(1:end - 1) + bin_w / 2;
counts = zeros(1, n_bins);
overlap_sum = zeros(1, n_bins);
n_spikes = numel(spikes);
cursor = 1;
for w = 1:numel(starts)
    duration = ends(w) - starts(w);
    if ~isfinite(duration) || duration <= 0
        continue
    end
    while cursor <= n_spikes && spikes(cursor) < starts(w)
        cursor = cursor + 1;
    end
    stop = cursor;
    while stop <= n_spikes && spikes(stop) < ends(w)
        stop = stop + 1;
    end
    if stop > cursor
        counts = counts + pair_acg_counts(spikes(cursor:stop - 1), bin_edges);
    end
    overlap = duration - abs(centers);
    overlap(overlap < bin_w | ~isfinite(overlap)) = 0;
    overlap_sum = overlap_sum + overlap;
    cursor = stop;
end
ok = overlap_sum > 0;
rate(ok) = counts(ok) ./ overlap_sum(ok);
end
