function result = group_pair_ccg(reference_spikes, target_spikes, windows, bin_edges)
%GROUP_PAIR_CCG  Rate excess of one pooled train around another.
%
% windows are [start, end] in seconds. Intervals that touch are joined, so
% a spike on a shared boundary is not used twice. Reference spikes are
% inset from each end by the longest lag. Every such spike can see the
% full lag axis inside its window, so the histogram is not divided by
% (T - |lag|). Target spikes are every spike in the window.
%
% conditional_rate is target spikes per second around one reference spike.
% target_rate is the target's own rate in the windows that supplied those
% references, weighted by the reference count. excess_rate is conditional
% rate minus target rate. A negative excess means the target fired less
% than its own rate when the reference group fired. rate_ratio is the
% conditional rate divided by that same target rate. One is independence.
% A window shorter than two lag margins adds nothing.

bin_edges = double(bin_edges(:)');
n_bins = numel(bin_edges) - 1;
result = empty_result(bin_edges, n_bins);
if n_bins < 1
    return
end
gaps = diff(bin_edges);
bin_w = gaps(1);
if any(abs(gaps - bin_w) > 1e-9) || ~(bin_w > 0)
    error('group_ccg:bins', 'CCG bin edges must be evenly spaced and increasing.');
end
windows = merge_time_intervals(windows);
if isempty(windows)
    return
end
reference_spikes = finite_sorted(reference_spikes);
target_spikes = finite_sorted(target_spikes);
margin = max(abs(bin_edges(1)), abs(bin_edges(end)));
counts = zeros(1, n_bins);
n_reference = 0;
n_target = 0;
duration = 0;
expected_counts = 0;
n_windows = 0;

for w = 1:size(windows, 1)
    t0 = windows(w, 1);
    t1 = windows(w, 2);
    window_duration = t1 - t0;
    if ~(window_duration > 2 * margin)
        continue
    end
    reference = reference_spikes( ...
        reference_spikes >= (t0 + margin) & reference_spikes <= (t1 - margin));
    if isempty(reference)
        continue
    end
    target = target_spikes(target_spikes >= t0 & target_spikes <= t1);
    these_counts = offset_histogram(reference, target, bin_edges);
    this_rate = numel(target) / window_duration;
    counts = counts + these_counts;
    n_reference = n_reference + numel(reference);
    n_target = n_target + numel(target);
    duration = duration + window_duration;
    expected_counts = expected_counts + numel(reference) * this_rate * bin_w;
    n_windows = n_windows + 1;
end

result.counts = counts;
result.n_reference = n_reference;
result.n_target = n_target;
result.duration = duration;
result.n_windows = n_windows;
if n_reference < 1 || ~(duration > 0)
    result.counts = nan(1, n_bins);
    return
end
result.conditional_rate = counts / (n_reference * bin_w);
result.target_rate = expected_counts / (n_reference * bin_w);
result.excess_rate = result.conditional_rate - result.target_rate;
if result.target_rate > 0
    result.rate_ratio = result.conditional_rate / result.target_rate;
end
end


function result = empty_result(bin_edges, n_bins)
result.counts = nan(1, n_bins);
result.n_reference = 0;
result.n_target = 0;
result.duration = 0;
result.n_windows = 0;
result.conditional_rate = nan(1, n_bins);
result.target_rate = nan;
result.excess_rate = nan(1, n_bins);
result.rate_ratio = nan(1, n_bins);
result.centers = bin_edges(1:end - 1) + diff(bin_edges) / 2;
if n_bins < 1
    result.centers = zeros(1, 0);
end
end


function spikes = finite_sorted(spikes)
spikes = double(spikes(:));
spikes = sort(spikes(isfinite(spikes)));
end


function counts = offset_histogram(reference, target, bin_edges)
% Half-open bins, and the last bin includes its right edge, as histcounts does.
n_bins = numel(bin_edges) - 1;
counts = zeros(1, n_bins);
if isempty(reference) || isempty(target)
    return
end
bin_w = bin_edges(2) - bin_edges(1);
lo = bin_edges(1);
hi = bin_edges(end);
n_target = numel(target);
left = 1;
right = 1;
for i = 1:numel(reference)
    t = reference(i);
    while left <= n_target && target(left) < (t + lo)
        left = left + 1;
    end
    if right < left
        right = left;
    end
    while right <= n_target && target(right) <= (t + hi)
        right = right + 1;
    end
    if right <= left
        continue
    end
    lags = target(left:right - 1) - t;
    bin = floor((lags - lo) / bin_w + 1e-8) + 1;
    bin(bin == n_bins + 1) = n_bins;
    keep = bin >= 1 & bin <= n_bins;
    if any(keep)
        counts = counts + accumarray(bin(keep), 1, [n_bins, 1])';
    end
end
end
