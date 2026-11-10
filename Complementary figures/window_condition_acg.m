function block = window_condition_acg(spike_times, epochs, bin_edges, power)
%WINDOW_CONDITION_ACG  Overlap-weighted autocorrelogram of spikes in epochs.
%
% epochs is N x 2, [start, end). A spike at the end of one epoch belongs
% to the next epoch, so tiled epochs do not count it twice. Power is the
% linear delta power of each epoch and may be omitted. The autocorrelogram
% is pair counts divided by overlap, pooled across epochs the same way as
% the locking-period analysis. Silent epochs stay in the denominator.

if nargin < 4
    power = [];
end
spike_times = double(spike_times(:));
spike_times = spike_times(isfinite(spike_times));
spike_times = sort(spike_times);
bin_edges = double(bin_edges(:)');
n_bins = max(numel(bin_edges) - 1, 0);
block = empty_window_condition(n_bins);
if isempty(epochs) || n_bins < 1
    return
end
epochs = double(epochs(:, 1:2));
keep = all(isfinite(epochs), 2) & epochs(:, 2) > epochs(:, 1);
epochs = epochs(keep, :);
if ~isempty(power)
    power = double(power(:));
    power = power(keep);
end
n_epochs = size(epochs, 1);
if n_epochs < 1
    return
end
range_limits = [bin_edges(1), bin_edges(end)];
acg = nan(n_epochs, n_bins);
duration = epochs(:, 2) - epochs(:, 1);
n_spikes = 0;
for i = 1:n_epochs
    i0 = first_at_or_after(spike_times, epochs(i, 1));
    i1 = first_at_or_after(spike_times, epochs(i, 2));
    in_epoch = spike_times(i0:i1 - 1);
    n_spikes = n_spikes + numel(in_epoch);
    lags = autocorrelogram(in_epoch, range_limits);
    counts = lag_bin_counts(lags, bin_edges);
    acg(i, :) = edge_correct_acg(counts, bin_edges, duration(i));
end
block.acg = pool_selected_window_acgs(acg, bin_edges, duration, true);
block.n_spikes = n_spikes;
block.duration = sum(duration);
block.n_epochs = n_epochs;
block.rate = block.n_spikes / block.duration;
if ~isempty(power)
    finite_power = power(isfinite(power));
    if ~isempty(finite_power)
        block.mean_power = mean(finite_power);
    end
    positive = power(isfinite(power) & power > 0);
    if ~isempty(positive)
        block.mean_power_db = mean(10 * log10(positive));
    end
end
end


function block = empty_window_condition(n_bins)
block.acg = nan(1, n_bins);
block.n_spikes = 0;
block.duration = 0;
block.n_epochs = 0;
block.mean_power = NaN;
block.mean_power_db = NaN;
block.rate = NaN;
block.median_pair_gap_db = NaN;
end


function counts = lag_bin_counts(lags, edges)
if exist('histcounts', 'builtin') == 5 || exist('histcounts', 'file') == 2
    counts = histcounts(lags, edges);
    return
end
counts = zeros(1, numel(edges) - 1);
lags = double(lags(:));
lags = lags(isfinite(lags));
for b = 1:numel(counts)
    if b < numel(counts)
        in_bin = lags >= edges(b) & lags < edges(b + 1);
    else
        in_bin = lags >= edges(b) & lags <= edges(b + 1);
    end
    counts(b) = sum(in_bin);
end
end


function idx = first_at_or_after(x, value)
lo = 1;
hi = numel(x) + 1;
while lo < hi
    mid = floor((lo + hi) / 2);
    if x(mid) < value
        lo = mid + 1;
    else
        hi = mid;
    end
end
idx = lo;
end
