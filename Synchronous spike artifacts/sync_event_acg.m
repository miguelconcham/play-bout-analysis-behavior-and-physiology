function [centers, counts] = sync_event_acg(spike_times, acg_window, acg_bin)
%SYNC_EVENT_ACG  Autocorrelogram counts for one set of spike times.
%
% spike_times are in seconds. The zero lag is removed by autocorrelogram.
% centers are bin centers in seconds and counts are the number of spike
% pairs in each bin.

if nargin < 2 || isempty(acg_window)
    acg_window = [-0.05 0.05];
end
if nargin < 3 || isempty(acg_bin)
    acg_bin = 0.001;
end
spike_times = sort(double(spike_times(:)));
spike_times = spike_times(isfinite(spike_times));
edges = acg_window(1):acg_bin:acg_window(2);
if numel(edges) < 2
    edges = [acg_window(1) acg_window(2)];
end
lags = [];
if numel(spike_times) >= 2
    lags = autocorrelogram(spike_times, acg_window);
end
counts = bin_counts(lags, edges);
centers = (edges(1:end - 1) + edges(2:end)) / 2;
end

function counts = bin_counts(values, edges)
n_bins = numel(edges) - 1;
counts = zeros(1, n_bins);
if isempty(values) || n_bins < 1
    return
end
values = double(values(:)');
width = edges(2) - edges(1);
if ~(width > 0)
    return
end
index = floor((values - edges(1)) / width) + 1;
index(values == edges(end)) = n_bins;
keep = isfinite(index) & index >= 1 & index <= n_bins;
index = index(keep);
for value_index = 1:numel(index)
    counts(index(value_index)) = counts(index(value_index)) + 1;
end
end
