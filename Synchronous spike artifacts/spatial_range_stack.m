function summary = spatial_range_stack(distances, ranges, centers)
%SPATIAL_RANGE_STACK  Stack peak-to-peak range against distance on the probe.
%
% distances and ranges are cell arrays, one neuron each. Distance is the
% signed Euclidean distance from that neuron's channel, in micrometers.
% Each neuron is divided by the range on its own channel, so the channel
% sits at 1. Samples are then averaged in bins centered on `centers`.
% Rows whose own channel has no positive range are left out. The remaining
% rows are ordered by the distance of their peak.
%
% mean_range is the mean across neurons in each bin, ignoring bins a
% neuron does not cover. ci_low and ci_high are the 95% confidence
% interval of that mean, mean ± 1.96 standard errors.

if nargin < 3 || isempty(centers)
    centers = -400:20:400;
end
centers = double(centers(:)');
summary = empty_summary(centers);
if isempty(distances) || numel(centers) < 2
    return
end
width = centers(2) - centers(1);
if ~(width > 0) || any(abs(diff(centers) - width) > 1e-6)
    summary.reason = 'distance centers must be evenly spaced';
    return
end
n_neurons = numel(distances);
n_bins = numel(centers);
binned = nan(n_neurons, n_bins);
for neuron_index = 1:n_neurons
    [distance, range] = average_sites(distances{neuron_index}, ranges{neuron_index});
    if isempty(distance)
        continue
    end
    [~, nearest] = min(abs(distance));
    reference = range(nearest);
    if abs(distance(nearest)) > width / 2 || ~isfinite(reference) || reference <= 0
        continue
    end
    binned(neuron_index, :) = bin_ranges(distance, range / reference, centers, width);
end
keep = any(isfinite(binned), 2) & isfinite(binned(:, center_bin(centers)));
binned = binned(keep, :);
if isempty(binned)
    return
end
peak_bin = zeros(size(binned, 1), 1);
for row_index = 1:size(binned, 1)
    row = binned(row_index, :);
    row(~isfinite(row)) = -inf;
    [~, peak_bin(row_index)] = max(row);
end
[~, order] = sort(peak_bin);
binned = binned(order, :);
summary.normalized = binned;
summary.n = size(binned, 1);
[summary.mean_range, summary.ci_low, summary.ci_high] = mean_ci(binned);
end

function values = bin_ranges(distance, range, centers, width)
values = nan(1, numel(centers));
for bin_index = 1:numel(centers)
    left = centers(bin_index) - width / 2;
    right = centers(bin_index) + width / 2;
    if bin_index < numel(centers)
        in_bin = distance >= left & distance < right;
    else
        in_bin = distance >= left & distance <= right;
    end
    samples = range(in_bin);
    samples = samples(isfinite(samples));
    if ~isempty(samples)
        values(bin_index) = mean(samples);
    end
end
end

function index = center_bin(centers)
[~, index] = min(abs(centers));
end

function [distance_out, range_out] = average_sites(distance, range)
distance = double(distance(:));
range = double(range(:));
if numel(distance) ~= numel(range)
    distance_out = zeros(0, 1);
    range_out = zeros(0, 1);
    return
end
keep = isfinite(distance) & isfinite(range);
distance = distance(keep);
range = range(keep);
if isempty(distance)
    distance_out = zeros(0, 1);
    range_out = zeros(0, 1);
    return
end
[distance_out, ~, group] = unique(distance);
range_out = zeros(size(distance_out));
for group_index = 1:numel(distance_out)
    range_out(group_index) = mean(range(group == group_index));
end
end

function [mean_range, ci_low, ci_high] = mean_ci(matrix)
n_bins = size(matrix, 2);
mean_range = nan(1, n_bins);
ci_low = nan(1, n_bins);
ci_high = nan(1, n_bins);
for bin_index = 1:n_bins
    values = matrix(:, bin_index);
    values = values(isfinite(values));
    n_values = numel(values);
    if n_values < 1
        continue
    end
    mean_range(bin_index) = mean(values);
    if n_values > 1
        standard_error = std(values, 0) / sqrt(n_values);
    else
        standard_error = 0;
    end
    ci_low(bin_index) = mean_range(bin_index) - 1.96 * standard_error;
    ci_high(bin_index) = mean_range(bin_index) + 1.96 * standard_error;
end
end

function summary = empty_summary(centers)
n_bins = numel(centers);
summary.centers = centers;
summary.normalized = zeros(0, n_bins);
summary.n = 0;
summary.mean_range = nan(1, n_bins);
summary.ci_low = nan(1, n_bins);
summary.ci_high = nan(1, n_bins);
summary.reason = '';
end
