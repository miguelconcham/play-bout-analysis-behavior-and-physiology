function [bin_center, median_height, lower_ci, upper_ci] = ...
        frequency_bin_median_interval( ...
        centers, heights, frequency_limits, bin_hz, confidence)
%FREQUENCY_BIN_MEDIAN_INTERVAL  Median height and its CI in frequency bins.
%
% Bins are bin_hz wide and aligned to multiples of bin_hz, covering
% frequency_limits. median_height is the median FOOOF height of the bumps
% in that bin. lower_ci and upper_ci are the order-statistic confidence
% interval for that median: with z the normal quantile of confidence, the
% ranks are n/2 +/- z*sqrt(n)/2. A bin with too few bumps to support that
% interval has NaN bounds and still keeps its median. Empty bins are omitted.

if nargin < 5 || isempty(confidence)
    confidence = 0.95;
end
if ~isscalar(bin_hz) || ~isfinite(bin_hz) || bin_hz <= 0
    error('regionGaussianCenters:bin', ...
        'histogram_bin_hz must be a positive bin width.');
end
if ~isscalar(confidence) || ~isfinite(confidence) || ...
        confidence <= 0 || confidence >= 1
    error('regionGaussianCenters:confidence', ...
        'amplitude_confidence_level must lie in (0, 1).');
end

centers = double(centers(:));
heights = double(heights(:));
if numel(centers) ~= numel(heights)
    error('regionGaussianCenters:bin', ...
        'Centers and heights must have the same length.');
end
finite = isfinite(centers) & isfinite(heights);
centers = centers(finite);
heights = heights(finite);

frequency_limits = double(frequency_limits(:).');
if numel(frequency_limits) < 2 || ~all(isfinite(frequency_limits(1:2)))
    error('regionGaussianCenters:frequency', ...
        'frequency_limits must be [low high].');
end
span = frequency_limits(1:2);
if ~isempty(centers)
    span(1) = min(span(1), min(centers));
    span(2) = max(span(2), max(centers));
end
edges = frequency_bin_edges(span, bin_hz);
[~, ~, bin_index] = histcounts(centers, edges);
occupied = unique(bin_index(bin_index > 0));
occupied = occupied(:);
n_bins = numel(occupied);
bin_center = zeros(n_bins, 1);
median_height = zeros(n_bins, 1);
lower_ci = nan(n_bins, 1);
upper_ci = nan(n_bins, 1);
for bin_n = 1:n_bins
    in_bin = bin_index == occupied(bin_n);
    values = sort(heights(in_bin));
    bin_center(bin_n) = (edges(occupied(bin_n)) + ...
        edges(occupied(bin_n) + 1)) / 2;
    median_height(bin_n) = median(values);
    [lower_ci(bin_n), upper_ci(bin_n)] = median_order_interval( ...
        values, confidence);
end
end

function edges = frequency_bin_edges(span, bin_hz)
first_edge = floor(span(1) / bin_hz) * bin_hz;
last_edge = ceil(span(2) / bin_hz) * bin_hz;
if last_edge <= first_edge
    last_edge = first_edge + bin_hz;
end
n_bins = round((last_edge - first_edge) / bin_hz);
edges = first_edge + (0:n_bins).' * bin_hz;
end

function [lower_bound, upper_bound] = median_order_interval( ...
        sorted_values, confidence)
lower_bound = NaN;
upper_bound = NaN;
n_values = numel(sorted_values);
if n_values < 2
    return
end
z = sqrt(2) * erfinv(confidence);
lower_rank = floor((n_values + 1) / 2 - z * sqrt(n_values) / 2);
upper_rank = n_values - lower_rank + 1;
if lower_rank < 1 || upper_rank > n_values || lower_rank > upper_rank
    return
end
lower_bound = sorted_values(lower_rank);
upper_bound = sorted_values(upper_rank);
end
