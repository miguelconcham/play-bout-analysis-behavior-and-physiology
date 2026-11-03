function bins = speed_quantile_bins(speed, fraction_limits)
%SPEED_QUANTILE_BINS  Equal-count bins of a session's speed samples.
%
% fraction_limits is [min max], default [0.05 0.10]. Finite samples are
% sorted and split into the smallest number of bins that keeps every bin
% at or under the upper fraction, when that fraction contains at least one
% sample. The upper end (10%) is therefore used whenever the count divides
% evenly. A short vector can force a bin above the upper fraction or below
% the lower one; both flags are returned and every finite sample is still
% assigned.
%
% Equal speed values keep their original order and can fall in neighbouring
% bins so that one speed value does not swallow more than the upper fraction.

if nargin < 2 || isempty(fraction_limits)
    fraction_limits = [0.05 0.10];
end
fraction_limits = double(fraction_limits(:).');
if numel(fraction_limits) ~= 2 || fraction_limits(1) <= 0 || ...
        fraction_limits(2) < fraction_limits(1) || fraction_limits(2) >= 1
    error('speed_quantile_bins:limits', ...
        'Speed-bin fractions must be [min max] with 0 < min <= max < 1.');
end

speed = double(speed(:));
bin_id = nan(size(speed));
finite = isfinite(speed);
n = nnz(finite);
[counts, within_upper, within_lower] = bin_counts( ...
    n, fraction_limits(1), fraction_limits(2));
n_bins = numel(counts);

bins = struct();
bins.requested_fraction_limits = fraction_limits;
bins.rule = ['Equal-count split of finite speed samples. ' ...
    'The bin count is the smallest that keeps each bin at or under ' ...
    'the upper fraction when that fraction contains at least one sample.'];
bins.n_samples = n;
bins.n_bins = n_bins;
bins.within_upper_fraction = within_upper;
bins.within_lower_fraction = within_lower;
bins.within_requested_fractions = within_upper && within_lower;
bins.bin_id = bin_id;
bins.count = counts;
bins.fraction = zeros(n_bins, 1);
bins.speed_min = zeros(n_bins, 1);
bins.speed_max = zeros(n_bins, 1);
bins.speed_median = zeros(n_bins, 1);
bins.speed_mean = zeros(n_bins, 1);
bins.edges = zeros(n_bins + 1, 1);

if n_bins == 0
    bins.fraction = zeros(0, 1);
    bins.speed_min = zeros(0, 1);
    bins.speed_max = zeros(0, 1);
    bins.speed_median = zeros(0, 1);
    bins.speed_mean = zeros(0, 1);
    bins.edges = zeros(0, 1);
    return
end

values = speed(finite);
[~, order] = sortrows([values, (1:n).']);
finite_id = zeros(n, 1);
cursor = 1;
for bin_n = 1:n_bins
    finite_id(order(cursor:cursor + counts(bin_n) - 1)) = bin_n;
    cursor = cursor + counts(bin_n);
end
bin_id(finite) = finite_id;
bins.bin_id = bin_id;
bins.fraction = counts / n;

sorted_values = values(order);
edges = nan(n_bins + 1, 1);
edges(1) = sorted_values(1);
cursor = 1;
for bin_n = 1:n_bins - 1
    cursor = cursor + counts(bin_n);
    edges(bin_n + 1) = (sorted_values(cursor - 1) + sorted_values(cursor)) / 2;
end
edges(end) = sorted_values(end);
bins.edges = edges;

for bin_n = 1:n_bins
    this_speed = values(finite_id == bin_n);
    bins.speed_min(bin_n) = min(this_speed);
    bins.speed_max(bin_n) = max(this_speed);
    bins.speed_median(bin_n) = median(this_speed);
    bins.speed_mean(bin_n) = mean(this_speed);
end
end

function [counts, within_upper, within_lower] = bin_counts(n, min_fraction, max_fraction)
if n < 1
    counts = zeros(0, 1);
    within_upper = false;
    within_lower = false;
    return
end

max_count = floor(max_fraction * n + 1e-12);
forced_above_upper = max_count < 1;
if forced_above_upper
    max_count = 1;
end
n_bins = ceil(n / max_count);
base = floor(n / n_bins);
extra = n - base * n_bins;
counts = base * ones(n_bins, 1);
counts(1:extra) = counts(1:extra) + 1;

fractions = counts / n;
within_upper = ~forced_above_upper && all(fractions <= max_fraction + 1e-9);
within_lower = all(fractions >= min_fraction - 1e-12);
end
