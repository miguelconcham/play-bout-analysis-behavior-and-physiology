function summary = normalized_acg_stack(counts)
%NORMALIZED_ACG_STACK  Normalize autocorrelograms and average them.
%
% counts is one row per neuron and one column per lag bin. Each row is
% divided by its sum, so the value in a bin is the fraction of that
% neuron's spike pairs. Rows with no pairs are left out. The remaining
% rows are ordered by the lag of their peak.
%
% mean_acg is the mean across neurons. ci_low and ci_high are the 95%
% confidence interval of that mean, mean ± 1.96 standard errors.

summary = empty_summary(size(counts, 2));
if isempty(counts)
    return
end
counts = double(counts);
if isvector(counts)
    counts = counts(:)';
end
n_bins = size(counts, 2);
normalized = nan(size(counts, 1), n_bins);
for row_index = 1:size(counts, 1)
    row = counts(row_index, :);
    row(~isfinite(row)) = 0;
    total = sum(row);
    if total > 0
        normalized(row_index, :) = row / total;
    end
end
keep = all(isfinite(normalized), 2);
normalized = normalized(keep, :);
if isempty(normalized)
    summary = empty_summary(n_bins);
    return
end
[~, peak_bin] = max(normalized, [], 2);
[~, order] = sort(peak_bin);
normalized = normalized(order, :);
summary.normalized = normalized;
summary.n = size(normalized, 1);
summary.mean_acg = mean(normalized, 1);
if summary.n > 1
    standard_error = std(normalized, 0, 1) / sqrt(summary.n);
else
    standard_error = zeros(1, n_bins);
end
summary.ci_low = summary.mean_acg - 1.96 * standard_error;
summary.ci_high = summary.mean_acg + 1.96 * standard_error;
end

function summary = empty_summary(n_bins)
if nargin < 1 || isempty(n_bins) || n_bins < 1
    n_bins = 0;
end
summary.normalized = zeros(0, n_bins);
summary.n = 0;
summary.mean_acg = nan(1, n_bins);
summary.ci_low = nan(1, n_bins);
summary.ci_high = nan(1, n_bins);
end
