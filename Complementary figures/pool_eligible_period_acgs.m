function pooled = pool_eligible_period_acgs(partner_counts, n_eligible, run_length, score, low_percent, high_percent)
%POOL_ELIGIBLE_PERIOD_ACGS  Mean partner count of spikes that can see the bin.
%
% Low and high tails use the same 25/75 rule as pool_locking_period_acgs.
% In each lag bin the pooled value is the sum of partner counts divided
% by the number of spikes that could see that bin. A bin no spike can see
% is missing. It is not recorded as zero, and it is not divided by
% T-|lag|.

if nargin < 5 || isempty(low_percent)
    low_percent = 25;
end
if nargin < 6 || isempty(high_percent)
    high_percent = 75;
end

n_bins = size(partner_counts, 2);
pooled = empty_pooled(n_bins);
if n_bins < 1 || isempty(partner_counts)
    return
end
partner_counts = double(partner_counts);
n_eligible = double(n_eligible);
score = double(score(:));
run_length = double(run_length(:));
if size(n_eligible, 1) ~= size(partner_counts, 1) ...
        || size(n_eligible, 2) ~= n_bins
    error('pool_eligible_period_acgs:width', ...
        'Partner counts and eligible-spike counts must be the same size.');
end
n_runs = size(partner_counts, 1);
if numel(score) ~= n_runs || numel(run_length) ~= n_runs
    error('pool_eligible_period_acgs:rows', ...
        'Score and run length must have one entry per window.');
end

usable = isfinite(run_length) & run_length > 0;
score_ok = usable & isfinite(score);
low_idx = false(n_runs, 1);
high_idx = false(n_runs, 1);
if sum(score_ok) >= 4
    pooled.threshold_low = linear_percentile(score(score_ok), low_percent);
    pooled.threshold_high = linear_percentile(score(score_ok), high_percent);
end
if isfinite(pooled.threshold_low) && isfinite(pooled.threshold_high) ...
        && pooled.threshold_high > pooled.threshold_low
    low_idx = score_ok & score <= pooled.threshold_low;
    high_idx = score_ok & score >= pooled.threshold_high;
end
all_idx = usable & any(isfinite(partner_counts), 2);

[pooled.low, pooled.n_low] = pool_selected(partner_counts, n_eligible, low_idx);
[pooled.high, pooled.n_high] = pool_selected(partner_counts, n_eligible, high_idx);
[pooled.all, pooled.n_all] = pool_selected(partner_counts, n_eligible, all_idx);
pooled.all_index = all_idx;
pooled.score_low = mean_finite(score(low_idx));
pooled.score_high = mean_finite(score(high_idx));
pooled.score_all = mean_finite(score(score_ok));
pooled.low_index = low_idx;
pooled.high_index = high_idx;
end


function [mu, n_used] = pool_selected(partner_counts, n_eligible, idx)
n_bins = size(partner_counts, 2);
mu = nan(1, n_bins);
n_used = sum(idx);
if n_used < 1
    return
end
counts = partner_counts(idx, :);
nel = n_eligible(idx, :);
bad = ~isfinite(nel) | nel < 0;
nel(bad) = 0;
counts(~isfinite(counts) | bad) = 0;
den = sum(nel, 1);
num = sum(counts, 1);
ok = den > 0;
mu(ok) = num(ok) ./ den(ok);
end


function m = mean_finite(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end


function pooled = empty_pooled(n_bins)
pooled.low = nan(1, n_bins);
pooled.high = nan(1, n_bins);
pooled.all = nan(1, n_bins);
pooled.n_low = 0;
pooled.n_high = 0;
pooled.n_all = 0;
pooled.threshold_low = NaN;
pooled.threshold_high = NaN;
pooled.score_low = NaN;
pooled.score_high = NaN;
pooled.score_all = NaN;
pooled.low_index = false(0, 1);
pooled.high_index = false(0, 1);
pooled.all_index = false(0, 1);
end
