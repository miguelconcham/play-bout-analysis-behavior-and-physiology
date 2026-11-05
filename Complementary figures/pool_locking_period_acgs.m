function pooled = pool_locking_period_acgs(acg, bin_edges, run_length, ppc, already_corrected, low_percent, high_percent)
%POOL_LOCKING_PERIOD_ACGS  Overlap-weighted ACGs for low, high, and all runs.
%
% Each run has one autocorrelogram row and one PPC. Low runs are those at
% or below the low_percent percentile of that neuron's finite PPC values
% (default 25). High runs are at or above high_percent (default 75).
% Percentiles use linear_percentile. The tails are left empty unless at
% least four runs have a finite PPC and the high percentile is strictly
% above the low one. Fewer runs, or a neuron whose runs all share one
% PPC, have no separated low and high periods.
%
% The pooled value in a lag bin is the sum of pair counts divided by the
% sum of overlap (T - |lag|). Averaging per-run rates would let a short
% run count as much as a long one. Bins whose overlap is shorter than one
% bin are ignored. When already_corrected is true, rows are rates
% (counts / overlap) and the counts are recovered before summing.

if nargin < 6 || isempty(low_percent)
    low_percent = 25;
end
if nargin < 7 || isempty(high_percent)
    high_percent = 75;
end
if nargin < 5 || isempty(already_corrected)
    already_corrected = false;
end

bin_edges = bin_edges(:)';
n_bins = max(numel(bin_edges) - 1, 0);
pooled = empty_pooled(n_bins);
if n_bins < 1 || isempty(acg)
    return
end

acg = double(acg);
ppc = double(ppc(:));
run_length = double(run_length(:));
if size(acg, 2) ~= n_bins
    error('pool_locking_period_acgs:width', ...
        'ACG columns (%d) do not match the bin edges.', size(acg, 2));
end
n_runs = size(acg, 1);
if numel(ppc) ~= n_runs || numel(run_length) ~= n_runs
    error('pool_locking_period_acgs:rows', ...
        'PPC and run length must have one entry per ACG row.');
end

usable = isfinite(run_length) & run_length > 0;
ppc_ok = usable & isfinite(ppc);
low_idx = false(n_runs, 1);
high_idx = false(n_runs, 1);
if sum(ppc_ok) >= 4
    pooled.threshold_low = linear_percentile(ppc(ppc_ok), low_percent);
    pooled.threshold_high = linear_percentile(ppc(ppc_ok), high_percent);
end
if isfinite(pooled.threshold_low) && isfinite(pooled.threshold_high) ...
        && pooled.threshold_high > pooled.threshold_low
    low_idx = ppc_ok & ppc <= pooled.threshold_low;
    high_idx = ppc_ok & ppc >= pooled.threshold_high;
end
all_idx = usable & any(isfinite(acg), 2);

[pooled.low, pooled.n_low] = pool_selected(acg, bin_edges, run_length, low_idx, already_corrected);
[pooled.high, pooled.n_high] = pool_selected(acg, bin_edges, run_length, high_idx, already_corrected);
[pooled.all, pooled.n_all] = pool_selected(acg, bin_edges, run_length, all_idx, already_corrected);
pooled.low_index = low_idx;
pooled.high_index = high_idx;
end


function [mu, n_used] = pool_selected(acg, bin_edges, run_length, idx, already_corrected)
n_bins = size(acg, 2);
mu = nan(1, n_bins);
n_used = sum(idx);
if n_used < 1
    return
end
bin_w = mean(diff(bin_edges));
centers = bin_edges(1:end - 1) + 0.5 * bin_w;
overlap = run_length(idx) - abs(centers);
overlap(overlap < bin_w | ~isfinite(overlap)) = NaN;
counts = acg(idx, :);
if already_corrected
    counts = counts .* overlap;
end
counts(~isfinite(overlap)) = NaN;
den = sum(overlap, 1, 'omitnan');
num = sum(counts, 1, 'omitnan');
n_finite = sum(isfinite(counts) & isfinite(overlap), 1);
ok = den > 0 & n_finite > 0;
mu(ok) = num(ok) ./ den(ok);
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
pooled.low_index = false(0, 1);
pooled.high_index = false(0, 1);
end
