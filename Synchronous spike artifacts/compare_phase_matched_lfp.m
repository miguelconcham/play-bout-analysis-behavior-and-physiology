function comparison = compare_phase_matched_lfp(sync_power, sync_phase, pool_power, pool_phase, n_repeats, n_bins)
%COMPARE_PHASE_MATCHED_LFP  Synchronous LFP power versus same-phase samples.
%
% sync_power and sync_phase are one sample per synchronous event.
% pool_power and pool_phase are the LFP samples those events may be
% compared with (play-bout windows, with the events themselves removed).
% Phase is the analytic angle of the 1-5 Hz delta, in radians, without the
% circular-uniformize warp and without peak centering. Angle 0 is the
% positive peak. Bins are equal steps of 2*pi/n_bins (10 degrees when
% n_bins is 36). +pi and -pi share the first bin.
%
% Each null draw takes one pool sample from the same bin as each event,
% with replacement inside the bin, and averages those draws. percentile is
% the fraction of null means strictly below the synchronous mean. A high
% percentile means the synchronous events sit on larger raw-LFP power than
% other moments at the same delta phase. mean_phase_matched is the average,
% over events, of the mean pool power in that event's bin.

if nargin < 6 || isempty(n_bins)
    n_bins = 36;
end
comparison = empty_comparison(n_bins);
if ~isscalar(n_bins) || ~isfinite(n_bins) || n_bins < 1 || n_bins ~= round(n_bins)
    comparison.reason = 'phase bin count must be a positive integer';
    return
end
comparison.phase_bin_width = 2 * pi / n_bins;
if nargin < 5 || isempty(n_repeats) || ~isscalar(n_repeats) || n_repeats < 1
    comparison.reason = 'n_repeats must be a positive scalar';
    return
end
n_repeats = round(n_repeats);

sync_power = sync_power(:);
sync_phase = sync_phase(:);
pool_power = pool_power(:);
pool_phase = pool_phase(:);
if numel(sync_power) ~= numel(sync_phase) || numel(pool_power) ~= numel(pool_phase)
    comparison.reason = 'power and phase must have the same number of samples';
    return
end

sync_keep = isfinite(sync_power) & isfinite(sync_phase);
sync_power = sync_power(sync_keep);
sync_phase = sync_phase(sync_keep);
pool_keep = isfinite(pool_power) & isfinite(pool_phase);
pool_power = pool_power(pool_keep);
pool_phase = pool_phase(pool_keep);
comparison.n_sync = numel(sync_power);
comparison.n_pool = numel(pool_power);
if comparison.n_sync < 1
    comparison.reason = 'no synchronous sample with finite power and phase';
    return
end

sync_bin = phase_bin(sync_phase, n_bins);
pool_bin = phase_bin(pool_phase, n_bins);
candidates = cell(n_bins, 1);
for bin_index = 1:n_bins
    candidates{bin_index} = find(pool_bin == bin_index);
end
counts = zeros(comparison.n_sync, 1);
for event_index = 1:comparison.n_sync
    counts(event_index) = numel(candidates{sync_bin(event_index)});
end
comparison.min_candidates = min(counts);
if comparison.min_candidates < 1
    comparison.reason = 'a synchronous sample has no LFP in the same delta-phase bin';
    return
end

bin_mean = zeros(comparison.n_sync, 1);
for event_index = 1:comparison.n_sync
    bin_mean(event_index) = mean(pool_power(candidates{sync_bin(event_index)}));
end
comparison.mean_phase_matched = mean(bin_mean);
comparison.mean_sync = mean(sync_power);

null_means = zeros(n_repeats, 1);
for repeat_index = 1:n_repeats
    drawn = zeros(comparison.n_sync, 1);
    for event_index = 1:comparison.n_sync
        pool_index = candidates{sync_bin(event_index)};
        drawn(event_index) = pool_power(pool_index(randi(numel(pool_index))));
    end
    null_means(repeat_index) = mean(drawn);
end
comparison.null_means = null_means;
comparison.percentile = mean(null_means < comparison.mean_sync);
comparison.ok = true;
end

function bin = phase_bin(phase, n_bins)
phase = phase(:);
wrapped = mod(phase + pi, 2 * pi) - pi;
bin = floor((wrapped + pi) / (2 * pi / n_bins)) + 1;
bin(bin < 1) = 1;
bin(bin > n_bins) = n_bins;
bin(~isfinite(phase)) = 0;
end

function comparison = empty_comparison(n_bins)
if nargin < 1 || isempty(n_bins) || ~isfinite(n_bins) || n_bins < 1
    n_bins = 36;
end
comparison.ok = false;
comparison.reason = '';
comparison.n_sync = 0;
comparison.n_pool = 0;
comparison.mean_sync = NaN;
comparison.mean_phase_matched = NaN;
comparison.null_means = zeros(0, 1);
comparison.percentile = NaN;
comparison.min_candidates = 0;
comparison.phase_bin_width = 2 * pi / n_bins;
end
