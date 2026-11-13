function summary = pool_group_ccg_sessions(sessions, central_lag)
%POOL_GROUP_CCG_SESSIONS  Equal-weight mean of the per-session correlograms.
%
% Each saved session is one recording. A session that did not have both
% groups, or that did not reach the reference-spike minimum, stays out of
% that pair. The central value of a session is the mean excess, and the
% mean rate ratio, on lags whose absolute value is at most central_lag.
% The t-test compares those session means with zero. Holm-Sidak is applied
% across the directed pairs that have a finite p-value.

if nargin < 2 || isempty(central_lag)
    central_lag = 0.05;
end
if ~(isscalar(central_lag) && central_lag > 0)
    error('group_ccg:pool', 'central_lag must be a positive lag in seconds.');
end
names = group_ccg_names();
if isempty(sessions)
    error('group_ccg:pool', 'There are no session correlograms to pool.');
end
centers = double(sessions(1).centers(:)');
n_bins = numel(centers);
n_pairs = numel(names.pair_reference);
n_sessions = numel(sessions);
excess = nan(n_sessions, n_bins, n_pairs);
ratio = nan(n_sessions, n_bins, n_pairs);
session_names = cell(n_sessions, 1);
n_neurons_reference = nan(n_sessions, n_pairs);
n_neurons_target = nan(n_sessions, n_pairs);

for s = 1:n_sessions
    session_names{s} = ccg_plain_text(sessions(s).animal_code);
    if numel(sessions(s).centers) ~= n_bins ...
            || max(abs(double(sessions(s).centers(:))' - centers)) > 1e-9
        error('group_ccg:pool', ...
            'Session %s is not on the same lag axis as %s.', ...
            session_names{s}, session_names{1});
    end
    if numel(sessions(s).pairs) ~= n_pairs
        error('group_ccg:pool', ...
            'Session %s does not have %d directed pairs.', session_names{s}, n_pairs);
    end
    for k = 1:n_pairs
        pair = sessions(s).pairs(k);
        if ~strcmp(pair.reference, names.group_id{names.pair_reference(k)}) ...
                || ~strcmp(pair.target, names.group_id{names.pair_target(k)})
            error('group_ccg:pool', ...
                'Session %s pair %d is not %s against %s.', session_names{s}, k, ...
                names.group_id{names.pair_reference(k)}, names.group_id{names.pair_target(k)});
        end
        if ~isfield(pair, 'usable') || ~pair.usable || ~isfield(pair, 'excess_rate')
            continue
        end
        excess(s, :, k) = double(pair.excess_rate(:))';
        ratio(s, :, k) = double(pair.rate_ratio(:))';
        n_neurons_reference(s, k) = pair.n_neurons_reference;
        n_neurons_target(s, k) = pair.n_neurons_target;
    end
end

in_center = abs(centers) <= central_lag + 1e-12;
central_excess = squeeze_center(excess, in_center);
central_ratio = squeeze_center(ratio, in_center);
[mean_excess, ci_low, ci_high, n_bins_used] = mean_ci(excess);
[mean_ratio, ratio_low, ratio_high, ~] = mean_ci(ratio);
raw_p = nan(n_pairs, 1);
t_stat = nan(n_pairs, 1);
n_tested = zeros(n_pairs, 1);
central_mean = nan(n_pairs, 1);
central_low = nan(n_pairs, 1);
central_high = nan(n_pairs, 1);
for k = 1:n_pairs
    values = central_excess(:, k);
    values = values(isfinite(values));
    n_tested(k) = numel(values);
    if numel(values) < 2
        continue
    end
    central_mean(k) = mean(values);
    [~, raw_p(k), ci, stats] = ttest(values, 0);
    t_stat(k) = stats.tstat;
    central_low(k) = ci(1);
    central_high(k) = ci(2);
end

summary.centers = centers;
summary.group_id = names.group_id;
summary.title = names.title;
summary.color = names.color;
summary.pair_reference = names.pair_reference;
summary.pair_target = names.pair_target;
summary.session_names = session_names;
summary.central_lag = central_lag;
summary.excess = excess;
summary.rate_ratio = ratio;
summary.mean_excess = mean_excess;
summary.ci_low = ci_low;
summary.ci_high = ci_high;
summary.mean_ratio = mean_ratio;
summary.ratio_low = ratio_low;
summary.ratio_high = ratio_high;
summary.n_bins_used = n_bins_used;
summary.central_excess = central_excess;
summary.central_ratio = central_ratio;
summary.central_mean = central_mean;
summary.central_low = central_low;
summary.central_high = central_high;
summary.t_stat = t_stat;
summary.raw_p = raw_p;
summary.holm_p = holm_sidak_adjust(raw_p);
summary.n_sessions = n_tested;
summary.n_neurons_reference = mean(n_neurons_reference, 1, 'omitnan')';
summary.n_neurons_target = mean(n_neurons_target, 1, 'omitnan')';
end


function values = squeeze_center(curves, in_center)
n_sessions = size(curves, 1);
n_pairs = size(curves, 3);
values = nan(n_sessions, n_pairs);
for s = 1:n_sessions
    for k = 1:n_pairs
        row = curves(s, in_center, k);
        if any(isfinite(row))
            values(s, k) = mean(row, 'omitnan');
        end
    end
end
end


function [mu, lo, hi, n] = mean_ci(curves)
n_bins = size(curves, 2);
n_pairs = size(curves, 3);
mu = nan(n_pairs, n_bins);
lo = nan(n_pairs, n_bins);
hi = nan(n_pairs, n_bins);
n = zeros(n_pairs, n_bins);
for k = 1:n_pairs
    sample = squeeze(curves(:, :, k));
    if size(curves, 1) == 1
        sample = sample(:)';
    end
    for b = 1:n_bins
        values = sample(:, b);
        values = values(isfinite(values));
        n(k, b) = numel(values);
        if isempty(values)
            continue
        end
        mu(k, b) = mean(values);
        if numel(values) < 2
            continue
        end
        half = tinv(0.975, numel(values) - 1) * std(values) / sqrt(numel(values));
        lo(k, b) = mu(k, b) - half;
        hi(k, b) = mu(k, b) + half;
    end
end
end
