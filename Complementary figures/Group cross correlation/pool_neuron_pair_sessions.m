function summary = pool_neuron_pair_sessions(sessions, central_lag)
%POOL_NEURON_PAIR_SESSIONS  Normalize every pair, then mean and CI by session.
%
% Each neuron pair is divided by its own target rate. The session curve is
% the equal-weight mean of that session's usable pairs. The mean and the
% 95% t interval are then taken across sessions. pair_excess{k} stacks the
% normalized pairs for the images. session_excess is one row per recording.

if nargin < 2 || isempty(central_lag)
    central_lag = 0.05;
end
if ~(isscalar(central_lag) && central_lag > 0)
    error('group_ccg:pool', 'central_lag must be a positive lag in seconds.');
end
names = group_ccg_names();
if isempty(sessions)
    error('group_ccg:pool', 'There are no neuron-pair sessions to pool.');
end
centers = double(sessions(1).centers(:)');
n_bins = numel(centers);
bin_w = median(diff(centers));
n_combinations = numel(names.pair_reference);
n_sessions = numel(sessions);
session_excess = nan(n_sessions, n_bins, n_combinations);
pair_excess = cell(n_combinations, 1);
session_names = cell(n_sessions, 1);
n_pairs = zeros(n_combinations, 1);
for k = 1:n_combinations
    pair_excess{k} = zeros(0, n_bins);
end

for s = 1:n_sessions
    session_names{s} = ccg_plain_text(sessions(s).animal_code);
    if numel(sessions(s).centers) ~= n_bins ...
            || max(abs(double(sessions(s).centers(:))' - centers)) > 1e-9
        error('group_ccg:pool', ...
            'Session %s is not on the same lag axis as %s.', ...
            session_names{s}, session_names{1});
    end
    if numel(sessions(s).combinations) ~= n_combinations
        error('group_ccg:pool', ...
            'Session %s does not have %d group combinations.', session_names{s}, n_combinations);
    end
    for k = 1:n_combinations
        combination = sessions(s).combinations(k);
        if ~strcmp(combination.reference, names.group_id{names.pair_reference(k)}) ...
                || ~strcmp(combination.target, names.group_id{names.pair_target(k)})
            error('group_ccg:pool', ...
                'Session %s combination %d is not %s against %s.', session_names{s}, k, ...
                names.group_id{names.pair_reference(k)}, names.group_id{names.pair_target(k)});
        end
        normalized = normalized_rows(combination, bin_w);
        if isempty(normalized)
            continue
        end
        pair_excess{k} = [pair_excess{k}; normalized]; %#ok<AGROW>
        n_pairs(k) = n_pairs(k) + size(normalized, 1);
        session_excess(s, :, k) = mean(normalized, 1, 'omitnan');
    end
end

in_center = abs(centers) <= central_lag + 1e-12;
central = nan(n_sessions, n_combinations);
for s = 1:n_sessions
    for k = 1:n_combinations
        row = session_excess(s, in_center, k);
        if any(isfinite(row))
            central(s, k) = mean(row, 'omitnan');
        end
    end
end
[mean_excess, ci_low, ci_high, n_bins_used] = mean_ci(session_excess);
raw_p = nan(n_combinations, 1);
t_stat = nan(n_combinations, 1);
n_tested = zeros(n_combinations, 1);
central_mean = nan(n_combinations, 1);
central_low = nan(n_combinations, 1);
central_high = nan(n_combinations, 1);
for k = 1:n_combinations
    values = central(:, k);
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
summary.pair_reference = names.pair_reference;
summary.pair_target = names.pair_target;
summary.session_names = session_names;
summary.central_lag = central_lag;
summary.session_excess = session_excess;
summary.pair_excess = pair_excess;
summary.mean_excess = mean_excess;
summary.ci_low = ci_low;
summary.ci_high = ci_high;
summary.n_bins_used = n_bins_used;
summary.central_mean = central_mean;
summary.central_low = central_low;
summary.central_high = central_high;
summary.t_stat = t_stat;
summary.raw_p = raw_p;
summary.holm_p = holm_sidak_adjust(raw_p);
summary.n_sessions = n_tested;
summary.n_pairs = n_pairs;
end


function normalized = normalized_rows(combination, bin_w)
normalized = zeros(0, size(combination.counts, 2));
if ~isfield(combination, 'usable') || ~any(combination.usable)
    return
end
rows = find(combination.usable);
kept = zeros(numel(rows), size(combination.counts, 2));
n_kept = 0;
for i = 1:numel(rows)
    row = rows(i);
    if ~(combination.n_reference(row) > 0)
        continue
    end
    conditional = double(combination.counts(row, :)) / (combination.n_reference(row) * bin_w);
    [~, this_normalized] = normalize_pair_curve(conditional, combination.target_rate(row));
    if ~any(isfinite(this_normalized))
        continue
    end
    n_kept = n_kept + 1;
    kept(n_kept, :) = this_normalized;
end
normalized = kept(1:n_kept, :);
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
