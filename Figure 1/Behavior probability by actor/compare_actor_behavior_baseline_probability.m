function result = compare_actor_behavior_baseline_probability( ...
    recordings, target_names, settings)
%COMPARE_ACTOR_BEHAVIOR_BASELINE_PROBABILITY Compare post-offset to baseline.
%
% Baseline is the mean probability in pre_average_window relative to target
% onset. Post bins are in real time after target offset. The target duration
% is time-wrapped only for display. Run this function independently for
% relation='self' and relation='other'.

if nargin < 3
    settings = struct();
end
settings = baseline_settings(settings);
target_names = normalize_names(target_names);
if isempty(recordings)
    error('compare_actor_behavior_baseline_probability:Empty', ...
        'No actor-relative recordings were supplied.');
end

behavior_types = recordings(1).behavior_types;
plot_time = recordings(1).wrapped_plot_time;
pre_index = recordings(1).pre_index;
wrap_index = recordings(1).wrap_index;
post_index = recordings(1).post_index;
pre_lag = recordings(1).pre_lag;
n_time = numel(plot_time);
n_behaviors = numel(behavior_types);
[group_names, group_behaviors, group_index] = ...
    behavior_groups(behavior_types);
n_groups = numel(group_names);
n_targets = numel(target_names);
[pre_window, pre_average_index] = resolve_pre_window( ...
    settings.pre_average_window, pre_lag, pre_index);

occupancy = false(0, n_behaviors, n_time);
in_session = false(0, n_time);
target_label = zeros(0, 1);
stratum_index = zeros(0, 1);
stratum_names = cell(0, 1);
stratum_counter = 0;

for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    if ~isequal(recording.behavior_types, behavior_types) || ...
            ~isequal(recording.wrapped_plot_time, plot_time)
        error('compare_actor_behavior_baseline_probability:Parameters', ...
            'All recordings must share behavior and wrapped-time settings.');
    end
    for session_n = 1:numel(recording.sessions)
        session = recording.sessions(session_n);
        for actor_n = 1:numel(session.actor_data)
            stratum_counter = stratum_counter + 1;
            actor = session.actor_data(actor_n);
            stratum_names{stratum_counter, 1} = sprintf('%s P%d %s', ...
                recording.recording_name, session.play_session_number, ...
                actor.actor);
            [selected, labels] = ismember( ...
                actor.target_behavior, target_names);
            rows = find(selected);
            if isempty(rows)
                continue
            end
            if strcmp(settings.relation, 'self')
                values = actor.self_behavior_wrapped(rows, :, :);
            else
                values = actor.other_behavior_wrapped(rows, :, :);
            end
            new_rows = size(occupancy, 1) + (1:numel(rows));
            occupancy(new_rows, :, :) = values;
            in_session(new_rows, :) = actor.wrapped_in_session(rows, :);
            target_label(new_rows, 1) = labels(rows);
            stratum_index(new_rows, 1) = stratum_counter;
        end
    end
end

if any(arrayfun(@(target_n) ~any(target_label == target_n), 1:n_targets))
    error('compare_actor_behavior_baseline_probability:MissingTarget', ...
        'Every requested target behavior must be present.');
end

complete = all(in_session(:, pre_index), 2) & ...
    all(in_session(:, wrap_index), 2) & all(in_session(:, post_index), 2);
group_values = false(size(occupancy, 1), n_groups, n_time);
for group_n = 1:n_groups
    group_values(:, group_n, :) = reshape(any( ...
        occupancy(:, group_index{group_n}, :), 2), ...
        size(occupancy, 1), 1, n_time);
end

n_post = numel(post_index);
probability = nan(n_groups, n_targets, n_time);
behavior_probability = nan(n_behaviors, n_targets, n_time);
pre_probability = nan(n_groups, n_targets);
post_probability = nan(n_groups, n_targets);
observed_difference = nan(n_groups, n_targets);
observed_by_bin = nan(n_groups, n_post, n_targets);
p_value = nan(n_groups, n_targets);
p_value_by_bin = nan(n_groups, n_post, n_targets);
p_value_by_bin_fwer_group = nan(n_groups, n_post, n_targets);
p_value_by_bin_fwer_target = nan(n_groups, n_post, n_targets);
post_probability_ci_lower = nan(n_groups, n_post, n_targets);
post_probability_ci_upper = nan(n_groups, n_post, n_targets);
null_difference = cell(1, n_targets);
n_events_available = zeros(1, n_targets);
n_events_test = zeros(1, n_targets);
n_strata = zeros(1, n_targets);
stratum_names_by_target = cell(1, n_targets);

previous_rng = rng;
rng_cleanup = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(settings.random_seed, 'twister');
for target_n = 1:n_targets
    target_rows = target_label == target_n;
    n_events_available(target_n) = sum(target_rows);
    included = target_rows & complete;
    n_events_test(target_n) = sum(included);
    if n_events_test(target_n) < 1
        error('compare_actor_behavior_baseline_probability:NoCompleteEvents', ...
            '%s has no complete pre/post events.', target_names{target_n});
    end
    included_strata = unique(stratum_index(included));
    n_strata(target_n) = numel(included_strata);
    stratum_names_by_target{target_n} = ...
        stratum_names(included_strata);
    included_stratum_index = stratum_index(included);

    target_groups = double(group_values(included, :, :));
    target_behaviors = double(occupancy(included, :, :));
    probability(:, target_n, :) = reshape(mean(target_groups, 1), ...
        n_groups, 1, n_time);
    behavior_probability(:, target_n, :) = reshape( ...
        mean(target_behaviors, 1), n_behaviors, 1, n_time);
    pre_event = mean(target_groups(:, :, pre_average_index), 3);
    post_event = mean(target_groups(:, :, post_index), 3);
    pre_probability(:, target_n) = mean(pre_event, 1).';
    post_probability(:, target_n) = mean(post_event, 1).';
    score_diff = post_event - pre_event;
    observed_difference(:, target_n) = mean(score_diff, 1).';
    postbin_diff = bsxfun(@minus, ...
        target_groups(:, :, post_index), ...
        reshape(pre_event, size(pre_event, 1), n_groups, 1));
    observed_by_bin(:, :, target_n) = reshape( ...
        mean(postbin_diff, 1), n_groups, n_post);

    null_difference{target_n} = nan(settings.n_permutations, n_groups);
    exceedances = zeros(n_groups, n_post);
    null_max_group = nan(settings.n_permutations, n_groups);
    null_max_target = nan(settings.n_permutations, 1);
    postbin_matrix = reshape(postbin_diff, n_events_test(target_n), []);
    for permutation_n = 1:settings.n_permutations
        signs = ones(n_events_test(target_n), 1);
        for stratum_n = 1:numel(included_strata)
            rows = included_stratum_index == included_strata(stratum_n);
            signs(rows) = 1 - 2 * (rand < 0.5);
        end
        null_difference{target_n}(permutation_n, :) = ...
            (signs.' * score_diff) / n_events_test(target_n);
        shuffled = reshape((signs.' * postbin_matrix) / ...
            n_events_test(target_n), n_groups, n_post);
        exceedances = exceedances + (abs(shuffled) >= ...
            abs(observed_by_bin(:, :, target_n)));
        null_max_group(permutation_n, :) = max(abs(shuffled), [], 2).';
        null_max_target(permutation_n) = max(abs(shuffled), [], 'all');
    end
    p_value(:, target_n) = (1 + sum(abs(null_difference{target_n}) >= ...
        abs(observed_difference(:, target_n).'), 1)).' / ...
        (settings.n_permutations + 1);
    p_value_by_bin(:, :, target_n) = ...
        (1 + exceedances) / (settings.n_permutations + 1);
    for group_n = 1:n_groups
        p_value_by_bin_fwer_group(group_n, :, target_n) = ...
            (1 + sum(bsxfun(@ge, null_max_group(:, group_n), ...
            abs(observed_by_bin(group_n, :, target_n))), 1)) / ...
            (settings.n_permutations + 1);
    end
    p_value_by_bin_fwer_target(:, :, target_n) = reshape( ...
        (1 + sum(bsxfun(@ge, null_max_target, ...
        abs(reshape(observed_by_bin(:, :, target_n), 1, []))), 1)) / ...
        (settings.n_permutations + 1), n_groups, n_post);

    % Resample complete actor-session clusters while retaining every event
    % in a sampled cluster. The resulting pointwise interval describes
    % uncertainty in the observed post-offset mean, not the permutation null.
    rng_after_permutations = rng;
    bootstrap_seed = mod(settings.bootstrap_seed + target_n - 1, 2^32);
    rng(bootstrap_seed, 'twister');
    [~, event_stratum] = ismember( ...
        included_stratum_index, included_strata);
    post_matrix = reshape(target_groups(:, :, post_index), ...
        n_events_test(target_n), n_groups * n_post);
    bootstrap_probability = nan(settings.n_bootstraps, ...
        n_groups * n_post);
    for bootstrap_n = 1:settings.n_bootstraps
        sampled_strata = randi(numel(included_strata), ...
            numel(included_strata), 1);
        stratum_counts = accumarray(sampled_strata, 1, ...
            [numel(included_strata), 1]);
        event_weights = stratum_counts(event_stratum);
        bootstrap_probability(bootstrap_n, :) = ...
            (event_weights.' * post_matrix) / sum(event_weights);
    end
    tail_probability = (1 - settings.confidence_level) / 2;
    post_probability_ci_lower(:, :, target_n) = reshape( ...
        empirical_quantile(bootstrap_probability, tail_probability), ...
        n_groups, n_post);
    post_probability_ci_upper(:, :, target_n) = reshape( ...
        empirical_quantile(bootstrap_probability, 1 - tail_probability), ...
        n_groups, n_post);
    rng(rng_after_permutations);
end

result = struct();
result.target_names = target_names;
result.relation = settings.relation;
result.group_names = group_names;
result.group_behaviors = group_behaviors;
result.behavior_types = behavior_types;
result.wrapped_plot_time = plot_time;
result.pre_index = pre_index;
result.wrap_index = wrap_index;
result.post_index = post_index;
result.post_time = recordings(1).post_lag;
result.pre_average_window = pre_window;
result.pre_average_index = pre_average_index;
result.probability = probability;
result.behavior_probability = behavior_probability;
result.pre_probability = pre_probability;
result.post_probability = post_probability;
result.observed_difference = observed_difference;
result.p_value = p_value;
result.p_value_holm = reshape(holm_adjust(p_value(:)), size(p_value));
result.observed_difference_by_bin = observed_by_bin;
result.p_value_by_bin = p_value_by_bin;
result.p_value_by_bin_fwer_group = p_value_by_bin_fwer_group;
result.p_value_by_bin_fwer_target = p_value_by_bin_fwer_target;
result.post_probability_ci_lower = post_probability_ci_lower;
result.post_probability_ci_upper = post_probability_ci_upper;
result.confidence_level = settings.confidence_level;
result.n_bootstraps = settings.n_bootstraps;
result.bootstrap_seed = settings.bootstrap_seed;
result.null_difference = null_difference;
result.n_events_available = n_events_available;
result.n_events_test = n_events_test;
result.n_actor_session_strata = n_strata;
result.stratum_names = stratum_names_by_target;
result.n_permutations = settings.n_permutations;
result.random_seed = settings.random_seed;
result.bin_significance_alpha = 0.05;
result.test_description = [ ...
    'Two-sided actor-session-cluster sign-flip test of post-offset minus baseline'];
result.confidence_interval_description = [ ...
    'Pointwise percentile interval from actor-session cluster bootstrap'];
end

function settings = baseline_settings(overrides)
settings = struct('relation', 'self', 'pre_average_window', [], ...
    'n_permutations', 10000, 'random_seed', 1, ...
    'n_bootstraps', 1000, 'bootstrap_seed', 2, ...
    'confidence_level', 0.95);
names = fieldnames(overrides);
for field_n = 1:numel(names)
    if ~isfield(settings, names{field_n})
        error('compare_actor_behavior_baseline_probability:Settings', ...
            'Unknown setting "%s".', names{field_n});
    end
    settings.(names{field_n}) = overrides.(names{field_n});
end
settings.relation = validatestring(settings.relation, {'self', 'other'});
if ~isempty(settings.pre_average_window) && ...
        (~isnumeric(settings.pre_average_window) || ...
        numel(settings.pre_average_window) ~= 2 || ...
        any(~isfinite(settings.pre_average_window)) || ...
        settings.pre_average_window(2) <= settings.pre_average_window(1) || ...
        settings.pre_average_window(2) > 0)
    error('compare_actor_behavior_baseline_probability:PreWindow', ...
        'pre_average_window must be empty or an ascending nonpositive pair.');
end
if ~isnumeric(settings.n_permutations) || ...
        ~isscalar(settings.n_permutations) || ...
        ~isfinite(settings.n_permutations) || settings.n_permutations < 1 || ...
        settings.n_permutations ~= round(settings.n_permutations)
    error('compare_actor_behavior_baseline_probability:Permutations', ...
        'n_permutations must be a positive integer.');
end
if ~isnumeric(settings.random_seed) || ...
        ~isscalar(settings.random_seed) || ...
        ~isfinite(settings.random_seed) || settings.random_seed < 0 || ...
        settings.random_seed > 2^32 - 1 || ...
        settings.random_seed ~= round(settings.random_seed)
    error('compare_actor_behavior_baseline_probability:Seed', ...
        'random_seed must be an integer from 0 through 2^32 - 1.');
end
if ~isnumeric(settings.n_bootstraps) || ...
        ~isscalar(settings.n_bootstraps) || ...
        ~isfinite(settings.n_bootstraps) || settings.n_bootstraps < 1 || ...
        settings.n_bootstraps ~= round(settings.n_bootstraps)
    error('compare_actor_behavior_baseline_probability:Bootstraps', ...
        'n_bootstraps must be a positive integer.');
end
if ~isnumeric(settings.bootstrap_seed) || ...
        ~isscalar(settings.bootstrap_seed) || ...
        ~isfinite(settings.bootstrap_seed) || settings.bootstrap_seed < 0 || ...
        settings.bootstrap_seed > 2^32 - 1 || ...
        settings.bootstrap_seed ~= round(settings.bootstrap_seed)
    error('compare_actor_behavior_baseline_probability:BootstrapSeed', ...
        'bootstrap_seed must be an integer from 0 through 2^32 - 1.');
end
if ~isnumeric(settings.confidence_level) || ...
        ~isscalar(settings.confidence_level) || ...
        ~isfinite(settings.confidence_level) || ...
        settings.confidence_level <= 0 || settings.confidence_level >= 1
    error('compare_actor_behavior_baseline_probability:ConfidenceLevel', ...
        'confidence_level must be strictly between zero and one.');
end
end

function output = empirical_quantile(values, probability)
values = sort(values, 1);
n = size(values, 1);
position = 1 + (n - 1) * probability;
lower_index = floor(position);
upper_index = ceil(position);
fraction = position - lower_index;
output = (1 - fraction) * values(lower_index, :) + ...
    fraction * values(upper_index, :);
end

function [window, index] = resolve_pre_window(requested, pre_lag, pre_index)
if isempty(requested)
    window = [pre_lag(1), 0];
    index = pre_index;
    return
end
window = double(requested(:).');
tolerance = eps(max(1, max(abs(pre_lag)))) * 16;
if numel(pre_lag) > 1
    tolerance = max(tolerance, min(abs(diff(pre_lag))) * 1e-9);
end
local = pre_lag >= window(1) - tolerance & ...
    pre_lag <= window(2) + tolerance;
index = pre_index(local);
if isempty(index)
    error('compare_actor_behavior_baseline_probability:PreWindowBins', ...
        'pre_average_window contains no sampled pre-onset bins.');
end
end

function names = normalize_names(values)
if ischar(values)
    names = {values};
elseif isstring(values)
    names = cellstr(values(:).');
elseif iscell(values)
    names = cellfun(@char, values, 'UniformOutput', false);
    names = names(:).';
else
    error('compare_actor_behavior_baseline_probability:Targets', ...
        'target_names must be text.');
end
end

function [names, behaviors, index] = behavior_groups(behavior_types)
names = {'Play', 'Aggressive', 'Exploratory'};
behaviors = { ...
    {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'}, ...
    {'PounceI', 'Bite'}, ...
    {'Grooming', 'Rearing', 'Sniffing', 'Scratching'}};
index = cell(size(behaviors));
for group_n = 1:numel(behaviors)
    missing = setdiff(behaviors{group_n}, behavior_types);
    if ~isempty(missing)
        error('compare_actor_behavior_baseline_probability:BehaviorTypes', ...
            'Recording behavior list is missing: %s.', strjoin(missing, ', '));
    end
    index{group_n} = find(ismember(behavior_types, behaviors{group_n}));
end
end

function adjusted = holm_adjust(p_values)
[sorted, order] = sort(p_values);
adjusted_sorted = nan(size(sorted));
running_max = 0;
for test_n = 1:numel(sorted)
    running_max = max(running_max, ...
        (numel(sorted) - test_n + 1) * sorted(test_n));
    adjusted_sorted(test_n) = min(1, running_max);
end
adjusted = nan(size(p_values));
adjusted(order) = adjusted_sorted;
end
