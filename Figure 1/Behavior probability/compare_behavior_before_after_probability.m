function result = compare_behavior_before_after_probability( ...
    recordings, target_names, settings)
%COMPARE_BEHAVIOR_BEFORE_AFTER_PROBABILITY Test pre-onset vs post-offset.
%
% Each target event is placed on a shared axis: real time from -T up to
% but not including onset, a fixed number of samples from onset through
% offset, then real time from just after offset through offset+T.
%
% For each requested target behavior and each behavior group, the
% statistic is the mean across events of (mean post-flank occupancy minus
% mean pre-flank occupancy). Under the null, the pre-onset and post-offset
% flanks of each event are exchangeable, so event-level signs of that
% paired difference are flipped. Flips are drawn independently inside each
% play session. Events whose flanks are not entirely inside the Partners
% session are excluded. The six target-by-group window p-values are
% Holm-corrected. Lag-matched bins use the same sign flips, with max-T
% correction across lags and groups within each target.
%
% A second, one-sided bin-wise test uses those same signs. Each post-offset
% bin is compared with the mean occupancy over the whole pre-onset flank,
% so "before" is the rate at any time in (-T, onset) rather than one
% matched lag. The alternative is that the post bin is higher than that
% pre-onset rate. Max-T correction is one-sided and within each target.

if nargin < 3
    settings = struct();
end
settings = before_after_settings(settings);
target_names = normalize_targets(target_names);
if isempty(recordings)
    error('compare_behavior_before_after_probability:Empty', ...
        'No behavior-probability recordings were supplied.');
end
if ~isfield(recordings, 'wrapped_plot_time') || ...
        ~isfield(recordings(1).sessions, 'self_behavior_wrapped')
    error('compare_behavior_before_after_probability:MissingWrap', ...
        ['Recordings have no time-wrapped matrices. Re-run ' ...
        'Estimate_behavior_probability.m before this test.']);
end

behavior_types = recordings(1).behavior_types;
plot_time = recordings(1).wrapped_plot_time;
pre_index = recordings(1).pre_index;
wrap_index = recordings(1).wrap_index;
post_index = recordings(1).post_index;
pre_lag = recordings(1).pre_lag;
post_lag = recordings(1).post_lag;
n_time = numel(plot_time);
n_behaviors = numel(behavior_types);
group_names = {'Play', 'Aggressive', 'Exploratory'};
group_behaviors = { ...
    {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'}, ...
    {'PounceI', 'Bite'}, ...
    {'Grooming', 'Rearing', 'Sniffing', 'Scratching'}};
n_groups = numel(group_names);
n_targets = numel(target_names);

for group_n = 1:n_groups
    missing = setdiff(group_behaviors{group_n}, behavior_types);
    if ~isempty(missing)
        error('compare_behavior_before_after_probability:BehaviorTypes', ...
            'Recording behavior list is missing: %s.', strjoin(missing, ', '));
    end
end
group_index = cell(1, n_groups);
for group_n = 1:n_groups
    group_index{group_n} = find(ismember(behavior_types, ...
        group_behaviors{group_n}));
end
if numel(pre_index) ~= numel(post_index) || isempty(pre_index)
    error('compare_behavior_before_after_probability:Flanks', ...
        'The pre-onset and post-offset flanks must contain the same bins.');
end

self_occupancy = false(0, n_behaviors, n_time);
partner_occupancy = false(0, n_behaviors, n_time);
in_session = false(0, n_time);
target_label = zeros(0, 1);
session_index = zeros(0, 1);
session_names = cell(0, 1);
session_counter = 0;

for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    if ~isequal(recording.behavior_types, behavior_types) || ...
            ~isequal(recording.wrapped_plot_time, plot_time) || ...
            ~isequal(recording.pre_index, pre_index) || ...
            ~isequal(recording.post_index, post_index)
        error('compare_behavior_before_after_probability:Parameters', ...
            ['All recordings must share the behavior list and the ' ...
            'time-wrapped axis.']);
    end
    for session_n = 1:numel(recording.sessions)
        session_counter = session_counter + 1;
        session = recording.sessions(session_n);
        [selected, labels] = ismember(session.target_behavior, target_names);
        rows = find(selected);
        session_names{session_counter, 1} = sprintf('%s P%d', ...
            recording.recording_name, session.play_session_number);
        if isempty(rows)
            continue
        end
        labels = labels(rows);
        new_rows = size(self_occupancy, 1) + (1:numel(rows));
        self_occupancy(new_rows, :, :) = ...
            session.self_behavior_wrapped(rows, :, :);
        partner_occupancy(new_rows, :, :) = ...
            session.partner_behavior_wrapped(rows, :, :);
        in_session(new_rows, :) = session.wrapped_in_session(rows, :);
        target_label(new_rows, 1) = labels;
        session_index(new_rows, 1) = session_counter;
    end
end

if any(arrayfun(@(x) ~any(target_label == x), 1:n_targets))
    error('compare_behavior_before_after_probability:MissingTarget', ...
        ['Every requested target must be present. Re-run Estimate with ' ...
        'target_behaviors containing %s.'], strjoin(target_names, ', '));
end

switch settings.actor_relation
    case 'self'
        event_matrix = self_occupancy;
    case 'partner'
        event_matrix = partner_occupancy;
    case 'either'
        event_matrix = self_occupancy | partner_occupancy;
end

complete = all(in_session(:, pre_index), 2) & ...
    all(in_session(:, wrap_index), 2) & ...
    all(in_session(:, post_index), 2);
group_values = group_occupancy(event_matrix, group_index, in_session);
self_values = masked_occupancy(self_occupancy, in_session);
partner_values = masked_occupancy(partner_occupancy, in_session);
actor_values = masked_occupancy(event_matrix, in_session);

n_lags = numel(post_index);
probability = nan(n_groups, n_targets, n_time);
behavior_probability = nan(n_behaviors, n_targets, n_time);
behavior_probability_self = nan(n_behaviors, n_targets, n_time);
behavior_probability_partner = nan(n_behaviors, n_targets, n_time);
pre_probability = nan(n_groups, n_targets);
post_probability = nan(n_groups, n_targets);
behavior_pre_probability = nan(n_behaviors, n_targets);
behavior_post_probability = nan(n_behaviors, n_targets);
observed_difference = nan(n_groups, n_targets);
observed_difference_by_lag = nan(n_groups, n_lags, n_targets);
observed_postbin_minus_pre = nan(n_groups, n_lags, n_targets);
p_value = nan(n_groups, n_targets);
p_value_postbin_vs_pre = nan(n_groups, n_lags, n_targets);
p_value_postbin_vs_pre_fwer_group = nan(n_groups, n_lags, n_targets);
p_value_postbin_vs_pre_fwer_target = nan(n_groups, n_lags, n_targets);
p_value_by_lag = nan(n_groups, n_lags, n_targets);
p_value_by_lag_fwer_group = nan(n_groups, n_lags, n_targets);
p_value_by_lag_fwer_target = nan(n_groups, n_lags, n_targets);
null_difference = cell(1, n_targets);
n_events_test = zeros(1, n_targets);
n_events_available = zeros(1, n_targets);
n_sessions = zeros(1, n_targets);
session_names_by_target = cell(1, n_targets);

previous_rng = rng;
rng_cleanup = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(settings.random_seed, 'twister');

for target_n = 1:n_targets
    target_rows = target_label == target_n;
    n_events_available(target_n) = sum(target_rows);
    included = target_rows & complete;
    n_events_test(target_n) = sum(included);
    if n_events_test(target_n) < 1
        error('compare_behavior_before_after_probability:NoCompleteEvents', ...
            ['%s has no events whose pre, wrapped, and post samples all ' ...
            'fall inside the Partners session.'], target_names{target_n});
    end
    included_sessions = unique(session_index(included));
    n_sessions(target_n) = numel(included_sessions);
    session_names_by_target{target_n} = session_names(included_sessions);

    target_groups = group_values(included, :, :);
    probability(:, target_n, :) = reshape(mean(target_groups, 1), ...
        n_groups, 1, n_time);
    pre_probability(:, target_n) = reshape(mean(mean( ...
        target_groups(:, :, pre_index), 3), 1), n_groups, 1);
    post_probability(:, target_n) = reshape(mean(mean( ...
        target_groups(:, :, post_index), 3), 1), n_groups, 1);
    behavior_probability(:, target_n, :) = reshape(mean( ...
        actor_values(included, :, :), 1), n_behaviors, 1, n_time);
    behavior_probability_self(:, target_n, :) = reshape(mean( ...
        self_values(included, :, :), 1), n_behaviors, 1, n_time);
    behavior_probability_partner(:, target_n, :) = reshape(mean( ...
        partner_values(included, :, :), 1), n_behaviors, 1, n_time);
    behavior_pre_probability(:, target_n) = reshape(mean(mean( ...
        actor_values(included, :, pre_index), 3), 1), n_behaviors, 1);
    behavior_post_probability(:, target_n) = reshape(mean(mean( ...
        actor_values(included, :, post_index), 3), 1), n_behaviors, 1);

    score_diff = reshape(mean(target_groups(:, :, post_index), 3) - ...
        mean(target_groups(:, :, pre_index), 3), [], n_groups);
    % Pair the sample T seconds before onset with the sample T seconds
    % after offset, and the last pre-onset sample with the first post sample.
    lag_diff = reshape(target_groups(:, :, post_index) - ...
        target_groups(:, :, pre_index(end:-1:1)), ...
        n_events_test(target_n), n_groups, n_lags);
    observed_difference(:, target_n) = mean(score_diff, 1).';
    observed_difference_by_lag(:, :, target_n) = reshape( ...
        mean(lag_diff, 1), n_groups, n_lags);
    pre_mean_event = mean(target_groups(:, :, pre_index), 3);
    postbin_minus_pre = bsxfun(@minus, target_groups(:, :, post_index), ...
        reshape(pre_mean_event, size(pre_mean_event, 1), n_groups, 1));
    observed_postbin_minus_pre(:, :, target_n) = reshape( ...
        mean(postbin_minus_pre, 1), n_groups, n_lags);

    included_session_index = session_index(included);
    null_difference{target_n} = nan(settings.n_permutations, n_groups);
    lag_exceedances = zeros(n_groups, n_lags);
    null_max_by_group = nan(settings.n_permutations, n_groups);
    null_max_target = nan(settings.n_permutations, 1);
    postbin_exceedances = zeros(n_groups, n_lags);
    null_max_postbin_group = nan(settings.n_permutations, n_groups);
    null_max_postbin_target = nan(settings.n_permutations, 1);
    lag_matrix = reshape(lag_diff, n_events_test(target_n), []);
    observed_lag = reshape(observed_difference_by_lag(:, :, target_n), 1, []);
    postbin_matrix = reshape(postbin_minus_pre, n_events_test(target_n), []);
    observed_postbin = reshape( ...
        observed_postbin_minus_pre(:, :, target_n), 1, []);
    for permutation_n = 1:settings.n_permutations
        signs = ones(n_events_test(target_n), 1);
        for session_n = 1:numel(included_sessions)
            rows = find(included_session_index == included_sessions(session_n));
            signs(rows) = 1 - 2 * (rand(numel(rows), 1) < 0.5);
        end
        null_difference{target_n}(permutation_n, :) = ...
            (signs.' * score_diff) / n_events_test(target_n);
        shuffled_lag = reshape((signs.' * lag_matrix) / ...
            n_events_test(target_n), n_groups, n_lags);
        lag_exceedances = lag_exceedances + ...
            (abs(shuffled_lag) >= abs(observed_difference_by_lag(:, :, target_n)));
        null_max_by_group(permutation_n, :) = max(abs(shuffled_lag), [], 2).';
        null_max_target(permutation_n) = max(abs(shuffled_lag), [], 'all');
        shuffled_postbin = reshape((signs.' * postbin_matrix) / ...
            n_events_test(target_n), n_groups, n_lags);
        postbin_exceedances = postbin_exceedances + (shuffled_postbin >= ...
            observed_postbin_minus_pre(:, :, target_n));
        null_max_postbin_group(permutation_n, :) = ...
            max(shuffled_postbin, [], 2).';
        null_max_postbin_target(permutation_n) = ...
            max(shuffled_postbin, [], 'all');
    end
    p_value(:, target_n) = (1 + sum(abs(null_difference{target_n}) >= ...
        abs(observed_difference(:, target_n).'), 1)).' / ...
        (settings.n_permutations + 1);
    p_value_by_lag(:, :, target_n) = (1 + lag_exceedances) / ...
        (settings.n_permutations + 1);
    for group_n = 1:n_groups
        p_value_by_lag_fwer_group(group_n, :, target_n) = ...
            (1 + sum(bsxfun(@ge, null_max_by_group(:, group_n), ...
            abs(observed_difference_by_lag(group_n, :, target_n))), 1)) / ...
            (settings.n_permutations + 1);
    end
    target_fwer = (1 + sum(bsxfun(@ge, null_max_target, ...
        abs(observed_lag)), 1)) / (settings.n_permutations + 1);
    p_value_by_lag_fwer_target(:, :, target_n) = reshape( ...
        target_fwer, n_groups, n_lags);
    p_value_postbin_vs_pre(:, :, target_n) = ...
        (1 + postbin_exceedances) / (settings.n_permutations + 1);
    for group_n = 1:n_groups
        p_value_postbin_vs_pre_fwer_group(group_n, :, target_n) = ...
            (1 + sum(bsxfun(@ge, null_max_postbin_group(:, group_n), ...
            observed_postbin_minus_pre(group_n, :, target_n)), 1)) / ...
            (settings.n_permutations + 1);
    end
    postbin_fwer = (1 + sum(bsxfun(@ge, null_max_postbin_target, ...
        observed_postbin), 1)) / (settings.n_permutations + 1);
    p_value_postbin_vs_pre_fwer_target(:, :, target_n) = reshape( ...
        postbin_fwer, n_groups, n_lags);
end

p_value_holm = reshape(holm_adjust(p_value(:)), size(p_value));

result = struct();
result.target_names = target_names;
result.actor_relation = settings.actor_relation;
result.group_names = group_names;
result.group_behaviors = group_behaviors;
result.behavior_types = behavior_types;
result.wrapped_plot_time = plot_time;
result.pre_index = pre_index;
result.wrap_index = wrap_index;
result.post_index = post_index;
result.pre_lag = pre_lag;
result.post_lag = post_lag;
result.lag_seconds = post_lag;
result.flank_seconds = recordings(1).flank_seconds;
result.wrapped_bins = recordings(1).wrapped_bins;
result.probability = probability;
result.pre_probability = pre_probability;
result.post_probability = post_probability;
result.behavior_probability = behavior_probability;
result.behavior_probability_self = behavior_probability_self;
result.behavior_probability_partner = behavior_probability_partner;
result.behavior_pre_probability = behavior_pre_probability;
result.behavior_post_probability = behavior_post_probability;
result.observed_difference = observed_difference;
result.observed_difference_by_lag = observed_difference_by_lag;
result.p_value = p_value;
result.p_value_holm = p_value_holm;
result.p_value_by_lag = p_value_by_lag;
result.p_value_by_lag_fwer_group = p_value_by_lag_fwer_group;
result.p_value_by_lag_fwer_target = p_value_by_lag_fwer_target;
result.observed_postbin_minus_pre = observed_postbin_minus_pre;
result.p_value_postbin_vs_pre = p_value_postbin_vs_pre;
result.p_value_postbin_vs_pre_fwer_group = p_value_postbin_vs_pre_fwer_group;
result.p_value_postbin_vs_pre_fwer_target = ...
    p_value_postbin_vs_pre_fwer_target;
result.null_difference = null_difference;
result.n_events_test = n_events_test;
result.n_events_available = n_events_available;
result.n_sessions = n_sessions;
result.session_names_by_target = session_names_by_target;
result.n_permutations = settings.n_permutations;
result.random_seed = settings.random_seed;
result.significance_alpha = 0.05;
result.test_description = [ ...
    'Two-sided paired sign-flip of (post-offset minus pre-onset) mean ' ...
    'occupancy. Signs flip independently for each event within a play ' ...
    'session. Difference is post minus pre. Window p-values are Holm ' ...
    'corrected across targets and behavior groups. Lag-wise p-values use ' ...
    'the same signs, with max-T correction within each target.'];
result.postbin_test_description = [ ...
    'One-sided paired sign-flip. Each post-offset bin is compared with ' ...
    'the mean occupancy across the whole pre-onset flank. The counted ' ...
    'tail is post higher than that before rate. The same event signs as ' ...
    'the window test are reused. Max-T correction is within group across ' ...
    'post bins, and across groups and post bins within each target.'];
end

function values = group_occupancy(event_matrix, group_index, in_session)
n_events = size(event_matrix, 1);
n_time = size(event_matrix, 3);
n_groups = numel(group_index);
values = nan(n_events, n_groups, n_time);
missing = ~in_session;
for group_n = 1:n_groups
    occupied = reshape(any(event_matrix(:, group_index{group_n}, :), 2), ...
        n_events, n_time);
    occupied = double(occupied);
    occupied(missing) = NaN;
    values(:, group_n, :) = reshape(occupied, n_events, 1, n_time);
end
end

function values = masked_occupancy(event_matrix, in_session)
values = double(event_matrix);
n_behaviors = size(event_matrix, 2);
missing = ~in_session;
for behavior_n = 1:n_behaviors
    slice = reshape(values(:, behavior_n, :), size(in_session, 1), []);
    slice(missing) = NaN;
    values(:, behavior_n, :) = reshape(slice, size(slice, 1), 1, []);
end
end

function settings = before_after_settings(overrides)
if ~isstruct(overrides) || ~isscalar(overrides)
    error('compare_behavior_before_after_probability:Settings', ...
        'settings must be a scalar structure.');
end
settings = struct( ...
    'actor_relation', 'either', ...
    'n_permutations', 10000, ...
    'random_seed', 1);
names = fieldnames(overrides);
for field_n = 1:numel(names)
    if ~isfield(settings, names{field_n})
        error('compare_behavior_before_after_probability:Settings', ...
            'Unknown setting "%s".', names{field_n});
    end
    settings.(names{field_n}) = overrides.(names{field_n});
end
settings.actor_relation = validatestring(settings.actor_relation, ...
    {'self', 'partner', 'either'});
if ~isnumeric(settings.n_permutations) || ...
        ~isscalar(settings.n_permutations) || ...
        ~isfinite(settings.n_permutations) || ...
        settings.n_permutations < 1 || ...
        settings.n_permutations ~= round(settings.n_permutations)
    error('compare_behavior_before_after_probability:Permutations', ...
        'n_permutations must be a positive integer.');
end
if ~isnumeric(settings.random_seed) || ...
        ~isscalar(settings.random_seed) || ...
        ~isfinite(settings.random_seed) || ...
        settings.random_seed < 0 || ...
        settings.random_seed > 2^32 - 1 || ...
        settings.random_seed ~= round(settings.random_seed)
    error('compare_behavior_before_after_probability:Seed', ...
        'random_seed must be an integer from 0 through 2^32 - 1.');
end
end

function names = normalize_targets(values)
if isstring(values)
    names = cellstr(values(:).');
elseif iscell(values)
    names = cellfun(@char, values, 'UniformOutput', false);
    names = names(:).';
else
    error('compare_behavior_before_after_probability:Targets', ...
        'target_names must be a cell or string array.');
end
if isempty(names) || any(cellfun(@isempty, names)) || ...
        numel(unique(names)) ~= numel(names)
    error('compare_behavior_before_after_probability:Targets', ...
        'Select one or more different, non-empty target behaviors.');
end
end

function adjusted = holm_adjust(p_values)
[sorted, order] = sort(p_values);
adjusted_sorted = nan(size(sorted));
running_max = 0;
n_tests = numel(sorted);
for test_n = 1:n_tests
    running_max = max(running_max, ...
        (n_tests - test_n + 1) * sorted(test_n));
    adjusted_sorted(test_n) = min(1, running_max);
end
adjusted = nan(size(p_values));
adjusted(order) = adjusted_sorted;
end
