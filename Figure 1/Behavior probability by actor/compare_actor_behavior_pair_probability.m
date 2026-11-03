function comparison = compare_actor_behavior_pair_probability( ...
    recordings, behavior_pair, settings)
%COMPARE_ACTOR_BEHAVIOR_PAIR_PROBABILITY Compare targets within actor strata.
%
% Target labels are shuffled only among events performed by the same animal
% in the same Partners session. Run separately with relation='self' and
% relation='other' to preserve the two actor-relative categories.

if nargin < 3
    settings = struct();
end
settings = comparison_settings(settings);
behavior_pair = normalize_pair(behavior_pair);
if isempty(recordings)
    error('compare_actor_behavior_pair_probability:Empty', ...
        'No actor-relative recordings were supplied.');
end

behavior_types = recordings(1).behavior_types;
time = recordings(1).time;
n_time = numel(time);
[group_names, group_behaviors, group_index] = ...
    behavior_groups(behavior_types);
n_groups = numel(group_names);

presence = false(0, n_groups, n_time);
target_label = zeros(0, 1);
stratum_index = zeros(0, 1);
stratum_names = cell(0, 1);
stratum_counter = 0;

for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    if ~isequal(recording.behavior_types, behavior_types) || ...
            ~isequal(recording.time, time)
        error('compare_actor_behavior_pair_probability:Parameters', ...
            'All recordings must share behavior types and time bins.');
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
                actor.target_behavior, behavior_pair);
            rows = find(selected);
            if isempty(rows)
                continue
            end
            labels = labels(rows);
            if strcmp(settings.relation, 'self')
                event_matrix = actor.self_behavior(rows, :, :);
            else
                event_matrix = actor.other_behavior(rows, :, :);
            end
            this_presence = false(numel(rows), n_groups, n_time);
            for group_n = 1:n_groups
                this_presence(:, group_n, :) = reshape(any( ...
                    event_matrix(:, group_index{group_n}, :), 2), ...
                    numel(rows), 1, n_time);
            end
            new_rows = size(presence, 1) + (1:numel(rows));
            presence(new_rows, :, :) = this_presence;
            target_label(new_rows, 1) = labels;
            stratum_index(new_rows, 1) = stratum_counter;
        end
    end
end

if any(arrayfun(@(target_n) ~any(target_label == target_n), 1:2))
    error('compare_actor_behavior_pair_probability:MissingTarget', ...
        'Both selected target behaviors must be present.');
end

probability_all_strata = nan(n_groups, 2, n_time);
n_events_all_strata = zeros(1, 2);
for target_n = 1:2
    rows = target_label == target_n;
    n_events_all_strata(target_n) = sum(rows);
    probability_all_strata(:, target_n, :) = reshape(mean( ...
        double(presence(rows, :, :)), 1), n_groups, 1, n_time);
end

all_strata = unique(stratum_index);
eligible = false(size(all_strata));
for stratum_n = 1:numel(all_strata)
    labels = unique(target_label(stratum_index == all_strata(stratum_n)));
    eligible(stratum_n) = all(ismember(1:2, labels));
end
eligible_ids = all_strata(eligible);
test_rows = ismember(stratum_index, eligible_ids);
test_presence = presence(test_rows, :, :);
test_labels = target_label(test_rows);
test_strata = stratum_index(test_rows);
if isempty(test_labels)
    error('compare_actor_behavior_pair_probability:NoPairedStratum', ...
        ['No target actor in a session performed both selected target ' ...
        'behaviors.']);
end

time_tolerance = max(eps(max(1, max(abs(time)))) * 16, ...
    min(abs(diff(time))) * 1e-9);
window_index = time >= settings.test_window(1) - time_tolerance & ...
    time <= settings.test_window(2) + time_tolerance & ...
    time > time_tolerance;
if ~any(window_index)
    error('compare_actor_behavior_pair_probability:TestWindow', ...
        'test_window contains no strictly post-offset sample.');
end
n_window_bins = sum(window_index);
event_score = reshape(mean(double( ...
    test_presence(:, :, window_index)), 3), [], n_groups);
bin_score = reshape(double(test_presence(:, :, window_index)), ...
    numel(test_labels), n_groups * n_window_bins);
n_by_target = [sum(test_labels == 1), sum(test_labels == 2)];
observed_weight = double(test_labels == 1) / n_by_target(1) - ...
    double(test_labels == 2) / n_by_target(2);
observed_difference = observed_weight.' * event_score;
observed_by_bin = reshape(observed_weight.' * bin_score, ...
    n_groups, n_window_bins);

previous_rng = rng;
rng_cleanup = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(settings.random_seed, 'twister');
null_difference = nan(settings.n_permutations, n_groups);
bin_exceedances = zeros(n_groups, n_window_bins);
null_max_group = nan(settings.n_permutations, n_groups);
null_max_all = nan(settings.n_permutations, 1);
for permutation_n = 1:settings.n_permutations
    shuffled = test_labels;
    for stratum_n = 1:numel(eligible_ids)
        rows = find(test_strata == eligible_ids(stratum_n));
        shuffled(rows) = shuffled(rows(randperm(numel(rows))));
    end
    shuffled_weight = double(shuffled == 1) / n_by_target(1) - ...
        double(shuffled == 2) / n_by_target(2);
    null_difference(permutation_n, :) = shuffled_weight.' * event_score;
    shuffled_by_bin = reshape(shuffled_weight.' * bin_score, ...
        n_groups, n_window_bins);
    bin_exceedances = bin_exceedances + ...
        (abs(shuffled_by_bin) >= abs(observed_by_bin));
    null_max_group(permutation_n, :) = ...
        max(abs(shuffled_by_bin), [], 2).';
    null_max_all(permutation_n) = max(abs(shuffled_by_bin), [], 'all');
end

p_value = (1 + sum(abs(null_difference) >= ...
    abs(observed_difference), 1)) / (settings.n_permutations + 1);
p_value_holm = holm_adjust(p_value);
p_bin = (1 + bin_exceedances) / (settings.n_permutations + 1);
p_bin_group = nan(n_groups, n_window_bins);
for group_n = 1:n_groups
    p_bin_group(group_n, :) = (1 + sum(bsxfun(@ge, ...
        null_max_group(:, group_n), abs(observed_by_bin(group_n, :))), 1)) / ...
        (settings.n_permutations + 1);
end
p_bin_all = reshape((1 + sum(bsxfun(@ge, null_max_all, ...
    abs(observed_by_bin(:).')), 1)) / (settings.n_permutations + 1), ...
    n_groups, n_window_bins);

probability = nan(n_groups, 2, n_time);
for target_n = 1:2
    rows = test_labels == target_n;
    probability(:, target_n, :) = reshape(mean( ...
        double(test_presence(rows, :, :)), 1), n_groups, 1, n_time);
end

comparison = struct();
comparison.behavior_pair = behavior_pair;
comparison.relation = settings.relation;
comparison.group_names = group_names;
comparison.group_behaviors = group_behaviors;
comparison.time = time;
comparison.probability = probability;
comparison.probability_all_strata = probability_all_strata;
comparison.n_events_all_strata = n_events_all_strata;
comparison.test_window = settings.test_window;
comparison.test_time = time(window_index);
comparison.n_permutations = settings.n_permutations;
comparison.random_seed = settings.random_seed;
comparison.eligible_stratum_ids = eligible_ids;
comparison.eligible_stratum_names = stratum_names(eligible_ids);
comparison.n_eligible_strata = numel(eligible_ids);
comparison.n_events_test = n_by_target;
comparison.observed_difference = observed_difference;
comparison.p_value = p_value;
comparison.p_value_holm = p_value_holm;
comparison.null_difference = null_difference;
comparison.observed_difference_by_bin = observed_by_bin;
comparison.p_value_by_bin = p_bin;
comparison.p_value_by_bin_fwer_group = p_bin_group;
comparison.p_value_by_bin_fwer_all = p_bin_all;
comparison.bin_significance_alpha = 0.05;
comparison.test_description = [ ...
    'Two-sided target-label permutation within target actor and session'];
end

function settings = comparison_settings(overrides)
settings = struct('relation', 'self', 'test_window', [0 2], ...
    'n_permutations', 10000, 'random_seed', 1);
names = fieldnames(overrides);
for field_n = 1:numel(names)
    if ~isfield(settings, names{field_n})
        error('compare_actor_behavior_pair_probability:Settings', ...
            'Unknown setting "%s".', names{field_n});
    end
    settings.(names{field_n}) = overrides.(names{field_n});
end
settings.relation = validatestring(settings.relation, {'self', 'other'});
if ~isnumeric(settings.test_window) || numel(settings.test_window) ~= 2 || ...
        any(~isfinite(settings.test_window)) || ...
        settings.test_window(1) < 0 || ...
        settings.test_window(2) <= settings.test_window(1)
    error('compare_actor_behavior_pair_probability:TestWindow', ...
        'test_window must be an ascending post-offset range.');
end
settings.test_window = double(settings.test_window(:).');
if ~isnumeric(settings.n_permutations) || ...
        ~isscalar(settings.n_permutations) || ...
        ~isfinite(settings.n_permutations) || ...
        settings.n_permutations < 1 || ...
        settings.n_permutations ~= round(settings.n_permutations)
    error('compare_actor_behavior_pair_probability:Permutations', ...
        'n_permutations must be a positive integer.');
end
if ~isnumeric(settings.random_seed) || ...
        ~isscalar(settings.random_seed) || ...
        ~isfinite(settings.random_seed) || settings.random_seed < 0 || ...
        settings.random_seed > 2^32 - 1 || ...
        settings.random_seed ~= round(settings.random_seed)
    error('compare_actor_behavior_pair_probability:Seed', ...
        'random_seed must be an integer from 0 through 2^32 - 1.');
end
end

function pair = normalize_pair(values)
if isstring(values)
    pair = cellstr(values(:).');
elseif iscell(values)
    pair = cellfun(@char, values, 'UniformOutput', false);
    pair = pair(:).';
else
    error('compare_actor_behavior_pair_probability:Pair', ...
        'behavior_pair must contain two names.');
end
if numel(pair) ~= 2 || strcmp(pair{1}, pair{2})
    error('compare_actor_behavior_pair_probability:Pair', ...
        'Select two different target behaviors.');
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
        error('compare_actor_behavior_pair_probability:BehaviorTypes', ...
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
