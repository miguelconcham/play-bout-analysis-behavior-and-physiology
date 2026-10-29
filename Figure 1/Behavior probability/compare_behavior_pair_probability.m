function comparison = compare_behavior_pair_probability( ...
    recordings, behavior_pair, settings)
%COMPARE_BEHAVIOR_PAIR_PROBABILITY Compare grouped behavior probabilities.
%
% Target events are pooled across recordings. For inference, target labels
% are shuffled within each play session, preserving each session's event
% counts and controlling session-level differences. A session contributes
% to the permutation test only when both target behaviors occur in it.

if nargin < 3
    settings = struct();
end
settings = comparison_settings(settings);
behavior_pair = normalize_pair(behavior_pair);
if isempty(recordings)
    error('compare_behavior_pair_probability:Empty', ...
        'No behavior-probability recordings were supplied.');
end

behavior_types = recordings(1).behavior_types;
time = recordings(1).time;
n_time = numel(time);
group_names = {'Play', 'Aggressive', 'Exploratory'};
group_behaviors = { ...
    {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'}, ...
    {'PounceI', 'Bite'}, ...
    {'Grooming', 'Rearing', 'Sniffing', 'Scratching'}};
n_groups = numel(group_names);

for group_n = 1:n_groups
    missing = setdiff(group_behaviors{group_n}, behavior_types);
    if ~isempty(missing)
        error('compare_behavior_pair_probability:BehaviorTypes', ...
            'Recording behavior list is missing: %s.', strjoin(missing, ', '));
    end
end

presence = false(0, n_groups, n_time);
target_label = zeros(0, 1);
session_index = zeros(0, 1);
session_names = cell(0, 1);
session_counter = 0;

for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    if ~isequal(recording.behavior_types, behavior_types) || ...
            ~isequal(recording.time, time)
        error('compare_behavior_pair_probability:Parameters', ...
            'All recordings must use the same behavior list and time bins.');
    end
    for session_n = 1:numel(recording.sessions)
        session_counter = session_counter + 1;
        session = recording.sessions(session_n);
        [selected, labels] = ismember(session.target_behavior, behavior_pair);
        rows = find(selected);
        if isempty(rows)
            continue
        end
        labels = labels(rows);
        switch settings.actor_relation
            case 'self'
                event_matrix = session.self_behavior(rows, :, :);
            case 'partner'
                event_matrix = session.partner_behavior(rows, :, :);
            case 'either'
                event_matrix = session.self_behavior(rows, :, :) | ...
                    session.partner_behavior(rows, :, :);
        end

        this_presence = false(numel(rows), n_groups, n_time);
        for group_n = 1:n_groups
            group_index = ismember(behavior_types, ...
                group_behaviors{group_n});
            this_presence(:, group_n, :) = reshape(any( ...
                event_matrix(:, group_index, :), 2), ...
                numel(rows), 1, n_time);
        end
        new_rows = size(presence, 1) + (1:numel(rows));
        presence(new_rows, :, :) = this_presence;
        target_label(new_rows, 1) = labels;
        session_index(new_rows, 1) = session_counter;
        session_names{session_counter, 1} = sprintf('%s P%d', ...
            recording.recording_name, session.play_session_number);
    end
end

if any(arrayfun(@(x) ~any(target_label == x), 1:2))
    error('compare_behavior_pair_probability:MissingTarget', ...
        ['Both selected target behaviors must be present. Re-run Estimate ' ...
        'with target_behaviors containing %s and %s.'], ...
        behavior_pair{1}, behavior_pair{2});
end

probability_all_sessions = nan(n_groups, 2, n_time);
n_events_all_sessions = zeros(1, 2);
for target_n = 1:2
    target_rows = target_label == target_n;
    n_events_all_sessions(target_n) = sum(target_rows);
    probability_all_sessions(:, target_n, :) = reshape(mean( ...
        double(presence(target_rows, :, :)), 1), ...
        n_groups, 1, n_time);
end

all_sessions = unique(session_index);
eligible_sessions = false(size(all_sessions));
for session_n = 1:numel(all_sessions)
    labels = unique(target_label(session_index == all_sessions(session_n)));
    eligible_sessions(session_n) = all(ismember(1:2, labels));
end
eligible_session_ids = all_sessions(eligible_sessions);
test_rows = ismember(session_index, eligible_session_ids);
test_presence = presence(test_rows, :, :);
test_labels = target_label(test_rows);
test_sessions = session_index(test_rows);
if isempty(test_labels)
    error('compare_behavior_pair_probability:NoPairedSession', ...
        'No play session contains both selected target behaviors.');
end

if numel(time) > 1
    time_resolution = min(abs(diff(time)));
else
    time_resolution = 1;
end
time_tolerance = max(eps(max(1, max(abs(time)))) * 16, ...
    time_resolution * 1e-9);
window_index = time >= settings.test_window(1) - time_tolerance & ...
    time <= settings.test_window(2) + time_tolerance & ...
    time > time_tolerance;
if ~any(window_index)
    error('compare_behavior_pair_probability:TestWindow', ...
        'test_window must contain at least one strictly post-offset time bin.');
end
event_score = reshape(mean(double( ...
    test_presence(:, :, window_index)), 3), [], n_groups);
n_window_bins = sum(window_index);
bin_score = reshape(double(test_presence(:, :, window_index)), ...
    numel(test_labels), n_groups * n_window_bins);
n_target_1 = sum(test_labels == 1);
n_target_2 = sum(test_labels == 2);
observed_weight = double(test_labels == 1) / n_target_1 - ...
    double(test_labels == 2) / n_target_2;
observed_difference = observed_weight.' * event_score;
observed_difference_by_bin = reshape( ...
    observed_weight.' * bin_score, n_groups, n_window_bins);

previous_rng = rng;
rng_cleanup = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(settings.random_seed, 'twister');
null_difference = nan(settings.n_permutations, n_groups);
bin_exceedances = zeros(n_groups, n_window_bins);
null_max_by_group = nan(settings.n_permutations, n_groups);
null_max_all = nan(settings.n_permutations, 1);
for permutation_n = 1:settings.n_permutations
    shuffled_labels = test_labels;
    for session_n = 1:numel(eligible_session_ids)
        rows = find(test_sessions == eligible_session_ids(session_n));
        shuffled_labels(rows) = shuffled_labels(rows(randperm(numel(rows))));
    end
    shuffled_weight = double(shuffled_labels == 1) / n_target_1 - ...
        double(shuffled_labels == 2) / n_target_2;
    null_difference(permutation_n, :) = ...
        shuffled_weight.' * event_score;
    shuffled_bin_difference = reshape( ...
        shuffled_weight.' * bin_score, n_groups, n_window_bins);
    bin_exceedances = bin_exceedances + ...
        (abs(shuffled_bin_difference) >= ...
        abs(observed_difference_by_bin));
    null_max_by_group(permutation_n, :) = ...
        max(abs(shuffled_bin_difference), [], 2).';
    null_max_all(permutation_n) = ...
        max(abs(shuffled_bin_difference), [], 'all');
end
p_value = (1 + sum(abs(null_difference) >= ...
    abs(observed_difference), 1)) / (settings.n_permutations + 1);
p_value_holm = holm_adjust(p_value);
p_value_by_bin = (1 + bin_exceedances) / ...
    (settings.n_permutations + 1);
p_value_by_bin_fwer_group = nan(n_groups, n_window_bins);
for group_n = 1:n_groups
    p_value_by_bin_fwer_group(group_n, :) = (1 + sum(bsxfun(@ge, ...
        null_max_by_group(:, group_n), ...
        abs(observed_difference_by_bin(group_n, :))), 1)) / ...
        (settings.n_permutations + 1);
end
p_value_by_bin_fwer_all = (1 + sum(bsxfun(@ge, null_max_all, ...
    abs(observed_difference_by_bin(:).')), 1)) / ...
    (settings.n_permutations + 1);
p_value_by_bin_fwer_all = reshape( ...
    p_value_by_bin_fwer_all, n_groups, n_window_bins);

n_events_test = [sum(test_labels == 1), sum(test_labels == 2)];
probability = nan(n_groups, 2, n_time);
for target_n = 1:2
    target_rows = test_labels == target_n;
    probability(:, target_n, :) = reshape(mean( ...
        double(test_presence(target_rows, :, :)), 1), ...
        n_groups, 1, n_time);
end
comparison = struct();
comparison.behavior_pair = behavior_pair;
comparison.actor_relation = settings.actor_relation;
comparison.group_names = group_names;
comparison.group_behaviors = group_behaviors;
comparison.time = time;
comparison.probability = probability;
comparison.n_events_plot = n_events_test;
comparison.probability_all_sessions = probability_all_sessions;
comparison.n_events_all_sessions = n_events_all_sessions;
comparison.test_window = settings.test_window;
comparison.test_time = time(window_index);
comparison.n_permutations = settings.n_permutations;
comparison.random_seed = settings.random_seed;
comparison.eligible_session_ids = eligible_session_ids;
comparison.eligible_session_names = session_names(eligible_session_ids);
comparison.n_eligible_sessions = numel(eligible_session_ids);
comparison.n_events_test = n_events_test;
comparison.observed_difference = observed_difference;
comparison.p_value = p_value;
comparison.p_value_holm = p_value_holm;
comparison.null_difference = null_difference;
comparison.observed_difference_by_bin = observed_difference_by_bin;
comparison.p_value_by_bin = p_value_by_bin;
comparison.p_value_by_bin_fwer_group = p_value_by_bin_fwer_group;
comparison.p_value_by_bin_fwer_all = p_value_by_bin_fwer_all;
comparison.bin_significance_alpha = 0.05;
comparison.test_description = [ ...
    'Two-sided post-offset event-label permutation stratified by play ' ...
    'session; time zero excluded; three behavior-group p-values Holm corrected'];
comparison.bin_test_description = [ ...
    'Bin-wise two-sided stratified permutation; raw, within-group max-T, ' ...
    'and all-group/all-bin max-T family-wise corrected p-values'];
end

function settings = comparison_settings(overrides)
if ~isstruct(overrides) || ~isscalar(overrides)
    error('compare_behavior_pair_probability:Settings', ...
        'settings must be a scalar structure.');
end
settings = struct( ...
    'actor_relation', 'either', ...
    'test_window', [0 2], ...
    'n_permutations', 10000, ...
    'random_seed', 1);
names = fieldnames(overrides);
for field_n = 1:numel(names)
    if ~isfield(settings, names{field_n})
        error('compare_behavior_pair_probability:Settings', ...
            'Unknown setting "%s".', names{field_n});
    end
    settings.(names{field_n}) = overrides.(names{field_n});
end
settings.actor_relation = validatestring(settings.actor_relation, ...
    {'self', 'partner', 'either'});
if ~isnumeric(settings.test_window) || numel(settings.test_window) ~= 2 || ...
        any(~isfinite(settings.test_window)) || ...
        settings.test_window(1) < 0 || ...
        settings.test_window(2) <= settings.test_window(1)
    error('compare_behavior_pair_probability:TestWindow', ...
        ['test_window must be an ascending post-offset range. Time zero ' ...
        'is excluded from inference because it contains the target itself.']);
end
settings.test_window = double(settings.test_window(:).');
if ~isnumeric(settings.n_permutations) || ...
        ~isscalar(settings.n_permutations) || ...
        ~isfinite(settings.n_permutations) || ...
        settings.n_permutations < 1 || ...
        settings.n_permutations ~= round(settings.n_permutations)
    error('compare_behavior_pair_probability:Permutations', ...
        'n_permutations must be a positive integer.');
end
if ~isnumeric(settings.random_seed) || ...
        ~isscalar(settings.random_seed) || ...
        ~isfinite(settings.random_seed) || ...
        settings.random_seed < 0 || ...
        settings.random_seed > 2^32 - 1 || ...
        settings.random_seed ~= round(settings.random_seed)
    error('compare_behavior_pair_probability:Seed', ...
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
    error('compare_behavior_pair_probability:Pair', ...
        'behavior_pair must be a two-element cell or string array.');
end
if numel(pair) ~= 2 || any(cellfun(@isempty, pair)) || ...
        strcmp(pair{1}, pair{2})
    error('compare_behavior_pair_probability:Pair', ...
        'Select two different, non-empty target behaviors.');
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
