function recording = GENERATE_ACTOR_BEHAVIOR_PROBABILITY( ...
    behavior_source, options)
%GENERATE_ACTOR_BEHAVIOR_PROBABILITY Build matrices by target-event actor.
%
% Every Partners session contains one actor_data entry per animal observed
% in that session. For actor A, self_behavior contains behavior performed
% by A and other_behavior contains the logical OR across every animal
% other than A. Arrays are target events x behavior types x time bins.

if nargin < 2 || isempty(options)
    options = actor_behavior_probability_options();
else
    options = actor_behavior_probability_options(options);
end

[behavior, recording_name] = load_behavior(behavior_source);
behavior.Animal = text_column(behavior.Animal);
behavior.Type = text_column(behavior.Type);
behavior.Type2 = normalize_type(behavior.Type);

if any(~isfinite(behavior.Start)) || any(~isfinite(behavior.End)) || ...
        any(behavior.End < behavior.Start)
    error('GENERATE_ACTOR_BEHAVIOR_PROBABILITY:Times', ...
        '%s contains invalid behavior intervals.', recording_name);
end

session_mask = strcmp(behavior.Animal, 'Session_structure') & ...
    strcmp(behavior.Type, 'Partners session');
session_rows = find(session_mask);
session_ranges = behavior{session_mask, {'Start', 'End'}};
if isempty(session_ranges)
    error('GENERATE_ACTOR_BEHAVIOR_PROBABILITY:Sessions', ...
        '%s has no Session_structure / Partners session rows.', ...
        recording_name);
end
[session_ranges, order] = sortrows(session_ranges, 1);
session_rows = session_rows(order);
if any(session_ranges(:, 2) <= session_ranges(:, 1)) || ...
        any(session_ranges(2:end, 1) < session_ranges(1:end-1, 2))
    error('GENERATE_ACTOR_BEHAVIOR_PROBABILITY:Sessions', ...
        '%s contains invalid or overlapping partner sessions.', ...
        recording_name);
end

n_behaviors = numel(options.all_behaviors);
n_time = numel(options.time);
n_wrapped = numel(options.wrapped_plot_time);
empty_actor = struct( ...
    'actor', '', ...
    'other_animals', {{}}, ...
    'target_behavior', {{}}, ...
    'target_start', zeros(0, 1), ...
    'target_offset', zeros(0, 1), ...
    'target_source_row', zeros(0, 1), ...
    'self_behavior', false(0, n_behaviors, n_time), ...
    'other_behavior', false(0, n_behaviors, n_time), ...
    'self_behavior_wrapped', false(0, n_behaviors, n_wrapped), ...
    'other_behavior_wrapped', false(0, n_behaviors, n_wrapped), ...
    'wrapped_in_session', false(0, n_wrapped));
empty_session = struct( ...
    'play_session_number', [], ...
    'play_session', [], ...
    'play_session_table_row', [], ...
    'participant_animals', {{}}, ...
    'actor_data', {repmat(empty_actor, 0, 1)});
sessions = repmat(empty_session, size(session_ranges, 1), 1);

real_animal_row = ~ismember(behavior.Animal, ...
    {'', 'Session_structure', 'Reversal'});
measured_row = real_animal_row & ...
    ismember(behavior.Type2, options.all_behaviors);

for session_n = 1:size(session_ranges, 1)
    session_start = session_ranges(session_n, 1);
    session_end = session_ranges(session_n, 2);
    overlaps = behavior.End > session_start;
    if session_n < size(session_ranges, 1)
        overlaps = overlaps & behavior.Start < session_end;
    else
        overlaps = overlaps & behavior.Start <= session_end;
    end
    participant_rows = real_animal_row & overlaps;
    participants = unique(behavior.Animal(participant_rows), 'stable');
    participants(cellfun(@isempty, participants)) = [];
    in_session = measured_row & overlaps;

    actor_data = repmat(empty_actor, numel(participants), 1);
    for actor_n = 1:numel(participants)
        actor = participants{actor_n};
        others = participants(~strcmp(participants, actor));
        target_mask = measured_row & ...
            ismember(behavior.Type2, options.target_behaviors) & ...
            strcmp(behavior.Animal, actor) & ...
            behavior.End > session_start & behavior.End <= session_end;
        target_rows = find(target_mask);
        [~, target_order] = sortrows( ...
            [behavior.End(target_rows), behavior.Start(target_rows)], [1 2]);
        target_rows = target_rows(target_order);
        n_targets = numel(target_rows);

        self_behavior = false(n_targets, n_behaviors, n_time);
        other_behavior = false(n_targets, n_behaviors, n_time);
        self_wrapped = false(n_targets, n_behaviors, n_wrapped);
        other_wrapped = false(n_targets, n_behaviors, n_wrapped);
        wrapped_in_session = false(n_targets, n_wrapped);

        for target_n = 1:n_targets
            target_row = target_rows(target_n);
            onset = behavior.Start(target_row);
            offset = behavior.End(target_row);
            absolute_time = offset + options.time;
            sampled_in_session = absolute_time >= session_start & ...
                absolute_time <= session_end;
            wrapped_time = [onset + options.pre_lag, ...
                onset + (offset - onset) * options.wrap_fraction, ...
                offset + options.post_lag];
            wrapped_mask = wrapped_time >= session_start & ...
                wrapped_time <= session_end;
            wrapped_in_session(target_n, :) = wrapped_mask;

            for behavior_n = 1:n_behaviors
                type_mask = strcmp(behavior.Type2, ...
                    options.all_behaviors{behavior_n});
                self_rows = find(in_session & type_mask & ...
                    strcmp(behavior.Animal, actor));
                other_rows = find(in_session & type_mask & ...
                    ismember(behavior.Animal, others));
                self_behavior(target_n, behavior_n, :) = reshape( ...
                    interval_occupancy(absolute_time, behavior, self_rows) & ...
                    sampled_in_session, 1, 1, n_time);
                other_behavior(target_n, behavior_n, :) = reshape( ...
                    interval_occupancy(absolute_time, behavior, other_rows) & ...
                    sampled_in_session, 1, 1, n_time);
                self_wrapped(target_n, behavior_n, :) = reshape( ...
                    interval_occupancy(wrapped_time, behavior, self_rows) & ...
                    wrapped_mask, 1, 1, n_wrapped);
                other_wrapped(target_n, behavior_n, :) = reshape( ...
                    interval_occupancy(wrapped_time, behavior, other_rows) & ...
                    wrapped_mask, 1, 1, n_wrapped);
            end
        end

        actor_data(actor_n).actor = actor;
        actor_data(actor_n).other_animals = others(:).';
        actor_data(actor_n).target_behavior = behavior.Type2(target_rows);
        actor_data(actor_n).target_start = behavior.Start(target_rows);
        actor_data(actor_n).target_offset = behavior.End(target_rows);
        actor_data(actor_n).target_source_row = target_rows;
        actor_data(actor_n).self_behavior = self_behavior;
        actor_data(actor_n).other_behavior = other_behavior;
        actor_data(actor_n).self_behavior_wrapped = self_wrapped;
        actor_data(actor_n).other_behavior_wrapped = other_wrapped;
        actor_data(actor_n).wrapped_in_session = wrapped_in_session;
    end

    sessions(session_n).play_session_number = session_n;
    sessions(session_n).play_session = session_ranges(session_n, :);
    sessions(session_n).play_session_table_row = session_rows(session_n);
    sessions(session_n).participant_animals = participants(:).';
    sessions(session_n).actor_data = actor_data;
end

recording = struct();
recording.recording_name = recording_name;
recording.target_behaviors = options.target_behaviors;
recording.behavior_types = options.all_behaviors;
recording.time = options.time;
recording.time_window = options.time_window;
recording.bin_size = options.bin_size;
recording.comparison_window = options.comparison_window;
recording.flank_seconds = options.flank_seconds;
recording.wrapped_bins = options.wrapped_bins;
recording.pre_lag = options.pre_lag;
recording.post_lag = options.post_lag;
recording.wrap_fraction = options.wrap_fraction;
recording.wrapped_plot_time = options.wrapped_plot_time;
recording.pre_index = 1:numel(options.pre_lag);
recording.wrap_index = numel(options.pre_lag) + (1:options.wrapped_bins);
recording.post_index = numel(options.pre_lag) + options.wrapped_bins + ...
    (1:numel(options.post_lag));
recording.alignment = options.alignment;
recording.self_definition = options.self_definition;
recording.other_definition = options.other_definition;
recording.sessions = sessions;
end

function occupancy = interval_occupancy(time, behavior, rows)
if isempty(rows)
    occupancy = false(1, numel(time));
    return
end
occupancy = any(bsxfun(@ge, time, behavior.Start(rows)) & ...
    bsxfun(@le, time, behavior.End(rows)), 1);
end

function values = normalize_type(values)
values = text_column(values);
values(ismember(values, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
values(ismember(values, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
values(ismember(values, {'Scratch', 'Scratching'})) = {'Scratching'};
values(cellfun(@isempty, values)) = {'Other'};
end

function [behavior, recording_name] = load_behavior(source)
if istable(source)
    behavior = source;
    recording_name = 'table_input';
elseif ischar(source) || (isstring(source) && isscalar(source))
    source = char(source);
    if ~isfile(source)
        error('GENERATE_ACTOR_BEHAVIOR_PROBABILITY:File', ...
            'Behavior file not found: %s', source);
    end
    behavior = readtable(source);
    [~, recording_name] = fileparts(source);
else
    error('GENERATE_ACTOR_BEHAVIOR_PROBABILITY:Input', ...
        'behavior_source must be a table or behavior-file path.');
end

required = {'Animal', 'Start', 'End', 'Length', 'Type'};
if ~all(ismember(required, behavior.Properties.VariableNames))
    if width(behavior) == 6
        behavior(:, 2) = [];
    end
    if width(behavior) ~= 5
        error('GENERATE_ACTOR_BEHAVIOR_PROBABILITY:Columns', ...
            'Behavior input must contain Animal, Start, End, Length, and Type.');
    end
    behavior.Properties.VariableNames = required;
end
behavior = behavior(:, required);
end

function values = text_column(values)
values = cellstr(string(values));
values = strtrim(values);
values(ismissing(string(values))) = {''};
end
