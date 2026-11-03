function recording = GENERATE_BEHAVIOR_PROBABILITY( ...
    behavior_source, focal_animal, options)
%GENERATE_BEHAVIOR_PROBABILITY Build binary behavior matrices by play session.
%
% Each target-behavior occurrence contributes one row, aligned to its
% offset. self_behavior describes the target performer; partner_behavior
% describes the other animal(s) in that Partners session. Both arrays are
% N target events x N behavior types x N time bins.
%
% self_behavior_wrapped and partner_behavior_wrapped use a second axis:
% real-time samples from -T to just before onset, then wrapped_bins
% samples evenly spaced from onset through offset, then real-time samples
% from just after offset through offset+T. wrapped_in_session marks which
% of those samples fall inside the Partners session.

if nargin < 2
    focal_animal = '';
end
if nargin < 3 || isempty(options)
    options = behavior_probability_options();
else
    options = behavior_probability_options(options);
end

[behavior, recording_name] = load_behavior(behavior_source);
if isempty(focal_animal)
    focal_animal = infer_focal_animal(recording_name);
end
focal_animal = char(string(focal_animal));

behavior.Animal = text_column(behavior.Animal);
behavior.Type = text_column(behavior.Type);
behavior.Type2 = behavior.Type;
behavior.Type2(ismember(behavior.Type2, {'Pounce_A', 'Pounce_B'})) = ...
    {'Pounce'};
behavior.Type2(ismember(behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = ...
    {'PounceI'};
behavior.Type2(ismember(behavior.Type2, {'Scratch', 'Scratching'})) = ...
    {'Scratching'};
behavior.Type2(cellfun(@isempty, behavior.Type2)) = {'Other'};

if any(~isfinite(behavior.Start)) || any(~isfinite(behavior.End)) || ...
        any(behavior.End < behavior.Start)
    error('GENERATE_BEHAVIOR_PROBABILITY:Times', ...
        '%s contains invalid behavior intervals.', recording_name);
end

session_mask = ismember(behavior.Animal, 'Session_structure') & ...
    ismember(behavior.Type, 'Partners session');
session_rows = find(session_mask);
session_ranges = behavior{session_mask, {'Start', 'End'}};
if isempty(session_ranges)
    error('GENERATE_BEHAVIOR_PROBABILITY:Sessions', ...
        '%s has no Session_structure / Partners session rows.', ...
        recording_name);
end
[session_ranges, session_order] = sortrows(session_ranges, 1);
session_rows = session_rows(session_order);
if any(session_ranges(:, 2) <= session_ranges(:, 1))
    error('GENERATE_BEHAVIOR_PROBABILITY:Sessions', ...
        '%s contains a partner-play session with non-positive duration.', ...
        recording_name);
end
if any(session_ranges(2:end, 1) < session_ranges(1:end-1, 2))
    error('GENERATE_BEHAVIOR_PROBABILITY:Sessions', ...
        '%s contains overlapping partner-play sessions.', recording_name);
end

valid_event = ~ismember(behavior.Animal, ...
    {'Session_structure', 'Reversal'}) & ...
    ismember(behavior.Type2, options.all_behaviors);
n_behaviors = numel(options.all_behaviors);
n_time = numel(options.time);
n_pre = numel(options.pre_lag);
n_post = numel(options.post_lag);
n_wrapped = n_pre + options.wrapped_bins + n_post;

empty_session = struct( ...
    'play_session_number', [], ...
    'play_session', [], ...
    'play_session_table_row', [], ...
    'participant_animals', {{}}, ...
    'self_behavior', false(0, n_behaviors, n_time), ...
    'partner_behavior', false(0, n_behaviors, n_time), ...
    'self_behavior_wrapped', false(0, n_behaviors, n_wrapped), ...
    'partner_behavior_wrapped', false(0, n_behaviors, n_wrapped), ...
    'wrapped_in_session', false(0, n_wrapped), ...
    'target_behavior', {{}}, ...
    'target_actor', {{}}, ...
    'target_actor_is_focal', false(0, 1), ...
    'target_start', zeros(0, 1), ...
    'target_offset', zeros(0, 1), ...
    'partner_actors', {{}});
sessions = repmat(empty_session, size(session_ranges, 1), 1);

for session_n = 1:size(session_ranges, 1)
    session_start = session_ranges(session_n, 1);
    session_end = session_ranges(session_n, 2);
    in_session = valid_event & behavior.End > session_start;
    if session_n < size(session_ranges, 1)
        in_session = in_session & behavior.Start < session_end;
    else
        in_session = in_session & behavior.Start <= session_end;
    end
    participants = unique([{focal_animal}; ...
        behavior.Animal(in_session)], 'stable');
    participants(cellfun(@isempty, participants)) = [];

    target_mask = valid_event & ...
        ismember(behavior.Type2, options.target_behaviors) & ...
        ismember(behavior.Animal, participants) & ...
        behavior.End > session_start & behavior.End <= session_end;
    target_rows = find(target_mask);
    [~, target_order] = sortrows( ...
        [behavior.End(target_rows), behavior.Start(target_rows)], [1 2]);
    target_rows = target_rows(target_order);

    n_targets = numel(target_rows);
    self_behavior = false(n_targets, n_behaviors, n_time);
    partner_behavior = false(n_targets, n_behaviors, n_time);
    self_behavior_wrapped = false(n_targets, n_behaviors, n_wrapped);
    partner_behavior_wrapped = false(n_targets, n_behaviors, n_wrapped);
    wrapped_in_session = false(n_targets, n_wrapped);
    partner_actors = cell(n_targets, 1);

    for target_n = 1:n_targets
        target_row = target_rows(target_n);
        target_actor = behavior.Animal{target_row};
        other_actors = participants(~ismember(participants, target_actor));
        partner_actors{target_n} = other_actors(:).';
        onset = behavior.Start(target_row);
        offset = behavior.End(target_row);
        absolute_time = offset + options.time;
        sampled_in_session = absolute_time >= session_start & ...
            absolute_time <= session_end;
        wrapped_time = [onset + options.pre_lag, ...
            onset + (offset - onset) * options.wrap_fraction, ...
            offset + options.post_lag];
        wrapped_in_session(target_n, :) = ...
            wrapped_time >= session_start & wrapped_time <= session_end;

        for behavior_n = 1:n_behaviors
            type_mask = ismember(behavior.Type2, ...
                options.all_behaviors{behavior_n});
            self_rows = find(in_session & type_mask & ...
                ismember(behavior.Animal, target_actor));
            other_rows = find(in_session & type_mask & ...
                ismember(behavior.Animal, other_actors));
            self_occupancy = interval_occupancy( ...
                absolute_time, behavior, self_rows) & sampled_in_session;
            partner_occupancy = interval_occupancy( ...
                absolute_time, behavior, other_rows) & sampled_in_session;
            self_behavior(target_n, behavior_n, :) = ...
                reshape(self_occupancy, 1, 1, n_time);
            partner_behavior(target_n, behavior_n, :) = ...
                reshape(partner_occupancy, 1, 1, n_time);
            self_wrapped = interval_occupancy( ...
                wrapped_time, behavior, self_rows);
            partner_wrapped = interval_occupancy( ...
                wrapped_time, behavior, other_rows);
            self_behavior_wrapped(target_n, behavior_n, :) = ...
                reshape(self_wrapped, 1, 1, n_wrapped);
            partner_behavior_wrapped(target_n, behavior_n, :) = ...
                reshape(partner_wrapped, 1, 1, n_wrapped);
        end
    end

    sessions(session_n).play_session_number = session_n;
    sessions(session_n).play_session = session_ranges(session_n, :);
    sessions(session_n).play_session_table_row = session_rows(session_n);
    sessions(session_n).participant_animals = participants(:).';
    sessions(session_n).self_behavior = self_behavior;
    sessions(session_n).partner_behavior = partner_behavior;
    sessions(session_n).self_behavior_wrapped = self_behavior_wrapped;
    sessions(session_n).partner_behavior_wrapped = partner_behavior_wrapped;
    sessions(session_n).wrapped_in_session = wrapped_in_session;
    sessions(session_n).target_behavior = behavior.Type2(target_rows);
    sessions(session_n).target_actor = behavior.Animal(target_rows);
    sessions(session_n).target_actor_is_focal = ...
        ismember(behavior.Animal(target_rows), focal_animal);
    sessions(session_n).target_start = behavior.Start(target_rows);
    sessions(session_n).target_offset = behavior.End(target_rows);
    sessions(session_n).partner_actors = partner_actors;
end

recording = struct();
recording.recording_name = recording_name;
recording.focal_animal = focal_animal;
recording.target_behaviors = options.target_behaviors;
recording.behavior_types = options.all_behaviors;
recording.time = options.time;
recording.time_window = options.time_window;
recording.bin_size = options.bin_size;
recording.flank_seconds = options.flank_seconds;
recording.wrapped_bins = options.wrapped_bins;
recording.pre_lag = options.pre_lag;
recording.post_lag = options.post_lag;
recording.wrap_fraction = options.wrap_fraction;
recording.wrapped_plot_time = options.wrapped_plot_time;
recording.pre_index = 1:n_pre;
recording.wrap_index = n_pre + (1:options.wrapped_bins);
recording.post_index = n_pre + options.wrapped_bins + (1:n_post);
recording.alignment = options.alignment;
recording.self_definition = options.self_definition;
recording.partner_definition = options.partner_definition;
recording.sessions = sessions;
end

function occupancy = interval_occupancy(time, behavior, rows)
if isempty(rows)
    occupancy = false(1, numel(time));
    return
end
starts = behavior.Start(rows);
ends = behavior.End(rows);
occupancy = any(bsxfun(@ge, time, starts) & ...
    bsxfun(@le, time, ends), 1);
end

function [behavior, recording_name] = load_behavior(source)
if istable(source)
    behavior = source;
    recording_name = 'table_input';
elseif ischar(source) || (isstring(source) && isscalar(source))
    source = char(source);
    if ~isfile(source)
        error('GENERATE_BEHAVIOR_PROBABILITY:File', ...
            'Behavior file not found: %s', source);
    end
    behavior = readtable(source);
    [~, recording_name] = fileparts(source);
else
    error('GENERATE_BEHAVIOR_PROBABILITY:Input', ...
        'behavior_source must be a table or behavior-file path.');
end

required = {'Animal', 'Start', 'End', 'Length', 'Type'};
if ~all(ismember(required, behavior.Properties.VariableNames))
    if width(behavior) == 6
        behavior(:, 2) = [];
    end
    if width(behavior) ~= 5
        error('GENERATE_BEHAVIOR_PROBABILITY:Columns', ...
            'Behavior input must contain Animal, Start, End, Length, and Type.');
    end
    behavior.Properties.VariableNames = required;
end
behavior = behavior(:, required);
end

function focal_animal = infer_focal_animal(recording_name)
parts = strsplit(recording_name, ' ');
if numel(parts) < 3 || isempty(parts{3})
    error('GENERATE_BEHAVIOR_PROBABILITY:FocalAnimal', ...
        ['Could not infer the focal animal from "%s". Pass focal_animal ' ...
        'explicitly.'], recording_name);
end
focal_animal = parts{3};
end

function values = text_column(values)
values = cellstr(string(values));
values = strtrim(values);
missing = ismissing(string(values));
values(missing) = {''};
end
