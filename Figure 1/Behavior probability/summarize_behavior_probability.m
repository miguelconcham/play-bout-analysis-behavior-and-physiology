function analysis = summarize_behavior_probability(recordings)
%SUMMARIZE_BEHAVIOR_PROBABILITY Estimate probabilities within each session.
%
% The primary output is analysis.sessions: one entry per retained play
% session. Grand means are also returned in two forms. mean_session_* gives
% every session equal weight; event_weighted_* gives every target event
% equal weight.

if isempty(recordings)
    error('summarize_behavior_probability:Empty', ...
        'No behavior-probability recordings were supplied.');
end

behavior_types = recordings(1).behavior_types;
time = recordings(1).time;
n_behaviors = numel(behavior_types);
n_time = numel(time);

empty_summary = struct( ...
    'recording_name', '', ...
    'focal_animal', '', ...
    'play_session_number', [], ...
    'play_session', [], ...
    'participant_animals', {{}}, ...
    'n_target_events', 0, ...
    'target_behavior', {{}}, ...
    'target_actor', {{}}, ...
    'self_probability', nan(n_behaviors, n_time), ...
    'partner_probability', nan(n_behaviors, n_time));
session_summaries = repmat(empty_summary, 0, 1);
session_self = nan(0, n_behaviors, n_time);
session_partner = nan(0, n_behaviors, n_time);
all_self = false(0, n_behaviors, n_time);
all_partner = false(0, n_behaviors, n_time);

for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    if ~isequal(recording.behavior_types, behavior_types) || ...
            ~isequal(recording.time, time)
        error('summarize_behavior_probability:Parameters', ...
            'All recordings must use the same behavior list and time bins.');
    end
    for session_n = 1:numel(recording.sessions)
        source = recording.sessions(session_n);
        n_events = size(source.self_behavior, 1);
        summary = empty_summary;
        summary.recording_name = recording.recording_name;
        summary.focal_animal = recording.focal_animal;
        summary.play_session_number = source.play_session_number;
        summary.play_session = source.play_session;
        summary.participant_animals = source.participant_animals;
        summary.n_target_events = n_events;
        summary.target_behavior = source.target_behavior;
        summary.target_actor = source.target_actor;
        if n_events > 0
            summary.self_probability = reshape( ...
                mean(double(source.self_behavior), 1), ...
                n_behaviors, n_time);
            summary.partner_probability = reshape( ...
                mean(double(source.partner_behavior), 1), ...
                n_behaviors, n_time);
            all_self(end + (1:n_events), :, :) = ...
                source.self_behavior; %#ok<AGROW>
            all_partner(end + (1:n_events), :, :) = ...
                source.partner_behavior; %#ok<AGROW>
        end
        session_summaries(end + 1, 1) = summary; %#ok<AGROW>
        session_self(end + 1, :, :) = summary.self_probability; %#ok<AGROW>
        session_partner(end + 1, :, :) = ...
            summary.partner_probability; %#ok<AGROW>
    end
end

analysis = struct();
analysis.target_behaviors = recordings(1).target_behaviors;
analysis.behavior_types = behavior_types;
analysis.time = time;
analysis.bin_size = recordings(1).bin_size;
analysis.alignment = recordings(1).alignment;
analysis.self_definition = recordings(1).self_definition;
analysis.partner_definition = recordings(1).partner_definition;
analysis.sessions = session_summaries;
analysis.n_sessions_total = numel(session_summaries);
analysis.n_sessions_with_targets = ...
    sum([session_summaries.n_target_events] > 0);
analysis.total_target_events = size(all_self, 1);
analysis.mean_session_self = collapse_mean(session_self, ...
    n_behaviors, n_time);
analysis.mean_session_partner = collapse_mean(session_partner, ...
    n_behaviors, n_time);
analysis.sem_session_self = collapse_sem(session_self, ...
    n_behaviors, n_time);
analysis.sem_session_partner = collapse_sem(session_partner, ...
    n_behaviors, n_time);
analysis.event_weighted_self = collapse_mean(all_self, ...
    n_behaviors, n_time);
analysis.event_weighted_partner = collapse_mean(all_partner, ...
    n_behaviors, n_time);
end

function output = collapse_mean(values, n_behaviors, n_time)
if isempty(values)
    output = nan(n_behaviors, n_time);
else
    output = reshape(mean(double(values), 1, 'omitmissing'), ...
        n_behaviors, n_time);
end
end

function output = collapse_sem(values, n_behaviors, n_time)
if isempty(values)
    output = nan(n_behaviors, n_time);
    return
end
values = double(values);
n = sum(~isnan(values), 1);
sd = std(values, 0, 1, 'omitmissing');
output = reshape(sd ./ sqrt(n), n_behaviors, n_time);
output(reshape(n, n_behaviors, n_time) < 2) = NaN;
end
