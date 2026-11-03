function analysis = summarize_actor_behavior_probability(recordings)
%SUMMARIZE_ACTOR_BEHAVIOR_PROBABILITY Summarize each actor-session stratum.
%
% The primary sampling unit is one target actor in one Partners session.
% Equal-stratum maps therefore prevent an actor with many target events
% from dominating the descriptive average. Event-weighted maps are also
% retained for analyses where event pooling is intentional.

if isempty(recordings)
    error('summarize_actor_behavior_probability:Empty', ...
        'No actor-relative recordings were supplied.');
end

behavior_types = recordings(1).behavior_types;
target_behaviors = recordings(1).target_behaviors;
time = recordings(1).time;
n_behaviors = numel(behavior_types);
n_targets = numel(target_behaviors);
n_time = numel(time);

empty_stratum = struct( ...
    'recording_name', '', ...
    'play_session_number', [], ...
    'play_session', [], ...
    'participant_animals', {{}}, ...
    'target_actor', '', ...
    'other_animals', {{}}, ...
    'target_behavior', '', ...
    'n_target_events', 0, ...
    'self_probability', nan(n_behaviors, n_time), ...
    'other_probability', nan(n_behaviors, n_time));
strata = repmat(empty_stratum, 0, 1);
stratum_self = cell(1, n_targets);
stratum_other = cell(1, n_targets);
event_self = cell(1, n_targets);
event_other = cell(1, n_targets);
for target_n = 1:n_targets
    stratum_self{target_n} = nan(0, n_behaviors, n_time);
    stratum_other{target_n} = nan(0, n_behaviors, n_time);
    event_self{target_n} = false(0, n_behaviors, n_time);
    event_other{target_n} = false(0, n_behaviors, n_time);
end

for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    check_parameters(recording, behavior_types, target_behaviors, time);
    for session_n = 1:numel(recording.sessions)
        session = recording.sessions(session_n);
        for actor_n = 1:numel(session.actor_data)
            actor = session.actor_data(actor_n);
            for target_n = 1:n_targets
                rows = strcmp(actor.target_behavior, ...
                    target_behaviors{target_n});
                n_events = sum(rows);
                summary = empty_stratum;
                summary.recording_name = recording.recording_name;
                summary.play_session_number = session.play_session_number;
                summary.play_session = session.play_session;
                summary.participant_animals = session.participant_animals;
                summary.target_actor = actor.actor;
                summary.other_animals = actor.other_animals;
                summary.target_behavior = target_behaviors{target_n};
                summary.n_target_events = n_events;
                if n_events > 0
                    self_values = actor.self_behavior(rows, :, :);
                    other_values = actor.other_behavior(rows, :, :);
                    summary.self_probability = reshape(mean( ...
                        double(self_values), 1), n_behaviors, n_time);
                    summary.other_probability = reshape(mean( ...
                        double(other_values), 1), n_behaviors, n_time);
                    event_self{target_n}(end + (1:n_events), :, :) = ...
                        self_values; %#ok<AGROW>
                    event_other{target_n}(end + (1:n_events), :, :) = ...
                        other_values; %#ok<AGROW>
                end
                strata(end + 1, 1) = summary; %#ok<AGROW>
                stratum_self{target_n}(end + 1, :, :) = ...
                    summary.self_probability; %#ok<AGROW>
                stratum_other{target_n}(end + 1, :, :) = ...
                    summary.other_probability; %#ok<AGROW>
            end
        end
    end
end

analysis = struct();
analysis.target_behaviors = target_behaviors;
analysis.behavior_types = behavior_types;
analysis.time = time;
analysis.bin_size = recordings(1).bin_size;
analysis.alignment = recordings(1).alignment;
analysis.self_definition = recordings(1).self_definition;
analysis.other_definition = recordings(1).other_definition;
analysis.strata = strata;
analysis.mean_stratum_self = nan(n_behaviors, n_targets, n_time);
analysis.mean_stratum_other = nan(n_behaviors, n_targets, n_time);
analysis.sem_stratum_self = nan(n_behaviors, n_targets, n_time);
analysis.sem_stratum_other = nan(n_behaviors, n_targets, n_time);
analysis.event_weighted_self = nan(n_behaviors, n_targets, n_time);
analysis.event_weighted_other = nan(n_behaviors, n_targets, n_time);
analysis.n_strata_with_target = zeros(1, n_targets);
analysis.n_events = zeros(1, n_targets);

for target_n = 1:n_targets
    analysis.mean_stratum_self(:, target_n, :) = reshape( ...
        mean(stratum_self{target_n}, 1, 'omitmissing'), ...
        n_behaviors, 1, n_time);
    analysis.mean_stratum_other(:, target_n, :) = reshape( ...
        mean(stratum_other{target_n}, 1, 'omitmissing'), ...
        n_behaviors, 1, n_time);
    analysis.sem_stratum_self(:, target_n, :) = reshape( ...
        sem(stratum_self{target_n}), n_behaviors, 1, n_time);
    analysis.sem_stratum_other(:, target_n, :) = reshape( ...
        sem(stratum_other{target_n}), n_behaviors, 1, n_time);
    analysis.event_weighted_self(:, target_n, :) = reshape( ...
        mean(double(event_self{target_n}), 1, 'omitmissing'), ...
        n_behaviors, 1, n_time);
    analysis.event_weighted_other(:, target_n, :) = reshape( ...
        mean(double(event_other{target_n}), 1, 'omitmissing'), ...
        n_behaviors, 1, n_time);
    analysis.n_strata_with_target(target_n) = ...
        sum([strata.n_target_events] > 0 & ...
        strcmp({strata.target_behavior}, target_behaviors{target_n}));
    analysis.n_events(target_n) = size(event_self{target_n}, 1);
end
analysis.n_actor_session_strata = numel(strata) / n_targets;
analysis.total_target_events = sum(analysis.n_events);
end

function output = sem(values)
if isempty(values)
    output = nan(size(values, 2), size(values, 3));
    return
end
n = sum(~isnan(values), 1);
output = reshape(std(values, 0, 1, 'omitmissing') ./ sqrt(n), ...
    size(values, 2), size(values, 3));
output(reshape(n, size(values, 2), size(values, 3)) < 2) = NaN;
end

function check_parameters(recording, behaviors, targets, time)
if ~isequal(recording.behavior_types, behaviors) || ...
        ~isequal(recording.target_behaviors, targets) || ...
        ~isequal(recording.time, time)
    error('summarize_actor_behavior_probability:Parameters', ...
        'All recordings must share behavior, target, and time settings.');
end
end
