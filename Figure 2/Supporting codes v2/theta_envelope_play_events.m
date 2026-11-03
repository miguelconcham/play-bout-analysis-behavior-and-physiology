function events = theta_envelope_play_events(behavior, repeated_animal, play_behaviors, clean_pre_seconds)
%THETA_ENVELOPE_PLAY_EVENTS  Recorded-animal play onsets on the NPX clock.
%
% Type labels follow GENERATE_PSTH_PLAY_BOUT: Pounce_A/B become Pounce and
% Pounce_Ai/Bi become PounceI. clean_pre is true when no play onset, from
% either animal, falls in the open interval (onset - clean_pre_seconds, onset).

[animal, start_s, type_label] = behavior_columns(behavior);
type2 = normalize_type(type_label);
play = ismember(type2, play_behaviors) & ...
    ~strcmp(animal, 'Reversal') & ~strcmp(animal, 'Session_structure');
play = play & isfinite(start_s);
all_onsets = start_s(play);
recorded = play & strcmp(animal, repeated_animal);
events.onsets = start_s(recorded);
events.behavior_type = type2(recorded);
events.clean_pre = false(numel(events.onsets), 1);
events.row_index = find(recorded);
for event_n = 1:numel(events.onsets)
    onset = events.onsets(event_n);
    earlier = all_onsets > (onset - clean_pre_seconds) & all_onsets < (onset - 1e-6);
    events.clean_pre(event_n) = ~any(earlier);
end
events.onsets = events.onsets(:);
events.behavior_type = events.behavior_type(:);
events.clean_pre = events.clean_pre(:);
events.row_index = events.row_index(:);
end

function [animal, start_s, type_label] = behavior_columns(behavior)
animal = {};
start_s = [];
type_label = {};
if isempty(behavior)
    return
end
if isstruct(behavior)
    animal = as_text(behavior.Animal);
    start_s = double(behavior.Start(:));
    type_label = as_text(behavior.Type);
    return
end
if is_table(behavior)
    animal = as_text(behavior.Animal);
    start_s = double(behavior.Start(:));
    type_label = as_text(behavior.Type);
end
end

function type2 = normalize_type(type_label)
type2 = type_label;
if isempty(type2)
    return
end
type2(ismember(type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
type2(ismember(type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
blank = strcmp(type2, '');
type2(blank) = {'Other'};
end

function text = as_text(value)
if iscell(value)
    text = cell(numel(value), 1);
    for k = 1:numel(value)
        text{k} = char(value{k});
    end
elseif ischar(value)
    text = cellstr(value);
elseif exist('isstring', 'builtin') == 5 && isstring(value)
    text = cellstr(value);
else
    text = cellstr(string(value));
end
text = text(:);
end

function tf = is_table(behavior)
tf = false;
if exist('istable', 'builtin') == 5 || exist('istable', 'file') == 2
    tf = istable(behavior);
end
end
