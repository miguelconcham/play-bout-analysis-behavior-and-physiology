function catalog = theta_envelope_behavior_events(behavior, repeated_animal, play_behaviors, other_behaviors)
%THETA_ENVELOPE_BEHAVIOR_EVENTS  One onset list per behavior and per actor.
%
% Alignment is the start of each annotated behavior. Play types are not
% pooled, and a play bout is not an event. self_onsets belong to the
% implanted animal. partner_onsets belong to every other animal in the
% file except Reversal and Session_structure.
%
% Type labels follow GENERATE_PSTH_ALL_BEHAVIORS: Pounce_A/B become Pounce,
% Pounce_Ai/Bi become PounceI, and Scratch becomes Scratching.

[animal, start_s, type_label] = behavior_columns(behavior);
type2 = normalize_type(type_label);
usable = ~strcmp(animal, 'Reversal') & ~strcmp(animal, 'Session_structure');
usable = usable & isfinite(start_s);
names = [play_behaviors(:).', other_behaviors(:).'];
is_play = [true(1, numel(play_behaviors)), false(1, numel(other_behaviors))];
catalog = repmat(empty_entry(), numel(names), 1);
self = usable & strcmp(animal, repeated_animal);
partner = usable & ~strcmp(animal, repeated_animal);
for behavior_n = 1:numel(names)
    match = strcmp(type2, names{behavior_n});
    catalog(behavior_n).name = names{behavior_n};
    catalog(behavior_n).is_play = is_play(behavior_n);
    catalog(behavior_n).self_onsets = start_s(self & match);
    catalog(behavior_n).partner_onsets = start_s(partner & match);
end
end

function entry = empty_entry()
entry = struct('name', '', 'is_play', false, ...
    'self_onsets', zeros(0, 1), 'partner_onsets', zeros(0, 1));
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
type2(ismember(type2, {'Scratch', 'Scratching'})) = {'Scratching'};
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
