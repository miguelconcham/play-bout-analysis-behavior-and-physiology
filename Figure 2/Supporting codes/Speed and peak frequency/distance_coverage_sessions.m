function sessions = distance_coverage_sessions(psth_structure, animal_names)
%DISTANCE_COVERAGE_SESSIONS  Electrode-1 distance, power, and play mask.
%
% Distance is the saved within-session z-score of inter-animal separation.
% Positive values are farther than that session's own mean. Speed is ignored.

if size(animal_names, 1) ~= numel(psth_structure)
    error('distance_coverage_sessions:NameMismatch', ...
        'animal_names rows (%d) do not match psth_structure (%d).', ...
        size(animal_names, 1), numel(psth_structure));
end

electrode_id = cellfun(@double, animal_names(:, 3));
session_rows = find(electrode_id == 1);
if isempty(session_rows)
    error('distance_coverage_sessions:NoElectrode', ...
        'animal_names has no electrode-1 rows.');
end

n_sess = numel(session_rows);
sessions = repmat(struct( ...
    'animal', '', ...
    'partner', NaN, ...
    'distance_z', [], ...
    'power', [], ...
    'is_play', [], ...
    'hop_s', NaN, ...
    'wind_s', NaN), n_sess, 1);

for s = 1:n_sess
    row = session_rows(s);
    entry = psth_structure(row);
    mdl = entry.model_data;
    [distance_z, found_distance] = model_column(mdl, 'Distance');
    [power, found_power] = model_column(mdl, 'Power');
    [play_col, found_play] = model_column(mdl, 'Play');
    if ~found_distance
        error('distance_coverage_sessions:NoDistance', ...
            'Session %d has no Distance column.', row);
    end
    if ~found_power
        error('distance_coverage_sessions:NoPower', ...
            'Session %d has no Power column.', row);
    end
    if ~found_play
        error('distance_coverage_sessions:NoPlay', ...
            'Session %d has no Play column.', row);
    end

    distance_z = double(distance_z(:));
    power = double(power(:));
    is_play = play_indicator(play_col);
    if numel(distance_z) ~= numel(power) || numel(is_play) ~= numel(distance_z)
        error('distance_coverage_sessions:LengthMismatch', ...
            'Distance, Power, and Play differ in length in session %d.', row);
    end

    [hop_s, wind_s] = frame_timing(entry);
    name = animal_names{row, 1};
    if exist('isstring', 'builtin') == 5 && isstring(name)
        name = char(name);
    end
    sessions(s).animal = char(name);
    sessions(s).partner = double(animal_names{row, 2});
    sessions(s).distance_z = distance_z;
    sessions(s).power = power;
    sessions(s).is_play = is_play;
    sessions(s).hop_s = hop_s;
    sessions(s).wind_s = wind_s;
end
end

function [values, found] = model_column(mdl, name)
if isstruct(mdl)
    found = isfield(mdl, name);
else
    found = any(strcmp(mdl.Properties.VariableNames, name));
end
values = [];
if found
    values = mdl.(name);
end
end

function is_play = play_indicator(play_col)
play_col = play_col(:);
if isempty(play_col)
    is_play = false(0, 1);
    return
end
if islogical(play_col)
    is_play = play_col;
    return
end
if isnumeric(play_col)
    is_play = play_col ~= 0 & isfinite(play_col);
    return
end
if exist('isstring', 'builtin') == 5 && isstring(play_col)
    play_col = cellstr(play_col);
end
if iscell(play_col)
    tokens = play_col;
else
    tokens = cellstr(play_col);
end
is_play = false(numel(tokens), 1);
for k = 1:numel(tokens)
    token = lower(strtrim(char(tokens{k})));
    is_play(k) = strcmp(token, 'true') || strcmp(token, '1');
end
end

function [hop_s, wind_s] = frame_timing(entry)
% Same hop as Analyze_speed_freq_relation: spectrogram frames at 2500 Hz.
fs = 2500;
wind_s = 0.25;
if isfield(entry, 'wind_length') && ~isempty(entry.wind_length)
    wind_s = entry.wind_length;
end
if isfield(entry, 'wind_length') && isfield(entry, 'wind_overlap') ...
        && ~isempty(entry.wind_length) && ~isempty(entry.wind_overlap)
    hop_samples = entry.wind_length * fs - floor(entry.wind_overlap * fs);
    hop_s = hop_samples / fs;
else
    hop_s = 0.0012;
end
if ~isfinite(hop_s) || hop_s <= 0
    hop_s = 0.0012;
end
end
