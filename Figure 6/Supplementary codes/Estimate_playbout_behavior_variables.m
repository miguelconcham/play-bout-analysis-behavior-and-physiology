%% Estimate_playbout_behavior_variables
% For each TWO ANIMALS session:
%   1) Load ELAN behavior, tracking, USV Stats, video2NPX + audio2NPX
%   2) Map tracking animal_pos = focal (focal txt); partner_pos = other
%   3) Transform behavior/tracking -> NPX (video2NPX), calls -> NPX (audio2NPX)
%   4) Estimate play bouts with canonical play_bout() (same params as MI GENERATE)
%   5) Build peri-onset matrices (-20..+20 s, bin 10 ms) for each animal:
%        vars: speed, acceleration, inter-animal distance, vocalization
%        size: NP x (NVAR+1) x T
%      Distance and vocalization are identical in both animals' matrices.
%   6) Save one .mat per session with clear animal-code field names.
%
% Plot the saved matrices across animals with
% Plot_playbout_behavior_variable_stacks.m.
% Does not delete/replace Test_* plot scripts.

%% Paths / options
experiment_folder = '\\experimentfs.home.bccn-berlin.de\experiment';
repo_root = [experiment_folder, ...
    '\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology'];
datasets_two = [experiment_folder, ...
    '\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\TWO ANIMALS'];
saving_folder = fullfile(repo_root, 'Data', 'Analysis results', 'playbout behavior regressors two animals');
if exist(saving_folder, 'dir') ~= 7
    mkdir(saving_folder);
end

% Ensure canonical play_bout is on path
addpath(fullfile(repo_root, 'Custom functions'));
addpath(fullfile(repo_root, 'Figure 6'));

% Play-bout params (match GENERATE_STRUCUTRE_animal_synch)
play_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'};
beh_bin        = 0.01;
conv_length    = 1;
behavior_window = 0;

% Peri-onset alignment
hist_range = [-20 20];
bin_size   = 0.01;   % 10 ms
variable_names = {'speed', 'acceleration', 'distance', 'vocalization'};  % NVAR+1 = 4
% indices: 1 speed, 2 accel (per animal); 3 distance, 4 voc (shared)

% Which sessions: [] = all
sessions_to_run = [];

%% Focal map
focal_file = fullfile(datasets_two, 'focal = enters 1st in the arena or.txt');
focal_map = local_load_focal_map(focal_file);

session_dirs = dir(datasets_two);
session_dirs = session_dirs([session_dirs.isdir] & ~ismember({session_dirs.name}, {'.','..'}));

if isempty(sessions_to_run)
    session_idx = 1:numel(session_dirs);
elseif isnumeric(sessions_to_run)
    session_idx = sessions_to_run;
else
    names = {session_dirs.name};
    session_idx = find(ismember(names, sessions_to_run));
end

psth_time = hist_range(1):bin_size:hist_range(2);
T = numel(psth_time);
n_var = numel(variable_names);

all_session_names = {};

for si = session_idx
    session_name = session_dirs(si).name;
    session_path = fullfile(datasets_two, session_name);
    fprintf('\n==== %s ====\n', session_name);

    tokens = strsplit(strtrim(session_name), ' ');
    if numel(tokens) < 2
        warning('Cannot parse animal names from folder "%s". Skipping.', session_name);
        continue
    end
    animal_A = char(tokens{1});  % folder / MI animal 1
    animal_B = char(tokens{2});  % folder / MI animal 2

    focal_name = local_resolve_focal(focal_map, session_name);
    if isempty(focal_name)
        warning('No focal entry for "%s". Skipping.', session_name);
        continue
    end
    if strcmp(focal_name, animal_A)
        other_name = animal_B;
    elseif strcmp(focal_name, animal_B)
        other_name = animal_A;
    else
        warning('Focal "%s" not in {%s,%s}. Skipping.', focal_name, animal_A, animal_B);
        continue
    end
    fprintf('folder animals: %s, %s | tracking animal_pos=%s | partner_pos=%s\n', ...
        animal_A, animal_B, focal_name, other_name);

    %% Sync
    video_sync_file = fullfile(session_path, 'synch_model_video2NPX.mat');
    audio_sync_file = fullfile(session_path, 'synch_model_audio2NPX.mat');
    if exist(video_sync_file, 'file') ~= 2 || exist(audio_sync_file, 'file') ~= 2
        warning('Missing sync models in %s. Skipping.', session_name);
        continue
    end
    V = load(video_sync_file);
    A = load(audio_sync_file);
    synch_model_video2NPX = V.synch_model_video2NPX;
    synch_model_audio2NPX = A.synch_model_audio2NPX;

    %% Behavior -> NPX, then play bouts
    Behavior = local_load_behavior(session_path);
    Behavior.Start = predict(synch_model_video2NPX, Behavior.Start);
    Behavior.End   = predict(synch_model_video2NPX, Behavior.End);

    animal_types = unique(Behavior.Animal);
    animal_types = animal_types(~ismember(animal_types, {'Session_structure','Reversal',''}));
    animal_types = cellstr(animal_types(:));

    % repeated_animal = folder animal 1 (same convention as MI GENERATE)
    config = struct();
    config.Behavior         = Behavior;
    config.repeated_animal  = animal_A;
    config.animal_types     = animal_types;
    config.play_behaviors   = play_behaviors;
    config.beh_bin          = beh_bin;
    config.conv_length      = conv_length;
    config.behavior_window  = behavior_window;

    [play_bouts_table, play_bout_behaviors] = play_bout(config);
    if isempty(play_bouts_table)
        warning('play_bout returned empty for %s. Skipping.', session_name);
        continue
    end
    NP = size(play_bouts_table, 1);
    fprintf('Play bouts: %d\n', NP);

    %% Tracking -> NPX kinematics
    tr_file = dir(fullfile(session_path, '*traking_structure*.mat'));
    if isempty(tr_file)
        warning('No tracking in %s. Skipping.', session_name);
        continue
    end
    S = load(fullfile(session_path, tr_file(1).name));
    if isfield(S, 'traking_structure')
        TS = S.traking_structure;
    else
        ff = fieldnames(S);
        TS = S.(ff{1});
    end

    if isfield(TS, 'time') && ~isempty(TS.time)
        track_time_video = TS.time(:)';
    elseif isfield(TS, 'frames2stract') && ~isempty(TS.frames2stract)
        warning('%s: TS.time missing; fallback frames2stract/30', session_name);
        track_time_video = (TS.frames2stract(:))' / 30;
    else
        error('No usable tracking time in %s', tr_file(1).name);
    end
    if numel(track_time_video) ~= size(TS.animal_pos, 1)
        error('track time length mismatch in %s', tr_file(1).name);
    end
    track_time = predict(synch_model_video2NPX, track_time_video(:))';

    % Identify tracking slots with folder animal codes (consistent with focal)
    % animal_pos = focal, partner_pos = other
    pos = struct();
    pos.(focal_name) = TS.animal_pos;
    pos.(other_name) = TS.partner_pos;

    dt = diff(track_time);
    dt(dt <= 0) = NaN;
    mid_time = track_time(1:end-1) + dt(:)' / 2;

    speed = struct();
    accel = struct();
    for an = {animal_A, animal_B}
        name = an{1};
        dpos = diff(pos.(name));
        sp = sqrt(sum(dpos.^2, 2)) ./ dt(:);
        % acceleration magnitude of speed change (px/s^2)
        dt_sp = diff(mid_time);
        dt_sp(dt_sp <= 0) = NaN;
        ac = abs(diff(sp)) ./ dt_sp(:);
        ac_time = mid_time(1:end-1) + dt_sp(:)' / 2;
        % Interpolate accel onto speed mid_time for a shared kinematic clock
        ac_on_speed_time = interp1(ac_time, ac, mid_time, 'linear', NaN);
        speed.(name) = sp;
        accel.(name) = ac_on_speed_time;
    end

    % Shared inter-animal distance on track samples, then midpoints for interp
    dist_track = sqrt(sum((TS.animal_pos - TS.partner_pos).^2, 2));
    dist_mid = interp1(track_time, dist_track, mid_time, 'linear', NaN);

    %% Vocalizations -> NPX
    CallStats = local_load_calls(session_path);
    if isempty(CallStats)
        warning('No calls in %s; vocalization channel = 0.', session_name);
        call_begin = [];
        call_end = [];
    else
        call_begin = predict(synch_model_audio2NPX, CallStats.BeginTimes);
        call_end   = predict(synch_model_audio2NPX, CallStats.EndTimes);
        fprintf('Calls: %d\n', numel(call_begin));
    end

    %% Build NP x (NVAR+1) x T matrices per animal
    vars_A = nan(NP, n_var, T);
    vars_B = nan(NP, n_var, T);

    for pb = 1:NP
        onset = play_bouts_table(pb, 1);
        t_abs = psth_time + onset;

        % Per-animal speed / accel
        vars_A(pb, 1, :) = interp1(mid_time, speed.(animal_A), t_abs, 'linear', NaN);
        vars_A(pb, 2, :) = interp1(mid_time, accel.(animal_A), t_abs, 'linear', NaN);
        vars_B(pb, 1, :) = interp1(mid_time, speed.(animal_B), t_abs, 'linear', NaN);
        vars_B(pb, 2, :) = interp1(mid_time, accel.(animal_B), t_abs, 'linear', NaN);

        % Shared distance
        dist_trace = interp1(mid_time, dist_mid, t_abs, 'linear', NaN);
        vars_A(pb, 3, :) = dist_trace;
        vars_B(pb, 3, :) = dist_trace;

        % Shared vocalization occupancy
        if isempty(call_begin)
            voc = zeros(1, T);
        else
            voc = zeros(1, T);
            for c = 1:numel(call_begin)
                voc = voc | (t_abs >= call_begin(c) & t_abs <= call_end(c));
            end
            voc = double(voc);
        end
        vars_A(pb, 4, :) = voc;
        vars_B(pb, 4, :) = voc;
    end

    %% Package + save (animal-code field names)
    playbout_behavior = struct();
    playbout_behavior.session_name       = session_name;
    playbout_behavior.animal_names       = {animal_A, animal_B};  % folder order
    playbout_behavior.tracking_animal_pos = focal_name;           % = TS animal_pos identity
    playbout_behavior.tracking_partner_pos = other_name;          % = TS partner_pos identity
    playbout_behavior.play_behaviors     = play_behaviors;
    playbout_behavior.play_bouts_table   = play_bouts_table;      % NP x 2 [onset offset] NPX
    playbout_behavior.psth_time          = psth_time;             % 1 x T
    playbout_behavior.hist_range         = hist_range;
    playbout_behavior.bin_size           = bin_size;
    playbout_behavior.variable_names     = variable_names;        % 1 x (NVAR+1)
    playbout_behavior.play_bout_behaviors = play_bout_behaviors;

    % Matrices named by animal code for easy identification
    playbout_behavior.(['vars_' local_safe_field(animal_A)]) = vars_A;  % NP x 4 x T
    playbout_behavior.(['vars_' local_safe_field(animal_B)]) = vars_B;

    % Also store ordered cell for programmatic access
    playbout_behavior.vars_by_animal = {vars_A, vars_B};

    out_name = sprintf('playbout_behavior_vars_%s.mat', local_safe_filename(session_name));
    out_path = fullfile(saving_folder, out_name);
    save(out_path, 'playbout_behavior', '-v7.3');
    fprintf('Saved %s (NP=%d, T=%d, n_var=%d)\n', out_name, NP, T, n_var);

    all_session_names{end+1} = session_name; %#ok<AGROW>
end

% Index file of sessions processed
session_list = all_session_names;
save(fullfile(saving_folder, 'playbout_behavior_vars_session_list.mat'), 'session_list');
fprintf('\nDone. Outputs in:\n%s\n', saving_folder);

%% ========================= local helpers =========================

function focal_map = local_load_focal_map(focal_file)
focal_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
if exist(focal_file, 'file') ~= 2
    error('Focal file not found: %s', focal_file);
end
lines = readlines(focal_file);
for i = 1:numel(lines)
    line = strtrim(lines(i));
    if line == "" || ~contains(line, ',')
        continue
    end
    parts = split(line, ',');
    if numel(parts) < 2
        continue
    end
    sess = char(strtrim(parts(1)));
    rest = char(strtrim(parts(2)));
    tok = regexp(rest, '^(\S+)\s+is\s+focal', 'tokens', 'once', 'ignorecase');
    if isempty(tok)
        continue
    end
    focal_map(sess) = tok{1};
end
end

function focal_name = local_resolve_focal(focal_map, session_name)
focal_name = '';
if isKey(focal_map, session_name)
    focal_name = focal_map(session_name);
    return
end
key_alt = regexprep(session_name, '\s+', ' ');
ks = keys(focal_map);
for i = 1:numel(ks)
    if strcmp(regexprep(ks{i}, '\s+', ' '), key_alt)
        focal_name = focal_map(ks{i});
        return
    end
end
end

function Behavior = local_load_behavior(session_path)
elan = dir(fullfile(session_path, 'ELAN*'));
if isempty(elan)
    error('No ELAN* behavior file in %s', session_path);
end
raw = readtable(fullfile(session_path, elan(1).name), ...
    'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);
if width(raw) >= 6
    raw(:, 2) = [];
end
if width(raw) < 5
    error('Unexpected ELAN format in %s (%d columns)', elan(1).name, width(raw));
end
raw = raw(:, 1:5);
raw.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
Behavior = raw;
Behavior.Animal = string(Behavior.Animal);
Behavior.Type = string(Behavior.Type);
Behavior(Behavior.Animal == "Reversal", :) = [];
Behavior.Type2 = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A','Pounce_B'})) = "Pounce";
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai','Pounce_Bi'})) = "PounceI";
Behavior.Type2(Behavior.Type2 == "") = "Other";
% play_bout expects cellstr-like Animal for ismember with animal_types cell
Behavior.Animal = cellstr(Behavior.Animal);
Behavior.Type = cellstr(Behavior.Type);
Behavior.Type2 = cellstr(Behavior.Type2);
end

function CallStats = local_load_calls(session_path)
CallStats = [];
xf = dir(fullfile(session_path, '*_Stats.xlsx'));
if isempty(xf)
    xf = dir(fullfile(session_path, '*.xlsx'));
end
if isempty(xf)
    return
end
CallStats = readtable(fullfile(session_path, xf(1).name));
CallStats.Properties.VariableNames = cellfun(@(x) strrep(x, '_', ''), ...
    CallStats.Properties.VariableNames, 'UniformOutput', false);
if ~ismember('BeginTimes', CallStats.Properties.VariableNames) || ...
        ~ismember('EndTimes', CallStats.Properties.VariableNames)
    error('Call table in %s lacks BeginTimes/EndTimes.', xf(1).name);
end
end

function s = local_safe_field(name)
% MATLAB struct field from animal code (e.g. Single2 -> Single2)
s = regexprep(char(name), '[^A-Za-z0-9_]', '_');
if isempty(s)
    s = 'animal';
end
if s(1) >= '0' && s(1) <= '9'
    s = ['A_' s];
end
end

function s = local_safe_filename(name)
s = regexprep(strtrim(char(name)), '\s+', '_');
s = regexprep(s, '[^A-Za-z0-9_\-]', '');
end
