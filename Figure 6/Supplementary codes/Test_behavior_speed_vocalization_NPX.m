%% Test_behavior_speed_vocalization_NPX
% Two-animal peri-event speed + vocalization, all clocks transformed to NPX.
%
% For each session in DataSets\TWO ANIMALS:
%   1) Read ELAN behavior + tracking + USV Stats.xlsx
%   2) Load synch_model_video2NPX and synch_model_audio2NPX
%   3) Map tracking animal_pos = focal, partner_pos = other (focal txt)
%   4) Transform behavior Start/End and tracking time -> NPX (video2NPX)
%   5) Transform call Begin/End times -> NPX (audio2NPX)
%   6) Align speed of actor and partner, plus call rate, to behavior onsets
%
% Behavior selection:
%   behavior_labels = {'Escape'};           % single behavior
%   behavior_labels = {'Escape','Evasion'}; % any onset in the list (union)
%
% Figure (2 x 3 subplots):
%   Row 1 — when animal A performs selected behavior(s):
%     (1) A speed (actor)  (2) B speed (partner)  (3) vocalization rate
%   Row 2 — when animal B performs selected behavior(s):
%     (4) B speed (actor)  (5) A speed (partner)  (6) vocalization rate

%% Paths / options
experiment_folder = '\\experimentfs.home.bccn-berlin.de\experiment';
datasets_two = [experiment_folder, ...
    '\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\TWO ANIMALS'];

% Behavior label(s). Match Type or Type2 (merged Pounce etc.).
% For a list, every event whose Type/Type2 is in the list is used (onset of any).
behavior_labels = {'Escape', 'CC','CB','CD', 'Pounce', 'Evade'};

% Peri-event window (seconds, NPX time) and plot binning
hist_range = [-1 2];
bin_size   = 1/30;   % ~video frame rate

% Call-rate kernel (s): binary occupancy smoothed with this half-width moving mean
call_smooth_sec = 0.05;

% Which sessions to plot: [] = all, or indices / folder-name cellstr
sessions_to_plot = [];  % e.g. 1, or {'C1 C2 0620'}; [] = all seven

%% Focal map: animal_pos = focal animal
focal_file = fullfile(datasets_two, 'focal = enters 1st in the arena or.txt');
focal_map = local_load_focal_map(focal_file);

session_dirs = dir(datasets_two);
session_dirs = session_dirs([session_dirs.isdir] & ~ismember({session_dirs.name}, {'.','..'}));

if isempty(sessions_to_plot)
    session_idx = 1:numel(session_dirs);
elseif isnumeric(sessions_to_plot)
    session_idx = sessions_to_plot;
else
    names = {session_dirs.name};
    session_idx = find(ismember(names, sessions_to_plot));
end

psth_time = hist_range(1):bin_size:hist_range(2);
call_smooth_bins = max(1, round(call_smooth_sec / bin_size));

for si = session_idx
    session_name = session_dirs(si).name;
    session_path = fullfile(datasets_two, session_name);
    fprintf('\n==== %s ====\n', session_name);

    %% Playing animals from folder tokens (same convention as MI GENERATE)
    tokens = strsplit(strtrim(session_name), ' ');
    if numel(tokens) < 2
        warning('Cannot parse animal names from folder "%s". Skipping.', session_name);
        continue
    end
    animal_A = tokens{1};  % folder / delta animal 1
    animal_B = tokens{2};  % folder / delta animal 2

    %% Focal / tracking identity
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
        warning('Focal "%s" is not one of {%s, %s}. Skipping.', focal_name, animal_A, animal_B);
        continue
    end
    fprintf('Animals: %s, %s | focal(animal_pos)=%s | partner_pos=%s\n', ...
        animal_A, animal_B, focal_name, other_name);

    %% Sync models (required for NPX alignment)
    video_sync_file = fullfile(session_path, 'synch_model_video2NPX.mat');
    audio_sync_file = fullfile(session_path, 'synch_model_audio2NPX.mat');
    if exist(video_sync_file, 'file') ~= 2 || exist(audio_sync_file, 'file') ~= 2
        warning('Missing video2NPX and/or audio2NPX in %s. Skipping.', session_name);
        continue
    end
    V = load(video_sync_file);
    A = load(audio_sync_file);
    synch_model_video2NPX = V.synch_model_video2NPX;
    synch_model_audio2NPX = A.synch_model_audio2NPX;

    %% Load behavior (ELAN) and transform to NPX
    Behavior = local_load_behavior(session_path);
    Behavior.Start = predict(synch_model_video2NPX, Behavior.Start);
    Behavior.End   = predict(synch_model_video2NPX, Behavior.End);

    beh_animals = unique(Behavior.Animal);
    beh_animals = beh_animals(~ismember(beh_animals, {'Session_structure','Reversal',''}));
    fprintf('Behavior actors: %s\n', strjoin(string(beh_animals), ', '));

    is_beh = ismember(Behavior.Type, behavior_labels) | ismember(Behavior.Type2, behavior_labels);
    beh_A = Behavior(is_beh & strcmp(Behavior.Animal, animal_A), :);
    beh_B = Behavior(is_beh & strcmp(Behavior.Animal, animal_B), :);
    fprintf('Behavior (%s): %s=%d, %s=%d\n', strjoin(behavior_labels, '/'), ...
        animal_A, height(beh_A), animal_B, height(beh_B));

    %% Load tracking, transform time to NPX, compute speed
    tr_file = dir(fullfile(session_path, '*traking_structure*.mat'));
    if isempty(tr_file)
        warning('No tracking file in %s. Skipping.', session_name);
        continue
    end
    S = load(fullfile(session_path, tr_file(1).name));
    if isfield(S, 'traking_structure')
        TS = S.traking_structure;
    else
        ff = fieldnames(S);
        TS = S.(ff{1});
    end

    % Use TS.time (same video clock as ELAN), matching Test_darting_speed_alignment.
    % Do NOT prefer frames2stract/30: for C-batch sessions frames/30 is offset from
    % ELAN/TS.time (~2x), so video2NPX keeps that offset and peri-event speeds are all NaN.
    if isfield(TS, 'time') && ~isempty(TS.time)
        track_time_video = TS.time(:)';
    elseif isfield(TS, 'frames2stract') && ~isempty(TS.frames2stract)
        warning('%s: TS.time missing; falling back to frames2stract/30', session_name);
        track_time_video = (TS.frames2stract(:))' / 30;
    else
        error('No usable tracking time in %s', tr_file(1).name);
    end
    if numel(track_time_video) ~= size(TS.animal_pos, 1)
        error('track time length (%d) != animal_pos rows (%d) in %s', ...
            numel(track_time_video), size(TS.animal_pos, 1), tr_file(1).name);
    end
    track_time = predict(synch_model_video2NPX, track_time_video(:))';

    dt = diff(track_time);
    dt(dt <= 0) = NaN;
    focal_speed   = sqrt(sum(diff(TS.animal_pos).^2, 2)) ./ dt(:);
    partner_speed = sqrt(sum(diff(TS.partner_pos).^2, 2)) ./ dt(:);
    speed_time = track_time(1:end-1) + dt(:)' / 2;

    speed = struct();
    speed.(focal_name) = focal_speed;
    speed.(other_name) = partner_speed;

    %% Load vocalizations and transform to NPX
    CallStats = local_load_calls(session_path);
    if isempty(CallStats)
        warning('No *_Stats.xlsx in %s. Vocalization panels will be empty.', session_name);
        call_begin = [];
        call_end = [];
    else
        call_begin = predict(synch_model_audio2NPX, CallStats.BeginTimes);
        call_end   = predict(synch_model_audio2NPX, CallStats.EndTimes);
        fprintf('Calls: %d (audio -> NPX)\n', numel(call_begin));
    end

    %% Peri-event matrices (NPX-aligned)
    spd_A_when_A = local_align_speed(speed.(animal_A), speed_time, beh_A.Start, psth_time);
    spd_B_when_A = local_align_speed(speed.(animal_B), speed_time, beh_A.Start, psth_time);
    spd_B_when_B = local_align_speed(speed.(animal_B), speed_time, beh_B.Start, psth_time);
    spd_A_when_B = local_align_speed(speed.(animal_A), speed_time, beh_B.Start, psth_time);

    voc_when_A = local_align_call_rate(call_begin, call_end, beh_A.Start, psth_time, call_smooth_bins);
    voc_when_B = local_align_call_rate(call_begin, call_end, beh_B.Start, psth_time, call_smooth_bins);

    fprintf('Finite peri-event speed rows: %s=%d/%d, %s=%d/%d | track NPX [%.1f %.1f]\n', ...
        animal_A, sum(any(isfinite(spd_A_when_A), 2)), size(spd_A_when_A, 1), ...
        animal_B, sum(any(isfinite(spd_B_when_B), 2)), size(spd_B_when_B, 1), ...
        min(speed_time), max(speed_time));

    %% Figure: 2 x 3 (speed pair + vocalization per actor)
    fig = figure('Name', ['Behavior speed+USV NPX | ', session_name], 'Color', 'w');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    beh_str = strjoin(behavior_labels, '/');

    local_plot_mean_sem(nexttile, psth_time, spd_A_when_A, ...
        sprintf('When %s %s: %s speed (actor)  n=%d', animal_A, beh_str, animal_A, size(spd_A_when_A, 1)), ...
        'Speed (px/s)');
    local_plot_mean_sem(nexttile, psth_time, spd_B_when_A, ...
        sprintf('When %s %s: %s speed (partner)  n=%d', animal_A, beh_str, animal_B, size(spd_B_when_A, 1)), ...
        'Speed (px/s)');
    local_plot_mean_sem(nexttile, psth_time, voc_when_A, ...
        sprintf('When %s %s: vocalization rate  n=%d', animal_A, beh_str, size(voc_when_A, 1)), ...
        'Call occupancy');

    local_plot_mean_sem(nexttile, psth_time, spd_B_when_B, ...
        sprintf('When %s %s: %s speed (actor)  n=%d', animal_B, beh_str, animal_B, size(spd_B_when_B, 1)), ...
        'Speed (px/s)');
    local_plot_mean_sem(nexttile, psth_time, spd_A_when_B, ...
        sprintf('When %s %s: %s speed (partner)  n=%d', animal_B, beh_str, animal_A, size(spd_A_when_B, 1)), ...
        'Speed (px/s)');
    local_plot_mean_sem(nexttile, psth_time, voc_when_B, ...
        sprintf('When %s %s: vocalization rate  n=%d', animal_B, beh_str, size(voc_when_B, 1)), ...
        'Call occupancy');

    sgtitle(sprintf('%s | behavior={%s} | NPX-aligned | focal=%s', ...
        session_name, beh_str, focal_name), 'Interpreter', 'none');
end

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
% Expect BeginTimes / EndTimes (DeepSqueak style after underscore strip)
if ~ismember('BeginTimes', CallStats.Properties.VariableNames) || ...
        ~ismember('EndTimes', CallStats.Properties.VariableNames)
    error('Call table in %s lacks BeginTimes/EndTimes. Columns: %s', ...
        xf(1).name, strjoin(CallStats.Properties.VariableNames, ', '));
end
end

function M = local_align_speed(speed_trace, speed_time, onsets, psth_time)
n = numel(onsets);
M = nan(n, numel(psth_time));
if n == 0 || isempty(speed_trace)
    return
end
for k = 1:n
    t = psth_time + onsets(k);
    M(k, :) = interp1(speed_time, speed_trace, t, 'linear', NaN);
end
end

function M = local_align_call_rate(call_begin, call_end, onsets, psth_time, smooth_bins)
n = numel(onsets);
M = nan(n, numel(psth_time));
if n == 0
    return
end
if isempty(call_begin)
    M(:) = 0;
    return
end
for k = 1:n
    t = psth_time + onsets(k);
    occ = zeros(1, numel(t));
    for c = 1:numel(call_begin)
        occ = occ | (t >= call_begin(c) & t <= call_end(c));
    end
    if smooth_bins > 1
        occ = movmean(double(occ), smooth_bins);
    else
        occ = double(occ);
    end
    M(k, :) = occ;
end
end

function local_plot_mean_sem(ax, t, M, title_str, ylab)
axes(ax); %#ok<LAXES>
hold on
if isempty(M)
    text(0.5, 0.5, 'no events', 'HorizontalAlignment', 'center');
    axis([t(1) t(end) 0 1]);
elseif all(isnan(M(:)))
    text(0.5, 0.5, sprintf('n=%d but no overlap with tracking', size(M, 1)), ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none');
    axis([t(1) t(end) 0 1]);
else
    mu = mean(M, 1, 'omitnan');
    sem = std(M, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(M), 1));
    fill([t, fliplr(t)], [mu + sem, fliplr(mu - sem)], ...
        [0.7 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(t, mu, 'k-', 'LineWidth', 1.5);
    xline(0, '--', 'Color', [0.4 0.4 0.4]);
    xlim([t(1) t(end)]);
end
xlabel('Time from behavior onset (s, NPX)');
ylabel(ylab);
title(title_str, 'Interpreter', 'none');
box off
end
