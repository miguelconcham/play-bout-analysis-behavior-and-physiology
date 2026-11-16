%% Test_darting_speed_alignment
% Two-animal darting speed check (no sync models).
%
% For each session in DataSets\TWO ANIMALS:
%   1) Read ELAN behavior + tracking structure
%   2) Identify the two playing animals from behavior
%   3) Map tracking animal_pos = focal (from focal txt), partner_pos = other
%   4) Extract Escape ("darting") events per actor
%   5) Align speed of darter and partner to dart onset (video time)
%
% Figure (4 subplots):
%   (1) Animal A darts -> mean speed of A (darter)
%   (2) Animal A darts -> mean speed of B (partner)
%   (3) Animal B darts -> mean speed of B (darter)
%   (4) Animal B darts -> mean speed of A (partner)
%
% Uses video-clock alignment only (Behavior Start/End and traking_structure.time).

%% Paths / options
experiment_folder = '\\experimentfs.home.bccn-berlin.de\experiment';
datasets_two = [experiment_folder, ...
    '\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\TWO ANIMALS'];

% Darting label(s) in the ethogram. Change if needed (e.g. {'Escape','Evasion'}).
darting_labels = {'Escape'};

% Peri-event window (seconds, video time) and plot binning
hist_range = [-1 2];
bin_size   = 1/30;   % ~video frame rate

% Which sessions to plot: [] = all, or indices / folder-name cellstr
sessions_to_plot = [];  % e.g. 1, or {'C1 C2 0620'}

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
    if ~isKey(focal_map, session_name)
        % tolerate double-space variants in the focal txt
        key_alt = regexprep(session_name, '\s+', ' ');
        matched = '';
        for k = keys(focal_map)
            if strcmp(regexprep(k{1}, '\s+', ' '), key_alt)
                matched = k{1};
                break
            end
        end
        if isempty(matched)
            warning('No focal entry for "%s". Skipping.', session_name);
            continue
        end
        focal_name = focal_map(matched);
    else
        focal_name = focal_map(session_name);
    end

    if strcmp(focal_name, animal_A)
        other_name = animal_B;
    elseif strcmp(focal_name, animal_B)
        other_name = animal_A;
    else
        warning('Focal "%s" is not one of {%s, %s}. Skipping.', focal_name, animal_A, animal_B);
        continue
    end
    fprintf('Behavior/folder animals: %s, %s | focal(animal_pos)=%s | partner_pos=%s\n', ...
        animal_A, animal_B, focal_name, other_name);

    %% Load behavior (ELAN, video time)
    Behavior = local_load_behavior(session_path);
    beh_animals = unique(Behavior.Animal);
    beh_animals = beh_animals(~ismember(beh_animals, {'Session_structure','Reversal',''}));
    fprintf('Behavior actors: %s\n', strjoin(string(beh_animals), ', '));
    if ~all(ismember({animal_A, animal_B}, beh_animals))
        warning('Folder animals not both present in ELAN for %s.', session_name);
    end

    is_dart = ismember(Behavior.Type, darting_labels) | ismember(Behavior.Type2, darting_labels);
    dart_A = Behavior(is_dart & strcmp(Behavior.Animal, animal_A), :);
    dart_B = Behavior(is_dart & strcmp(Behavior.Animal, animal_B), :);
    fprintf('Darting (%s): %s=%d, %s=%d\n', strjoin(darting_labels, '/'), ...
        animal_A, height(dart_A), animal_B, height(dart_B));

    %% Load tracking and compute speed (video time via .time)
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

    track_time = TS.time(:)';
    if numel(track_time) ~= size(TS.animal_pos, 1)
        error('time length (%d) != animal_pos rows (%d) in %s', ...
            numel(track_time), size(TS.animal_pos, 1), tr_file(1).name);
    end

    % Speed at midpoints between consecutive samples (same idea as GENERATE)
    dt = diff(track_time);
    dt(dt <= 0) = NaN;
    focal_speed  = sqrt(sum(diff(TS.animal_pos).^2, 2)) ./ dt(:);
    partner_speed = sqrt(sum(diff(TS.partner_pos).^2, 2)) ./ dt(:);
    speed_time = track_time(1:end-1) + dt(:)' / 2;

    % Named animal -> speed series
    speed = struct();
    speed.(focal_name) = focal_speed;
    speed.(other_name) = partner_speed;

    %% Peri-event speed matrices
    spd_A_when_A = local_align_speed(speed.(animal_A), speed_time, dart_A.Start, psth_time);
    spd_B_when_A = local_align_speed(speed.(animal_B), speed_time, dart_A.Start, psth_time);
    spd_B_when_B = local_align_speed(speed.(animal_B), speed_time, dart_B.Start, psth_time);
    spd_A_when_B = local_align_speed(speed.(animal_A), speed_time, dart_B.Start, psth_time);

    %% Four-panel figure
    fig = figure('Name', ['Darting speed | ', session_name], 'Color', 'w');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    local_plot_mean_sem(nexttile, psth_time, spd_A_when_A, ...
        sprintf('When %s darts: %s speed (darter)  n=%d', animal_A, animal_A, size(spd_A_when_A, 1)));
    local_plot_mean_sem(nexttile, psth_time, spd_B_when_A, ...
        sprintf('When %s darts: %s speed (partner)  n=%d', animal_A, animal_B, size(spd_B_when_A, 1)));
    local_plot_mean_sem(nexttile, psth_time, spd_B_when_B, ...
        sprintf('When %s darts: %s speed (darter)  n=%d', animal_B, animal_B, size(spd_B_when_B, 1)));
    local_plot_mean_sem(nexttile, psth_time, spd_A_when_B, ...
        sprintf('When %s darts: %s speed (partner)  n=%d', animal_B, animal_A, size(spd_A_when_B, 1)));

    sgtitle(sprintf('%s | darting=%s | focal(animal\\_pos)=%s', ...
        session_name, strjoin(darting_labels, '/'), focal_name), 'Interpreter', 'tex');
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
    % e.g. "C1 C2  0620, C1 is focal"
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

function Behavior = local_load_behavior(session_path)
elan = dir(fullfile(session_path, 'ELAN*'));
if isempty(elan)
    error('No ELAN* behavior file in %s', session_path);
end
raw = readtable(fullfile(session_path, elan(1).name), ...
    'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);

% Same convention as GENERATE_STRUCUTRE_animal_synch: drop 2nd column
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

function local_plot_mean_sem(ax, t, M, title_str)
axes(ax); %#ok<LAXES>
hold on
if isempty(M) || all(isnan(M(:)))
    text(0.5, 0.5, 'no events', 'HorizontalAlignment', 'center');
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
xlabel('Time from dart onset (s)');
ylabel('Speed (px/s)');
title(title_str, 'Interpreter', 'none');
box off
end
