%% Test_behavior_speed_vocalization_NPX_timewarp
% Same loading/alignment as Test_behavior_speed_vocalization_NPX.m, but
% peri-event traces are TIMEWARPED:
%   - 10 s before onset           (real time)
%   - 10 equal bins onset->offset (warped across event duration)
%   - 10 s after offset           (real time)
%
% Keeps the fixed-window NPX script untouched.
% Optionally plots a pooled figure (all animals/sessions) via plot_pooled.
%
% Behavior selection:
%   behavior_labels = {'Escape'};
%   behavior_labels = {'Escape','Evasion'};  % union of onsets (any in list)

%% Paths / options
experiment_folder = '\\experimentfs.home.bccn-berlin.de\experiment';
datasets_two = [experiment_folder, ...
    '\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\TWO ANIMALS'];

behavior_labels = {'Escape', 'CC','CB','CD', 'Pounce', 'Evasion'};
behavior_labels = {'Pounce'};

% Time-warp layout
pre_sec       = 10;     % seconds before onset
n_warp_bins   = 500;     % equal bins spanning onset -> offset
post_sec      = 10;     % seconds after offset
bin_size      = 1/30;   % sampling for pre/post (real time)
call_smooth_sec = 0.05;

% Which sessions to plot: [] = all, or indices / folder-name cellstr
sessions_to_plot = [];

% Plotting mode
plot_per_session = true;   % keep individual session figures
plot_pooled      = true;   % one figure pooling all animals/sessions

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

pre_t  = (-pre_sec):bin_size:0;
post_t = 0:bin_size:post_sec;
warp_frac = ((1:n_warp_bins) - 0.5) / n_warp_bins;  % bin centers in [0 1]

% Display x-axis: three EQUAL-width segments (pre | warped bout | post).
% Sampling still uses real pre/post seconds and 10 warp bins; only plot spacing changes.
x_pre  = linspace(0, 1, numel(pre_t));
x_warp = linspace(1, 2, numel(warp_frac));
x_post = linspace(2, 3, numel(post_t));
x_all  = [x_pre, x_warp, x_post];
n_pre  = numel(pre_t);
call_smooth_bins = max(1, round(call_smooth_sec / bin_size));

% Accumulators for pooled (all animals, all sessions)
pool_actor_speed   = [];
pool_partner_speed = [];
pool_voc           = [];

for si = session_idx
    session_name = session_dirs(si).name;
    session_path = fullfile(datasets_two, session_name);
    fprintf('\n==== %s ====\n', session_name);

    tokens = strsplit(strtrim(session_name), ' ');
    if numel(tokens) < 2
        warning('Cannot parse animal names from folder "%s". Skipping.', session_name);
        continue
    end
    animal_A = tokens{1};
    animal_B = tokens{2};

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
        warning('Focal "%s" is not one of {%s,%s}. Skipping.', focal_name, animal_A, animal_B);
        continue
    end
    fprintf('Animals: %s, %s | focal(animal_pos)=%s | partner_pos=%s\n', ...
        animal_A, animal_B, focal_name, other_name);

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

    %% Timewarped peri-event matrices
    spd_A_when_A = local_align_timewarp(speed.(animal_A), speed_time, beh_A.Start, beh_A.End, pre_t, warp_frac, post_t);
    spd_B_when_A = local_align_timewarp(speed.(animal_B), speed_time, beh_A.Start, beh_A.End, pre_t, warp_frac, post_t);
    spd_B_when_B = local_align_timewarp(speed.(animal_B), speed_time, beh_B.Start, beh_B.End, pre_t, warp_frac, post_t);
    spd_A_when_B = local_align_timewarp(speed.(animal_A), speed_time, beh_B.Start, beh_B.End, pre_t, warp_frac, post_t);

    voc_when_A = local_align_call_timewarp(call_begin, call_end, beh_A.Start, beh_A.End, ...
        pre_t, warp_frac, post_t, call_smooth_bins, n_pre);
    voc_when_B = local_align_call_timewarp(call_begin, call_end, beh_B.Start, beh_B.End, ...
        pre_t, warp_frac, post_t, call_smooth_bins, n_pre);

    fprintf('Finite timewarp speed rows: %s=%d/%d, %s=%d/%d\n', ...
        animal_A, sum(any(isfinite(spd_A_when_A), 2)), size(spd_A_when_A, 1), ...
        animal_B, sum(any(isfinite(spd_B_when_B), 2)), size(spd_B_when_B, 1));

    % Pool actor/partner roles across both animals in this session
    pool_actor_speed   = [pool_actor_speed;   spd_A_when_A; spd_B_when_B];
    pool_partner_speed = [pool_partner_speed; spd_B_when_A; spd_A_when_B];
    pool_voc           = [pool_voc;           voc_when_A;   voc_when_B];

    %% Figure: 2 x 3 (per session)
    if plot_per_session
        figure('Name', ['Behavior speed+USV NPX timewarp | ', session_name], 'Color', 'w');
        tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        beh_str = strjoin(behavior_labels, '/');

        local_plot_timewarp(nexttile, x_all, spd_A_when_A, ...
            sprintf('When %s %s: %s speed (actor)  n=%d', animal_A, beh_str, animal_A, size(spd_A_when_A, 1)), 'Speed (px/s)');
        local_plot_timewarp(nexttile, x_all, spd_B_when_A, ...
            sprintf('When %s %s: %s speed (partner)  n=%d', animal_A, beh_str, animal_B, size(spd_B_when_A, 1)), 'Speed (px/s)');
        local_plot_timewarp(nexttile, x_all, voc_when_A, ...
            sprintf('When %s %s: vocalization rate  n=%d', animal_A, beh_str, size(voc_when_A, 1)), 'Call occupancy');

        local_plot_timewarp(nexttile, x_all, spd_B_when_B, ...
            sprintf('When %s %s: %s speed (actor)  n=%d', animal_B, beh_str, animal_B, size(spd_B_when_B, 1)), 'Speed (px/s)');
        local_plot_timewarp(nexttile, x_all, spd_A_when_B, ...
            sprintf('When %s %s: %s speed (partner)  n=%d', animal_B, beh_str, animal_A, size(spd_A_when_B, 1)), 'Speed (px/s)');
        local_plot_timewarp(nexttile, x_all, voc_when_B, ...
            sprintf('When %s %s: vocalization rate  n=%d', animal_B, beh_str, size(voc_when_B, 1)), 'Call occupancy');

        sgtitle(sprintf('%s | {%s} | NPX timewarp: -%gs | %d bins | +%gs | focal=%s', ...
            session_name, beh_str, pre_sec, n_warp_bins, post_sec, focal_name), 'Interpreter', 'none');
    end
end

%% Pooled figure: all animals / all sessions
if plot_pooled
    beh_str = strjoin(behavior_labels, '/');
    fprintf('\n==== POOLED all animals/sessions ====\n');
    fprintf('Pooled events: actor=%d partner=%d voc=%d\n', ...
        size(pool_actor_speed, 1), size(pool_partner_speed, 1), size(pool_voc, 1));

    figure('Name', ['Behavior speed+USV NPX timewarp | POOLED'], 'Color', 'w');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    local_plot_timewarp(nexttile, x_all, pool_actor_speed, ...
        sprintf('Pooled: actor speed  n=%d', size(pool_actor_speed, 1)), 'Speed (px/s)');
    local_plot_timewarp(nexttile, x_all, pool_partner_speed, ...
        sprintf('Pooled: partner speed  n=%d', size(pool_partner_speed, 1)), 'Speed (px/s)');
    local_plot_timewarp(nexttile, x_all, pool_voc, ...
        sprintf('Pooled: vocalization rate  n=%d', size(pool_voc, 1)), 'Call occupancy');

    sgtitle(sprintf('ALL animals/sessions | {%s} | NPX timewarp: -%gs | %d bins | +%gs', ...
        beh_str, pre_sec, n_warp_bins, post_sec), 'Interpreter', 'none');
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
if ~ismember('BeginTimes', CallStats.Properties.VariableNames) || ...
        ~ismember('EndTimes', CallStats.Properties.VariableNames)
    error('Call table in %s lacks BeginTimes/EndTimes. Columns: %s', ...
        xf(1).name, strjoin(CallStats.Properties.VariableNames, ', '));
end
end

function M = local_align_timewarp(trace, trace_time, onsets, offsets, pre_t, warp_frac, post_t)
n = numel(onsets);
M = nan(n, numel(pre_t) + numel(warp_frac) + numel(post_t));
if n == 0 || isempty(trace)
    return
end
for k = 1:n
    ons = onsets(k);
    off = offsets(k);
    dur = off - ons;
    if ~(isfinite(ons) && isfinite(off)) || dur <= 0
        continue
    end
    t = [ons + pre_t, ons + dur * warp_frac, off + post_t];
    M(k, :) = interp1(trace_time, trace, t, 'linear', NaN);
end
end

function M = local_align_call_timewarp(call_begin, call_end, onsets, offsets, ...
    pre_t, warp_frac, post_t, smooth_bins, n_pre)
n = numel(onsets);
n_cols = numel(pre_t) + numel(warp_frac) + numel(post_t);
M = nan(n, n_cols);
if n == 0
    return
end
if isempty(call_begin)
    M(:) = 0;
    return
end
for k = 1:n
    ons = onsets(k);
    off = offsets(k);
    dur = off - ons;
    if ~(isfinite(ons) && isfinite(off)) || dur <= 0
        continue
    end
    t = [ons + pre_t, ons + dur * warp_frac, off + post_t];
    occ = zeros(1, numel(t));
    for c = 1:numel(call_begin)
        occ = occ | (t >= call_begin(c) & t <= call_end(c));
    end
    occ = double(occ);
    if smooth_bins > 1
        n_warp_local = numel(warp_frac);
        occ(1:n_pre) = movmean(occ(1:n_pre), smooth_bins);
        occ(n_pre + n_warp_local + 1:end) = movmean(occ(n_pre + n_warp_local + 1:end), smooth_bins);
    end
    M(k, :) = occ;
end
end

function local_plot_timewarp(ax, x_all, M, title_str, ylab)
axes(ax); %#ok<LAXES>
hold on
if isempty(M)
    text(mean(x_all([1 end])), 0.5, 'no events', 'HorizontalAlignment', 'center');
    ylim([0 1]);
elseif all(isnan(M(:)))
    text(mean(x_all([1 end])), 0.5, sprintf('n=%d but no overlap with tracking', size(M, 1)), ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none');
    ylim([0 1]);
else
    mu = mean(M, 1, 'omitnan');
    mu = movmean(mu,25);
    sem = std(M, 0, 1, 'omitnan') ./ sqrt(max(1, sum(~isnan(M), 1)));
    fill([x_all, fliplr(x_all)], [mu + sem, fliplr(mu - sem)], ...
        [0.7 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(x_all, mu, 'k-', 'LineWidth', 1.5);
end
xline(1, '--', 'Color', [0.2 0.2 0.8]);  % onset
xline(2, '--', 'Color', [0.8 0.2 0.2]);  % offset (end of warp)
xlim([x_all(1) x_all(end)]);
xticks([0 1 2 3]);
xticklabels({'-10 s', 'onset', 'offset', '+10 s'});
xlabel({'equal x segments: pre | timewarped bout (10 bins) | post'; 'blue dashed = onset, red dashed = offset'});
ylabel(ylab);
title(title_str, 'Interpreter', 'none');
box off
end
