function summary = GENERATE_BEHAVIOR_ALIGNED_CSD_LFP(npx_data_dir, saving_root, t1, t2, recommended_min_events)
%GENERATE_BEHAVIOR_ALIGNED_CSD_LFP Export paired raw LFP around behavior events.
%
% The function writes one onset/offset pair at a time to avoid retaining all
% grouped LFP stacks in memory. Raw, unfiltered, unreferenced LFP is exported.
% It also saves an evoked LFP PSTH for each alignment after subtracting the
% event's pre-onset baseline from BOTH its onset and offset snippets.
%
% Analyses written:
%   behavior/<type>/...       each of 14 behaviors, self + partners pooled
%   category/<name>/...       play, aggression, exploration (roles pooled)
%   play_role/self_play/...   all implanted-animal play behaviors
%   play_role/other_play/...  all partner-animal play behaviors
%   play_bout/all/...         bouts returned by play_bout.m
%
% Only events whose complete onset and offset windows are inside a
% 'Partners session' interval are retained. Events are also excluded if any
% genuine behavior annotation from either animal overlaps the onset baseline.
%
% Every onset/offset leaf contains matching event_ids.npy. The Python paired
% loader refuses mismatched IDs and applies the same onset baseline again
% before kCSD.

if nargin < 5 || isempty(recommended_min_events)
    recommended_min_events = 5;
end
if t1 >= 0 || t2 <= 0
    error('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:Window', 'Need t1 < 0 < t2');
end

data_root        = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
synch_directory  = fullfile(data_root, 'Synch data');
area_limit_table = fullfile(data_root, 'Area_limits_GoodLooking.xlsx');
behavior_data    = fullfile(data_root, 'Behavior backups');
area2analyse     = 'PAG';
sr_LFP           = 2500;
hard_x_npx2      = [8 40; 258 290; 508 540; 758 790];

all_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', ...
    'CB', 'CD', 'Grooming', 'Rearing', 'Sniffing', 'Scratching', ...
    'PounceI', 'Bite'};
play_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', ...
    'CB', 'CD'};
aggression_behaviors = {'PounceI', 'Bite'};
exploration_behaviors = {'Grooming', 'Rearing', 'Sniffing', 'Scratching'};

animal_code = local_last_path_part(npx_data_dir);
animal_code_params = strsplit(animal_code, ' ');
if numel(animal_code_params) < 3
    error('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:AnimalName', ...
        'Cannot parse implanted animal from folder name: %s', animal_code);
end
animal_batch    = animal_code_params{1};
repeated_animal = animal_code_params{3};

%% Behavior annotations and video -> NPX synchronization
synch_file = fullfile(synch_directory, animal_code, 'synch_model_video2NPX.mat');
load(synch_file, 'synch_model_video2NPX');

behavior_file = fullfile(behavior_data, [animal_code, '.txt']);
Behavior = readtable(behavior_file);
Behavior(:, 2) = [];
Behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
Behavior.Animal = local_to_cellstr(Behavior.Animal);
Behavior.Type = local_to_cellstr(Behavior.Type);
Behavior.OriginalRow = (1:height(Behavior))';
Behavior.Type2 = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
Behavior.Type2(ismember(Behavior.Type2, {'Scratch', 'Scratching'})) = {'Scratching'};
Behavior.Type2(ismember(Behavior.Type2, '')) = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'), :) = [];
Behavior.Start = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End   = predict(synch_model_video2NPX, Behavior.End);

partner_sessions = Behavior{ismember(Behavior.Type, 'Partners session'), {'Start', 'End'}};
if isempty(partner_sessions)
    warning('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:NoPartnerSession', ...
        '%s: no Partners session annotation; skipping', animal_code);
    summary = local_empty_summary();
    return
end

% Baseline contamination is checked against every genuine annotation, not
% only the 14 response labels. Session-structure rows are metadata.
baseline_annotations = Behavior( ...
    ~ismember(Behavior.Animal, 'Session_structure') & ...
    ~ismember(Behavior.Type, 'Partners session') & ...
    isfinite(Behavior.Start) & isfinite(Behavior.End) & ...
    Behavior.End > Behavior.Start, :);

% Remove exact duplicate ethogram rows once, before building pooled groups.
ethogram_mask = ismember(Behavior.Type2, all_behaviors);
ethogram = Behavior(ethogram_mask, :);
[~, unique_rows] = unique(ethogram(:, {'Animal', 'Start', 'End', 'Type2'}), ...
    'rows', 'stable');
n_duplicate_annotations = height(ethogram) - numel(unique_rows);
ethogram = ethogram(sort(unique_rows), :);

animal_types = unique(ethogram.Animal);
partner_names = animal_types;
partner_names(ismember(partner_names, repeated_animal)) = [];

%% Define individual, pooled-category, role-specific, and bout datasets
definitions = local_empty_definitions();

for b = 1:numel(all_behaviors)
    rows = ethogram(ismember(ethogram.Type2, all_behaviors{b}), :);
    definitions = local_append_rows(definitions, 'behavior', all_behaviors{b}, ...
        rows, repeated_animal);
end

definitions = local_append_rows(definitions, 'category', 'play', ...
    ethogram(ismember(ethogram.Type2, play_behaviors), :), repeated_animal);
definitions = local_append_rows(definitions, 'category', 'aggression', ...
    ethogram(ismember(ethogram.Type2, aggression_behaviors), :), repeated_animal);
definitions = local_append_rows(definitions, 'category', 'exploration', ...
    ethogram(ismember(ethogram.Type2, exploration_behaviors), :), repeated_animal);

self_play = ethogram(ismember(ethogram.Type2, play_behaviors) & ...
    ismember(ethogram.Animal, repeated_animal), :);
other_play = ethogram(ismember(ethogram.Type2, play_behaviors) & ...
    ismember(ethogram.Animal, partner_names), :);
definitions = local_append_rows(definitions, 'play_role', 'self_play', ...
    self_play, repeated_animal);
definitions = local_append_rows(definitions, 'play_role', 'other_play', ...
    other_play, repeated_animal);

% play_bout.m pools play from all animals and snaps each bout to its first
% and last play annotation.
config.Behavior        = Behavior;
config.repeated_animal = {repeated_animal};
config.animal_types    = animal_types;
config.play_behaviors  = play_behaviors;
config.beh_bin         = 0.01;
config.conv_length     = 1;
config.behavior_window = 0;
try
    play_bouts_table = play_bout(config);
catch ME
    warning('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:PlayBout', ...
        '%s: play_bout failed (%s); behavior exports will continue', ...
        animal_code, ME.message);
    play_bouts_table = zeros(0, 2);
end
if isempty(play_bouts_table)
    play_bouts_table = zeros(0, 2);
elseif size(play_bouts_table, 2) ~= 2
    error('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:PlayBoutShape', ...
        'play_bout must return an N-by-2 onset/offset matrix');
end
n_bouts = size(play_bouts_table, 1);
bout_rows = table( ...
    repmat({'Both'}, n_bouts, 1), ...
    play_bouts_table(:, 1), play_bouts_table(:, 2), ...
    play_bouts_table(:, 2) - play_bouts_table(:, 1), ...
    repmat({'PlayBout'}, n_bouts, 1), nan(n_bouts, 1), ...
    'VariableNames', {'Animal', 'Start', 'End', 'Length', 'Type2', 'OriginalRow'});
definitions = local_append_rows(definitions, 'play_bout', 'all', ...
    bout_rows, repeated_animal);

%% Channel map, LPAG shanks, and raw LFP
load(fullfile(npx_data_dir, ['chann_map_', area2analyse, '.mat']), ...
    'chanMap', 'xcoords', 'ycoords');
if any(ismember(xcoords, hard_x_npx2))
    NPX_Type = 2;
else
    NPX_Type = 1;
    [chanMap, xcoords, ycoords, map_ok] = ...
        local_restore_npx1_channel_192(chanMap, xcoords, ycoords);
    if ~map_ok
        warning('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:ChanMap', ...
            '%s: inconsistent NPX1 channel map; skipping', animal_code);
        summary = local_empty_summary();
        return
    end
end

area_limit = readtable(area_limit_table);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);
lpag_rows = area_limit(ismember(area_limit.area, {'LPAG'}), :);
if isempty(lpag_rows)
    warning('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:NoLPAG', ...
        '%s: no LPAG in area limits; skipping', animal_code);
    summary = local_empty_summary();
    return
end

jobs = local_build_shank_jobs(NPX_Type, lpag_rows, xcoords, ycoords, ...
    chanMap, hard_x_npx2);
if isempty(jobs)
    warning('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:NoLPAGChannels', ...
        '%s: LPAG listed but no matching channels; skipping', animal_code);
    summary = local_empty_summary();
    return
end

disp(['LOADING RAW LFP: ', animal_code])
if NPX_Type == 2
    load(fullfile(npx_data_dir, ['LFP_', area2analyse, '.mat']), 'LFP');
else
    fid = fopen(fullfile(npx_data_dir, ['LFP_', area2analyse, '.dat']), 'r');
    if fid < 0
        error('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:LFP', ...
            'Could not open LFP file for %s', animal_code);
    end
    cleaner = onCleanup(@() fclose(fid));
    LFP = fread(fid, '*int16');
    clear cleaner
    LFP = reshape(LFP, 384, numel(LFP) / 384);
end
n_time_total = size(LFP, 2);
disp('RAW LFP LOADED')

pre_samp   = round(abs(t1) * sr_LFP);
post_samp  = round(t2 * sr_LFP);
n_time_win = pre_samp + post_samp;
T0         = pre_samp; % zero-based alignment sample for NumPy
time_ms    = ((0:n_time_win - 1) - T0) / sr_LFP * 1000;

summary = local_empty_summary();
for j = 1:numel(jobs)
    job = jobs(j);
    [y_sh, order] = sort(job.y_coords(:));
    x_sh = job.x_coords(order);
    shank_rows = job.lfp_rows(order);
    valid_rows = isfinite(shank_rows) & shank_rows == round(shank_rows) & ...
        shank_rows >= 1 & shank_rows <= size(LFP, 1);
    if ~all(valid_rows)
        warning('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:InvalidChannelRows', ...
            '%s %s: dropping %d invalid channel-map entries', ...
            animal_code, job.shank_label, sum(~valid_rows));
        shank_rows = shank_rows(valid_rows);
        x_sh = x_sh(valid_rows);
        y_sh = y_sh(valid_rows);
    end
    if isempty(shank_rows)
        warning('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:EmptyShank', ...
            '%s %s: no valid LFP channel rows; skipping shank', ...
            animal_code, job.shank_label);
        continue
    end

    disp(sprintf('  %s: n_chan=%d', job.shank_label, numel(shank_rows)));
    for d = 1:numel(definitions)
        definition = definitions(d);
        events_all = local_qc_events(definition, baseline_annotations, partner_sessions, ...
            t1, t2, sr_LFP, n_time_total);
        included = events_all.Included;
        events = events_all(included, :);
        n_saved = height(events);

        group_type = local_safe_path(definition.group_type);
        group_name = local_safe_path(definition.group_name);
        if NPX_Type == 1
            session_tail = {animal_code};
        else
            session_tail = {animal_code, job.shank_label};
        end
        onset_dir = fullfile(saving_root, group_type, group_name, 'onset', ...
            session_tail{:});
        offset_dir = fullfile(saving_root, group_type, group_name, 'offset', ...
            session_tail{:});

        entry = local_summary_entry(animal_code, job.shank_label, ...
            definition.group_type, definition.group_name, height(events_all), ...
            n_saved, sum(events_all.BaselineContaminated), ...
            sum(ismember(events_all.ExclusionReason, 'outside_recording')), ...
            sum(ismember(events_all.ExclusionReason, 'outside_partner_session')), ...
            n_saved < recommended_min_events, n_duplicate_annotations);

        if n_saved == 0
            % Invalidate old exports so a rerun cannot leave a stale pair
            % discoverable by Python.
            local_invalidate_alignment(onset_dir);
            local_invalidate_alignment(offset_dir);
            entry.status = 'NO_VALID_EVENTS';
            summary(end + 1) = entry; %#ok<AGROW>
            continue
        end

        onset_idx = events.OnsetSample;
        offset_idx = events.OffsetSample;
        lfp_onset = local_cut_snippets(LFP, shank_rows, onset_idx, ...
            pre_samp, post_samp, n_time_win);
        lfp_offset = local_cut_snippets(LFP, shank_rows, offset_idx, ...
            pre_samp, post_samp, n_time_win);

        % LFP PSTHs use exactly the baseline requested for paired kCSD:
        % one pre-onset baseline per event/channel, also applied to offset.
        onset_baseline = mean(lfp_onset(:, :, 1:pre_samp), 3);
        mean_baseline = reshape(mean(onset_baseline, 1), numel(shank_rows), 1);
        evoked_onset = bsxfun(@minus, ...
            reshape(mean(lfp_onset, 1), numel(shank_rows), n_time_win), ...
            mean_baseline);
        evoked_offset = bsxfun(@minus, ...
            reshape(mean(lfp_offset, 1), numel(shank_rows), n_time_win), ...
            mean_baseline);

        common_meta = struct();
        common_meta.animal_code = animal_code;
        common_meta.implanted_animal = repeated_animal;
        common_meta.npx_type = NPX_Type;
        common_meta.probe_num = job.probe_num;
        common_meta.shank_label = job.shank_label;
        common_meta.channel_ids = shank_rows(:);
        common_meta.group_type = definition.group_type;
        common_meta.group_name = definition.group_name;
        common_meta.fs = sr_LFP;
        common_meta.t1 = t1;
        common_meta.t2 = t2;
        common_meta.T0 = T0;
        common_meta.n_events = n_saved;
        common_meta.n_candidates = height(events_all);
        common_meta.n_channels = numel(shank_rows);
        common_meta.n_time = n_time_win;
        common_meta.raw_lfp_filtered = false;
        common_meta.rereferenced = false;
        common_meta.raw_lfp_baseline_corrected = false;
        common_meta.evoked_lfp_baseline_source = 'onset';
        common_meta.onset_baseline_window_s = [t1 0];
        common_meta.partner_session_only = true;
        common_meta.behavior_free_baseline_required = true;
        common_meta.recommended_min_events = recommended_min_events;
        common_meta.eligible_for_csd = n_saved >= recommended_min_events;
        common_meta.play_behaviors = play_behaviors;
        common_meta.aggression_behaviors = aggression_behaviors;
        common_meta.exploration_behaviors = exploration_behaviors;
        common_meta.partner_sessions = partner_sessions;

        local_save_alignment(onset_dir, lfp_onset, evoked_onset, ...
            x_sh, y_sh, T0, sr_LFP, time_ms, events, events_all, ...
            common_meta, 'onset');
        local_save_alignment(offset_dir, lfp_offset, evoked_offset, ...
            x_sh, y_sh, T0, sr_LFP, time_ms, events, events_all, ...
            common_meta, 'offset');

        disp(sprintf('    saved %s/%s: n=%d%s', ...
            definition.group_type, definition.group_name, n_saved, ...
            local_low_n_text(n_saved, recommended_min_events)));
        entry.status = 'OK';
        summary(end + 1) = entry; %#ok<AGROW>

        clear lfp_onset lfp_offset onset_baseline mean_baseline
    end
end
end

%% -------------------------------------------------------------------------
function definitions = local_empty_definitions()
definitions = struct('group_type', {}, 'group_name', {}, 'animal', {}, ...
    'behavior', {}, 'role', {}, 'onset', {}, 'offset', {}, 'source_row', {});
end

function definitions = local_append_rows(definitions, group_type, group_name, rows, repeated_animal)
n = height(rows);
entry.group_type = group_type;
entry.group_name = group_name;
entry.animal = local_to_cellstr(rows.Animal);
entry.behavior = local_to_cellstr(rows.Type2);
entry.role = repmat({'other'}, n, 1);
entry.role(ismember(rows.Animal, repeated_animal)) = {'self'};
if strcmp(group_type, 'play_bout')
    entry.role(:) = {'pooled'};
end
entry.onset = rows.Start;
entry.offset = rows.End;
entry.source_row = rows.OriginalRow;
definitions(end + 1) = entry;
end

function events = local_qc_events(definition, baseline_annotations, partner_sessions, t1, t2, fs, n_time)
n = numel(definition.onset);
group_type = repmat({definition.group_type}, n, 1);
group_name = repmat({definition.group_name}, n, 1);
animal = definition.animal(:);
role = definition.role(:);
behavior = definition.behavior(:);
onset_s = definition.onset(:);
offset_s = definition.offset(:);
duration_s = offset_s - onset_s;
source_row = definition.source_row(:);

onset_sample = round(onset_s * fs);
offset_sample = round(offset_s * fs);
% Pack synchronized onset and offset sample numbers into one stable int64.
% This detects stale onset/offset pairs across reruns, unlike row positions.
event_id = bitshift(int64(onset_sample), 32) + int64(offset_sample);
partner_session_index = zeros(n, 1);
baseline_contaminated = false(n, 1);
other_behavior_onset_window = false(n, 1);
other_behavior_offset_window = false(n, 1);
onset_offset_windows_overlap = duration_s < (t2 - t1);
included = true(n, 1);
exclusion_reason = repmat({''}, n, 1);

for e = 1:n
    if ~isfinite(onset_s(e)) || ~isfinite(offset_s(e)) || offset_s(e) <= onset_s(e)
        included(e) = false;
        exclusion_reason{e} = 'invalid_duration';
        continue
    end

    recording_ok = onset_sample(e) > round(abs(t1) * fs) && ...
        onset_sample(e) + round(t2 * fs) - 1 <= n_time && ...
        offset_sample(e) > round(abs(t1) * fs) && ...
        offset_sample(e) + round(t2 * fs) - 1 <= n_time;
    if ~recording_ok
        included(e) = false;
        exclusion_reason{e} = 'outside_recording';
        continue
    end

    for s = 1:size(partner_sessions, 1)
        if onset_s(e) + t1 >= partner_sessions(s, 1) && ...
                offset_s(e) + t2 <= partner_sessions(s, 2)
            partner_session_index(e) = s;
            break
        end
    end
    if partner_session_index(e) == 0
        included(e) = false;
        exclusion_reason{e} = 'outside_partner_session';
        continue
    end

    % Strict interval overlap with [onset+t1, onset). Events beginning
    % exactly at onset do not contaminate the baseline.
    baseline_overlap = baseline_annotations.Start < onset_s(e) & ...
        baseline_annotations.End > onset_s(e) + t1;
    if isfinite(source_row(e))
        baseline_overlap = baseline_overlap & ...
            baseline_annotations.OriginalRow ~= source_row(e);
    end
    baseline_contaminated(e) = any(baseline_overlap);
    if baseline_contaminated(e)
        included(e) = false;
        exclusion_reason{e} = 'baseline_behavior_overlap';
    end

    onset_window = baseline_annotations.Start < onset_s(e) + t2 & ...
        baseline_annotations.End > onset_s(e) + t1;
    offset_window = baseline_annotations.Start < offset_s(e) + t2 & ...
        baseline_annotations.End > offset_s(e) + t1;
    if isfinite(source_row(e))
        not_source = baseline_annotations.OriginalRow ~= source_row(e);
        onset_window = onset_window & not_source;
        offset_window = offset_window & not_source;
    end
    other_behavior_onset_window(e) = any(onset_window);
    other_behavior_offset_window(e) = any(offset_window);
end

events = table(event_id, group_type, group_name, animal, role, behavior, ...
    onset_s, offset_s, duration_s, source_row, onset_sample, offset_sample, ...
    partner_session_index, baseline_contaminated, ...
    other_behavior_onset_window, other_behavior_offset_window, ...
    onset_offset_windows_overlap, included, exclusion_reason, ...
    'VariableNames', {'EventID', 'GroupType', 'GroupName', 'Animal', 'Role', ...
    'Behavior', 'Onset_s', 'Offset_s', 'Duration_s', 'OriginalRow', ...
    'OnsetSample', 'OffsetSample', 'PartnerSessionIndex', ...
    'BaselineContaminated', 'OtherBehaviorInOnsetWindow', ...
    'OtherBehaviorInOffsetWindow', 'OnsetOffsetWindowsOverlap', ...
    'Included', 'ExclusionReason'});
end

function local_save_alignment(folder, lfp_stack, evoked_lfp, x, y, T0, fs, ...
        time_ms, events, events_all, common_meta, alignment)
if ~exist(folder, 'dir')
    mkdir(folder);
end
invalid_file = fullfile(folder, 'INVALID_NO_VALID_EVENTS.txt');
if exist(invalid_file, 'file')
    delete(invalid_file);
end
write_csd_npy(fullfile(folder, 'lfp.npy'), single(lfp_stack));
write_csd_npy(fullfile(folder, 'evoked_lfp.npy'), single(evoked_lfp));
write_csd_npy(fullfile(folder, 'x_coords.npy'), double(x(:)));
write_csd_npy(fullfile(folder, 'y_coords.npy'), double(y(:)));
write_csd_npy(fullfile(folder, 'T0.npy'), double(T0));
write_csd_npy(fullfile(folder, 'fs.npy'), double(fs));
write_csd_npy(fullfile(folder, 'time_ms.npy'), double(time_ms(:)));
write_csd_npy(fullfile(folder, 'event_ids.npy'), int64(events.EventID));

meta = common_meta;
meta.alignment = alignment;
meta.event_ids = events.EventID;
save(fullfile(folder, 'meta.mat'), 'meta', 'events', 'events_all');
end

function stack = local_cut_snippets(LFP, channel_rows, event_idx, pre_samp, post_samp, n_time_win)
n_events = numel(event_idx);
n_channels = numel(channel_rows);
stack = zeros(n_events, n_channels, n_time_win, 'single');
for e = 1:n_events
    first_sample = event_idx(e) - pre_samp;
    last_sample = event_idx(e) + post_samp - 1;
    stack(e, :, :) = LFP(channel_rows, first_sample:last_sample);
end
end

function jobs = local_build_shank_jobs(NPX_Type, lpag_rows, xcoords, ycoords, chanMap, hard_x)
jobs = struct('probe_num', {}, 'shank_label', {}, 'lfp_rows', {}, ...
    'x_coords', {}, 'y_coords', {});

if NPX_Type == 1
    range = lpag_rows{1, {'ProbeNum', 'depth_start', 'depth_end'}};
    lpag_mask = ycoords >= range(2) & ycoords <= range(3);
    if ~any(lpag_mask)
        return
    end
    valid = ~isnan(xcoords) & ~isnan(ycoords) & ~isnan(chanMap);
    jobs(1).probe_num = 1;
    jobs(1).shank_label = 'NPX1';
    jobs(1).lfp_rows = chanMap(valid);
    jobs(1).x_coords = xcoords(valid);
    jobs(1).y_coords = ycoords(valid);
    return
end

probe_list = unique(lpag_rows.ProbeNum);
for p = 1:numel(probe_list)
    probe_n = probe_list(p);
    if probe_n < 1 || probe_n > size(hard_x, 1)
        continue
    end
    rows = lpag_rows(lpag_rows.ProbeNum == probe_n, :);
    shank_mask = ismember(xcoords, hard_x(probe_n, :)) & ...
        isfinite(xcoords) & isfinite(ycoords) & isfinite(chanMap);
    lpag_mask = shank_mask & ycoords >= min(rows.depth_start) & ...
        ycoords <= max(rows.depth_end);
    if ~any(lpag_mask)
        continue
    end
    jobs(end + 1).probe_num = probe_n; %#ok<AGROW>
    jobs(end).shank_label = sprintf('probe%d', probe_n);
    jobs(end).lfp_rows = chanMap(shank_mask);
    jobs(end).x_coords = xcoords(shank_mask);
    jobs(end).y_coords = ycoords(shank_mask);
end
end

function [chanMap, xcoords, ycoords, ok] = local_restore_npx1_channel_192(chanMap, xcoords, ycoords)
ok = true;
if ismember(192, chanMap)
    return
end
pos_191 = find(chanMap == 191);
pos_193 = find(chanMap == 193);
if isempty(pos_191) || isempty(pos_193) || pos_193 ~= pos_191 + 1
    ok = false;
    return
end
xcoords = [xcoords; NaN];
xcoords(pos_193 + 1:end) = xcoords(pos_193:end - 1);
xcoords(pos_193) = 43;
ycoords = [ycoords; NaN];
ycoords(pos_193 + 1:end) = ycoords(pos_193:end - 1);
ycoords(pos_193) = 1900;
chanMap = [chanMap; NaN];
chanMap(pos_193 + 1:end) = chanMap(pos_193:end - 1);
chanMap(pos_193) = 192;
end

function summary = local_empty_summary()
summary = struct('animal_code', {}, 'shank_label', {}, 'group_type', {}, ...
    'group_name', {}, 'n_candidates', {}, 'n_saved', {}, ...
    'n_baseline_contaminated', {}, 'n_outside_recording', {}, ...
    'n_outside_partner_session', {}, 'low_event_count', {}, ...
    'n_duplicate_annotations_removed', {}, 'status', {});
end

function entry = local_summary_entry(animal, shank, group_type, group_name, ...
        n_candidates, n_saved, n_baseline, n_recording, n_session, low_count, n_duplicates)
entry.animal_code = animal;
entry.shank_label = shank;
entry.group_type = group_type;
entry.group_name = group_name;
entry.n_candidates = n_candidates;
entry.n_saved = n_saved;
entry.n_baseline_contaminated = n_baseline;
entry.n_outside_recording = n_recording;
entry.n_outside_partner_session = n_session;
entry.low_event_count = low_count;
entry.n_duplicate_annotations_removed = n_duplicates;
entry.status = '';
end

function text = local_low_n_text(n_events, recommended_min_events)
if n_events < recommended_min_events
    text = sprintf(' [LOW N < %d]', recommended_min_events);
else
    text = '';
end
end

function local_invalidate_alignment(folder)
% Removing lfp.npy is enough to make the Python discovery ignore this leaf.
lfp_file = fullfile(folder, 'lfp.npy');
if exist(lfp_file, 'file')
    delete(lfp_file);
end
invalid_file = fullfile(folder, 'INVALID_NO_VALID_EVENTS.txt');
if exist(folder, 'dir')
    fid = fopen(invalid_file, 'w');
    if fid >= 0
        fprintf(fid, 'No valid paired events in the latest MATLAB export.\n');
        fclose(fid);
    end
end
end

function values = local_to_cellstr(values)
if iscell(values)
    values = cellfun(@char, values, 'UniformOutput', false);
elseif isstring(values) || iscategorical(values)
    values = cellstr(values);
elseif ischar(values)
    values = cellstr(values);
else
    error('GENERATE_BEHAVIOR_ALIGNED_CSD_LFP:TextType', ...
        'Expected text annotation column; got %s', class(values));
end
values = values(:);
end

function value = local_safe_path(value)
value = strrep(value, '/', '_');
value = strrep(value, '\', '_');
value = strrep(value, ' ', '_');
end

function name = local_last_path_part(path_value)
parts = strsplit(path_value, {filesep, '/', '\'});
parts = parts(~cellfun(@isempty, parts));
name = parts{end};
end
