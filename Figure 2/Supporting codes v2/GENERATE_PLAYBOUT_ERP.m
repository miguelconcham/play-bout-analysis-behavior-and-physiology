function erp_struct = GENERATE_PLAYBOUT_ERP(current_dir, window_s)
%GENERATE_PLAYBOUT_ERP  Raw mid-LPAG LFP around each play-bout onset.
%
% This is the play-bout onset cut used for delta power, with the band-power
% summary replaced by the unfiltered LFP waveform. The default window is
% -2 s to +2 s. One struct is returned per mid-LPAG site. lfp_trials is
% trials by time, in LFP_PAG sample values. erp is the mean of every trial
% whose window lies entirely inside the recording, saturated trials included.
% Analyze_playbout_erp rejects those saturated trials before the grand average.

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(fileparts(this_dir)), 'Custom functions'));
end

options = playbout_erp_defaults();
if nargin < 2 || isempty(window_s)
    window_s = options.window_s;
end
sampling_rate = options.sampling_rate;

[~, animal_code] = fileparts(current_dir);
name_parts = regexp(strtrim(animal_code), '\s+', 'split');
name_parts = name_parts(~cellfun('isempty', name_parts));
if numel(name_parts) < 3
    error('GENERATE_PLAYBOUT_ERP:AnimalName', ...
        'Folder name "%s" is not "{animal} {date} {Single|Dual|Single2}".', ...
        animal_code);
end
animal_batch = name_parts{1};
repeated_animal = name_parts{3};
animal_id = animal_batch;

data_root = recording_data_root(current_dir);
[Behavior, animal_types] = load_play_behavior( ...
    current_dir, animal_code, data_root);
play_bouts = detect_play_bouts(Behavior, animal_types, repeated_animal, ...
    options.play_behaviors);

[LFP, npx_type] = load_lfp_pag(current_dir);
selected = select_lpag_sites(current_dir, data_root, animal_batch, ...
    repeated_animal, npx_type, options.area, size(LFP, 1));

erp_struct = repmat(new_playbout_erp_session(), numel(selected), 1);
for site_n = 1:numel(selected)
    [trials, time_s] = cut_event_lfp( ...
        LFP(selected(site_n).lfp_row, :), play_bouts(:, 1), window_s, sampling_rate);

    session = new_playbout_erp_session();
    session.animal_name = animal_code;
    session.animal_id = animal_id;
    session.session_id = animal_code;
    session.neuropixels_directory = current_dir;
    session.npx_type = npx_type;
    session.electrode_index = site_n;
    session.lpag_channel = selected(site_n).channel;
    session.lpag_lfp_row = selected(site_n).lfp_row;
    session.lpag_x_um = selected(site_n).x_um;
    session.lpag_y_um = selected(site_n).y_um;
    session.sampling_rate = sampling_rate;
    session.window_s = window_s(:).';
    session.time_s = time_s;
    session.lfp_trials = trials;
    session.play_bout_onset_s = play_bouts(:, 1);
    session.play_bout_offset_s = play_bouts(:, 2);
    session.erp = mean_event_erp(trials);
    session.n_trials = size(play_bouts, 1);
    session.n_complete_trials = sum(all(isfinite(trials), 2));
    session.lfp_units = options.lfp_units;
    erp_struct(site_n) = session;

    fprintf(['%s  electrode %d  channel %g  %d play bouts  ' ...
        '%d complete LFP epochs\n'], ...
        animal_code, site_n, session.lpag_channel, session.n_trials, ...
        session.n_complete_trials);
end
clear LFP
end

function data_root = recording_data_root(current_dir)
this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
repository_data_root = fullfile(repo_root, 'Data');
candidate = fileparts(fileparts(fileparts(current_dir)));
if isfile(fullfile(candidate, 'Area_limits_GoodLooking.xlsx'))
    data_root = candidate;
else
    data_root = repository_data_root;
end
end

function [Behavior, animal_types] = load_play_behavior(current_dir, animal_code, data_root)
local_synch_file = fullfile(current_dir, 'synch_model_video2NPX.mat');
data_synch_file = fullfile(data_root, 'Synch data', animal_code, ...
    'synch_model_video2NPX.mat');
if isfile(local_synch_file)
    synch_file = local_synch_file;
elseif isfile(data_synch_file)
    synch_file = data_synch_file;
else
    error('GENERATE_PLAYBOUT_ERP:synch', ...
        'Synchronization model not found for %s.', animal_code);
end
loaded = load(synch_file, 'synch_model_video2NPX');

local_behavior_files = dir(fullfile(current_dir, 'ELAN*'));
data_behavior_file = fullfile(data_root, 'Behavior backups', [animal_code, '.txt']);
if ~isempty(local_behavior_files)
    behavior_file = fullfile(current_dir, local_behavior_files(1).name);
elseif isfile(data_behavior_file)
    behavior_file = data_behavior_file;
else
    error('GENERATE_PLAYBOUT_ERP:behavior', ...
        'Behavior file not found for %s.', animal_code);
end

Behavior = readtable(behavior_file);
Behavior(:, 2) = [];
Behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
Behavior.Type2 = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
Behavior.Type2(ismember(Behavior.Type2, {''})) = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'), :) = [];
Behavior.Start = predict(loaded.synch_model_video2NPX, Behavior.Start);
Behavior.End = predict(loaded.synch_model_video2NPX, Behavior.End);

animal_types = unique(Behavior.Animal);
animal_types(ismember(animal_types, {'Session_structure'})) = [];
end

function play_bouts = detect_play_bouts(Behavior, animal_types, repeated_animal, play_behaviors)
config.Behavior = Behavior;
config.repeated_animal = repeated_animal;
config.animal_types = animal_types;
config.play_behaviors = play_behaviors;
config.beh_bin = 0.01;
config.conv_length = 1;
config.behavior_window = 0;
play_bouts = play_bout(config);
play_bouts = play_bout_matrix(play_bouts);
if isempty(play_bouts)
    play_bouts = zeros(0, 2);
    return
end
[~, order] = sort(play_bouts(:, 1));
play_bouts = play_bouts(order, :);
end

function play_bouts = play_bout_matrix(play_bouts)
if isempty(play_bouts)
    play_bouts = zeros(0, 2);
    return
end
if isa(play_bouts, 'table')
        play_bouts = play_bouts{:, 1:min(2, size(play_bouts, 2))};
    end
play_bouts = double(play_bouts);
if size(play_bouts, 2) < 2
    error('GENERATE_PLAYBOUT_ERP:PlayBouts', ...
        'play_bout must return onset and offset times.');
end
play_bouts = play_bouts(:, 1:2);
play_bouts = play_bouts(all(isfinite(play_bouts), 2), :);
end

function [LFP, npx_type] = load_lfp_pag(current_dir)
lfp_mat_file = fullfile(current_dir, 'LFP_PAG.mat');
lfp_dat_file = fullfile(current_dir, 'LFP_PAG.dat');
if isfile(lfp_mat_file)
    npx_type = 2;
    loaded = load(lfp_mat_file, 'LFP');
    if ~isfield(loaded, 'LFP')
        error('GENERATE_PLAYBOUT_ERP:lfp', ...
            'LFP_PAG.mat does not contain a variable named LFP.');
    end
    LFP = loaded.LFP;
elseif isfile(lfp_dat_file)
    npx_type = 1;
    file_pointer = fopen(lfp_dat_file, 'r');
    if file_pointer < 0
        error('GENERATE_PLAYBOUT_ERP:lfp', 'Could not open %s.', lfp_dat_file);
    end
    cleaner = onCleanup(@() fclose(file_pointer));
    LFP = fread(file_pointer, 'int16=>int16');
    clear cleaner
    if mod(numel(LFP), 384) ~= 0
        error('GENERATE_PLAYBOUT_ERP:lfp', ...
            'LFP_PAG.dat does not contain a whole number of 384-channel samples.');
    end
    LFP = reshape(LFP, 384, []);
else
    error('GENERATE_PLAYBOUT_ERP:lfp', ...
        'Neither LFP_PAG.mat nor LFP_PAG.dat exists in %s.', current_dir);
end
if size(LFP, 1) ~= 384
    error('GENERATE_PLAYBOUT_ERP:lfp', ...
        'Expected 384 LFP channels; found %d.', size(LFP, 1));
end
end

function selected = select_lpag_sites(current_dir, data_root, animal_batch, ...
        repeated_animal, npx_type, area_name, n_lfp_rows)
channel_map_file = fullfile(current_dir, ['chann_map_', area_name, '.mat']);
if ~isfile(channel_map_file)
    error('GENERATE_PLAYBOUT_ERP:map', ...
        'Channel map not found: %s', channel_map_file);
end
loaded = load(channel_map_file, 'xcoords', 'ycoords', 'chanMap');

area_limit_file = fullfile(data_root, 'Area_limits_GoodLooking.xlsx');
if ~isfile(area_limit_file)
    error('GENERATE_PLAYBOUT_ERP:areas', ...
        'Area-limit table not found: %s', area_limit_file);
end
area_limit = readtable(area_limit_file);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    if numel(animal_batch) < 4
        error('GENERATE_PLAYBOUT_ERP:AnimalName', ...
            'Cannot derive the area-table animal name from %s.', animal_batch);
    end
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);
lpag_rows = ismember(area_limit.area, 'LPAG');
lpag_limits = area_limit{lpag_rows, {'ProbeNum', 'depth_start', 'depth_end'}};
if isempty(lpag_limits)
    error('GENERATE_PLAYBOUT_ERP:areas', ...
        'No LPAG limits found for %s.', this_animal);
end

% chanMap is the LFP row, the same index GENERATE_PSTH_PLAY_BOUT uses.
selected = select_mid_lpag_channels( ...
    loaded.xcoords, loaded.ycoords, loaded.chanMap, lpag_limits, npx_type);
for site_n = 1:numel(selected)
    if selected(site_n).lfp_row > n_lfp_rows
        error('GENERATE_PLAYBOUT_ERP:channel', ...
            'Selected LPAG channel %g is outside the LFP matrix.', ...
            selected(site_n).lfp_row);
    end
end
end
