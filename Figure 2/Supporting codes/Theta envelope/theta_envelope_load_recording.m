function loaded = theta_envelope_load_recording(current_dir, options)
%THETA_ENVELOPE_LOAD_RECORDING  Mid-LPAG LFP, channel map, and synced behavior.
%
% Behavior times are converted with synch_model_video2NPX, the same clock as
% sample k at time k/sampling_rate. Partner-session rows stay in the table;
% partner_ranges are only the Session_structure / Partners session intervals.

if nargin < 2 || isempty(options)
    options = theta_envelope_delta_options();
end
[animal_code, table_name, repeated_animal] = recording_names(current_dir);
data_root = resolve_data_root(current_dir, options);
[xcoords, ycoords, chanMap] = load_channel_map(current_dir);
npx_type = neuropixels_type(xcoords);
area_limit = readtable(fullfile(data_root, 'Area_limits_GoodLooking.xlsx'));
area_limit = area_limit(ismember(area_limit.AnimalName, table_name), :);
lpag_limit = area_limit(ismember(area_limit.area, options.area_name), :);

loaded = struct();
loaded.animal_code = animal_code;
loaded.table_name = table_name;
loaded.repeated_animal = repeated_animal;
loaded.directory = current_dir;
loaded.npx_type = npx_type;
loaded.selected = empty_selection();
loaded.lfp = [];
loaded.behavior = [];
loaded.partner_ranges = zeros(0, 2);
loaded.sampling_rate = options.sampling_rate;
if isempty(lpag_limit)
    warning('theta_envelope_load_recording:NoLPAG', ...
        '%s: no %s rows for %s.', animal_code, options.area_name, table_name);
    return
end

limit_matrix = lpag_limit{:, {'ProbeNum', 'depth_start', 'depth_end'}};
loaded.selected = select_mid_lpag_channels( ...
    xcoords, ycoords, chanMap, limit_matrix, npx_type);
[loaded.lfp, ~] = read_selected_lfp(current_dir, [loaded.selected.lfp_row], ...
    options.n_lfp_channels);
[loaded.behavior, loaded.partner_ranges] = load_behavior( ...
    current_dir, data_root, animal_code);
if strcmp(options.interval_mode, 'full_recording')
    n_samples = size(loaded.lfp, 2);
    loaded.partner_ranges = [1, n_samples] / options.sampling_rate;
end
end

function [behavior_table, ranges] = load_behavior(current_dir, data_root, animal_code)
local_synch = fullfile(current_dir, 'synch_model_video2NPX.mat');
data_synch = fullfile(data_root, 'Synch data', animal_code, 'synch_model_video2NPX.mat');
if isfile(local_synch)
    synch_file = local_synch;
elseif isfile(data_synch)
    synch_file = data_synch;
else
    error('theta_envelope_load_recording:synch', ...
        'Synchronization model not found for %s.', animal_code);
end
loaded_synch = load(synch_file, 'synch_model_video2NPX');
synch_model = loaded_synch.synch_model_video2NPX;

local_behavior = dir(fullfile(current_dir, 'ELAN*'));
data_behavior = fullfile(data_root, 'Behavior backups', [animal_code, '.txt']);
if ~isempty(local_behavior)
    behavior_file = fullfile(current_dir, local_behavior(1).name);
elseif isfile(data_behavior)
    behavior_file = data_behavior;
else
    error('theta_envelope_load_recording:behavior', ...
        'Behavior file not found for %s.', animal_code);
end

behavior_table = readtable(behavior_file);
behavior_table(:, 2) = [];
behavior_table.Properties.VariableNames = ...
    {'Animal', 'Start', 'End', 'Length', 'Type'};
behavior_table.Start = predict(synch_model, behavior_table.Start);
behavior_table.End = predict(synch_model, behavior_table.End);
ranges = partner_windows_from_behavior(behavior_table);
end

function [channel_lfp, n_samples] = read_selected_lfp(current_dir, rows, n_lfp_channels)
rows = rows(:);
mat_file = fullfile(current_dir, 'LFP_PAG.mat');
dat_file = fullfile(current_dir, 'LFP_PAG.dat');
if isfile(mat_file)
    loaded = load(mat_file, 'LFP');
    if ~isfield(loaded, 'LFP')
        error('theta_envelope_load_recording:lfp', 'LFP_PAG.mat does not contain LFP.');
    end
    if any(rows < 1 | rows > size(loaded.LFP, 1) | rows ~= round(rows))
        error('theta_envelope_load_recording:lfp', ...
            'A selected LPAG row is outside LFP_PAG.mat.');
    end
    channel_lfp = double(loaded.LFP(rows, :));
    n_samples = size(channel_lfp, 2);
    return
end
if ~isfile(dat_file)
    error('theta_envelope_load_recording:lfp', ...
        'Neither LFP_PAG.mat nor LFP_PAG.dat exists in %s.', current_dir);
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (n_lfp_channels * 2);
if n_samples ~= round(n_samples) || n_samples < 1
    error('theta_envelope_load_recording:lfp', ...
        'LFP_PAG.dat is not a multiple of %d int16 channels.', n_lfp_channels);
end
if any(rows < 1 | rows > n_lfp_channels | rows ~= round(rows))
    error('theta_envelope_load_recording:lfp', ...
        'A selected LPAG row is outside the %d-channel LFP file.', n_lfp_channels);
end
channel_lfp = zeros(numel(rows), n_samples);
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('theta_envelope_load_recording:lfp', 'Could not open %s.', dat_file);
end
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
for row_n = 1:numel(rows)
    if fseek(file_id, 2 * (rows(row_n) - 1), 'bof') ~= 0
        error('theta_envelope_load_recording:lfp', ...
            'Could not seek to LFP row %d.', rows(row_n));
    end
    [samples, count] = fread(file_id, n_samples, 'int16=>double', ...
        2 * (n_lfp_channels - 1));
    if count ~= n_samples
        error('theta_envelope_load_recording:lfp', ...
            'Read %d of %d samples on row %d.', count, n_samples, rows(row_n));
    end
    channel_lfp(row_n, :) = samples(:).';
end
end

function [xcoords, ycoords, chanMap] = load_channel_map(current_dir)
map_file = fullfile(current_dir, 'chann_map_PAG.mat');
if ~isfile(map_file)
    error('theta_envelope_load_recording:map', 'Channel map not found: %s', map_file);
end
loaded = load(map_file, 'xcoords', 'ycoords', 'chanMap');
xcoords = loaded.xcoords;
ycoords = loaded.ycoords;
chanMap = loaded.chanMap;
end

function npx_type = neuropixels_type(xcoords)
hard_coded_x = [8 40; 258 290; 508 540; 758 790];
if any(ismember(xcoords(:), hard_coded_x(:)))
    npx_type = 2;
else
    npx_type = 1;
end
end

function data_root = resolve_data_root(current_dir, options)
recording_data_root = fileparts(fileparts(fileparts(current_dir)));
if isfile(fullfile(recording_data_root, 'Area_limits_GoodLooking.xlsx'))
    data_root = recording_data_root;
else
    data_root = options.data_root;
end
end

function [animal_code, table_name, repeated_animal] = recording_names(current_dir)
[~, animal_code] = fileparts(current_dir);
parts = strsplit(strtrim(animal_code), ' ');
parts = parts(~cellfun(@isempty, parts));
if numel(parts) < 3
    error('theta_envelope_load_recording:name', ...
        'Unexpected Neuropixels folder name: %s', animal_code);
end
animal_batch = parts{1};
repeated_animal = parts{3};
if strcmp(repeated_animal, 'Single2')
    table_name = ['Batch', animal_batch(2), repeated_animal];
else
    if numel(animal_batch) < 4
        error('theta_envelope_load_recording:name', ...
            'Cannot derive the area-table animal name from %s.', animal_code);
    end
    table_name = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
end

function selected = empty_selection()
selected = struct('channel', {}, 'lfp_row', {}, 'map_index', {}, ...
    'probe_number', {}, 'x_um', {}, 'y_um', {}, 'depth_range', {});
end
