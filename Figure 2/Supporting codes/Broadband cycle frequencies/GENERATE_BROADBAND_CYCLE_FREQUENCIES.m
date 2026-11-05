function recordings = GENERATE_BROADBAND_CYCLE_FREQUENCIES(current_dir, options)
%GENERATE_BROADBAND_CYCLE_FREQUENCIES  Low-pass cycles for one recording.

if nargin < 2 || isempty(options)
    options = struct();
end
options = broadband_cycle_options(options);
if isempty(options.lowpass_coefficients)
    options.lowpass_coefficients = broadband_lowpass_coefficients(options);
end

[animal_code, table_name] = recording_names(current_dir);
data_root = resolve_data_root(current_dir, options);
sampling_rate = options.sampling_rate;
[xcoords, ycoords, chanMap] = load_channel_map(current_dir);
npx_type = neuropixels_type(xcoords);
area_limit = readtable(fullfile(data_root, 'Area_limits_GoodLooking.xlsx'));
area_limit = area_limit(ismember(area_limit.AnimalName, table_name), :);
lpag_limit = area_limit(ismember(area_limit.area, options.area_name), :);
if isempty(lpag_limit)
    warning('GENERATE_BROADBAND_CYCLE_FREQUENCIES:NoLPAG', ...
        '%s: no %s rows for %s. Skipping.', animal_code, ...
        options.area_name, table_name);
    recordings = [];
    return
end
limit_matrix = lpag_limit{:, {'ProbeNum', 'depth_start', 'depth_end'}};
selected = select_mid_lpag_channels(xcoords, ycoords, chanMap, limit_matrix, npx_type);
[channel_lfp, n_samples] = read_selected_lfp(current_dir, [selected.lfp_row], ...
    options.n_lfp_channels);
sample_intervals = interval_sample_ranges(current_dir, data_root, animal_code, ...
    n_samples, sampling_rate, options);
if isempty(sample_intervals)
    warning('GENERATE_BROADBAND_CYCLE_FREQUENCIES:NoInterval', ...
        '%s: no continuous interval to analyze. Skipping.', animal_code);
    recordings = [];
    return
end

recordings = repmat(empty_recording(), numel(selected), 1);
for channel_n = 1:numel(selected)
    recordings(channel_n) = detect_one_channel( ...
        channel_lfp(channel_n, :), selected(channel_n), sample_intervals, ...
        animal_code, table_name, npx_type, sampling_rate, options);
end
end

function recording = detect_one_channel(lfp_row, selected, sample_intervals, ...
        animal_code, table_name, npx_type, sampling_rate, options)
recording = empty_recording();
recording.animal_name = animal_code;
recording.animal_id = table_name;
recording.npx_type = npx_type;
recording.lpag_channel = selected.channel;
recording.lpag_lfp_row = selected.lfp_row;
recording.lpag_probe_number = selected.probe_number;
recording.lpag_y_um = selected.y_um;
recording.parameters = options;

frequency = zeros(0, 1);
amplitude = zeros(0, 1);
rdsym = zeros(0, 1);
monotonicity = zeros(0, 1);
accepted = false(0, 1);
best_amplitude = -Inf;
example = empty_example();
analyzed_seconds = 0;

for interval_n = 1:numel(sample_intervals)
    bounds = sample_intervals(interval_n);
    if bounds.end_sample < bounds.start_sample
        continue
    end
    segment = lfp_row(bounds.start_sample:bounds.end_sample);
    try
        detected = detect_broadband_cycles(segment, sampling_rate, options);
    catch
        continue
    end
    analyzed_seconds = analyzed_seconds + numel(segment) / sampling_rate;
    frequency = [frequency; detected.frequency_hz]; %#ok<AGROW>
    amplitude = [amplitude; detected.amplitude]; %#ok<AGROW>
    rdsym = [rdsym; detected.rdsym]; %#ok<AGROW>
    monotonicity = [monotonicity; detected.monotonicity]; %#ok<AGROW>
    accepted = [accepted; detected.accepted]; %#ok<AGROW>
    accepted_index = find(detected.accepted);
    if isempty(accepted_index)
        continue
    end
    [interval_amplitude, relative] = max(detected.amplitude(accepted_index));
    if interval_amplitude > best_amplitude
        best_amplitude = interval_amplitude;
        cycle_index = accepted_index(relative);
        example = example_around_cycle(detected, cycle_index, ...
            bounds.start_sample, sampling_rate, options.example_half_seconds);
    end
end

recording.analyzed_seconds = analyzed_seconds;
recording.n_candidate_cycles = numel(frequency);
recording.n_accepted_cycles = nnz(accepted);
recording.frequency_hz = frequency;
recording.amplitude = amplitude;
recording.rdsym = rdsym;
recording.monotonicity = monotonicity;
recording.accepted = accepted;
recording.example = example;
end

function example = example_around_cycle(detected, cycle_index, ...
        interval_start, sampling_rate, half_seconds)
center = detected.peak_sample(cycle_index);
half_samples = round(half_seconds * sampling_rate);
first_sample = max(1, center - half_samples);
last_sample = min(numel(detected.broad), center + half_samples);
example = empty_example();
example.time_s = ((interval_start + first_sample - 2): ...
    (interval_start + last_sample - 2)) / sampling_rate;
example.broad = detected.broad(first_sample:last_sample);
example.peak_time_s = (interval_start + detected.peak_sample(cycle_index) - 2) / ...
    sampling_rate;
example.trough_time_s = (interval_start + ...
    [detected.trough_start(cycle_index), detected.trough_end(cycle_index)] - 2) / ...
    sampling_rate;
example.frequency_hz = detected.frequency_hz(cycle_index);
end

function recording = empty_recording()
recording = struct( ...
    'animal_name', '', ...
    'animal_id', '', ...
    'npx_type', NaN, ...
    'lpag_channel', NaN, ...
    'lpag_lfp_row', NaN, ...
    'lpag_probe_number', NaN, ...
    'lpag_y_um', NaN, ...
    'parameters', struct(), ...
    'analyzed_seconds', 0, ...
    'n_candidate_cycles', 0, ...
    'n_accepted_cycles', 0, ...
    'frequency_hz', zeros(0, 1), ...
    'amplitude', zeros(0, 1), ...
    'rdsym', zeros(0, 1), ...
    'monotonicity', zeros(0, 1), ...
    'accepted', false(0, 1), ...
    'example', empty_example());
end

function example = empty_example()
example = struct('time_s', [], 'broad', [], 'peak_time_s', NaN, ...
    'trough_time_s', [NaN NaN], 'frequency_hz', NaN);
end

function ranges = interval_sample_ranges(current_dir, data_root, animal_code, ...
        n_samples, sampling_rate, options)
if strcmp(options.interval_mode, 'full_recording')
    ranges = struct('start_sample', 1, 'end_sample', n_samples);
    return
end
behavior_ranges = load_partner_sessions(current_dir, data_root, animal_code);
if isempty(behavior_ranges)
    ranges = struct('start_sample', {}, 'end_sample', {});
    return
end
ranges = repmat(struct('start_sample', 1, 'end_sample', 1), ...
    size(behavior_ranges, 1), 1);
for interval_n = 1:size(behavior_ranges, 1)
    ranges(interval_n).start_sample = max(1, ...
        ceil(behavior_ranges(interval_n, 1) * sampling_rate));
    ranges(interval_n).end_sample = min(n_samples, ...
        floor(behavior_ranges(interval_n, 2) * sampling_rate));
end
end

function ranges = load_partner_sessions(current_dir, data_root, animal_code)
local_synch = fullfile(current_dir, 'synch_model_video2NPX.mat');
data_synch = fullfile(data_root, 'Synch data', animal_code, ...
    'synch_model_video2NPX.mat');
if isfile(local_synch)
    synch_file = local_synch;
elseif isfile(data_synch)
    synch_file = data_synch;
else
    error('GENERATE_BROADBAND_CYCLE_FREQUENCIES:synch', ...
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
    error('GENERATE_BROADBAND_CYCLE_FREQUENCIES:behavior', ...
        'Behavior file not found for %s.', animal_code);
end
behavior_table = readtable(behavior_file);
behavior_table(:, 2) = [];
behavior_table.Properties.VariableNames = ...
    {'Animal', 'Start', 'End', 'Length', 'Type'};
behavior_table.Start = predict(synch_model, behavior_table.Start);
behavior_table.End = predict(synch_model, behavior_table.End);
session_mask = ismember(behavior_table.Animal, 'Session_structure') & ...
    ismember(behavior_table.Type, 'Partners session');
behavior_table = behavior_table(session_mask, :);
if isempty(behavior_table)
    ranges = zeros(0, 2);
    return
end
ranges = behavior_table{:, {'Start', 'End'}};
ranges = sortrows(ranges, 1);
end

function [channel_lfp, n_samples] = read_selected_lfp(current_dir, rows, n_lfp_channels)
rows = rows(:);
mat_file = fullfile(current_dir, 'LFP_PAG.mat');
dat_file = fullfile(current_dir, 'LFP_PAG.dat');
if isfile(mat_file)
    loaded = load(mat_file, 'LFP');
    channel_lfp = double(loaded.LFP(rows, :));
    n_samples = size(channel_lfp, 2);
    return
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (n_lfp_channels * 2);
channel_lfp = zeros(numel(rows), n_samples);
file_id = fopen(dat_file, 'r', 'ieee-le');
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
for row_n = 1:numel(rows)
    fseek(file_id, 2 * (rows(row_n) - 1), 'bof');
    channel_lfp(row_n, :) = fread(file_id, n_samples, 'int16=>double', ...
        2 * (n_lfp_channels - 1)).';
end
end

function [xcoords, ycoords, chanMap] = load_channel_map(current_dir)
loaded = load(fullfile(current_dir, 'chann_map_PAG.mat'), ...
    'xcoords', 'ycoords', 'chanMap');
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

function [animal_code, table_name] = recording_names(current_dir)
[~, animal_code] = fileparts(current_dir);
parts = strsplit(strtrim(animal_code), ' ');
parts = parts(~cellfun(@isempty, parts));
animal_batch = parts{1};
repeated_animal = parts{3};
if strcmp(repeated_animal, 'Single2')
    table_name = ['Batch', animal_batch(2), repeated_animal];
else
    table_name = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
end
