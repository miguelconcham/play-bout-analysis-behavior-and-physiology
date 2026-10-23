function recordings = GENERATE_DELTA_OSCILLATORY_EPISODES(current_dir, options)
%GENERATE_DELTA_OSCILLATORY_EPISODES  Delta episodes on the mid-LPAG channel.
%
% One recording folder in Data/NPX data/NPX raw data. For each LPAG depth
% range, the channel nearest the middle of that range is read from
% LFP_PAG. Partner-session blocks (Session_structure / Partners session)
% are filtered with the existing 1-5 Hz Hamming FIR. Candidate episodes
% are elevated Hilbert-amplitude runs with successive positive peaks.
%
% Set options.interval_mode to 'full_recording' to use the whole file
% instead of the partner-session rows. This does not cut play bouts.
%
% The saved episodes are candidates. A narrow bandpass can ring around a
% transient, so the minimally processed LFP is kept on the example snippet.

if nargin < 2 || isempty(options)
    options = struct();
end
options = delta_episode_options(options);
if isempty(options.filter_coefficients)
    options.filter_coefficients = delta_bandpass_coefficients(options);
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
    warning('GENERATE_DELTA_OSCILLATORY_EPISODES:NoLPAG', ...
        '%s: no %s rows for %s. Skipping.', animal_code, options.area_name, table_name);
    recordings = [];
    return
end
limit_matrix = lpag_limit{:, {'ProbeNum', 'depth_start', 'depth_end'}};
selected = select_mid_lpag_channels(xcoords, ycoords, chanMap, limit_matrix, npx_type);

[channel_lfp, n_samples] = read_selected_lfp(current_dir, [selected.lfp_row], ...
    options.n_lfp_channels);
if size(channel_lfp, 1) ~= numel(selected)
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
        'Expected %d LFP rows and read %d.', numel(selected), size(channel_lfp, 1));
end

sample_intervals = interval_sample_ranges(current_dir, data_root, animal_code, ...
    n_samples, sampling_rate, options);
if isempty(sample_intervals)
    warning('GENERATE_DELTA_OSCILLATORY_EPISODES:NoInterval', ...
        '%s: no continuous interval to analyze. Skipping.', animal_code);
    recordings = [];
    return
end

recordings = repmat(empty_recording(), numel(selected), 1);
for channel_n = 1:numel(selected)
    fprintf('%s  probe %g  LFP row %g\n', animal_code, ...
        selected(channel_n).probe_number, selected(channel_n).lfp_row);
    recordings(channel_n) = detect_one_channel( ...
        channel_lfp(channel_n, :), selected(channel_n), sample_intervals, ...
        animal_code, table_name, current_dir, npx_type, sampling_rate, options);
end
end

function recording = detect_one_channel(lfp_row, selected, sample_intervals, ...
        animal_code, table_name, current_dir, npx_type, sampling_rate, options)
recording = empty_recording();
recording.animal_name = animal_code;
recording.animal_area_table_name = table_name;
recording.neuropixels_directory = current_dir;
recording.npx_type = npx_type;
recording.lpag_channel = selected.channel;
recording.lpag_lfp_row = selected.lfp_row;
recording.lpag_map_index = selected.map_index;
recording.lpag_probe_number = selected.probe_number;
recording.lpag_x_um = selected.x_um;
recording.lpag_y_um = selected.y_um;
recording.lpag_depth_range = selected.depth_range;
recording.interval_mode = options.interval_mode;
recording.parameters = options;

n_intervals = numel(sample_intervals);
intervals = repmat(empty_interval(), n_intervals, 1);
episodes = empty_saved_episodes();
best_example = [];

min_filter_samples = 3 * (numel(options.filter_coefficients) - 1);
for interval_n = 1:n_intervals
    bounds = sample_intervals(interval_n);
    n_in_interval = bounds.end_sample - bounds.start_sample + 1;
    intervals(interval_n) = empty_interval();
    intervals(interval_n).index = interval_n;
    intervals(interval_n).behavior_start_s = bounds.behavior_start_s;
    intervals(interval_n).behavior_end_s = bounds.behavior_end_s;
    intervals(interval_n).start_sample = bounds.start_sample;
    intervals(interval_n).end_sample = bounds.end_sample;
    intervals(interval_n).start_time_s = bounds.start_sample / sampling_rate;
    intervals(interval_n).end_time_s = bounds.end_sample / sampling_rate;
    if n_in_interval > 0
        intervals(interval_n).duration_s = n_in_interval / sampling_rate;
    else
        intervals(interval_n).duration_s = 0;
    end

    if n_in_interval <= min_filter_samples
        intervals(interval_n).skipped = true;
        intervals(interval_n).skip_reason = 'shorter than the delta filter';
        warning('GENERATE_DELTA_OSCILLATORY_EPISODES:ShortInterval', ...
            '%s interval %d is too short to filter (%d samples).', ...
            animal_code, interval_n, max(0, bounds.end_sample - bounds.start_sample + 1));
        continue
    end

    segment = double(lfp_row(bounds.start_sample:bounds.end_sample));
    segment = segment - mean(segment);
    detected = detect_delta_episodes(segment, sampling_rate, options);
    intervals(interval_n).amplitude_threshold = detected.threshold;
    intervals(interval_n).amplitude_mean = detected.amplitude_mean;
    intervals(interval_n).amplitude_std = detected.amplitude_std;
    intervals(interval_n).amplitude_median = detected.amplitude_median;
    intervals(interval_n).n_threshold_runs = detected.n_threshold_runs;
    intervals(interval_n).n_episodes = detected.n_episodes;
    intervals(interval_n).edge_margin_seconds = detected.edge_margin_seconds;

    saved = globalize_episodes(detected.episodes, intervals(interval_n), sampling_rate);
    saved = saved(:);
    if isempty(episodes)
        episodes = saved;
    elseif ~isempty(saved)
        episodes = [episodes(:); saved]; %#ok<AGROW>
    end

    snippet = make_snippet(detected, intervals(interval_n), sampling_rate, ...
        options, best_saved_episode(saved));
    if isempty(best_example) || snippet.score > best_example.score
        best_example = snippet;
    end
end

recording.intervals = intervals;
recording.episodes = episodes;
if isempty(best_example)
    best_example = empty_example();
end
best_example = rmfield(best_example, 'score');
recording.example = best_example;
recording.summary = summarize_delta_episode_recording(recording);
end

function snippet = make_snippet(detected, interval, sampling_rate, options, episode)
n_samples = numel(detected.lfp);
if isempty(episode)
    local1 = min(n_samples, detected.edge_margin_samples + 1);
    local2 = min(n_samples, local1 + round(10 * sampling_rate) - 1);
    if local2 < local1
        local1 = 1;
        local2 = n_samples;
    end
    has_episode = false;
    episode_range = zeros(0, 2);
    peak_index = zeros(0, 1);
    duration_s = NaN;
    n_cycles = 0;
    median_period = NaN;
    period_iqr = NaN;
    score = -inf;
else
    pad = round(options.example_pad_seconds * sampling_rate);
    half_cap = round(options.example_max_seconds / 2 * sampling_rate);
    center = round((episode.start_sample + episode.end_sample) / 2);
    half_episode = ceil((episode.end_sample - episode.start_sample + 1) / 2);
    half = min(half_episode, half_cap);
    view1 = max(interval.start_sample, center - half - pad);
    view2 = min(interval.end_sample, center + half + pad);
    local1 = view1 - interval.start_sample + 1;
    local2 = view2 - interval.start_sample + 1;
    has_episode = true;
    snippet_length = local2 - local1 + 1;
    episode_range = [ ...
        max(1, episode.start_sample - view1 + 1), ...
        min(snippet_length, episode.end_sample - view1 + 1)];
    in_view = episode.peak_samples >= view1 & episode.peak_samples <= view2;
    peak_index = episode.peak_samples(in_view) - view1 + 1;
    duration_s = episode.duration_s;
    n_cycles = episode.n_cycles;
    median_period = episode.median_cycle_period_s;
    period_iqr = episode.cycle_period_iqr_s;
    score = n_cycles + 1e-3 * duration_s;
end

view_samples = (interval.start_sample + local1 - 1):(interval.start_sample + local2 - 1);
snippet = empty_example();
snippet.has_episode = has_episode;
snippet.interval_index = interval.index;
snippet.time_s = view_samples / sampling_rate;
snippet.lfp_minimal = single(detected.lfp(local1:local2));
snippet.lfp_delta = single(detected.filtered(local1:local2));
snippet.amplitude = single(detected.amplitude(local1:local2));
snippet.threshold = detected.threshold;
snippet.episode_sample_range = episode_range;
snippet.peak_index = peak_index(:);
snippet.episode_duration_s = duration_s;
snippet.n_cycles = n_cycles;
snippet.median_cycle_period_s = median_period;
snippet.cycle_period_iqr_s = period_iqr;
snippet.score = score;
end

function episode = best_saved_episode(episodes)
episode = [];
if isempty(episodes)
    return
end
score = [episodes.n_cycles] + 1e-3 * [episodes.duration_s];
[~, index] = max(score);
episode = episodes(index);
end

function saved = globalize_episodes(local_episodes, interval, sampling_rate)
saved = empty_saved_episodes();
for episode_n = 1:numel(local_episodes)
    local = local_episodes(episode_n);
    peak_samples = local.peak_samples + interval.start_sample - 1;
    saved(episode_n).interval_index = interval.index; %#ok<AGROW>
    saved(episode_n).start_sample = local.start_sample + interval.start_sample - 1;
    saved(episode_n).end_sample = local.end_sample + interval.start_sample - 1;
    saved(episode_n).start_time_s = saved(episode_n).start_sample / sampling_rate;
    saved(episode_n).end_time_s = saved(episode_n).end_sample / sampling_rate;
    saved(episode_n).duration_s = local.duration_s;
    saved(episode_n).n_peaks = local.n_peaks;
    saved(episode_n).n_cycles = local.n_cycles;
    saved(episode_n).peak_samples = peak_samples(:);
    saved(episode_n).peak_times_s = peak_samples(:) / sampling_rate;
    saved(episode_n).cycle_periods_s = local.cycle_periods_s(:);
    saved(episode_n).median_cycle_period_s = local.median_cycle_period_s;
    saved(episode_n).cycle_period_iqr_s = local.cycle_period_iqr_s;
    saved(episode_n).cycle_period_std_s = local.cycle_period_std_s;
    saved(episode_n).cycle_period_cv = local.cycle_period_cv;
    saved(episode_n).fraction_cycles_in_delta_band = local.fraction_cycles_in_delta_band;
end
saved = saved(:);
end

function ranges = interval_sample_ranges(current_dir, data_root, animal_code, ...
        n_samples, sampling_rate, options)
if strcmp(options.interval_mode, 'full_recording')
    ranges = struct('behavior_start_s', NaN, 'behavior_end_s', NaN, ...
        'start_sample', 1, 'end_sample', n_samples);
    return
end

[behavior_ranges, ~] = load_partner_sessions(current_dir, data_root, animal_code);
if isempty(behavior_ranges)
    ranges = struct('behavior_start_s', {}, 'behavior_end_s', {}, ...
        'start_sample', {}, 'end_sample', {});
    return
end
if size(behavior_ranges, 1) > 1 && ...
        any(behavior_ranges(2:end, 1) < behavior_ranges(1:end - 1, 2))
    warning('GENERATE_DELTA_OSCILLATORY_EPISODES:Overlap', ...
        '%s: partner-session intervals overlap. Episodes are detected in each row.', ...
        animal_code);
end

ranges = repmat(struct('behavior_start_s', NaN, 'behavior_end_s', NaN, ...
    'start_sample', 1, 'end_sample', 1), size(behavior_ranges, 1), 1);
for interval_n = 1:size(behavior_ranges, 1)
    ranges(interval_n).behavior_start_s = behavior_ranges(interval_n, 1);
    ranges(interval_n).behavior_end_s = behavior_ranges(interval_n, 2);
    ranges(interval_n).start_sample = max(1, ...
        ceil(behavior_ranges(interval_n, 1) * sampling_rate));
    ranges(interval_n).end_sample = min(n_samples, ...
        floor(behavior_ranges(interval_n, 2) * sampling_rate));
end
end

function [ranges, behavior_table] = load_partner_sessions(current_dir, data_root, animal_code)
local_synch = fullfile(current_dir, 'synch_model_video2NPX.mat');
data_synch = fullfile(data_root, 'Synch data', animal_code, 'synch_model_video2NPX.mat');
if isfile(local_synch)
    synch_file = local_synch;
elseif isfile(data_synch)
    synch_file = data_synch;
else
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:synch', ...
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
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:behavior', ...
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
[ranges, order] = sortrows(ranges, 1);
behavior_table = behavior_table(order, :);
if any(~isfinite(ranges), 'all') || any(ranges(:, 2) <= ranges(:, 1))
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:behavior', ...
        'Invalid partner-session times for %s.', animal_code);
end
end

function [channel_lfp, n_samples] = read_selected_lfp(current_dir, rows, n_lfp_channels)
rows = rows(:);
mat_file = fullfile(current_dir, 'LFP_PAG.mat');
dat_file = fullfile(current_dir, 'LFP_PAG.dat');
if isfile(mat_file)
    loaded = load(mat_file, 'LFP');
    if ~isfield(loaded, 'LFP')
        error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
            'LFP_PAG.mat does not contain LFP.');
    end
    if any(rows < 1 | rows > size(loaded.LFP, 1) | rows ~= round(rows))
        error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
            'A selected LPAG row is outside LFP_PAG.mat.');
    end
    channel_lfp = double(loaded.LFP(rows, :));
    n_samples = size(channel_lfp, 2);
    return
end
if ~isfile(dat_file)
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
        'Neither LFP_PAG.mat nor LFP_PAG.dat exists in %s.', current_dir);
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (n_lfp_channels * 2);
if n_samples ~= round(n_samples) || n_samples < 1
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
        'LFP_PAG.dat is not a multiple of %d int16 channels.', n_lfp_channels);
end
if any(rows < 1 | rows > n_lfp_channels | rows ~= round(rows))
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
        'A selected LPAG row is outside the %d-channel LFP file.', n_lfp_channels);
end
channel_lfp = zeros(numel(rows), n_samples);
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', 'Could not open %s.', dat_file);
end
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
for row_n = 1:numel(rows)
    if fseek(file_id, 2 * (rows(row_n) - 1), 'bof') ~= 0
        error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
            'Could not seek to LFP row %d.', rows(row_n));
    end
    [samples, count] = fread(file_id, n_samples, 'int16=>double', ...
        2 * (n_lfp_channels - 1));
    if count ~= n_samples
        error('GENERATE_DELTA_OSCILLATORY_EPISODES:lfp', ...
            'Read %d of %d samples on row %d.', count, n_samples, rows(row_n));
    end
    channel_lfp(row_n, :) = samples(:).';
end
end

function [xcoords, ycoords, chanMap] = load_channel_map(current_dir)
map_file = fullfile(current_dir, 'chann_map_PAG.mat');
if ~isfile(map_file)
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:map', ...
        'Channel map not found: %s', map_file);
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

function [animal_code, table_name] = recording_names(current_dir)
[~, animal_code] = fileparts(current_dir);
parts = strsplit(strtrim(animal_code), ' ');
parts = parts(~cellfun(@isempty, parts));
if numel(parts) < 3
    error('GENERATE_DELTA_OSCILLATORY_EPISODES:name', ...
        'Unexpected Neuropixels folder name: %s', animal_code);
end
animal_batch = parts{1};
repeated_animal = parts{3};
if strcmp(repeated_animal, 'Single2')
    table_name = ['Batch', animal_batch(2), repeated_animal];
else
    if numel(animal_batch) < 4
        error('GENERATE_DELTA_OSCILLATORY_EPISODES:name', ...
            'Cannot derive the area-table animal name from %s.', animal_code);
    end
    table_name = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
end

function recording = empty_recording()
recording = struct( ...
    'animal_name', '', ...
    'animal_area_table_name', '', ...
    'neuropixels_directory', '', ...
    'npx_type', NaN, ...
    'lpag_channel', NaN, ...
    'lpag_lfp_row', NaN, ...
    'lpag_map_index', NaN, ...
    'lpag_probe_number', NaN, ...
    'lpag_x_um', NaN, ...
    'lpag_y_um', NaN, ...
    'lpag_depth_range', [NaN NaN], ...
    'interval_mode', '', ...
    'parameters', struct(), ...
    'intervals', empty_interval(), ...
    'episodes', empty_saved_episodes(), ...
    'example', empty_example(), ...
    'summary', struct());
end

function interval = empty_interval()
interval = struct( ...
    'index', NaN, ...
    'behavior_start_s', NaN, ...
    'behavior_end_s', NaN, ...
    'start_sample', NaN, ...
    'end_sample', NaN, ...
    'start_time_s', NaN, ...
    'end_time_s', NaN, ...
    'duration_s', NaN, ...
    'amplitude_threshold', NaN, ...
    'amplitude_mean', NaN, ...
    'amplitude_std', NaN, ...
    'amplitude_median', NaN, ...
    'n_threshold_runs', 0, ...
    'n_episodes', 0, ...
    'edge_margin_seconds', NaN, ...
    'skipped', false, ...
    'skip_reason', '');
end

function episodes = empty_saved_episodes()
episodes = struct( ...
    'interval_index', {}, ...
    'start_sample', {}, ...
    'end_sample', {}, ...
    'start_time_s', {}, ...
    'end_time_s', {}, ...
    'duration_s', {}, ...
    'n_peaks', {}, ...
    'n_cycles', {}, ...
    'peak_samples', {}, ...
    'peak_times_s', {}, ...
    'cycle_periods_s', {}, ...
    'median_cycle_period_s', {}, ...
    'cycle_period_iqr_s', {}, ...
    'cycle_period_std_s', {}, ...
    'cycle_period_cv', {}, ...
    'fraction_cycles_in_delta_band', {});
end

function example = empty_example()
example = struct( ...
    'has_episode', false, ...
    'interval_index', NaN, ...
    'time_s', zeros(1, 0), ...
    'lfp_minimal', zeros(1, 0, 'single'), ...
    'lfp_delta', zeros(1, 0, 'single'), ...
    'amplitude', zeros(1, 0, 'single'), ...
    'threshold', NaN, ...
    'episode_sample_range', zeros(0, 2), ...
    'peak_index', zeros(0, 1), ...
    'episode_duration_s', NaN, ...
    'n_cycles', 0, ...
    'median_cycle_period_s', NaN, ...
    'cycle_period_iqr_s', NaN, ...
    'score', -inf);
end
