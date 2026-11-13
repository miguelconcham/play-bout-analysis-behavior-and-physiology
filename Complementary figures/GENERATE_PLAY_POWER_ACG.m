function session = GENERATE_PLAY_POWER_ACG(npx_data_dir, Hd_delta, options)
%GENERATE_PLAY_POWER_ACG  Play, matched non-play, and power-tail autocorrelograms.
%
% One-second epochs are cut from play bouts and from Partners-session time
% outside those bouts. Delta power is the mean squared 1-5 Hz envelope on
% the neuron's channel, the same quantity as the phase-locking power
% control. Play epochs are split at the 25th and 75th percentiles of that
% power. Outside epochs are matched one-to-one to play epochs on the
% decibel scale. Each neuron then gets an overlap-weighted autocorrelogram
% in the low-power play tail, the high-power play tail, and in the
% power-matched outside epochs.
%
% Hd_delta is the bandpass filter used by Estimate_play_power_acg.
% options.epoch_duration, acg_lag, acg_bin, low_percent, high_percent,
% caliper_sd, and areas_to_keep override the defaults.

if nargin < 1 || isempty(npx_data_dir)
    error('GENERATE_PLAY_POWER_ACG:folder', 'Provide the NPX session folder.');
end
if nargin < 2 || isempty(Hd_delta)
    error('GENERATE_PLAY_POWER_ACG:filter', ...
        'Provide the delta bandpass filter used to measure epoch power.');
end
if nargin < 3 || isempty(options)
    options = struct();
end
options = fill_generator_options(options);

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
synch_directory = fullfile(data_root, 'Synch data');
behavior_data = fullfile(data_root, 'Behavior backups');
area_limit_file = fullfile(data_root, 'Area_limits_GoodLooking.xlsx');
play_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'};
area_name = 'PAG';

[~, animal_code] = fileparts(npx_data_dir);
animal_parts = strsplit(animal_code, ' ');
if numel(animal_parts) < 3
    error('GENERATE_PLAY_POWER_ACG:name', ...
        'Session folder "%s" does not have batch, date, and animal tokens.', animal_code);
end
animal_batch = animal_parts{1};
repeated_animal = animal_parts{3};

behavior = load_synced_behavior(synch_directory, behavior_data, animal_code);
animal_types = unique(behavior.Animal);
animal_types(ismember(animal_types, {'Session_structure'})) = [];
config.Behavior = behavior;
config.repeated_animal = repeated_animal;
config.animal_types = animal_types;
config.play_behaviors = play_behaviors;
config.beh_bin = 0.01;
config.conv_length = 1;
config.behavior_window = 0;
play_bouts = play_bout_intervals(play_bout(config));

load(fullfile(npx_data_dir, ['chann_map_', area_name, '.mat']), ...
    'chanMap', 'xcoords', 'ycoords');
[areas_by_channel, ~, npx_type] = map_channels( ...
    chanMap, xcoords, ycoords, area_limit_file, animal_batch, repeated_animal);

spike_times = double(readNPY(fullfile(npx_data_dir, ...
    ['spike_times_', area_name, '.npy']))) / 30000;
spike_clusters = double(readNPY(fullfile(npx_data_dir, ...
    ['spike_clusters_', area_name, '.npy'])));
cluster_info = readtable(fullfile(npx_data_dir, ...
    ['cluster_info_', area_name, '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
if ismember('group', cluster_info.Properties.VariableNames)
    cluster_info = cluster_info(ismember(cluster_info.group, {'mua', 'good'}), :);
end
id_col = cluster_id_column(cluster_info);
cluster_info.area = areas_by_channel(cluster_info.ch + 1);
cluster_info.ch = cluster_info.ch + 1;

[lfp, n_lfp_samples, lfp_file] = open_lfp(npx_data_dir, area_name, npx_type);
lfp_time_limit = [1 / options.sample_rate, n_lfp_samples / options.sample_rate];
sessions = behavior(strcmp(behavior.Type2, options.session_type), :);
if isempty(sessions)
    error('GENERATE_PLAY_POWER_ACG:sessions', ...
        'No "%s" intervals in %s.', options.session_type, animal_code);
end
session_windows = [sessions.Start, sessions.End];
epochs = play_and_outside_epochs(session_windows, play_bouts, ...
    options.epoch_duration, lfp_time_limit);
if isempty(epochs.play) || isempty(epochs.outside)
    error('GENERATE_PLAY_POWER_ACG:epochs', ...
        'Play epochs (%d) or outside epochs (%d) are missing for %s.', ...
        size(epochs.play, 1), size(epochs.outside, 1), animal_code);
end

bin_edges = (-options.acg_lag):options.acg_bin:options.acg_lag;
centers = bin_edges(1:end - 1) + options.acg_bin / 2;
names = play_power_condition_names();
kept = false(height(cluster_info), 1);
kept_area = cell(height(cluster_info), 1);
for nn = 1:height(cluster_info)
    kept_area{nn} = canonical_kept_area(cluster_info.area{nn}, options.areas_to_keep);
    kept(nn) = ~isempty(kept_area{nn});
end
if ~any(kept)
    error('GENERATE_PLAY_POWER_ACG:areas', ...
        'No good or mua neurons in %s for %s.', strjoin(options.areas_to_keep, ', '), animal_code);
end

channels_with_spikes = unique(cluster_info.ch(kept))';
channel_store = repmat(empty_channel(), numel(channels_with_spikes), 1);
neuron_store = repmat(empty_neuron(names), sum(kept), 1);
neuron_count = 0;

for cc = 1:numel(channels_with_spikes)
    channel = channels_with_spikes(cc);
    fprintf('Play-power ACG: %s channel %d (%d/%d)\n', ...
        animal_code, channel, cc, numel(channels_with_spikes));
    if npx_type == 2
        raw_lfp = double(lfp(channel, :));
    else
        raw_lfp = double(read_npx1_lfp_channel(lfp_file, channel, n_lfp_samples));
    end
    filtered = filtfilt(Hd_delta.Coefficients, 1, raw_lfp);
    amplitude = abs(hilbert(filtered));
    clear raw_lfp filtered
    play_power = epoch_mean_square(amplitude, options.sample_rate, epochs.play);
    outside_power = epoch_mean_square(amplitude, options.sample_rate, epochs.outside);
    clear amplitude
    assigned = assign_play_power_conditions(play_power, outside_power, options);
    channel_store(cc).channel = channel;
    channel_store(cc).play_power = play_power;
    channel_store(cc).outside_power = outside_power;
    channel_store(cc).conditions = assigned;

    on_channel = find(kept & cluster_info.ch == channel);
    gap_of = pair_gap_by_condition(assigned);
    for k = 1:numel(on_channel)
        nn = on_channel(k);
        this_id = double(id_col(nn));
        unit_spikes = spike_times(spike_clusters == this_id);
        blocks = struct();
        for c = 1:numel(names)
            name = names{c};
            idx = assigned.(name);
            if condition_uses_outside(name)
                these_epochs = epochs.outside(idx, :);
                these_power = outside_power(idx);
            else
                these_epochs = epochs.play(idx, :);
                these_power = play_power(idx);
            end
            block = condition_window_acg(unit_spikes, these_epochs, bin_edges, these_power, ...
                options.acg_method, options.acg_count);
            block.median_pair_gap_db = gap_of.(name);
            blocks.(name) = block;
        end
        neuron_count = neuron_count + 1;
        neuron_store(neuron_count).cluster_id = this_id;
        neuron_store(neuron_count).channel = channel;
        neuron_store(neuron_count).area = kept_area{nn};
        neuron_store(neuron_count).conditions = blocks;
    end
end
clear lfp

session.animal_code = animal_code;
session.epoch_duration = options.epoch_duration;
session.acg_lag = options.acg_lag;
session.acg_bin = options.acg_bin;
session.bin_edges = bin_edges;
session.centers = centers;
session.sample_rate = options.sample_rate;
session.delta_band = options.delta_band;
session.low_percent = options.low_percent;
session.high_percent = options.high_percent;
session.caliper_sd = options.caliper_sd;
session.acg_method = options.acg_method;
session.acg_count = options.acg_count;
session.areas_to_keep = options.areas_to_keep;
session.session_type = options.session_type;
session.play_epochs = epochs.play;
session.outside_epochs = epochs.outside;
session.channels = channel_store;
session.neurons = neuron_store(1:neuron_count);
end


function options = fill_generator_options(options)
if ~isfield(options, 'epoch_duration') || isempty(options.epoch_duration)
    options.epoch_duration = 1;
end
if ~isfield(options, 'acg_lag') || isempty(options.acg_lag)
    options.acg_lag = 1;
end
if ~isfield(options, 'acg_bin') || isempty(options.acg_bin)
    options.acg_bin = 0.01;
end
if ~isfield(options, 'sample_rate') || isempty(options.sample_rate)
    options.sample_rate = 2500;
end
if ~isfield(options, 'delta_band') || isempty(options.delta_band)
    options.delta_band = [1 5];
end
if ~isfield(options, 'low_percent') || isempty(options.low_percent)
    options.low_percent = 25;
end
if ~isfield(options, 'high_percent') || isempty(options.high_percent)
    options.high_percent = 75;
end
if ~isfield(options, 'caliper_sd') || isempty(options.caliper_sd)
    options.caliper_sd = 0.25;
end
if ~isfield(options, 'areas_to_keep') || isempty(options.areas_to_keep)
    options.areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
end
if ~isfield(options, 'session_type') || isempty(options.session_type)
    options.session_type = 'Partners session';
end
if ~isfield(options, 'acg_method') || isempty(options.acg_method)
    options.acg_method = 'overlap';
end
options.acg_method = canonical_acg_method(options.acg_method);
if ~isfield(options, 'acg_count') || isempty(options.acg_count)
    options.acg_count = 'pairs';
end
options.acg_count = canonical_count_method(options.acg_count);
if strcmp(options.acg_method, 'eligible') && ~strcmp(options.acg_count, 'pairs')
    error('GENERATE_PLAY_POWER_ACG:count', ...
        ['acg_count is %s and acg_method is eligible. Eligible-spike ' ...
        'autocorrelograms keep their own count.'], options.acg_count);
end
end


function behavior = load_synced_behavior(synch_directory, behavior_data, animal_code)
load(fullfile(synch_directory, animal_code, 'synch_model_video2NPX.mat'), ...
    'synch_model_video2NPX');
behavior = readtable(fullfile(behavior_data, [animal_code, '.txt']));
behavior(:, 2) = [];
behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
behavior.Type2 = behavior.Type;
behavior.Type2(ismember(behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
behavior.Type2(ismember(behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
behavior.Type2(ismember(behavior.Type2, {''})) = {'Other'};
behavior(ismember(behavior.Animal, {'Reversal'}), :) = [];
behavior.Start = predict(synch_model_video2NPX, behavior.Start);
behavior.End = predict(synch_model_video2NPX, behavior.End);
end


function intervals = play_bout_intervals(play_bouts)
if istable(play_bouts)
    intervals = play_bouts{:, 1:2};
elseif isempty(play_bouts)
    intervals = zeros(0, 2);
else
    intervals = play_bouts(:, 1:2);
end
intervals = sortrows(double(intervals), 1);
end


function id_col = cluster_id_column(ci)
if ismember('cluster_id', ci.Properties.VariableNames)
    id_col = ci.cluster_id;
elseif ismember('id', ci.Properties.VariableNames)
    id_col = ci.id;
else
    error('GENERATE_PLAY_POWER_ACG:id', 'Cluster table has neither cluster_id nor id.');
end
end


function tf = condition_uses_outside(name)
tf = strncmp(name, 'outside_', 8);
end


function gap_of = pair_gap_by_condition(assigned)
gap_of.play_low = NaN;
gap_of.play_high = NaN;
gap_of.play_matched = assigned.gap_db.play_matched;
gap_of.outside_matched = assigned.gap_db.play_matched;
gap_of.play_low_matched = assigned.gap_db.play_low_matched;
gap_of.outside_low_matched = assigned.gap_db.play_low_matched;
gap_of.play_high_matched = assigned.gap_db.play_high_matched;
gap_of.outside_high_matched = assigned.gap_db.play_high_matched;
end


function channel = empty_channel()
channel = struct( ...
    'channel', [], ...
    'play_power', [], ...
    'outside_power', [], ...
    'conditions', []);
end


function neuron = empty_neuron(names)
blocks = struct();
for c = 1:numel(names)
    blocks.(names{c}) = [];
end
neuron = struct( ...
    'cluster_id', [], ...
    'channel', [], ...
    'area', '', ...
    'conditions', blocks);
end


function [lfp, n_samples, dat_file] = open_lfp(npx_data_dir, area_name, npx_type)
mat_file = fullfile(npx_data_dir, ['LFP_', area_name, '.mat']);
dat_file = fullfile(npx_data_dir, ['LFP_', area_name, '.dat']);
lfp = [];
if npx_type == 2
    if exist(mat_file, 'file') ~= 2
        error('GENERATE_PLAY_POWER_ACG:lfp', ...
            'Session was identified as NPX2, but %s is missing.', mat_file);
    end
    loaded = load(mat_file, 'LFP');
    lfp = loaded.LFP;
    clear loaded
    if isempty(lfp) || size(lfp, 1) ~= 384
        error('GENERATE_PLAY_POWER_ACG:lfp', ...
            'NPX2 LFP must be a nonempty 384-by-N array in %s.', mat_file);
    end
    n_samples = size(lfp, 2);
    dat_file = '';
    return
end
if exist(dat_file, 'file') ~= 2
    error('GENERATE_PLAY_POWER_ACG:lfp', ...
        'Session was identified as NPX1, but %s is missing.', dat_file);
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (384 * 2);
if n_samples < 1 || n_samples ~= floor(n_samples)
    error('GENERATE_PLAY_POWER_ACG:lfp', ...
        'NPX1 LFP binary is empty or not divisible by 384 int16 channels.');
end
end


function [areas_by_channel, channel_map, npx_type] = map_channels( ...
    chan_map, xcoords, ycoords, area_limit_file, animal_batch, repeated_animal)
hard_coded_x = [8 40; 258 290; 508 540; 758 790];
chan_map = chan_map(:);
xcoords = xcoords(:);
ycoords = ycoords(:);
npx_type = 1 + any(ismember(xcoords, hard_coded_x(:)));
if npx_type == 1 && ~ismember(192, chan_map)
    pos_191 = find(chan_map == 191, 1);
    pos_193 = find(chan_map == 193, 1);
    if ~isempty(pos_191) && ~isempty(pos_193) && pos_193 == pos_191 + 1
        xcoords = [xcoords; NaN];
        xcoords(pos_193 + 1:end) = xcoords(pos_193:end - 1);
        xcoords(pos_193) = 43;
        ycoords = [ycoords; NaN];
        ycoords(pos_193 + 1:end) = ycoords(pos_193:end - 1);
        ycoords(pos_193) = 1900;
        chan_map = [chan_map; NaN];
        chan_map(pos_193 + 1:end) = chan_map(pos_193:end - 1);
        chan_map(pos_193) = 192;
    else
        error('GENERATE_PLAY_POWER_ACG:map', ...
            'Inconsistent NPX1 channel map: channel 192 is missing.');
    end
end
areas_by_channel = repmat({''}, 384, 1);
channel_map = nan(384, 2);
area_limit = readtable(area_limit_file);
if strcmp(repeated_animal, 'Single2')
    animal_key = ['Batch', animal_batch(2), repeated_animal];
else
    animal_key = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, animal_key), :);

for ii = 1:numel(chan_map)
    channel = chan_map(ii);
    if ~isfinite(channel) || channel < 1 || channel > 384
        continue
    end
    channel_map(channel, :) = [xcoords(ii), ycoords(ii)];
    if npx_type == 1
        match = ycoords(ii) >= area_limit.depth_start & ...
            ycoords(ii) < area_limit.depth_end + 1 & ...
            ismember(area_limit.Probe_Area, {'PAG'});
    else
        probe = find(any(ismember(hard_coded_x, xcoords(ii)), 2), 1);
        match = ycoords(ii) >= area_limit.depth_start & ...
            ycoords(ii) < area_limit.depth_end + 1 & ...
            area_limit.ProbeNum == probe & ismember(area_limit.Probe_Area, {'PAG'});
    end
    if any(match)
        areas_by_channel{channel} = area_limit.area{find(match, 1)};
    end
end
end


function channel_data = read_npx1_lfp_channel(dat_file, channel, n_samples)
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('GENERATE_PLAY_POWER_ACG:lfp', 'Could not open NPX1 LFP binary: %s', dat_file);
end
cleanup = onCleanup(@() fclose(file_id)); %#ok<NASGU>
status = fseek(file_id, 2 * (channel - 1), 'bof');
if status ~= 0
    error('GENERATE_PLAY_POWER_ACG:lfp', 'Could not seek to channel %d in %s.', channel, dat_file);
end
[channel_data, count] = fread(file_id, n_samples, '*int16', 2 * (384 - 1));
if count ~= n_samples
    error('GENERATE_PLAY_POWER_ACG:lfp', ...
        'Channel %d in %s ended after %d of %d samples.', channel, dat_file, count, n_samples);
end
channel_data = double(channel_data(:))';
end
