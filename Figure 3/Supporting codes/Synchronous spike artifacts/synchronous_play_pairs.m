function selection = synchronous_play_pairs(npx_data_dir, options)
%SYNCHRONOUS_PLAY_PAIRS  Play-bout coincident pairs in one NPX session.
%
% Good units in DLPAG, LPAG, VLPAG, DR, and SupCol. Coincidence uses the
% same bout window and 5 ms tolerance as Estiamte_phase_coincidence.m.
% A neuron needs at least 1000 spikes before it can enter a pair. Pairs
% are added until each area has at least 90% of those neurons, and only
% when the pair brings in a neuron that is not yet observed. Every
% coincident spike is kept. Spike times in seconds are
% double(readNPY)/30000. spike_samples are those indices, shifted by +1
% when the file is 0-based, so sample 1 is the first frame of the
% continuous file.
%
% selection.pairs      one struct per kept pair, with sync_index into the
%                      sorted spike list of selection.neurons
% selection.neurons    only the neurons that appear in the kept pairs.
%                      interval_index lists spikes inside the bout windows.
% selection.coverage   how many good units in the areas were observed.
%
% options.spikes_only skips pair ranking and returns selection.all_neurons
% and selection.bout_starts. The waveform script uses that to build
% autocorrelograms without reading the continuous file again.

if nargin < 2 || isempty(options)
    options = struct();
end
options = sync_analysis_options(options);
[~, session_id] = fileparts(npx_data_dir);
selection.session_id = session_id;
selection.pairs = empty_pairs();
selection.neurons = empty_neurons();
selection.coverage = [];
selection.all_neurons = empty_neurons();
selection.bout_starts = [];
if options.max_pairs < 1
    [selection.pairs, selection.neurons, selection.coverage] = ...
        select_synchronous_pairs(empty_neurons(), [], options);
    return
end
if exist(npx_data_dir, 'dir') ~= 7
    error('synchronous_play_pairs:session', ...
        'Session folder was not found: %s', npx_data_dir);
end

name_parts = regexp(strtrim(session_id), '\s+', 'split');
name_parts = name_parts(~cellfun('isempty', name_parts));
if numel(name_parts) < 3
    error('synchronous_play_pairs:session', ...
        'Session name "%s" is not "{animal} {date} {Single|Dual|Single2}".', ...
        session_id);
end

bout_starts = play_bout_starts(session_id, name_parts, options);
if isempty(bout_starts)
    [selection.pairs, selection.neurons, selection.coverage] = ...
        select_synchronous_pairs(empty_neurons(), [], options);
    return
end
neurons = neurons_in_session(npx_data_dir, name_parts, options);
selection.all_neurons = neurons;
selection.bout_starts = bout_starts;
if isfield(options, 'spikes_only') && options.spikes_only
    return
end
[selection.pairs, selection.neurons, selection.coverage] = select_synchronous_pairs( ...
    neurons, bout_starts, options);
selection.session_id = session_id;
end

function starts = play_bout_starts(session_id, name_parts, options)
model_file = fullfile(options.data_root, 'Synch data', session_id, ...
    'synch_model_video2NPX.mat');
behavior_file = fullfile(options.data_root, 'Behavior backups', ...
    [session_id '.txt']);
if exist(model_file, 'file') ~= 2 || exist(behavior_file, 'file') ~= 2
    error('synchronous_play_pairs:behavior', ...
        'Missing the video-to-NPX model or the behavior file for %s.', session_id);
end
loaded = load(model_file, 'synch_model_video2NPX');
behavior = readtable(behavior_file);
behavior(:, 2) = [];
behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
behavior.Type2 = behavior.Type;
behavior.Type2(ismember(behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
behavior.Type2(ismember(behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
behavior.Type2(ismember(behavior.Type2, {''})) = {'Other'};
behavior(ismember(behavior.Animal, 'Reversal'), :) = [];
behavior.Start = predict(loaded.synch_model_video2NPX, behavior.Start);
behavior.End = predict(loaded.synch_model_video2NPX, behavior.End);
animal_types = unique(behavior.Animal);
animal_types(ismember(animal_types, {'Session_structure'})) = [];
config.Behavior = behavior;
config.repeated_animal = name_parts{3};
config.animal_types = animal_types;
config.play_behaviors = options.play_behaviors;
config.beh_bin = 0.01;
config.conv_length = 1;
config.behavior_window = 0;
play_bouts = play_bout(config);
if isempty(play_bouts)
    starts = [];
    return
end
if isa(play_bouts, 'table')
    play_bouts = play_bouts{:, 1:min(2, size(play_bouts, 2))};
end
starts = play_bouts(:, 1);
end

function neurons = neurons_in_session(session_dir, name_parts, options)
neurons = empty_neurons();
areas_by_channel = channel_areas(session_dir, name_parts, options);
[spike_times, spike_samples, spike_clusters, cluster_info] = ...
    load_sorting(session_dir);
ids = table_column(cluster_info, {'cluster_id', 'id'});
channels = table_column(cluster_info, {'ch', 'channel'});
groups = table_column(cluster_info, {'group', 'KSLabel'});
for unit_index = 1:numel(ids)
    if ~strcmp(group_label(groups, unit_index), 'good')
        continue
    end
    cluster_id = scalar_number(ids, unit_index);
    channel = scalar_number(channels, unit_index) + 1;
    if ~isfinite(cluster_id) || channel < 1 || channel > numel(areas_by_channel) ...
            || channel ~= round(channel)
        continue
    end
    area = areas_by_channel{channel};
    if ~any(strcmp(options.areas_to_keep, area))
        continue
    end
    keep = spike_clusters == cluster_id;
    times = spike_times(keep);
    samples = spike_samples(keep);
    finite_times = isfinite(times);
    times = times(finite_times);
    samples = samples(finite_times);
    if isempty(times)
        continue
    end
    neurons(end + 1).cluster_id = cluster_id; %#ok<AGROW>
    neurons(end).channel = channel;
    neurons(end).area = area;
    neurons(end).spike_times = times(:);
    neurons(end).spike_samples = samples(:);
end
end

function [spike_times, spike_samples, spike_clusters, cluster_info] = load_sorting(session_dir)
info_file = fullfile(session_dir, 'cluster_info_PAG.tsv');
suffix = '_PAG';
if exist(info_file, 'file') ~= 2
    info_file = fullfile(session_dir, 'cluster_info.tsv');
    suffix = '';
end
times_file = fullfile(session_dir, ['spike_times' suffix '.npy']);
cluster_file = fullfile(session_dir, ['spike_clusters' suffix '.npy']);
if exist(info_file, 'file') ~= 2 || exist(times_file, 'file') ~= 2 ...
        || exist(cluster_file, 'file') ~= 2
    error('synchronous_play_pairs:sorting', ...
        'Missing cluster_info, spike_times, or spike_clusters in %s.', session_dir);
end
cluster_info = readtable(info_file, 'FileType', 'text', 'Delimiter', '\t');
raw_times = double(readNPY(times_file));
raw_times = raw_times(:);
spike_samples = raw_times;
if ~isempty(raw_times) && min(raw_times) == 0
    spike_samples = raw_times + 1;
elseif any(raw_times < 0)
    error('synchronous_play_pairs:spikeTime', ...
        'spike_times must be non-negative sample indices.');
end
spike_times = raw_times / 30000;
spike_samples = round(spike_samples);
spike_clusters = double(readNPY(cluster_file));
spike_clusters = spike_clusters(:);
end

function areas_by_channel = channel_areas(session_dir, name_parts, options)
areas_by_channel = cell(384, 1);
for index = 1:numel(areas_by_channel)
    areas_by_channel{index} = '';
end
map_file = fullfile(session_dir, 'chann_map_PAG.mat');
if exist(map_file, 'file') ~= 2
    error('synchronous_play_pairs:channelMap', 'Missing chann_map_PAG.mat.');
end
loaded = load(map_file, 'chanMap', 'xcoords', 'ycoords');
chan_map = loaded.chanMap(:);
xcoords = loaded.xcoords(:);
ycoords = loaded.ycoords(:);
npx2_x = [8 40; 258 290; 508 540; 758 790];
npx_type = 1;
if any(ismember(xcoords, npx2_x))
    npx_type = 2;
elseif ~ismember(192, chan_map)
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
        error('synchronous_play_pairs:channelMap', 'Inconsistent channel map.');
    end
end

area_limit_file = fullfile(options.data_root, 'Area_limits_GoodLooking.xlsx');
area_limit = readtable(area_limit_file);
animal_batch = name_parts{1};
repeated_animal = name_parts{3};
if numel(animal_batch) < 2
    error('synchronous_play_pairs:animal', ...
        'Cannot read the batch digit from "%s".', animal_batch);
end
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    if numel(animal_batch) < 4
        error('synchronous_play_pairs:animal', ...
            'Cannot read the animal digit from "%s".', animal_batch);
    end
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
animal_names = table_column(area_limit, {'AnimalName', 'animalname'});
keep_animal = false(numel(animal_names), 1);
for row_index = 1:numel(animal_names)
    keep_animal(row_index) = strcmp(column_text(animal_names, row_index), this_animal);
end
area_limit = area_limit(keep_animal, :);
if isempty(area_limit)
    error('synchronous_play_pairs:area', 'No area limits for %s.', this_animal);
end
depth_start = double(table_column(area_limit, {'depth_start', 'depthstart'}));
depth_end = double(table_column(area_limit, {'depth_end', 'depthend'}));
area_names = table_column(area_limit, {'area', 'Area'});
probe_area = table_column(area_limit, {'Probe_Area', 'probearea'});
probe_num = double(table_column(area_limit, {'ProbeNum', 'probenum'}));
depth_start = depth_start(:);
depth_end = depth_end(:);
probe_num = probe_num(:);

for site_index = 1:numel(chan_map)
    channel = chan_map(site_index);
    if ~isfinite(channel) || channel < 1 || channel > 384 || channel ~= round(channel)
        continue
    end
    if npx_type == 2
        probe_n = find(any(ismember(npx2_x, xcoords(site_index)), 2), 1);
        if isempty(probe_n)
            probe_ok = false(size(probe_num));
        else
            probe_ok = probe_num == probe_n;
        end
    else
        probe_ok = true(size(probe_num));
    end
    on_probe = false(numel(probe_area), 1);
    for row_index = 1:numel(probe_area)
        on_probe(row_index) = probe_ok(row_index) ...
            && strcmpi(column_text(probe_area, row_index), 'PAG');
    end
    depth_ok = ycoords(site_index) >= depth_start ...
        & ycoords(site_index) < depth_end + 1 & on_probe;
    areas_by_channel{channel} = first_target_area(area_names, depth_ok);
end
end

function label = first_target_area(area_names, mask)
label = '';
targets = {'DLPAG', 'LPAG', 'VLPAG', 'DR', 'SupCol'};
hits = find(mask);
for index = 1:numel(hits)
    candidate = column_text(area_names, hits(index));
    if any(strcmp(targets, candidate))
        label = candidate;
        return
    end
end
end

function pairs = empty_pairs()
pairs = struct('neuron_index_1', {}, 'neuron_index_2', {}, ...
    'cluster_id_1', {}, 'cluster_id_2', {}, ...
    'channel_1', {}, 'channel_2', {}, ...
    'area_1', {}, 'area_2', {}, ...
    'n_spikes_1', {}, 'n_spikes_2', {}, ...
    'n_interval_1', {}, 'n_interval_2', {}, ...
    'n_sync_1', {}, 'n_sync_2', {}, ...
    'sync_index_1', {}, 'sync_index_2', {}, ...
    'rank_spikes', {});
end

function neurons = empty_neurons()
neurons = struct('cluster_id', {}, 'channel', {}, 'area', {}, ...
    'spike_times', {}, 'spike_samples', {}, 'interval_index', {});
end

function values = table_column(table_data, names)
vars = table_data.Properties.VariableNames;
keys = cell(1, numel(vars));
for index = 1:numel(vars)
    keys{index} = lower(regexprep(vars{index}, '[^A-Za-z0-9]', ''));
end
for name_index = 1:numel(names)
    key = lower(regexprep(names{name_index}, '[^A-Za-z0-9]', ''));
    hit = find(strcmp(keys, key), 1);
    if ~isempty(hit)
        values = table_data.(vars{hit});
        return
    end
end
error('synchronous_play_pairs:column', ...
    'Missing a column named %s.', strjoin(names, ' or '));
end

function text = column_text(values, index)
if iscell(values)
    value = values{index};
elseif exist('isstring', 'builtin') && isstring(values)
    value = char(values(index));
elseif exist('iscategorical', 'builtin') && iscategorical(values)
    value = char(values(index));
elseif ischar(values)
    value = strtrim(values(index, :));
else
    value = values(index);
end
if iscell(value)
    if numel(value) == 1
        value = value{1};
    else
        text = '';
        return
    end
end
if exist('isstring', 'builtin') && isstring(value)
    value = char(value);
end
if isnumeric(value)
    text = '';
    return
end
text = strtrim(char(value));
end

function label = group_label(groups, index)
text = column_text(groups, index);
label = lower(strtrim(text));
end

function number = scalar_number(values, index)
if iscell(values)
    number = double(values{index});
else
    number = double(values(index));
end
end
