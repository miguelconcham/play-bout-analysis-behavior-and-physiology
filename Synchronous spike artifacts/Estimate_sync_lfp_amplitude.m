% Estimate_sync_lfp_amplitude
%
% Second pipeline for the synchronous-spike check. The events are the same
% play-bout coincidences as Estiamte_phase_coincidence.m. For each neuron in
% a pair the measured value is the mean raw-LFP amplitude, not a waveform
% distance.
%
% Raw LFP is the unfiltered LFP_PAG trace on that neuron's channel (cluster
% ch + 1). Amplitude at a spike is the Hilbert envelope, abs(hilbert(raw)),
% interpolated to the spike time. The template context is the mean envelope
% at every spike of the neuron. The test itself compares two means of the
% same count:
%   observed  mean envelope at every synchronous spike of this pair
%   null      500 means, each from that many spikes of the same neuron that
%             are not in the synchronous set
% A high percentile means the synchronous spikes sit at a larger LFP
% amplitude than random spikes of that neuron.
%
% Run from the Editor, or type the file name. Section 2 only prints the pair
% list. Set compute_amplitudes to false to stop there.
%
% Caps, so one neuron is not reused as a new file read: at most 2 pairs per
% neuron, 8 pairs per session, 80 pairs in the run. A neuron is loaded once
% per session and reused if it is in a second pair. Neurons need at least
% 1000 spikes. Every synchronous spike is kept.

%% Parameters
npx_raw_data = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\NPX data\NPX raw data'];
synch_directory = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\Synch data'];
behavior_directory = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\Behavior backups'];
area_limit_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\Area_limits_GoodLooking.xlsx'];
results_dir = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\Analysis results\' ...
    'synchronous spike waveforms'];

areas_to_keep = {'DLPAG', 'LPAG', 'VLPAG', 'DR', 'SupCol'};
play_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'};
hist_range = [-1.5 1];
time_precision = 0.005;
min_spikes = 1000;
min_sync = 15;
n_repeats = 500;
max_pairs_per_neuron = 2;
max_pairs_per_session = 8;
max_pairs = 80;
lfp_sample_rate = 2500;
n_lfp_channels = 384;
compute_amplitudes = true;

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

%% 1. Choose pairs
rng('shuffle');
seed_state = rng;
fprintf('Random seed: %d\n', seed_state.Seed);

chosen_pairs = struct('session_id', {}, 'session_dir', {}, ...
    'cluster_id', {}, 'partner_id', {}, 'channel', {}, 'area', {}, ...
    'n_spikes', {}, 'n_sync', {}, 'sync_times', {});
neuron_pair_count = struct('session_id', {}, 'cluster_id', {}, 'n_pairs', {});

if exist(npx_raw_data, 'dir') ~= 7
    fprintf('NPX raw data folder was not found:\n  %s\n', npx_raw_data);
else
    listing = dir(npx_raw_data);
    for folder_index = 1:numel(listing)
        if numel(chosen_pairs) >= max_pairs
            break
        end
        if ~listing(folder_index).isdir
            continue
        end
        session_id = listing(folder_index).name;
        if strcmp(session_id, '.') || strcmp(session_id, '..')
            continue
        end
        session_dir = fullfile(npx_raw_data, session_id);
        try
            session_pairs = eligible_pairs_in_session(session_dir, session_id, ...
                synch_directory, behavior_directory, area_limit_file, ...
                areas_to_keep, play_behaviors, hist_range, time_precision, ...
                min_spikes, min_sync);
        catch session_error
            fprintf('%s  skipped: %s\n', session_id, session_error.message);
            continue
        end
        if isempty(session_pairs)
            fprintf('%s  no eligible pair\n', session_id);
            continue
        end
        taken = 0;
        for pair_index = 1:numel(session_pairs)
            if taken >= max_pairs_per_session || numel(chosen_pairs) >= max_pairs
                break
            end
            pair = session_pairs(pair_index);
            if neuron_is_full(neuron_pair_count, session_id, pair.cluster_id_1, max_pairs_per_neuron) ...
                    || neuron_is_full(neuron_pair_count, session_id, pair.cluster_id_2, max_pairs_per_neuron)
                continue
            end
            chosen_pairs = append_neuron(chosen_pairs, session_id, session_dir, ...
                pair.cluster_id_1, pair.cluster_id_2, pair.channel_1, pair.area_1, ...
                pair.n_spikes_1, pair.n_sync_1, pair.sync_times_1);
            chosen_pairs = append_neuron(chosen_pairs, session_id, session_dir, ...
                pair.cluster_id_2, pair.cluster_id_1, pair.channel_2, pair.area_2, ...
                pair.n_spikes_2, pair.n_sync_2, pair.sync_times_2);
            neuron_pair_count = note_neuron(neuron_pair_count, session_id, pair.cluster_id_1);
            neuron_pair_count = note_neuron(neuron_pair_count, session_id, pair.cluster_id_2);
            taken = taken + 1;
            fprintf(['%s  %s ch %g (%d spikes, %d sync)  with  ' ...
                '%s ch %g (%d spikes, %d sync)\n'], ...
                session_id, pair.area_1, pair.channel_1, pair.n_spikes_1, pair.n_sync_1, ...
                pair.area_2, pair.channel_2, pair.n_spikes_2, pair.n_sync_2);
        end
    end
end
fprintf('Chosen %d neuron rows (%d pairs).\n', numel(chosen_pairs), numel(chosen_pairs) / 2);

%% 2. Mean raw-LFP amplitude at those spikes
amplitude_results = struct('session_id', {}, 'cluster_id', {}, 'partner_id', {}, ...
    'channel', {}, 'area', {}, 'n_spikes', {}, 'n_sync', {}, ...
    'mean_amplitude_all', {}, 'mean_amplitude_sync', {}, ...
    'null_mean_amplitude', {}, 'percentile', {});

if ~compute_amplitudes || isempty(chosen_pairs)
    return
end

session_ids = unique_session_ids(chosen_pairs);
for session_index = 1:numel(session_ids)
    session_id = session_ids{session_index};
    session_rows = chosen_pairs(strcmp({chosen_pairs.session_id}, session_id));
    session_dir = session_rows(1).session_dir;
    try
        bank = load_neuron_amplitudes(session_dir, session_rows, lfp_sample_rate, n_lfp_channels);
    catch load_error
        fprintf('%s  LFP skipped: %s\n', session_id, load_error.message);
        continue
    end
    for row_index = 1:numel(session_rows)
        row = session_rows(row_index);
        neuron = find_neuron(bank, row.cluster_id);
        if isempty(neuron)
            fprintf('%s  cluster %g  no LFP amplitude\n', session_id, row.cluster_id);
            continue
        end
        comparison = compare_mean_amplitude(neuron.amplitude, neuron.spike_times, ...
            row.sync_times, n_repeats);
        if ~comparison.ok
            fprintf('%s  cluster %g  %s\n', session_id, row.cluster_id, comparison.reason);
            continue
        end
        amplitude_results(end + 1).session_id = session_id; %#ok<SAGROW>
        amplitude_results(end).cluster_id = row.cluster_id;
        amplitude_results(end).partner_id = row.partner_id;
        amplitude_results(end).channel = row.channel;
        amplitude_results(end).area = row.area;
        amplitude_results(end).n_spikes = row.n_spikes;
        amplitude_results(end).n_sync = comparison.n_sync;
        amplitude_results(end).mean_amplitude_all = comparison.mean_all;
        amplitude_results(end).mean_amplitude_sync = comparison.mean_sync;
        amplitude_results(end).null_mean_amplitude = comparison.null_means;
        amplitude_results(end).percentile = comparison.percentile;
        fprintf('%s  cluster %g  sync mean %.3g  all-spike mean %.3g  percentile %.1f\n', ...
            session_id, row.cluster_id, comparison.mean_sync, ...
            comparison.mean_all, 100 * comparison.percentile);
    end
    clear bank
    try
        save_amplitude_results(results_dir, chosen_pairs, amplitude_results, seed_state.Seed);
    catch save_error
        fprintf('Results are in amplitude_results. Save failed: %s\n', save_error.message);
    end
end

%% Local functions

function pairs = eligible_pairs_in_session(session_dir, session_id, ...
    synch_directory, behavior_directory, area_limit_file, areas_to_keep, ...
    play_behaviors, hist_range, time_precision, min_spikes, min_sync)
pairs = struct('cluster_id_1', {}, 'cluster_id_2', {}, ...
    'channel_1', {}, 'channel_2', {}, 'area_1', {}, 'area_2', {}, ...
    'n_spikes_1', {}, 'n_spikes_2', {}, 'n_sync_1', {}, 'n_sync_2', {}, ...
    'sync_times_1', {}, 'sync_times_2', {}, 'rank_spikes', {});
name_parts = regexp(strtrim(session_id), '\s+', 'split');
name_parts = name_parts(~cellfun('isempty', name_parts));
if numel(name_parts) < 3
    error('Estimate_sync_lfp_amplitude:session', ...
        'Session name "%s" is not "{animal} {date} {Single|Dual|Single2}".', session_id);
end
bout_starts = play_bout_starts(session_id, name_parts, synch_directory, ...
    behavior_directory, play_behaviors);
if isempty(bout_starts)
    return
end
neurons = neurons_in_session(session_dir, name_parts, area_limit_file, areas_to_keep);
if numel(neurons) < 2
    return
end
for first = 1:numel(neurons) - 1
    for second = first + 1:numel(neurons)
        sync_1 = synchronous_times(neurons(first).spike_times, ...
            neurons(second).spike_times, bout_starts, hist_range, time_precision);
        sync_2 = synchronous_times(neurons(second).spike_times, ...
            neurons(first).spike_times, bout_starts, hist_range, time_precision);
        n1 = numel(neurons(first).spike_times);
        n2 = numel(neurons(second).spike_times);
        if n1 < min_spikes || n2 < min_spikes
            continue
        end
        if numel(sync_1) < min_sync || numel(sync_2) < min_sync
            continue
        end
        if n1 - numel(sync_1) < numel(sync_1) || n2 - numel(sync_2) < numel(sync_2)
            continue
        end
        pairs(end + 1).cluster_id_1 = neurons(first).cluster_id; %#ok<AGROW>
        pairs(end).cluster_id_2 = neurons(second).cluster_id;
        pairs(end).channel_1 = neurons(first).channel;
        pairs(end).channel_2 = neurons(second).channel;
        pairs(end).area_1 = neurons(first).area;
        pairs(end).area_2 = neurons(second).area;
        pairs(end).n_spikes_1 = n1;
        pairs(end).n_spikes_2 = n2;
        pairs(end).n_sync_1 = numel(sync_1);
        pairs(end).n_sync_2 = numel(sync_2);
        pairs(end).sync_times_1 = sync_1;
        pairs(end).sync_times_2 = sync_2;
        pairs(end).rank_spikes = min(n1, n2);
    end
end
if isempty(pairs)
    return
end
[~, order] = sort([pairs.rank_spikes], 'descend');
pairs = pairs(order);
end

function starts = play_bout_starts(session_id, name_parts, synch_directory, ...
    behavior_directory, play_behaviors)
model_file = fullfile(synch_directory, session_id, 'synch_model_video2NPX.mat');
behavior_file = fullfile(behavior_directory, [session_id '.txt']);
if exist(model_file, 'file') ~= 2 || exist(behavior_file, 'file') ~= 2
    error('Estimate_sync_lfp_amplitude:behavior', ...
        'Missing the video-to-NPX model or the behavior file.');
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
config.play_behaviors = play_behaviors;
config.beh_bin = 0.01;
config.conv_length = 1;
config.behavior_window = 0;
play_bouts = play_bout(config);
if isempty(play_bouts)
    starts = [];
    return
end
starts = play_bouts(:, 1);
end

function neurons = neurons_in_session(session_dir, name_parts, area_limit_file, areas_to_keep)
neurons = struct('cluster_id', {}, 'channel', {}, 'area', {}, 'spike_times', {});
[areas_by_channel, channel_ok] = channel_areas(session_dir, name_parts, area_limit_file);
if ~channel_ok
    return
end
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
    error('Estimate_sync_lfp_amplitude:sorting', ...
        'Missing cluster_info, spike_times, or spike_clusters.');
end
cluster_info = readtable(info_file, 'FileType', 'text', 'Delimiter', '\t');
ids = table_column(cluster_info, {'cluster_id', 'id'});
channels = table_column(cluster_info, {'ch', 'channel'});
groups = table_column(cluster_info, {'group', 'KSLabel'});
spike_times = double(readNPY(times_file)) / 30000;
spike_clusters = double(readNPY(cluster_file));
spike_clusters = spike_clusters(:);
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
    if ~any(strcmp(areas_to_keep, area))
        continue
    end
    times = spike_times(spike_clusters == cluster_id);
    times = times(isfinite(times));
    times = sort(times(:));
    if isempty(times)
        continue
    end
    neurons(end + 1).cluster_id = cluster_id; %#ok<AGROW>
    neurons(end).channel = channel;
    neurons(end).area = area;
    neurons(end).spike_times = times;
end
end

function [areas_by_channel, ok] = channel_areas(session_dir, name_parts, area_limit_file)
areas_by_channel = cell(384, 1);
ok = false;
map_file = fullfile(session_dir, 'chann_map_PAG.mat');
if exist(map_file, 'file') ~= 2
    error('Estimate_sync_lfp_amplitude:channelMap', 'Missing chann_map_PAG.mat.');
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
        error('Estimate_sync_lfp_amplitude:channelMap', 'Inconsistent channel map.');
    end
end
area_limit = readtable(area_limit_file);
animal_batch = name_parts{1};
repeated_animal = name_parts{3};
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
animal_names = table_column(area_limit, {'AnimalName', 'animalname'});
keep_animal = false(numel(animal_names), 1);
for row_index = 1:numel(animal_names)
    keep_animal(row_index) = strcmp(column_text(animal_names, row_index), this_animal);
end
area_limit = area_limit(keep_animal, :);
if isempty(area_limit)
    error('Estimate_sync_lfp_amplitude:area', ...
        'No area limits for %s.', this_animal);
end
depth_start = double(table_column(area_limit, {'depth_start', 'depthstart'}));
depth_end = double(table_column(area_limit, {'depth_end', 'depthend'}));
area_names = table_column(area_limit, {'area', 'Area'});
probe_area = table_column(area_limit, {'Probe_Area', 'probearea'});
probe_num = double(table_column(area_limit, {'ProbeNum', 'probenum'}));
for site_index = 1:numel(chan_map)
    channel = chan_map(site_index);
    if ~isfinite(channel) || channel < 1 || channel > 384
        continue
    end
    if npx_type == 2
        probe_n = find(any(ismember(npx2_x, xcoords(site_index)), 2), 1);
        probe_ok = probe_num == probe_n;
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
ok = true;
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

function sync_times = synchronous_times(spikes, partner, bout_starts, hist_range, tolerance)
sync_times = spikes(synchronous_mask(spikes, partner, bout_starts, hist_range, tolerance));
end

function mask = synchronous_mask(spikes, partner, bout_starts, hist_range, tolerance)
mask = false(numel(spikes), 1);
for bout_index = 1:numel(bout_starts)
    t0 = bout_starts(bout_index);
    window = spikes >= t0 + hist_range(1) & spikes <= t0 + hist_range(2);
    partner_window = partner(partner >= t0 + hist_range(1) & partner <= t0 + hist_range(2));
    if ~any(window) || isempty(partner_window)
        continue
    end
    in_window = find(window);
    hit = mark_synchronous(spikes(in_window), partner_window, tolerance);
    mask(in_window(hit)) = true;
end
end

function hit = mark_synchronous(spikes, partner, tolerance)
n_spikes = numel(spikes);
n_partner = numel(partner);
hit = false(n_spikes, 1);
partner_start = 1;
for spike_index = 1:n_spikes
    while partner_start <= n_partner && partner(partner_start) < spikes(spike_index) - tolerance
        partner_start = partner_start + 1;
    end
    partner_index = partner_start;
    while partner_index <= n_partner && partner(partner_index) <= spikes(spike_index) + tolerance
        hit(spike_index) = true;
        partner_index = partner_index + 1;
    end
end
end

function bank = load_neuron_amplitudes(session_dir, rows, lfp_sample_rate, n_lfp_channels)
bank = struct('cluster_id', {}, 'spike_times', {}, 'amplitude', {});
spike_bank = full_spike_trains(session_dir, unique([rows.cluster_id]));
channels = unique([rows.channel]);
for channel_index = 1:numel(channels)
    channel = channels(channel_index);
    raw = read_lfp_channel(session_dir, channel, n_lfp_channels);
    envelope = abs(hilbert(double(raw(:))));
    clear raw
    lfp_time = (1:numel(envelope))' / lfp_sample_rate;
    neuron_ids = unique([rows([rows.channel] == channel).cluster_id]);
    for neuron_index = 1:numel(neuron_ids)
        cluster_id = neuron_ids(neuron_index);
        spike_times = spike_bank_times(spike_bank, cluster_id);
        if isempty(spike_times)
            continue
        end
        amplitude = interp1(lfp_time, envelope, spike_times, 'linear', NaN);
        bank(end + 1).cluster_id = cluster_id; %#ok<AGROW>
        bank(end).spike_times = spike_times;
        bank(end).amplitude = amplitude(:);
    end
    clear envelope lfp_time
end
end

function trains = full_spike_trains(session_dir, cluster_ids)
trains = struct('cluster_id', {}, 'spike_times', {});
info_file = fullfile(session_dir, 'cluster_info_PAG.tsv');
suffix = '_PAG';
if exist(info_file, 'file') ~= 2
    suffix = '';
end
times_file = fullfile(session_dir, ['spike_times' suffix '.npy']);
cluster_file = fullfile(session_dir, ['spike_clusters' suffix '.npy']);
spike_times = double(readNPY(times_file)) / 30000;
spike_clusters = double(readNPY(cluster_file));
for index = 1:numel(cluster_ids)
    times = spike_times(spike_clusters == cluster_ids(index));
    times = sort(times(isfinite(times)));
    trains(end + 1).cluster_id = cluster_ids(index); %#ok<AGROW>
    trains(end).spike_times = times(:);
end
end

function times = spike_bank_times(spike_bank, cluster_id)
times = [];
for index = 1:numel(spike_bank)
    if spike_bank(index).cluster_id == cluster_id
        times = spike_bank(index).spike_times;
        return
    end
end
end

function raw = read_lfp_channel(session_dir, channel, n_lfp_channels)
mat_file = fullfile(session_dir, 'LFP_PAG.mat');
dat_file = fullfile(session_dir, 'LFP_PAG.dat');
if exist(mat_file, 'file') == 2
    loaded = load(mat_file, 'LFP');
    raw = loaded.LFP(channel, :);
    return
end
if exist(dat_file, 'file') ~= 2
    error('Estimate_sync_lfp_amplitude:lfp', ...
        'Missing LFP_PAG.mat and LFP_PAG.dat.');
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (n_lfp_channels * 2);
if n_samples ~= round(n_samples) || n_samples < 1
    error('Estimate_sync_lfp_amplitude:lfp', ...
        'LFP_PAG.dat is not a multiple of %d int16 channels.', n_lfp_channels);
end
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('Estimate_sync_lfp_amplitude:lfp', 'Could not open %s.', dat_file);
end
closer = onCleanup(@() fclose(file_id)); %#ok<NASGU>
if fseek(file_id, 2 * (channel - 1), 'bof') ~= 0
    error('Estimate_sync_lfp_amplitude:lfp', ...
        'Could not seek to LFP channel %d.', channel);
end
[raw, count] = fread(file_id, n_samples, 'int16=>double', 2 * (n_lfp_channels - 1));
if count ~= n_samples
    error('Estimate_sync_lfp_amplitude:lfp', ...
        'Read %d of %d LFP samples on channel %d.', count, n_samples, channel);
end
end

function comparison = compare_mean_amplitude(amplitude, spike_times, sync_times, n_repeats)
comparison.ok = false;
comparison.reason = '';
comparison.n_sync = 0;
comparison.mean_all = NaN;
comparison.mean_sync = NaN;
comparison.null_means = nan(n_repeats, 1);
comparison.percentile = NaN;
sync_index = ismembertol_times(spike_times, sync_times);
sync_amplitude = amplitude(sync_index);
sync_amplitude = sync_amplitude(isfinite(sync_amplitude));
other = amplitude(~sync_index);
other = other(isfinite(other));
comparison.n_sync = numel(sync_amplitude);
comparison.mean_all = mean_finite(amplitude);
if isempty(sync_amplitude)
    comparison.reason = 'no finite LFP amplitude at the synchronous spikes';
    return
end
if numel(other) < numel(sync_amplitude)
    comparison.reason = 'fewer other spikes than synchronous spikes';
    return
end
comparison.mean_sync = mean(sync_amplitude);
n_sync = numel(sync_amplitude);
for repeat_index = 1:n_repeats
    draw = other(randperm(numel(other), n_sync));
    comparison.null_means(repeat_index) = mean(draw);
end
comparison.percentile = mean(comparison.null_means < comparison.mean_sync);
comparison.ok = true;
end

function mask = ismembertol_times(spike_times, sync_times)
mask = false(numel(spike_times), 1);
if isempty(sync_times) || isempty(spike_times)
    return
end
sync_times = sort(sync_times(:));
spike_index = 1;
for sync_index = 1:numel(sync_times)
    while spike_index <= numel(spike_times) ...
            && spike_times(spike_index) < sync_times(sync_index) - 1e-9
        spike_index = spike_index + 1;
    end
    if spike_index <= numel(spike_times) ...
            && abs(spike_times(spike_index) - sync_times(sync_index)) <= 1e-9
        mask(spike_index) = true;
    end
end
end

function value = mean_finite(values)
values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = mean(values);
end
end

function neuron = find_neuron(bank, cluster_id)
neuron = [];
for index = 1:numel(bank)
    if bank(index).cluster_id == cluster_id
        neuron = bank(index);
        return
    end
end
end

function rows = append_neuron(rows, session_id, session_dir, cluster_id, ...
    partner_id, channel, area, n_spikes, n_sync, sync_times)
rows(end + 1).session_id = session_id; %#ok<AGROW>
rows(end).session_dir = session_dir;
rows(end).cluster_id = cluster_id;
rows(end).partner_id = partner_id;
rows(end).channel = channel;
rows(end).area = area;
rows(end).n_spikes = n_spikes;
rows(end).n_sync = n_sync;
rows(end).sync_times = sync_times;
end

function counts = note_neuron(counts, session_id, cluster_id)
for index = 1:numel(counts)
    if strcmp(counts(index).session_id, session_id) && counts(index).cluster_id == cluster_id
        counts(index).n_pairs = counts(index).n_pairs + 1;
        return
    end
end
counts(end + 1).session_id = session_id; %#ok<AGROW>
counts(end).cluster_id = cluster_id;
counts(end).n_pairs = 1;
end

function full = neuron_is_full(counts, session_id, cluster_id, maximum)
full = false;
for index = 1:numel(counts)
    if strcmp(counts(index).session_id, session_id) && counts(index).cluster_id == cluster_id
        full = counts(index).n_pairs >= maximum;
        return
    end
end
end

function session_ids = unique_session_ids(rows)
session_ids = {};
for index = 1:numel(rows)
    if ~any(strcmp(session_ids, rows(index).session_id))
        session_ids{end + 1} = rows(index).session_id; %#ok<AGROW>
    end
end
end

function save_amplitude_results(results_dir, chosen_pairs, amplitude_results, seed)
if exist(results_dir, 'dir') ~= 7
    mkdir(results_dir);
end
save(fullfile(results_dir, 'sync_lfp_amplitude.mat'), ...
    'chosen_pairs', 'amplitude_results', 'seed');
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
error('Estimate_sync_lfp_amplitude:column', ...
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
if isempty(text) && ~iscell(groups)
    value = groups(index);
    if isnumeric(value)
        label = '';
        return
    end
end
label = lower(strtrim(text));
end

function number = scalar_number(values, index)
if iscell(values)
    number = double(values{index});
else
    number = double(values(index));
end
end
