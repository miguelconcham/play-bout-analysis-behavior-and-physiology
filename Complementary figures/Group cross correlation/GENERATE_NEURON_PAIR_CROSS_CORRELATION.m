function session = GENERATE_NEURON_PAIR_CROSS_CORRELATION(npx_data_dir, members, options)
%GENERATE_NEURON_PAIR_CROSS_CORRELATION  Neuron-pair correlograms for one recording.
%
% Spikes are not pooled across neurons. Every neuron in one of the six
% groups is cross-correlated with every neuron in each of the other
% groups. Lags come from crosscorrelogram and are counted in the estimate
% bins, 1 ms by default. Each group combination keeps a matrix with one
% row per neuron pair.

if nargin < 1 || isempty(npx_data_dir)
    error('group_ccg:folder', 'Provide the NPX session folder.');
end
if nargin < 2 || ~isstruct(members) || ~isfield(members, 'neuron_id')
    error('group_ccg:members', 'Provide the six-group membership.');
end
if nargin < 3 || isempty(options)
    options = struct();
end
settings = group_ccg_settings(options);
[~, animal_code] = fileparts(npx_data_dir);
names = group_ccg_names();
selected = session_name_matches(members.session, animal_code);
[spike_times, spike_clusters] = load_pag_spike_times(npx_data_dir);
neuron_id = members.neuron_id(selected);
group_index = members.group_index(selected);
raw_spikes = spikes_by_neuron(spike_times, spike_clusters, neuron_id);
if isempty(settings.session_type)
    windows = spike_span_window(raw_spikes);
else
    windows = partners_session_windows(animal_code, settings.session_type);
end
if isempty(windows)
    error('group_ccg:windows', ...
        'No "%s" intervals for %s.', settings.session_type, animal_code);
end
[group_trains, group_ids, missing_cluster_id, silent_cluster_id] = trains_by_group( ...
    raw_spikes, neuron_id, group_index, windows, names);
combinations = session_neuron_pair_ccgs(group_trains, group_ids, windows, ...
    settings.bin_edges, settings.min_reference_spikes);

session.version = settings.version;
session.animal_code = animal_code;
session.windows = windows;
session.bin_edges = settings.bin_edges;
session.centers = settings.centers;
session.areas_to_keep = settings.areas_to_keep;
session.session_type = settings.session_type;
session.min_reference_spikes = settings.min_reference_spikes;
session.acg_lag = settings.acg_lag;
session.acg_bin = settings.acg_bin;
session.group_ids = group_ids;
session.combinations = combinations;
session.n_members = sum(selected);
session.missing_cluster_id = missing_cluster_id;
session.silent_cluster_id = silent_cluster_id;
end


function raw_spikes = spikes_by_neuron(spike_times, spike_clusters, neuron_id)
raw_spikes = cell(numel(neuron_id), 1);
for i = 1:numel(neuron_id)
    spikes = spike_times(spike_clusters == neuron_id(i));
    raw_spikes{i} = sort(spikes(isfinite(spikes)));
end
end


function [group_trains, group_ids, missing_cluster_id, silent_cluster_id] = trains_by_group( ...
    raw_spikes, neuron_id, group_index, windows, names)
n_groups = numel(names.group_id);
group_trains = cell(1, n_groups);
group_ids = cell(1, n_groups);
for g = 1:n_groups
    group_trains{g} = {};
    group_ids{g} = zeros(0, 1);
end
missing_cluster_id = zeros(0, 1);
silent_cluster_id = zeros(0, 1);
for i = 1:numel(neuron_id)
    if isempty(raw_spikes{i})
        missing_cluster_id(end + 1, 1) = neuron_id(i); %#ok<AGROW>
        continue
    end
    spikes = spikes_inside_windows(raw_spikes{i}, windows);
    if isempty(spikes)
        silent_cluster_id(end + 1, 1) = neuron_id(i); %#ok<AGROW>
        continue
    end
    g = group_index(i);
    group_trains{g}{end + 1} = spikes; %#ok<AGROW>
    group_ids{g}(end + 1, 1) = neuron_id(i); %#ok<AGROW>
end
end


function windows = spike_span_window(raw_spikes)
first = inf;
last = -inf;
for i = 1:numel(raw_spikes)
    spikes = raw_spikes{i};
    if isempty(spikes)
        continue
    end
    first = min(first, spikes(1));
    last = max(last, spikes(end));
end
if isfinite(first) && last > first
    windows = [first, last];
else
    windows = zeros(0, 2);
end
end


function spikes = spikes_inside_windows(spikes, windows)
if isempty(spikes) || isempty(windows)
    spikes = zeros(0, 1);
    return
end
inside = false(numel(spikes), 1);
for w = 1:size(windows, 1)
    inside = inside | (spikes >= windows(w, 1) & spikes <= windows(w, 2));
end
spikes = spikes(inside);
end


function windows = partners_session_windows(animal_code, session_type)
data_root = group_ccg_data_root();
synch_file = fullfile(data_root, 'Synch data', animal_code, 'synch_model_video2NPX.mat');
behavior_file = fullfile(data_root, 'Behavior backups', [animal_code, '.txt']);
if exist(synch_file, 'file') ~= 2 || exist(behavior_file, 'file') ~= 2
    error('group_ccg:behavior', ...
        'Synchronized behavior for %s was not found.', animal_code);
end
loaded = load(synch_file, 'synch_model_video2NPX');
behavior = readtable(behavior_file);
behavior(:, 2) = [];
behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
behavior.Type2 = behavior.Type;
behavior.Type2(ismember(behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'};
behavior.Type2(ismember(behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
behavior.Type2(ismember(behavior.Type2, {''})) = {'Other'};
behavior(ismember(behavior.Animal, {'Reversal'}), :) = [];
behavior.Start = predict(loaded.synch_model_video2NPX, behavior.Start);
behavior.End = predict(loaded.synch_model_video2NPX, behavior.End);
kept = strcmp(behavior.Type2, session_type);
windows = merge_time_intervals([behavior.Start(kept), behavior.End(kept)]);
if isempty(windows)
    error('group_ccg:windows', ...
        'No "%s" intervals for %s.', session_type, animal_code);
end
end


function [spike_times, spike_clusters] = load_pag_spike_times(session_dir)
times_file = fullfile(session_dir, 'spike_times_PAG.npy');
clusters_file = fullfile(session_dir, 'spike_clusters_PAG.npy');
if exist(times_file, 'file') ~= 2 || exist(clusters_file, 'file') ~= 2
    error('group_ccg:spikes', ...
        'PAG spike files are missing in %s.', session_dir);
end
if exist('readNPY', 'file') ~= 2
    error('group_ccg:spikes', ...
        'readNPY is not on the MATLAB path. add_repo_paths does not install it.');
end
spike_times = double(readNPY(times_file)) / 30000;
spike_clusters = double(readNPY(clusters_file));
spike_times = spike_times(:);
spike_clusters = spike_clusters(:);
if numel(spike_times) ~= numel(spike_clusters)
    error('group_ccg:spikes', ...
        'Spike times and spike clusters in %s have different lengths.', session_dir);
end
end
