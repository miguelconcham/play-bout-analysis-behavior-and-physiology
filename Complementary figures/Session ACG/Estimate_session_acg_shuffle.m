%% Estimate_session_acg_shuffle
% Whole-session autocorrelograms, 10 s power windows, and their nulls.
%
% Run this file and leave it. It walks one session at a time. For every
% trough-locked, peak-locked, and non-entrained neuron it saves the
% autocorrelogram of all spikes, 500 uniform-time nulls, and the same pair
% for the low and high tails of FOOOF delta and of Hilbert delta amplitude.
% acg_method selects the autocorrelogram. 'edge' divides pair counts by
% (T-|lag|). 'crosscorr' uses crosscorrelogram, with reference spikes inset
% by acg_lag. On a 10 s window that is spikes in [1, 9] against spikes in
% the whole window.
% A stopped run continues from the neurons already saved when resume is
% true. After the last session, the grand mean of each lock group is built
% from those saved nulls. Delta-locked is trough-locked and peak-locked
% together. Nothing here is recomputed from spike times during that step.
%
% The morning file is grand_average_oscillation.txt, next to the neuron
% files. Each neuron file keeps the observed curve and all 500 null curves
% (about 3 MB for the five autocorrelograms). A neuron oscillates in delta
% when its autocorrelogram has a delta peak, scored as in acg_delta_peaks,
% and that peak's prominence is larger than its own null at
% oscillation_alpha. The same call is made on the grand mean.
%
% Set grand_only to true to rebuild the grand average without shuffling.
% only_sessions limits which folders are estimated on this run. Sessions
% already saved in the output folder still enter the grand average.
% resume false does not delete an existing folder.
%
% Check the helpers first, without the experiment files, by running
% test_session_acg_shuffle. This folder is Session ACG, next to the other
% complementary-figure scripts. add_repo_paths puts it on the path.

%% Parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_root = fullfile(data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(data_root, 'Analysis results', 'session acg shuffle');
phase_folder = fullfile(data_root, 'Analysis results', 'phase locking data');
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

n_shuffles = 500;
window_s = 10;
min_spikes = 200;
tail_percent = 25;
alpha = 0.01;
non_entrained_p = 0.1;
oscillation_alpha = 0.05;
acg_lag = 1;
bin_w = 0.01;
delta_band = [1 5];
smooth_half = 0.06;
acg_refractory_s = 0.005;
areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
fooof_max_peaks = 6;
fooof_fit_range = [0.5 40];
filter_order = 2000;
lfp_sample_rate = 2500;
% acg_method chooses the autocorrelogram.
%   'edge'       every pair in the interval, divided by (T-|lag|)
%   'crosscorr'  crosscorrelogram. Reference spikes are inset by acg_lag.
%                For a 10 s window and +/- 1 s lags, references are in
%                [1, 9] and partners are the spikes in the whole window.
acg_method = 'crosscorr';
resume = true;
grand_only = false;
only_sessions = {};
animal2exclude = {};

options = struct( ...
    'version', 1, ...
    'n_shuffles', n_shuffles, ...
    'window_s', window_s, ...
    'min_spikes', min_spikes, ...
    'tail_percent', tail_percent, ...
    'alpha', alpha, ...
    'non_entrained_p', non_entrained_p, ...
    'oscillation_alpha', oscillation_alpha, ...
    'acg_lag', acg_lag, ...
    'bin_w', bin_w, ...
    'delta_band', delta_band, ...
    'smooth_half', smooth_half, ...
    'acg_refractory_s', acg_refractory_s, ...
    'areas_to_keep', {areas_to_keep}, ...
    'fooof_max_peaks', fooof_max_peaks, ...
    'fooof_fit_range', fooof_fit_range, ...
    'filter_order', filter_order, ...
    'lfp_sample_rate', lfp_sample_rate, ...
    'acg_method', acg_method, ...
    'resume', resume);

if exist(saving_folder, 'dir') ~= 7
    mkdir(saving_folder);
end
log_file = fullfile(saving_folder, 'estimate_log.txt');
options.log_file = log_file;
settings = session_acg_settings(options);
options.bin_edges = settings.bin_edges;
options.smooth_n = settings.smooth_n;
param_file = fullfile(saving_folder, 'parameters.mat');
if exist(param_file, 'file') == 2
    saved = load(param_file, 'settings');
    if ~session_acg_parameters_match(saved.settings, settings)
        error('session_acg:parameters', ...
            ['%s was started with different settings. Move that folder ' ...
            'aside, or put the previous settings back, before running again.'], ...
            saving_folder);
    end
elseif grand_only
    error('session_acg:grand', ...
        'grand_only is set, but %s has no parameters.mat.', saving_folder);
else
    save(param_file, 'settings');
end

session_acg_log(log_file, sprintf( ...
    'Estimate folder %s. Method %s. Shuffles %d, window %g s, minimum spikes %d, tail %g%%, resume %d, grand only %d.', ...
    saving_folder, acg_method, n_shuffles, window_s, min_spikes, tail_percent, resume, grand_only));

%% Sessions
if ~grand_only
    delta_file = fullfile(phase_folder, 'delta_all_neurons_v2.mat');
    if exist(delta_file, 'file') ~= 2
        error('session_acg:delta', 'Delta-lock table was not found at %s.', delta_file);
    end
    loaded_delta = load(delta_file, 'all_neurons');
    all_neurons = loaded_delta.all_neurons;
    clear loaded_delta
    if ~(istable(all_neurons) && ismember('EntireSession', all_neurons.Properties.VariableNames))
        error('session_acg:delta', ...
            'delta_all_neurons_v2.mat has no EntireSession table.');
    end
    if ismember('area', all_neurons.Properties.VariableNames)
        all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
    end

    animal_list = dir(npx_root);
    animal_list = animal_list([animal_list.isdir]);
    animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));
    if ~isempty(animal2exclude)
        animal_list = animal_list(~ismember({animal_list.name}, animal2exclude));
    end
    if ~isempty(only_sessions)
        animal_list = animal_list(ismember({animal_list.name}, only_sessions));
    end
    [~, name_order] = sort({animal_list.name});
    animal_list = animal_list(name_order);
    session_acg_log(log_file, sprintf('%d session folders to visit.', numel(animal_list)));

    for session_index = 1:numel(animal_list)
        animal_name = animal_list(session_index).name;
        session_dir = fullfile(npx_root, animal_name);
        session_acg_log(log_file, sprintf('Session %d/%d: %s', ...
            session_index, numel(animal_list), animal_name));
        try
            recording = load_pag_recording(session_dir);
            if ~recording.ok
                session_acg_log(log_file, sprintf('Skip %s: %s', animal_name, recording.reason));
                clear recording
                continue
            end
            if recording.sample_rate ~= lfp_sample_rate
                session_acg_log(log_file, sprintf( ...
                    'Skip %s: LFP sample rate is %g, not %g.', ...
                    animal_name, recording.sample_rate, lfp_sample_rate));
                clear recording
                continue
            end
            neurons = neurons_for_session(all_neurons, animal_name, ...
                recording.cluster_info, areas_to_keep, alpha, non_entrained_p);
            if isempty(neurons)
                session_acg_log(log_file, sprintf('Skip %s: no trough, peak, or non-entrained neurons in the kept areas.', animal_name));
                clear recording neurons
                continue
            end
            session_acg_log(log_file, sprintf('%s: %d neurons.', animal_name, numel(neurons)));
            estimate_one_session_acg(animal_name, recording.spikes, recording.clusters, ...
                neurons, recording.lfp, recording.lfp_time, options, ...
                fullfile(saving_folder, animal_name));
            clear recording neurons
        catch failure
            if any(strcmp(failure.identifier, {'session_acg:parameters', 'session_acg:resume'}))
                rethrow(failure);
            end
            session_acg_log(log_file, sprintf('Session %s failed: %s', animal_name, failure.message));
            clear recording neurons
        end
    end
    clear all_neurons
end

%% Grand average from the saved nulls
session_acg_log(log_file, 'Building the grand average from the saved nulls.');
session_acg_grand_average(saving_folder);
session_acg_log(log_file, 'Grand average written.');


function neurons = neurons_for_session(all_neurons, animal_name, cluster_info, areas_to_keep, alpha, non_entrained_p)
id_all = table_cluster_id(all_neurons);
in_session = session_name_matches(all_neurons.session, animal_name);
rows = find(in_session);
cluster_ids = [];
channels = [];
if ~isempty(cluster_info) && height(cluster_info) > 0 ...
        && (ismember('cluster_id', cluster_info.Properties.VariableNames) ...
        || ismember('id', cluster_info.Properties.VariableNames))
    cluster_ids = table_cluster_id(cluster_info);
    if ismember('ch', cluster_info.Properties.VariableNames)
        channels = double(cluster_info.ch) + 1;
    else
        channels = nan(size(cluster_ids));
    end
end

n_rows = numel(rows);
keep = false(n_rows, 1);
cluster_id = nan(n_rows, 1);
channel = nan(n_rows, 1);
area = cell(n_rows, 1);
lock_group = cell(n_rows, 1);
for i = 1:n_rows
    row = rows(i);
    group_name = delta_lock_group( ...
        all_neurons.EntireSession.PPCPval(row), ...
        all_neurons.EntireSession.PreferedAngle(row), ...
        alpha, non_entrained_p);
    if isempty(group_name)
        continue
    end
    this_area = canonical_kept_area(table_text(all_neurons, 'area', row), areas_to_keep);
    if isempty(this_area)
        continue
    end
    keep(i) = true;
    cluster_id(i) = id_all(row);
    area{i} = this_area;
    lock_group{i} = group_name;
    if ~isempty(cluster_ids)
        hit = find(cluster_ids == id_all(row), 1);
        if ~isempty(hit)
            channel(i) = channels(hit);
        end
    end
end
if ~any(keep)
    neurons = struct('cluster_id', {}, 'channel', {}, 'area', {}, 'lock_group', {});
    return
end
cluster_id = cluster_id(keep);
channel = channel(keep);
area = area(keep);
lock_group = lock_group(keep);
[~, unique_index] = unique(cluster_id, 'stable');
neurons = struct( ...
    'cluster_id', num2cell(cluster_id(unique_index)), ...
    'channel', num2cell(channel(unique_index)), ...
    'area', area(unique_index), ...
    'lock_group', lock_group(unique_index));
end
