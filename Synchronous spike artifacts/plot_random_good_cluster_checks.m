function summary = plot_random_good_cluster_checks(varargin)
%PLOT_RANDOM_GOOD_CLUSTER_CHECKS  Waveform and ACG for one good unit per session.
%
%   summary = plot_random_good_cluster_checks
%
% Runs check_npx_raw_folder_composition first. Figures are drawn only when
% that comparison agrees. For each matched session it picks one random
% cluster whose group is good in the NPX raw-data folder, loads the
% waveform with load_waveforms_from_continuous, and plots that waveform on
% the channels around the unit together with the autocorrelogram.
%
% ChannelRange defaults to 5 (five channels on each side). A unit on
% channel 3 keeps channels 1:8; a unit in the middle of the probe keeps 11.
%
% Name-value pairs: ChannelRange, NWaveforms, ExperimentInfoFile,
% RawFoldersFile, NpxRawData, SampleRate, Report. Pass Report to reuse a
% comparison instead of reading the workbooks again.

options = parse_check_options(varargin);
if isempty(options.report)
    options.report = check_npx_raw_folder_composition( ...
        options.experiment_info_file, options.raw_folders_file);
end

summary = empty_summary();
if ~options.report.ok
    fprintf(['Folder composition does not match. ' ...
        'No cluster waveforms were plotted.\n']);
    return
end

pairs = matched_session_files(options.report);
if isempty(pairs)
    fprintf('No matched sessions to plot.\n');
    return
end

rng('shuffle');
random_state = rng;
fprintf('Random cluster draw, seed %d.\n', random_state.Seed);

for pair_index = 1:numel(pairs)
    session_name = pairs(pair_index).session_name;
    continuous_file = resolve_continuous_file(pairs(pair_index).continuous_file);
    session_dir = fullfile(options.npx_raw_data, session_name);
    fprintf('Session %d of %d: %s\n', pair_index, numel(pairs), session_name);

    entry = blank_summary(session_name, continuous_file);
    if isempty(continuous_file)
        entry.detail = 'Continuous file was not found.';
        summary(end + 1) = entry; %#ok<AGROW>
        fprintf('  %s\n', entry.detail);
        continue
    end
    if exist(session_dir, 'dir') ~= 7
        entry.detail = ['NPX raw-data folder was not found: ' session_dir];
        summary(end + 1) = entry; %#ok<AGROW>
        fprintf('  %s\n', entry.detail);
        continue
    end

    try
        clusters = list_good_clusters(session_dir);
        if isempty(clusters)
            entry.detail = 'No good cluster in this session folder.';
            summary(end + 1) = entry; %#ok<AGROW>
            fprintf('  %s\n', entry.detail);
            continue
        end
        cluster = clusters(randi(numel(clusters)));
        fprintf('  cluster %g, channel %d, %d spikes (%s)\n', ...
            cluster.id, cluster.channel, cluster.n_spikes, cluster.source);
        result = plot_cluster_waveform_and_acg( ...
            continuous_file, cluster.spike_samples, cluster.channel, ...
            'ChannelRange', options.channel_range, ...
            'NWaveforms', options.n_waveforms, ...
            'SampleRate', options.sample_rate, ...
            'SessionName', session_name, ...
            'ClusterID', cluster.id);
        entry.cluster_id = cluster.id;
        entry.channel = cluster.channel;
        entry.n_spikes = cluster.n_spikes;
        entry.source = cluster.source;
        entry.figure = result.figure;
        entry.detail = '';
    catch ME
        entry.detail = ME.message;
        fprintf('  %s\n', entry.detail);
    end
    summary(end + 1) = entry; %#ok<AGROW>
end
end

function options = parse_check_options(name_value)
options.channel_range = 5;
options.n_waveforms = 40;
options.sample_rate = 30000;
options.experiment_info_file = '';
options.raw_folders_file = '';
options.npx_raw_data = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\NPX data\NPX raw data'];
options.report = [];

if mod(numel(name_value), 2) ~= 0
    error('plot_random_good_cluster_checks:nameValue', ...
        'Optional arguments must be name-value pairs.');
end
for pair_index = 1:2:numel(name_value)
    name = name_value{pair_index};
    value = name_value{pair_index + 1};
    if exist('isstring', 'builtin') && isstring(name)
        name = char(name);
    end
    switch lower(strtrim(name))
        case 'channelrange'
            options.channel_range = value;
        case 'nwaveforms'
            options.n_waveforms = value;
        case 'samplerate'
            options.sample_rate = value;
        case 'experimentinfofile'
            options.experiment_info_file = char(value);
        case 'rawfoldersfile'
            options.raw_folders_file = char(value);
        case 'npxrawdata'
            options.npx_raw_data = char(value);
        case 'report'
            options.report = value;
        otherwise
            error('plot_random_good_cluster_checks:nameValue', ...
                'Unknown parameter "%s".', name);
    end
end
end

function pairs = matched_session_files(report)
pairs = struct('session_name', {}, 'continuous_file', {});
if ~isfield(report, 'rows') || isempty(report.rows)
    return
end
for row_index = 1:numel(report.rows)
    row = report.rows(row_index);
    if ~strcmp(row.status, 'ok') || isempty(row.session_name) || isempty(row.actual_path)
        continue
    end
    already = false;
    for pair_index = 1:numel(pairs)
        if strcmp(pairs(pair_index).session_name, row.session_name) ...
                && strcmp(pairs(pair_index).continuous_file, row.actual_path)
            already = true;
            break
        end
    end
    if already
        continue
    end
    pairs(end + 1).session_name = row.session_name; %#ok<AGROW>
    pairs(end).continuous_file = row.actual_path;
end
end

function clusters = list_good_clusters(session_dir)
info_files = dir(fullfile(session_dir, 'cluster_info*.tsv'));
clusters = struct('id', {}, 'channel', {}, 'spike_samples', {}, ...
    'n_spikes', {}, 'source', {});
for file_index = 1:numel(info_files)
    info_name = info_files(file_index).name;
    [~, info_base] = fileparts(info_name);
    suffix = '';
    prefix = 'cluster_info';
    if numel(info_base) > numel(prefix)
        suffix = info_base(numel(prefix) + 1:end);
    end
    times_file = fullfile(session_dir, ['spike_times' suffix '.npy']);
    cluster_file = fullfile(session_dir, ['spike_clusters' suffix '.npy']);
    if exist(times_file, 'file') ~= 2 || exist(cluster_file, 'file') ~= 2
        fprintf('  missing spike_times%s.npy or spike_clusters%s.npy\n', suffix, suffix);
        continue
    end
    cluster_info = readtable(fullfile(session_dir, info_name), ...
        'FileType', 'text', 'Delimiter', '\t');
    spike_times = readNPY(times_file);
    spike_clusters = readNPY(cluster_file);
    source_name = suffix;
    if ~isempty(source_name) && source_name(1) == '_'
        source_name = source_name(2:end);
    end
    if isempty(source_name)
        source_name = 'cluster_info';
    end
    these = good_clusters_from_sorting(cluster_info, spike_times, spike_clusters, source_name);
    if isempty(these)
        continue
    end
    if isempty(clusters)
        clusters = these;
    else
        clusters = [clusters, these]; %#ok<AGROW>
    end
end
end

function continuous_file = resolve_continuous_file(continuous_file)
if isempty(continuous_file)
    return
end
if exist(continuous_file, 'dir') == 7
    candidate = fullfile(continuous_file, 'continuous.dat');
    if exist(candidate, 'file') == 2
        continuous_file = candidate;
    else
        continuous_file = '';
    end
    return
end
if exist(continuous_file, 'file') ~= 2
    continuous_file = '';
end
end

function summary = empty_summary()
summary = blank_summary('', '');
summary(1) = [];
end

function entry = blank_summary(session_name, continuous_file)
entry.session_name = session_name;
entry.continuous_file = continuous_file;
entry.cluster_id = [];
entry.channel = [];
entry.n_spikes = [];
entry.source = '';
entry.figure = [];
entry.detail = '';
end
