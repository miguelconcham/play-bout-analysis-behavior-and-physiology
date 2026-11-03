% Check_NPX_folders_and_load_waveforms
%
% One script, three sections. Run this file from the Editor, or type its name
% in the command window. Later sections use the variables created above them.
% Supporting functions are at the bottom of this file.
%
% 1. Read the continuous-file locations in NPX_raw_data_folders.xlsx.
% 2. In each session folder under NPX raw data, keep the four good units with
%    the most spikes. selected_units has cluster_id, session_id, animal_name,
%    spike_times (samples), and channel (cluster_info.ch is 0-based, so this
%    is ch + 1).
% 3. Load a subset of those spikes from that session's continuous file and
%    draw one figure per session: four rows, waveform and autocorrelogram.
%    Each waveform row is scaled to its mean spike, after the pre-spike
%    baseline is removed, so a small spike is not drawn flat inside the
%    raw-count range.

%% Paths
raw_folders_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\NPX data\' ...
    'NPX_raw_data_folders.xlsx'];
npx_raw_data = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\NPX data\NPX raw data'];

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

%% 1. Continuous-file locations
if exist(raw_folders_file, 'file') ~= 2
    error('Check_NPX_folders_and_load_waveforms:file', ...
        'Could not find %s', raw_folders_file);
end
fprintf('Continuous-file locations: %s\n', raw_folders_file);

%% 2. Four good units per session, most spikes first
n_units_per_session = 4;
rng('shuffle');
seed_state = rng;
fprintf('Random seed: %d\n', seed_state.Seed);

selected_units = struct('cluster_id', {}, 'session_id', {}, ...
    'animal_name', {}, 'spike_times', {}, 'channel', {});

if exist(npx_raw_data, 'dir') ~= 7
    fprintf('NPX raw data folder was not found:\n  %s\n', npx_raw_data);
else
    listing = dir(npx_raw_data);
    for folder_index = 1:numel(listing)
        if ~listing(folder_index).isdir
            continue
        end
        session_id = listing(folder_index).name;
        if strcmp(session_id, '.') || strcmp(session_id, '..')
            continue
        end
        session_dir = fullfile(npx_raw_data, session_id);
        name_tokens = regexp(strtrim(session_id), '\s+', 'split');
        name_tokens = name_tokens(~cellfun('isempty', name_tokens));
        if isempty(name_tokens)
            animal_name = session_id;
        else
            animal_name = name_tokens{1};
        end
        units = good_units_in_session(session_dir);
        picked = units_with_most_spikes(units, n_units_per_session);
        if isempty(picked)
            fprintf('%s  no good unit\n', session_id);
            continue
        end
        for pick_index = 1:numel(picked)
            pick = picked(pick_index);
            selected_units(end + 1).cluster_id = pick.cluster_id; %#ok<SAGROW>
            selected_units(end).session_id = session_id;
            selected_units(end).animal_name = animal_name;
            selected_units(end).spike_times = pick.spike_times;
            selected_units(end).channel = pick.channel;
            fprintf('%s  cluster %g  channel %g  %d spikes\n', ...
                session_id, pick.cluster_id, pick.channel, numel(pick.spike_times));
        end
    end
end
fprintf('Selected %d units.\n', numel(selected_units));

%% 3. Load waveforms and plot four neurons per session
channel_range = 5;
n_waveforms = 40;
fs = 30000;
acg_window = [-0.05 0.05];
acg_bin = 0.001;
n_acg_spikes = 3000;

session_start = 1;
while session_start <= numel(selected_units)
    session_id = selected_units(session_start).session_id;
    session_stop = session_start;
    while session_stop < numel(selected_units) ...
            && strcmp(selected_units(session_stop + 1).session_id, session_id)
        session_stop = session_stop + 1;
    end
    try
        paths = read_npx_continuous_locations(raw_folders_file, session_id);
    catch location_error
        fprintf('%s  %s\n', session_id, location_error.message);
        session_start = session_stop + 1;
        continue
    end
    for path_index = 1:numel(paths)
        continuous_file = continuous_file_from_location(paths{path_index});
        loaded_units = struct('cluster_id', {}, 'session_id', {}, ...
            'animal_name', {}, 'spike_times', {}, 'channel', {});
        waveform_cells = {};
        channel_cells = {};
        time_samples = [];
        for unit_index = session_start:session_stop
            unit = selected_units(unit_index);
            n_draw = min(n_waveforms, numel(unit.spike_times));
            draw_index = randperm(numel(unit.spike_times), n_draw);
            waveform_samples = unit.spike_times(draw_index);
            fprintf('%s  cluster %g  %s\n', session_id, unit.cluster_id, continuous_file);
            try
                [waveforms, channels_kept, these_times] = load_waveforms_from_continuous( ...
                    continuous_file, waveform_samples, unit.channel, channel_range);
            catch load_error
                fprintf('  %s\n', load_error.message);
                continue
            end
            if ~isfield(selected_units, 'waveforms') || isempty(selected_units(unit_index).waveforms)
                selected_units(unit_index).continuous_file = continuous_file;
                selected_units(unit_index).waveform_samples = waveform_samples;
                selected_units(unit_index).waveforms = waveforms;
                selected_units(unit_index).channels_kept = channels_kept;
                selected_units(unit_index).time_samples = these_times;
            end
            loaded_units(end + 1).cluster_id = unit.cluster_id; %#ok<AGROW>
            loaded_units(end).session_id = unit.session_id;
            loaded_units(end).animal_name = unit.animal_name;
            loaded_units(end).spike_times = unit.spike_times;
            loaded_units(end).channel = unit.channel;
            waveform_cells{end + 1} = waveforms; %#ok<AGROW>
            channel_cells{end + 1} = channels_kept; %#ok<AGROW>
            time_samples = these_times;
        end
        if isempty(loaded_units)
            continue
        end
        try
            show_session_figure(loaded_units, continuous_file, waveform_cells, ...
                channel_cells, time_samples, fs, acg_window, acg_bin, n_acg_spikes);
        catch plot_error
            fprintf('  waveforms were loaded; the figure failed: %s\n', ...
                plot_error.message);
        end
    end
    session_start = session_stop + 1;
end

%% Local functions

function continuous_file = continuous_file_from_location(location)
continuous_file = location;
if exist(location, 'dir') == 7
    continuous_file = fullfile(location, 'continuous.dat');
end
end

function picked = units_with_most_spikes(units, n_keep)
if isempty(units) || n_keep < 1
    picked = units;
    return
end
counts = zeros(1, numel(units));
for unit_index = 1:numel(units)
    counts(unit_index) = numel(units(unit_index).spike_times);
end
[~, order] = sort(counts, 'descend');
n_keep = min(n_keep, numel(units));
picked = units(order(1:n_keep));
end

function units = good_units_in_session(session_dir)
units = struct('cluster_id', {}, 'channel', {}, 'spike_times', {});
info_files = dir(fullfile(session_dir, 'cluster_info*.tsv'));
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
    spike_times = spike_times_to_samples(readNPY(times_file));
    spike_clusters = double(readNPY(cluster_file));
    spike_clusters = spike_clusters(:);
    ids = table_column(cluster_info, {'id', 'cluster_id'});
    channels = table_column(cluster_info, {'ch', 'channel'});
    groups = table_column(cluster_info, {'group', 'KSLabel'});
    for unit_index = 1:numel(ids)
        if ~strcmp(group_label(groups, unit_index), 'good')
            continue
        end
        cluster_id = scalar_number(ids, unit_index);
        channel = scalar_number(channels, unit_index) + 1;
        if ~isfinite(cluster_id) || ~isfinite(channel) || channel ~= round(channel) || channel < 1
            continue
        end
        samples = spike_times(spike_clusters == cluster_id);
        samples = samples(isfinite(samples) & samples == round(samples) & samples >= 1);
        if isempty(samples)
            continue
        end
        units(end + 1).cluster_id = cluster_id; %#ok<AGROW>
        units(end).channel = channel;
        units(end).spike_times = samples(:);
    end
end
end

function spike_samples = spike_times_to_samples(spike_times)
spike_samples = double(spike_times(:));
if isempty(spike_samples)
    return
end
if any(~isfinite(spike_samples)) || any(spike_samples ~= round(spike_samples))
    error('Check_NPX_folders_and_load_waveforms:spikeTime', ...
        'spike_times must be integer sample indices.');
end
if min(spike_samples) == 0
    spike_samples = spike_samples + 1;
elseif any(spike_samples < 0)
    error('Check_NPX_folders_and_load_waveforms:spikeTime', ...
        'spike_times must be non-negative sample indices.');
end
spike_samples = round(spike_samples);
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
error('Check_NPX_folders_and_load_waveforms:column', ...
    'cluster_info needs a column named %s.', strjoin(names, ' or '));
end

function label = group_label(groups, index)
if iscell(groups)
    value = groups{index};
elseif ischar(groups)
    value = strtrim(groups(index, :));
elseif exist('isstring', 'builtin') && isstring(groups)
    value = char(groups(index));
elseif exist('iscategorical', 'builtin') && iscategorical(groups)
    value = char(groups(index));
else
    value = groups(index);
end
if exist('isstring', 'builtin') && isstring(value)
    value = char(value);
end
if isnumeric(value)
    label = '';
    return
end
label = lower(strtrim(value));
end

function number = scalar_number(values, index)
if iscell(values)
    number = double(values{index});
else
    number = double(values(index));
end
end

function show_session_figure(units, continuous_file, waveform_cells, channel_cells, ...
    time_samples, fs, acg_window, acg_bin, n_acg_spikes)
[parent_dir, file_name, file_ext] = fileparts(continuous_file);
[~, parent_name] = fileparts(parent_dir);
figure_name = sprintf('%s  %s', units(1).session_id, ...
    [parent_name filesep file_name file_ext]);
n_units = numel(units);
fig = figure('Color', 'w', 'Name', figure_name);
set(fig, 'Units', 'pixels', 'Position', [80 40 980 max(280, 230 * n_units)]);
for unit_index = 1:n_units
    subplot(n_units, 2, (unit_index - 1) * 2 + 1);
    draw_waveform_stack(waveform_cells{unit_index}, channel_cells{unit_index}, ...
        time_samples, units(unit_index), fs);
    subplot(n_units, 2, (unit_index - 1) * 2 + 2);
    draw_autocorrelogram(units(unit_index), fs, acg_window, acg_bin, n_acg_spikes);
end
end

function draw_waveform_stack(waveforms, channels_kept, time_samples, unit, fs)
n_channels = size(waveforms, 2);
n_times = size(waveforms, 3);
centered = cell(1, n_channels);
mean_wf = nan(n_channels, n_times);
for channel_index = 1:n_channels
    traces = reshape(waveforms(:, channel_index, :), size(waveforms, 1), n_times);
    traces = center_snippets(traces, time_samples);
    centered{channel_index} = traces;
    mean_wf(channel_index, :) = mean_omitnan(traces);
end
scale = mean_waveform_scale(mean_wf);
gap = 2.6 * scale;
time_ms = time_samples(:)' / fs * 1000;
hold on
y_ticks = (0:n_channels - 1) * gap;
for channel_index = 1:n_channels
    plot(time_ms, centered{channel_index}' + y_ticks(channel_index), ...
        'Color', [0.75 0.75 0.75]);
end
for channel_index = 1:n_channels
    mean_color = [0 0 0];
    mean_width = 1.5;
    if channels_kept(channel_index) == unit.channel
        mean_color = [0.85 0.1 0.1];
        mean_width = 2.5;
    end
    plot(time_ms, mean_wf(channel_index, :) + y_ticks(channel_index), ...
        'Color', mean_color, 'LineWidth', mean_width);
end
padding = 1.25 * scale;
ylim([y_ticks(1) - padding, y_ticks(end) + padding]);
xlim([time_ms(1) time_ms(end)]);
tick_labels = cell(1, n_channels);
for channel_index = 1:n_channels
    tick_labels{channel_index} = sprintf('%g', channels_kept(channel_index));
end
set(gca, 'YTick', y_ticks, 'YTickLabel', tick_labels);
xlabel('Time (ms)')
ylabel('Channel')
title(sprintf('cluster %g   ch %g   %d spikes', ...
    unit.cluster_id, unit.channel, numel(unit.spike_times)));
box off
end

function traces = center_snippets(traces, time_samples)
time_samples = time_samples(:)';
pre = time_samples <= -10;
if ~any(pre)
    pre = false(1, numel(time_samples));
    pre(1) = true;
end
baseline = zeros(size(traces, 1), 1);
for spike_index = 1:size(traces, 1)
    values = traces(spike_index, pre);
    values = values(isfinite(values));
    if isempty(values)
        baseline(spike_index) = 0;
    else
        baseline(spike_index) = mean(values);
    end
end
traces = bsxfun(@minus, traces, baseline);
end

function scale = mean_waveform_scale(mean_wf)
scale = max(abs(mean_wf(:)));
if isempty(scale) || ~isfinite(scale) || scale == 0
    scale = 1;
end
end

function draw_autocorrelogram(unit, fs, acg_window, acg_bin, n_acg_spikes)
spike_seconds = unit.spike_times / fs;
if numel(spike_seconds) > n_acg_spikes
    spike_seconds = spike_seconds(randperm(numel(spike_seconds), n_acg_spikes));
end
lags = autocorrelogram(spike_seconds, acg_window);
edges = acg_window(1):acg_bin:acg_window(2);
if numel(edges) < 2
    edges = [acg_window(1) acg_window(2)];
end
[counts, edges] = histcounts(lags, edges);
centers = (edges(1:end - 1) + edges(2:end)) / 2;
bar(centers, counts, 1, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none');
xlim(acg_window);
xlabel('Lag (s)')
ylabel('Count')
title(sprintf('cluster %g autocorrelogram', unit.cluster_id));
box off
end

function row = mean_omitnan(matrix)
n_times = size(matrix, 2);
row = nan(1, n_times);
for time_index = 1:n_times
    values = matrix(:, time_index);
    values = values(isfinite(values));
    if ~isempty(values)
        row(time_index) = mean(values);
    end
end
end
