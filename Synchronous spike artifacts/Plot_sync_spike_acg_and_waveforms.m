%% Waveform profile and autocorrelogram of the extracted synchronous spikes
% Run this after Estimate_sync_spike_waveforms. It does not repeat the
% 500 null draws and it does not read the continuous file again.
%
% The waveform is the saved synchronous mean, drawn on every channel that
% was kept, low channel at the bottom. The non-synchronous in-window mean
% is drawn on the same channels. Beside it, the peak-to-peak range of that
% synchronous mean is plotted against signed Euclidean distance from the
% neuron channel, in micrometers. Distance uses xcoords and ycoords from
% chann_map_PAG.mat. The neuron channel is at 0. Larger y is positive.
% The autocorrelogram uses the synchronous spike times of that pair,
% reloaded from the sorting and the play-bout windows. ±50 ms, 1 ms bins.
% The zero lag is omitted.
%
% One figure per session, six pairs per figure. Each pair is one row:
% neuron 1 waveform, neuron 1 spatial range, neuron 1 autocorrelogram,
% then the same three panels for neuron 2.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

acg_window = [-0.05 0.05];
acg_bin = 0.001;
pairs_per_figure = 6;

data_root = sync_analysis_options().data_root;
saving_folder = fullfile(data_root, 'Analysis results', 'synchronous spike waveforms');
npx_raw_data = fullfile(data_root, 'NPX data', 'NPX raw data');
session_files = dir(fullfile(saving_folder, '*_sync_spike_waveforms.mat'));
if isempty(session_files)
    error('Plot_sync_spike_acg_and_waveforms:files', ...
        'No per-session waveform files in %s.', saving_folder);
end

neuron_rows = struct([]);
for file_index = 1:numel(session_files)
    session_path = fullfile(session_files(file_index).folder, session_files(file_index).name);
    loaded = load(session_path, 'this_waveforms');
    if ~isfield(loaded, 'this_waveforms') || ~isfield(loaded.this_waveforms, 'rows') ...
            || isempty(loaded.this_waveforms.rows)
        continue
    end
    these_rows = loaded.this_waveforms.rows;
    if ~isfield(these_rows, 'session_id') || isempty(these_rows(1).session_id)
        [these_rows.session_id] = deal(loaded.this_waveforms.session_id);
    end
    these_rows = attach_plot_options(these_rows, loaded.this_waveforms);
    neuron_rows = append_struct_array(neuron_rows, these_rows);
end
if isempty(neuron_rows)
    error('Plot_sync_spike_acg_and_waveforms:rows', ...
        'The waveform files in %s have no neuron rows.', saving_folder);
end

session_ids = unique({neuron_rows.session_id}, 'stable');
for session_index = 1:numel(session_ids)
    session_id = session_ids{session_index};
    session_rows = neuron_rows(strcmp({neuron_rows.session_id}, session_id));
    [trains, bout_starts, time_options] = load_session_trains( ...
        npx_raw_data, session_id, session_rows(1));
    geometry = load_channel_geometry(npx_raw_data, session_id);
    pair_ids = unique([session_rows.pair_index], 'stable');
    n_figures = ceil(numel(pair_ids) / pairs_per_figure);
    for figure_index = 1:n_figures
        first_pair = (figure_index - 1) * pairs_per_figure + 1;
        last_pair = min(numel(pair_ids), figure_index * pairs_per_figure);
        shown = pair_ids(first_pair:last_pair);
        figure('Color', 'w', ...
            'Name', sprintf('%s synchronous spikes %d/%d', session_id, figure_index, n_figures), ...
            'Position', [40 40 1650 240 * numel(shown)]);
        for plot_index = 1:numel(shown)
            members = session_rows([session_rows.pair_index] == shown(plot_index));
            for side = 1:2
                panel = (plot_index - 1) * 6 + (side - 1) * 3;
                subplot(numel(shown), 6, panel + 1);
                if side <= numel(members)
                    draw_channel_profile(members(side), side == 1 && plot_index == 1);
                    subplot(numel(shown), 6, panel + 2);
                    draw_spatial_range(members(side), geometry);
                    subplot(numel(shown), 6, panel + 3);
                    draw_sync_acg(members(side), trains, bout_starts, time_options, ...
                        acg_window, acg_bin);
                else
                    title('No neuron');
                    axis off
                    subplot(numel(shown), 6, panel + 2);
                    axis off
                    subplot(numel(shown), 6, panel + 3);
                    axis off
                end
            end
        end
    end
end

function rows = attach_plot_options(rows, session_result)
if ~isfield(session_result, 'options')
    return
end
for index = 1:numel(rows)
    rows(index).plot_options = session_result.options;
end
end

function [trains, bout_starts, time_options] = load_session_trains(npx_raw_data, session_id, row)
trains = [];
bout_starts = [];
time_options = sync_analysis_options();
if isfield(row, 'plot_options') && isstruct(row.plot_options)
    time_options = sync_analysis_options(row.plot_options);
end
time_options.spikes_only = true;
time_options.max_pairs = inf;
session_dir = fullfile(npx_raw_data, session_id);
try
    selection = synchronous_play_pairs(session_dir, time_options);
    trains = selection.all_neurons;
    bout_starts = selection.bout_starts;
    fprintf('%s  reloaded %d neurons for autocorrelograms, no continuous-file read\n', ...
        session_id, numel(trains));
catch load_error
    warning('Plot_sync_spike_acg_and_waveforms:times', ...
        'Spike times for %s were not reloaded (%s). Waveforms are still drawn.', ...
        session_id, load_error.message);
end
end

function draw_channel_profile(row, show_legend)
sync_mean = [];
if isfield(row, 'sync_mean')
    sync_mean = row.sync_mean;
end
reference = reference_mean(row);
channels = [];
if isfield(row, 'channels_kept')
    channels = row.channels_kept;
end
if isempty(sync_mean) || isempty(channels)
    title(missing_title(row));
    axis off
    return
end
[channels, sync_mean, reference] = order_channels(channels, sync_mean, reference);
time_samples = row.time_samples(:)';
time_ms = time_samples / 30000 * 1000;
scale = max(abs(sync_mean(:)));
if ~isempty(reference)
    scale = max(scale, max(abs(reference(:))));
end
if isempty(scale) || ~isfinite(scale) || scale == 0
    scale = 1;
end
gap = 2.6 * scale;
n_channels = numel(channels);
hold on
sync_line = [];
reference_line = [];
for channel_index = 1:n_channels
    y_offset = (channel_index - 1) * gap;
    line_width = 1.2;
    if channels(channel_index) == row.channel
        line_width = 2.4;
    end
    if ~isempty(reference)
        reference_line = plot(time_ms, reference(channel_index, :) + y_offset, ...
            'Color', [0 0 0], 'LineWidth', line_width);
    end
    sync_line = plot(time_ms, sync_mean(channel_index, :) + y_offset, ...
        'Color', [0.85 0.1 0.1], 'LineWidth', line_width);
end
padding = 1.25 * scale;
ylim([-padding, (n_channels - 1) * gap + padding]);
if numel(time_ms) >= 2
    xlim([time_ms(1) time_ms(end)]);
end
y_ticks = (0:n_channels - 1) * gap;
tick_labels = cell(1, n_channels);
for channel_index = 1:n_channels
    tick_labels{channel_index} = sprintf('%g', channels(channel_index));
end
set(gca, 'YTick', y_ticks, 'YTickLabel', tick_labels);
xlabel('Time (ms)')
ylabel('Channel')
title(sprintf('%s   cluster %g   %s   ch %g   %d sync', ...
    row.session_id, row.cluster_id, row_area(row), row.channel, row_count(row, 'n_sync')), ...
    'Interpreter', 'none');
if show_legend && ~isempty(sync_line)
    if ~isempty(reference_line)
        legend([reference_line, sync_line], ...
            {'Non-synchronous mean', 'Synchronous mean'}, 'Location', 'northeast');
    else
        legend(sync_line, {'Synchronous mean'}, 'Location', 'northeast');
    end
end
box off
end

function geometry = load_channel_geometry(npx_raw_data, session_id)
geometry = [];
session_dir = fullfile(npx_raw_data, session_id);
try
    [chan_map, xcoords, ycoords] = read_npx_channel_geometry(session_dir);
    geometry.chan_map = chan_map;
    geometry.xcoords = xcoords;
    geometry.ycoords = ycoords;
    fprintf('%s  loaded %d probe sites for the spatial range\n', ...
        session_id, numel(chan_map));
catch load_error
    warning('Plot_sync_spike_acg_and_waveforms:geometry', ...
        'Channel positions for %s were not loaded (%s). Spatial range is skipped.', ...
        session_id, load_error.message);
end
end

function draw_spatial_range(row, geometry)
if isempty(geometry) || ~isfield(row, 'sync_mean') || isempty(row.sync_mean) ...
        || ~isfield(row, 'channels_kept') || isempty(row.channels_kept) ...
        || ~isfield(row, 'channel')
    title(missing_title(row));
    axis off
    return
end
profile = sync_waveform_spatial_range(row.channels_kept, row.sync_mean, row.channel, ...
    geometry.chan_map, geometry.xcoords, geometry.ycoords);
if isempty(profile.distance_um)
    title(sprintf('cluster %g   no probe position', row.cluster_id), 'Interpreter', 'none');
    axis off
    return
end
peak = max(profile.range);
if ~isfinite(peak) || peak <= 0
    peak = 1;
end
span = max(abs(profile.distance_um));
if ~isfinite(span) || span == 0
    span = 1;
end
hold on
plot([0 0], [0, peak * 1.15], 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
plot(profile.distance_um, profile.range, '-o', ...
    'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1], ...
    'LineWidth', 1.4, 'MarkerSize', 5);
neuron_point = profile.channel == row.channel;
if any(neuron_point)
    plot(profile.distance_um(neuron_point), profile.range(neuron_point), 'o', ...
        'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerSize', 8);
end
ylim([0, peak * 1.15]);
xlim(span * [-1.15 1.15]);
xlabel('Distance from neuron channel (\mum)')
ylabel('Peak-to-peak range')
title(sprintf('cluster %g   ch %g', row.cluster_id, row.channel), 'Interpreter', 'none');
box off
end

function draw_sync_acg(row, trains, bout_starts, time_options, acg_window, acg_bin)
spike_times = synchronous_times(trains, row, bout_starts, time_options);
[centers, counts] = sync_event_acg(spike_times, acg_window, acg_bin);
bar(centers, counts, 1, 'FaceColor', [0.85 0.1 0.1], 'EdgeColor', 'none');
xlim(acg_window);
xlabel('Lag (s)')
ylabel('Count')
if isempty(spike_times)
    title(sprintf('cluster %g   autocorrelogram unavailable', row.cluster_id), ...
        'Interpreter', 'none');
else
    title(sprintf('cluster %g   %d synchronous spikes', row.cluster_id, numel(spike_times)), ...
        'Interpreter', 'none');
end
box off
end

function spike_times = synchronous_times(trains, row, bout_starts, time_options)
spike_times = [];
if isempty(trains) || isempty(bout_starts) || ~isfield(row, 'partner_id')
    return
end
neuron = neuron_by_id(trains, row.cluster_id);
partner = neuron_by_id(trains, row.partner_id);
if isempty(neuron) || isempty(partner)
    return
end
times = sort(double(neuron.spike_times(:)));
times = times(isfinite(times));
partner_times = sort(double(partner.spike_times(:)));
partner_times = partner_times(isfinite(partner_times));
mask = sync_spike_mask(times, partner_times, bout_starts, ...
    time_options.hist_range, time_options.time_precision);
spike_times = times(mask);
end

function neuron = neuron_by_id(neurons, cluster_id)
neuron = [];
for index = 1:numel(neurons)
    if neurons(index).cluster_id == cluster_id
        neuron = neurons(index);
        return
    end
end
end

function [channels, sync_mean, reference] = order_channels(channels, sync_mean, reference)
[channels, order] = sort(double(channels(:)'));
sync_mean = sync_mean(order, :);
if ~isempty(reference) && size(reference, 1) == numel(order)
    reference = reference(order, :);
end
end

function waveform = reference_mean(row)
waveform = [];
if isfield(row, 'nonsync_mean') && ~isempty(row.nonsync_mean)
    waveform = row.nonsync_mean;
elseif isfield(row, 'template_mean') && ~isempty(row.template_mean)
    waveform = row.template_mean;
end
end

function label = missing_title(row)
if isfield(row, 'cluster_id')
    label = sprintf('cluster %g has no saved waveform', row.cluster_id);
else
    label = 'No saved waveform';
end
end

function label = row_area(row)
if isfield(row, 'area') && ~isempty(row.area)
    label = row.area;
else
    label = '';
end
end

function value = row_count(row, field_name)
value = 0;
if isfield(row, field_name) && ~isempty(row.(field_name)) && isfinite(row.(field_name))
    value = row.(field_name);
end
end
