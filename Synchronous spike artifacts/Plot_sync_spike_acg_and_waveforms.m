%% Waveform profile and autocorrelogram of the extracted synchronous spikes
% Run this after Estimate_sync_spike_waveforms. It does not repeat the
% 500 null draws. Both summary stacks keep only neurons that took part in a
% synchronous event in DLPAG, LPAG, VLPAG, DR, and SupCol. Those curves
% come from the saved analysis. Good units that were never in a synchronous
% pair are left out, and the continuous file is not read again.
%
% The waveform is the saved synchronous mean, drawn on every channel that
% was kept, low channel at the bottom. The non-synchronous in-window mean
% is drawn on the same channels. Beside it, the peak-to-peak range of that
% synchronous mean is plotted against signed Euclidean distance from the
% channel with the largest peak-to-peak range, in micrometers. Distance
% uses xcoords and ycoords from chann_map_PAG.mat. That peak channel is
% at 0. Larger y is positive.
% The autocorrelogram uses only spikes that belong to a synchronous event,
% reloaded from the sorting and the play-bout windows. ±50 ms, 1 ms bins.
% The zero lag is omitted. The summary figure stacks one normalized
% autocorrelogram for every such neuron. Each row is that neuron's
% synchronous-spike pair counts divided by their sum. A neuron that belongs
% to several pairs contributes the union of its synchronous spikes once.
% The upper panel is the stack. The lower panel is every neuron as a thin
% line, the mean across neurons, and the 95% confidence interval of that mean.
%
% max_example_neurons limits the detailed panels. Those panels show the
% first neurons in saved order, one row each: waveform, spatial range, and
% autocorrelogram. The waveform panel keeps example_channel_half_width
% channels on either side of the largest-range channel, even when the saved
% waveform covers a wider set. The other synchronous-event neurons stay in
% the two stacks and are not drawn as individual examples. Set the neuron
% cap to inf to show every neuron once.
% The spatial stack is a separate figure. Each curve is the peak-to-peak
% range of that neuron's saved synchronous mean, divided by the range on
% its largest-range channel. Distance 0 is that channel.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

acg_window = [-0.05 0.05];
acg_bin = 0.001;
max_example_neurons = 3;
example_channel_half_width = 3;
examples_per_figure = 6;
spatial_distance_um = -400:20:400;

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
neuron_spikes = empty_neuron_spikes();
spatial_profiles = empty_spatial_profiles();
examples = empty_examples();
for session_index = 1:numel(session_ids)
    session_id = session_ids{session_index};
    session_rows = neuron_rows(strcmp({neuron_rows.session_id}, session_id));
    [trains, bout_starts, time_options] = load_session_trains( ...
        npx_raw_data, session_id, session_rows(1));
    neuron_spikes = remember_analysis_neurons(neuron_spikes, session_id, ...
        session_rows, trains, bout_starts, time_options);
    geometry = load_channel_geometry(npx_raw_data, session_id);
    spatial_profiles = append_spatial_profiles(spatial_profiles, ...
        collect_session_profiles(session_id, session_rows, geometry, time_options));
    examples = remember_example_neurons(examples, max_example_neurons, ...
        session_rows, trains, bout_starts, time_options, geometry);
end
fprintf(['Example panels: %d neurons (cap %g). The stacks keep every ' ...
    'synchronous-event neuron in the analysis areas.\n'], ...
    numel(examples), max_example_neurons);
draw_example_figures(examples, examples_per_figure, acg_window, acg_bin, ...
    example_channel_half_width);
[acg_centers, acg_counts] = neuron_acg_counts(neuron_spikes, acg_window, acg_bin);
draw_normalized_acg_figure(acg_centers, acg_counts, acg_window);
draw_spatial_stack_figure(spatial_profiles, spatial_distance_um);

function examples = empty_examples()
examples = struct('row', {}, 'spike_times', {}, 'geometry', {});
end

function examples = remember_example_neurons(examples, max_neurons, session_rows, ...
    trains, bout_starts, time_options, geometry)
if ~isfinite(max_neurons)
    max_neurons = inf;
end
max_neurons = floor(double(max_neurons));
if max_neurons < 1 || numel(examples) >= max_neurons
    return
end
areas = analysis_areas(time_options);
saved_ids = [];
if ~isempty(session_rows) && isfield(session_rows, 'cluster_id')
    saved_ids = [session_rows.cluster_id];
end
for index = 1:numel(saved_ids)
    if numel(examples) >= max_neurons
        return
    end
    cluster_id = saved_ids(index);
    if index > 1 && any(saved_ids(1:index - 1) == cluster_id)
        continue
    end
    members = session_rows([session_rows.cluster_id] == cluster_id);
    if ~members_in_analysis_area(members, areas)
        continue
    end
    row = first_area_row(members, areas);
    if isempty(row)
        continue
    end
    examples(end + 1).row = row; %#ok<AGROW>
    examples(end).spike_times = synchronous_times(trains, row, bout_starts, time_options);
    examples(end).geometry = geometry;
end
end

function row = first_area_row(members, areas)
row = [];
for index = 1:numel(members)
    if row_in_analysis_area(members(index), areas)
        row = members(index);
        return
    end
end
end

function draw_example_figures(examples, per_figure, acg_window, acg_bin, channel_half_width)
if isempty(examples)
    return
end
if ~isfinite(per_figure) || per_figure < 1
    per_figure = numel(examples);
end
per_figure = floor(double(per_figure));
n_figures = ceil(numel(examples) / per_figure);
for figure_index = 1:n_figures
    first_neuron = (figure_index - 1) * per_figure + 1;
    last_neuron = min(numel(examples), figure_index * per_figure);
    shown = examples(first_neuron:last_neuron);
    figure('Color', 'w', ...
        'Name', sprintf('Synchronous-spike examples %d/%d', figure_index, n_figures), ...
        'Position', [40 40 1100 240 * numel(shown)]);
    for plot_index = 1:numel(shown)
        row = shown(plot_index).row;
        subplot(numel(shown), 3, (plot_index - 1) * 3 + 1);
        draw_channel_profile(row, figure_index == 1 && plot_index == 1, channel_half_width);
        subplot(numel(shown), 3, (plot_index - 1) * 3 + 2);
        draw_spatial_range(row, shown(plot_index).geometry);
        subplot(numel(shown), 3, (plot_index - 1) * 3 + 3);
        spike_times = shown(plot_index).spike_times;
        [centers, counts] = sync_event_acg(spike_times, acg_window, acg_bin);
        draw_sync_acg(row, centers, counts, spike_times, acg_window);
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
    fprintf('%s  reloaded %d spike trains for synchronous-event autocorrelograms\n', ...
        session_id, numel(trains));
catch load_error
    warning('Plot_sync_spike_acg_and_waveforms:times', ...
        'Spike times for %s were not reloaded (%s). Waveforms are still drawn.', ...
        session_id, load_error.message);
end
end

function draw_channel_profile(row, show_legend, channel_half_width)
if nargin < 3
    channel_half_width = inf;
end
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
[channels, sync_mean, reference] = limit_displayed_channels( ...
    channels, sync_mean, reference, row.channel, channel_half_width);
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
peak_channel = largest_range_channel(channels, sync_mean, row.channel);
hold on
sync_line = [];
reference_line = [];
for channel_index = 1:n_channels
    y_offset = (channel_index - 1) * gap;
    line_width = 1.2;
    if channels(channel_index) == peak_channel
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
title(sprintf('%s   cluster %g   %s   ch %g   with %s   %d sync', ...
    row.session_id, row.cluster_id, row_area(row), row.channel, ...
    partner_label(row), row_count(row, 'n_sync')), ...
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

function [channels, sync_mean, reference] = limit_displayed_channels( ...
    channels, sync_mean, reference, assigned_channel, half_width)
channels = double(channels(:).');
if ~isfinite(half_width) || half_width < 0 || isempty(channels)
    return
end
half_width = floor(double(half_width));
center = largest_range_channel(channels, sync_mean, assigned_channel);
keep = abs(channels - center) <= half_width;
if ~any(keep)
    keep = abs(channels - double(assigned_channel)) <= half_width;
end
if ~any(keep)
    return
end
n_channels = numel(channels);
channels = channels(keep);
sync_mean = sync_mean(keep, :);
if ~isempty(reference) && size(reference, 1) == n_channels
    reference = reference(keep, :);
end
end

function channel = largest_range_channel(channels, sync_mean, assigned_channel)
channel = double(assigned_channel);
ranges = nan(1, numel(channels));
for index = 1:numel(channels)
    samples = sync_mean(index, :);
    samples = samples(isfinite(samples));
    if isempty(samples)
        continue
    end
    ranges(index) = max(samples) - min(samples);
end
finite_range = find(isfinite(ranges));
if isempty(finite_range)
    return
end
[~, best] = max(ranges(finite_range));
channel = channels(finite_range(best));
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
peak_point = profile.channel == profile.center_channel;
if any(peak_point)
    plot(profile.distance_um(peak_point), profile.range(peak_point), 'o', ...
        'Color', [0.85 0.1 0.1], 'MarkerFaceColor', [0.85 0.1 0.1], 'MarkerSize', 8);
end
assigned_point = profile.channel == row.channel;
if any(assigned_point) && profile.center_channel ~= row.channel
    plot(profile.distance_um(assigned_point), profile.range(assigned_point), 'o', ...
        'Color', [0 0 0], 'MarkerSize', 8);
end
ylim([0, peak * 1.15]);
xlim(span * [-1.15 1.15]);
xlabel('Distance from peak channel (\mum)')
ylabel('Peak-to-peak range')
if profile.center_channel == row.channel
    range_title = sprintf('cluster %g   ch %g   with %s', ...
        row.cluster_id, row.channel, partner_label(row));
else
    range_title = sprintf('cluster %g   peak ch %g   assigned %g   with %s', ...
        row.cluster_id, profile.center_channel, row.channel, partner_label(row));
end
title(range_title, 'Interpreter', 'none');
box off
end

function entries = empty_neuron_spikes()
entries = struct('session_id', {}, 'cluster_id', {}, 'spike_times', {});
end

function entries = remember_analysis_neurons(entries, session_id, session_rows, ...
    trains, bout_starts, time_options)
n_before = numel(entries);
areas = analysis_areas(time_options);
saved_ids = [];
if ~isempty(session_rows) && isfield(session_rows, 'cluster_id')
    saved_ids = [session_rows.cluster_id];
end
n_analyzed = 0;
for index = 1:numel(saved_ids)
    cluster_id = saved_ids(index);
    if index > 1 && any(saved_ids(1:index - 1) == cluster_id)
        continue
    end
    members = session_rows([session_rows.cluster_id] == cluster_id);
    if ~members_in_analysis_area(members, areas)
        continue
    end
    n_analyzed = n_analyzed + 1;
    for member_index = 1:numel(members)
        row = members(member_index);
        if ~isfield(row, 'partner_id') || ~row_in_analysis_area(row, areas)
            continue
        end
        spike_times = synchronous_times(trains, row, bout_starts, time_options);
        entries = remember_neuron_spikes(entries, session_id, cluster_id, spike_times);
    end
end
n_neurons = numel(entries) - n_before;
fprintf(['%s  autocorrelogram stack  %d neurons with synchronous spikes, ' ...
    '%d neurons saved in the analysis areas\n'], session_id, n_neurons, n_analyzed);
end

function entries = remember_neuron_spikes(entries, session_id, cluster_id, spike_times)
if isempty(spike_times)
    return
end
for index = 1:numel(entries)
    if strcmp(entries(index).session_id, session_id) && entries(index).cluster_id == cluster_id
        entries(index).spike_times = unique([entries(index).spike_times(:); spike_times(:)]);
        return
    end
end
entries(end + 1).session_id = session_id; %#ok<AGROW>
entries(end).cluster_id = cluster_id;
entries(end).spike_times = unique(spike_times(:));
end

function [centers, counts] = neuron_acg_counts(entries, acg_window, acg_bin)
centers = [];
counts = [];
for index = 1:numel(entries)
    [these_centers, these_counts] = sync_event_acg(entries(index).spike_times, acg_window, acg_bin);
    if isempty(these_counts) || ~any(these_counts > 0)
        continue
    end
    if isempty(counts)
        centers = these_centers(:)';
        counts = these_counts(:)';
    elseif numel(these_counts) == size(counts, 2)
        counts(end + 1, :) = these_counts(:)'; %#ok<AGROW>
    end
end
end

function draw_normalized_acg_figure(centers, counts, acg_window)
if isempty(counts)
    warning('Plot_sync_spike_acg_and_waveforms:acg', ...
        ['No synchronous-spike autocorrelograms were available. The stack uses ' ...
         'only synchronous-event spikes of neurons in DLPAG, LPAG, VLPAG, DR, and SupCol.']);
    return
end
summary = normalized_acg_stack(counts);
if summary.n < 1
    return
end
figure('Color', 'w', 'Name', 'Synchronous-spike autocorrelograms', ...
    'Position', [80 40 760 860]);
subplot(2, 1, 1);
imagesc(centers, 1:summary.n, summary.normalized);
set(gca, 'YDir', 'reverse');
xlim(acg_window);
xlabel('Lag (s)')
ylabel('Neuron, sorted by peak lag')
title(sprintf('Synchronous spikes only, %d neurons in the analysis areas', summary.n));
color_bar = colorbar;
ylabel(color_bar, 'Fraction of spike pairs');
box off

subplot(2, 1, 2);
hold on
for row_index = 1:summary.n
    plot(centers, summary.normalized(row_index, :), ...
        'Color', [0.85 0.55 0.55], 'LineWidth', 0.4);
end
fill([centers, fliplr(centers)], [summary.ci_low, fliplr(summary.ci_high)], ...
    [0.85 0.1 0.1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(centers, summary.mean_acg, 'Color', [0.85 0.1 0.1], 'LineWidth', 2);
xlim(acg_window);
xlabel('Lag (s)')
ylabel('Fraction of spike pairs')
title(sprintf('Mean synchronous-spike autocorrelogram, 95%% CI, n = %d', summary.n));
box off
end

function profiles = empty_spatial_profiles()
profiles = struct('session_id', {}, 'cluster_id', {}, 'distance_um', {}, 'range', {});
end

function profiles = append_spatial_profiles(profiles, more)
if isempty(more)
    return
end
if isempty(profiles)
    profiles = more;
    return
end
profiles = [profiles, more];
end

function profiles = collect_session_profiles(session_id, session_rows, geometry, time_options)
profiles = empty_spatial_profiles();
if isempty(geometry)
    return
end
areas = analysis_areas(time_options);
saved_ids = [];
if ~isempty(session_rows) && isfield(session_rows, 'cluster_id')
    saved_ids = [session_rows.cluster_id];
end
n_kept = 0;
for index = 1:numel(saved_ids)
    cluster_id = saved_ids(index);
    if index > 1 && any(saved_ids(1:index - 1) == cluster_id)
        continue
    end
    members = session_rows([session_rows.cluster_id] == cluster_id);
    if ~members_in_analysis_area(members, areas)
        continue
    end
    [distance_um, range] = profile_from_saved(members, geometry);
    if isempty(distance_um)
        continue
    end
    profiles = add_spatial_profile(profiles, session_id, cluster_id, distance_um, range);
    n_kept = n_kept + 1;
end
fprintf(['%s  spatial stack  %d synchronous-event neurons in the analysis areas\n'], ...
    session_id, n_kept);
end

function [distance_um, range] = profile_from_saved(members, geometry)
distance_um = [];
range = [];
channels = [];
ranges = [];
neuron_channel = NaN;
for index = 1:numel(members)
    row = members(index);
    if ~isfield(row, 'sync_mean') || isempty(row.sync_mean) ...
            || ~isfield(row, 'channels_kept') || isempty(row.channels_kept)
        continue
    end
    profile = sync_waveform_spatial_range(row.channels_kept, row.sync_mean, ...
        row.channel, geometry.chan_map, geometry.xcoords, geometry.ycoords);
    channels = [channels; profile.channel(:)]; %#ok<AGROW>
    ranges = [ranges; profile.range(:)]; %#ok<AGROW>
    if isfield(row, 'channel')
        neuron_channel = row.channel;
    end
end
if isempty(channels)
    return
end
[channels, ranges] = average_channel_ranges(channels, ranges);
span = nan(numel(channels), 2);
finite_range = isfinite(ranges);
span(finite_range, 1) = 0;
span(finite_range, 2) = ranges(finite_range);
profile = sync_waveform_spatial_range(channels, span, neuron_channel, ...
    geometry.chan_map, geometry.xcoords, geometry.ycoords);
distance_um = profile.distance_um;
range = profile.range;
end

function [channels_out, ranges_out] = average_channel_ranges(channels, ranges)
channels = double(channels(:));
ranges = double(ranges(:));
keep = isfinite(channels) & isfinite(ranges);
channels = channels(keep);
ranges = ranges(keep);
if isempty(channels)
    channels_out = zeros(0, 1);
    ranges_out = zeros(0, 1);
    return
end
[channels_out, ~, group] = unique(channels);
ranges_out = zeros(size(channels_out));
for index = 1:numel(channels_out)
    ranges_out(index) = mean(ranges(group == index));
end
end

function profiles = add_spatial_profile(profiles, session_id, cluster_id, distance_um, range)
profiles(end + 1).session_id = session_id; %#ok<AGROW>
profiles(end).cluster_id = cluster_id;
profiles(end).distance_um = distance_um(:);
profiles(end).range = range(:);
end

function draw_spatial_stack_figure(profiles, centers)
if isempty(profiles)
    warning('Plot_sync_spike_acg_and_waveforms:spatial', ...
        ['No spatial profiles were available. The stack keeps synchronous-event ' ...
         'neurons in DLPAG, LPAG, VLPAG, DR, and SupCol.']);
    return
end
distances = cell(1, numel(profiles));
ranges = cell(1, numel(profiles));
for index = 1:numel(profiles)
    distances{index} = profiles(index).distance_um;
    ranges{index} = profiles(index).range;
end
summary = spatial_range_stack(distances, ranges, centers);
if summary.n < 1
    warning('Plot_sync_spike_acg_and_waveforms:spatial', ...
        'No neuron had a positive range on its peak channel.');
    return
end
used = any(isfinite(summary.normalized), 1) | any(isfinite(summary.mean_range));
if ~any(used)
    return
end
centers = summary.centers(used);
normalized = summary.normalized(:, used);
mean_range = summary.mean_range(used);
ci_low = summary.ci_low(used);
ci_high = summary.ci_high(used);
span = max(abs(centers));
if ~isfinite(span) || span == 0
    span = 1;
end
figure('Color', 'w', 'Name', 'Normalized spatial range', ...
    'Position', [120 40 760 860]);
subplot(2, 1, 1);
image_handle = imagesc(centers, 1:summary.n, normalized);
set(image_handle, 'AlphaData', double(isfinite(normalized)));
set(gca, 'YDir', 'reverse');
xlim(span * [-1.05 1.05]);
xlabel('Distance from peak channel (\mum)')
ylabel('Neuron')
title(sprintf('Synchronous-event neurons, range / peak channel, n = %d', summary.n));
color_bar = colorbar;
ylabel(color_bar, 'Normalized range');
color_peak = max(normalized(isfinite(normalized)));
if isempty(color_peak) || ~isfinite(color_peak) || color_peak < 1
    color_peak = 1;
elseif color_peak > 2
    color_peak = 2;
end
caxis([0 color_peak]);
box off

subplot(2, 1, 2);
hold on
for row_index = 1:summary.n
    plot(centers, normalized(row_index, :), 'Color', [0.85 0.55 0.55], 'LineWidth', 0.4);
end
fill([centers, fliplr(centers)], [ci_low, fliplr(ci_high)], ...
    [0.85 0.1 0.1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(centers, mean_range, 'Color', [0.85 0.1 0.1], 'LineWidth', 2);
xlim(span * [-1.05 1.05]);
xlabel('Distance from peak channel (\mum)')
ylabel('Range / range at the peak channel')
title(sprintf('Mean spatial range, centered on the peak channel, 95%% CI, n = %d', summary.n));
box off
end

function draw_sync_acg(row, centers, counts, spike_times, acg_window)
bar(centers, counts, 1, 'FaceColor', [0.85 0.1 0.1], 'EdgeColor', 'none');
xlim(acg_window);
xlabel('Lag (s)')
ylabel('Count')
if isempty(spike_times)
    title(sprintf('cluster %g   with %s   autocorrelogram unavailable', ...
        row.cluster_id, partner_label(row)), 'Interpreter', 'none');
else
    title(sprintf('cluster %g   %d synchronous spikes with %s', ...
        row.cluster_id, numel(spike_times), partner_label(row)), ...
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

function areas = analysis_areas(time_options)
areas = sync_analysis_options().areas_to_keep;
if isstruct(time_options) && isfield(time_options, 'areas_to_keep') ...
        && ~isempty(time_options.areas_to_keep)
    areas = time_options.areas_to_keep;
end
end

function keep = row_in_analysis_area(row, areas)
if ~isfield(row, 'area')
    keep = true;
    return
end
keep = area_label_kept(row.area, areas);
end

function keep = members_in_analysis_area(members, areas)
labeled = false;
for index = 1:numel(members)
    if ~isfield(members(index), 'area')
        continue
    end
    name = area_text(members(index).area);
    if isempty(name)
        continue
    end
    labeled = true;
    if area_matches(name, areas)
        keep = true;
        return
    end
end
keep = ~labeled;
end

function keep = area_label_kept(value, areas)
name = area_text(value);
if isempty(name)
    keep = true;
    return
end
keep = area_matches(name, areas);
end

function matches = area_matches(name, areas)
matches = false;
if ischar(areas)
    matches = strcmp(strtrim(areas), name);
    return
end
if exist('isstring', 'builtin') && isstring(areas)
    areas = cellstr(areas);
end
if iscell(areas)
    for index = 1:numel(areas)
        if strcmp(area_text(areas{index}), name)
            matches = true;
            return
        end
    end
end
end

function name = area_text(value)
name = '';
if exist('isstring', 'builtin') && isstring(value)
    if numel(value) < 1
        return
    end
    value = char(value(1));
end
if iscell(value)
    if isempty(value)
        return
    end
    value = value{1};
    if exist('isstring', 'builtin') && isstring(value)
        value = char(value);
    end
end
if ischar(value)
    name = strtrim(value(:).');
end
end

function label = partner_label(row)
partner_id = NaN;
if isfield(row, 'partner_id') && ~isempty(row.partner_id) && isnumeric(row.partner_id)
    partner_id = row.partner_id(1);
end
area = '';
if isfield(row, 'partner_area')
    area = area_text(row.partner_area);
end
partner_channel = NaN;
if isfield(row, 'partner_channel') && ~isempty(row.partner_channel) ...
        && isnumeric(row.partner_channel)
    partner_channel = row.partner_channel(1);
end
if isfinite(partner_id) && ~isempty(area) && isfinite(partner_channel)
    label = sprintf('partner %g %s ch %g', partner_id, area, partner_channel);
elseif isfinite(partner_id)
    label = sprintf('partner %g', partner_id);
else
    label = 'partner unknown';
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
