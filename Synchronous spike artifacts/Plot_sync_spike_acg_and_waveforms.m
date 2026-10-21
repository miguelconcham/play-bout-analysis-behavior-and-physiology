%% Waveform profile and autocorrelogram of the extracted synchronous spikes
% Run this after Estimate_sync_spike_waveforms. It does not repeat the
% 500 null draws. Paired neurons reuse the saved synchronous mean. Other
% good units are read from the continuous file, at most 100 spikes each.
%
% The waveform is the saved synchronous mean, drawn on every channel that
% was kept, low channel at the bottom. The non-synchronous in-window mean
% is drawn on the same channels. Beside it, the peak-to-peak range of that
% synchronous mean is plotted against signed Euclidean distance from the
% neuron channel, in micrometers. Distance uses xcoords and ycoords from
% chann_map_PAG.mat. The neuron channel is at 0. Larger y is positive.
% The autocorrelogram uses only spikes that belong to a synchronous event,
% reloaded from the sorting and the play-bout windows. ±50 ms, 1 ms bins.
% The zero lag is omitted. The summary figure stacks one normalized
% autocorrelogram for every neuron saved by the synchronous-spike analysis,
% not only the neurons drawn on one page. Each row is that neuron's
% synchronous-spike pair counts divided by their sum. A neuron that belongs
% to several pairs contributes the union of its synchronous spikes once.
% The upper panel is the stack. The lower panel is every neuron as a thin
% line, the mean across neurons, and the 95% confidence interval of that mean.
%
% One figure per session, six pairs per figure. Each pair is one row:
% neuron 1 waveform, neuron 1 spatial range, neuron 1 autocorrelogram,
% then the same three panels for neuron 2.
% The spatial stack is a separate figure and includes every good unit in
% the five areas, not only the neurons drawn in the pair panels.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

acg_window = [-0.05 0.05];
acg_bin = 0.001;
pairs_per_figure = 6;
template_max_spikes = 100;
template_min_in_window = 20;
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
for session_index = 1:numel(session_ids)
    session_id = session_ids{session_index};
    session_rows = neuron_rows(strcmp({neuron_rows.session_id}, session_id));
    [trains, bout_starts, time_options] = load_session_trains( ...
        npx_raw_data, session_id, session_rows(1));
    neuron_spikes = remember_analysis_neurons(neuron_spikes, session_id, ...
        session_rows, trains, bout_starts, time_options);
    geometry = load_channel_geometry(npx_raw_data, session_id);
    spatial_profiles = append_spatial_profiles(spatial_profiles, ...
        collect_session_profiles(session_id, session_rows, trains, bout_starts, ...
        geometry, time_options, template_max_spikes, template_min_in_window));
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
                    spike_times = synchronous_times(trains, members(side), bout_starts, time_options);
                    [centers, counts] = sync_event_acg(spike_times, acg_window, acg_bin);
                    draw_sync_acg(members(side), centers, counts, spike_times, acg_window);
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
[acg_centers, acg_counts] = neuron_acg_counts(neuron_spikes, acg_window, acg_bin);
draw_normalized_acg_figure(acg_centers, acg_counts, acg_window);
draw_spatial_stack_figure(spatial_profiles, spatial_distance_um);

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

function entries = empty_neuron_spikes()
entries = struct('session_id', {}, 'cluster_id', {}, 'spike_times', {});
end

function entries = remember_analysis_neurons(entries, session_id, session_rows, ...
    trains, bout_starts, time_options)
n_before = numel(entries);
for index = 1:numel(session_rows)
    row = session_rows(index);
    if ~isfield(row, 'cluster_id') || ~isfield(row, 'partner_id')
        continue
    end
    spike_times = synchronous_times(trains, row, bout_starts, time_options);
    entries = remember_neuron_spikes(entries, session_id, row.cluster_id, spike_times);
end
n_neurons = numel(entries) - n_before;
n_analyzed = 0;
if ~isempty(session_rows) && isfield(session_rows, 'cluster_id')
    n_analyzed = numel(unique([session_rows.cluster_id]));
end
fprintf(['%s  autocorrelogram stack  %d analyzed neurons with synchronous ' ...
    'spikes, %d neurons in the saved analysis\n'], session_id, n_neurons, n_analyzed);
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
         'only spikes from synchronous events, for neurons saved by the analysis.']);
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
title(sprintf('Synchronous spikes only, %d analyzed neurons', summary.n));
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

function profiles = collect_session_profiles(session_id, session_rows, trains, ...
    bout_starts, geometry, time_options, max_spikes, min_in_window)
profiles = empty_spatial_profiles();
if isempty(geometry)
    return
end
saved_ids = [];
if ~isempty(session_rows) && isfield(session_rows, 'cluster_id')
    saved_ids = [session_rows.cluster_id];
end
n_paired = 0;
for index = 1:numel(saved_ids)
    cluster_id = saved_ids(index);
    if index > 1 && any(saved_ids(1:index - 1) == cluster_id)
        continue
    end
    members = session_rows([session_rows.cluster_id] == cluster_id);
    [distance_um, range] = profile_from_saved(members, geometry);
    if isempty(distance_um)
        continue
    end
    profiles = add_spatial_profile(profiles, session_id, cluster_id, distance_um, range);
    n_paired = n_paired + 1;
end
[n_other, n_skipped] = add_other_neurons(profiles, session_id, session_rows, trains, ...
    saved_ids, bout_starts, geometry, time_options, max_spikes, min_in_window);
profiles = n_other.profiles;
fprintf(['%s  spatial stack  %d paired neurons  %d other good units  ' ...
    '%d skipped\n'], session_id, n_paired, n_other.n_added, n_skipped);
end

function [added, n_skipped] = add_other_neurons(profiles, session_id, session_rows, trains, ...
    saved_ids, bout_starts, geometry, time_options, max_spikes, min_in_window)
added.profiles = profiles;
added.n_added = 0;
n_skipped = 0;
if isempty(trains)
    return
end
missing = true(1, numel(trains));
for index = 1:numel(trains)
    if any(saved_ids == trains(index).cluster_id)
        missing(index) = false;
    end
end
if ~any(missing)
    return
end
continuous_file = session_continuous_file(session_rows, session_id, time_options, ...
    trains(find(missing, 1)));
if isempty(continuous_file)
    n_skipped = sum(missing);
    warning('Plot_sync_spike_acg_and_waveforms:continuous', ...
        ['%s  %d good units outside the pairs were left out of the spatial ' ...
         'stack because the continuous file was not available.'], ...
        session_id, n_skipped);
    return
end
for index = find(missing)
    try
        [distance_um, range] = profile_from_continuous(trains(index), bout_starts, ...
            geometry, time_options, continuous_file, max_spikes, min_in_window);
    catch read_error
        n_skipped = n_skipped + 1;
        warning('Plot_sync_spike_acg_and_waveforms:template', ...
            '%s cluster %g waveform was not read (%s).', ...
            session_id, trains(index).cluster_id, read_error.message);
        continue
    end
    if isempty(distance_um)
        n_skipped = n_skipped + 1;
        continue
    end
    added.profiles = add_spatial_profile(added.profiles, session_id, ...
        trains(index).cluster_id, distance_um, range);
    added.n_added = added.n_added + 1;
end
end

function [distance_um, range] = profile_from_saved(members, geometry)
distance_um = [];
range = [];
for index = 1:numel(members)
    row = members(index);
    if ~isfield(row, 'sync_mean') || isempty(row.sync_mean) ...
            || ~isfield(row, 'channels_kept') || isempty(row.channels_kept)
        continue
    end
    profile = sync_waveform_spatial_range(row.channels_kept, row.sync_mean, ...
        row.channel, geometry.chan_map, geometry.xcoords, geometry.ycoords);
    distance_um = [distance_um; profile.distance_um(:)]; %#ok<AGROW>
    range = [range; profile.range(:)]; %#ok<AGROW>
end
end

function [distance_um, range] = profile_from_continuous(neuron, bout_starts, ...
    geometry, time_options, continuous_file, max_spikes, min_in_window)
distance_um = [];
range = [];
samples = template_spike_samples(neuron.spike_times, neuron.spike_samples, ...
    bout_starts, time_options.hist_range, max_spikes, min_in_window);
if isempty(samples)
    return
end
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    continuous_file, samples, neuron.channel, time_options.channel_range, ...
    'NChannels', time_options.n_ap_channels, ...
    'SamplesBefore', time_options.samples_before, ...
    'SamplesAfter', time_options.samples_after);
waveforms = center_waveform_snippets(waveforms, time_samples);
profile = sync_waveform_spatial_range(channels_kept, mean_waveform_pages(waveforms), ...
    neuron.channel, geometry.chan_map, geometry.xcoords, geometry.ycoords);
distance_um = profile.distance_um(:);
range = profile.range(:);
end

function mean_wf = mean_waveform_pages(waveforms)
n_channels = size(waveforms, 2);
n_times = size(waveforms, 3);
mean_wf = nan(n_channels, n_times);
for channel_index = 1:n_channels
    traces = reshape(waveforms(:, channel_index, :), size(waveforms, 1), n_times);
    for time_index = 1:n_times
        values = traces(:, time_index);
        values = values(isfinite(values));
        if ~isempty(values)
            mean_wf(channel_index, time_index) = mean(values);
        end
    end
end
end

function text = path_text(value)
text = '';
if exist('isstring', 'builtin') && isstring(value)
    value = char(value);
end
if ischar(value)
    text = strtrim(value);
end
end

function path = session_continuous_file(session_rows, session_id, time_options, example_neuron)
path = '';
if ~isempty(session_rows) && isfield(session_rows, 'continuous_file')
    for index = 1:numel(session_rows)
        candidate = path_text(session_rows(index).continuous_file);
        if ~isempty(candidate) && exist(candidate, 'file') == 2
            path = candidate;
            return
        end
    end
end
try
    candidates = read_npx_continuous_locations(time_options.raw_folders_file, session_id);
    choice.min_scale = time_options.min_spike_scale;
    choice.n_channels = time_options.n_ap_channels;
    choice.n_check = time_options.n_check_spikes;
    choice.samples_before = time_options.samples_before;
    choice.samples_after = time_options.samples_after;
    path = choose_continuous_file(candidates, example_neuron.spike_samples, ...
        example_neuron.channel, time_options.channel_range, choice);
catch
    path = '';
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
        'No spatial profiles were available for the inverted-U stack.');
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
        'No neuron had a usable range on its own channel.');
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
xlabel('Distance from neuron channel (\mum)')
ylabel('Neuron, sorted by peak distance')
title(sprintf('Range / range at the neuron channel, %d neurons', summary.n));
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
xlabel('Distance from neuron channel (\mum)')
ylabel('Range / range at the neuron channel')
title(sprintf('Mean spatial range, 95%% CI, n = %d', summary.n));
box off
end

function draw_sync_acg(row, centers, counts, spike_times, acg_window)
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
