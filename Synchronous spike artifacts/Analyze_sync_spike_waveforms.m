%% Synchronous-spike waveforms inside play-bout windows
% Run Estimate_sync_spike_waveforms first. This loads the per-session mats,
% prints one row per neuron, and draws the distance-percentile histogram.
% The pairs with the highest percentiles are drawn with both neurons
% together: the non-synchronous mean, the synchronous mean, the null
% envelope on the neuron channel, and the time-by-time percentile.
%
% The comparison pool is every spike inside [bout onset - 1.5 s, bout
% onset + 1 s]. nonsync_mean is the mean of the spikes in that pool that
% are not synchronous for this pair. percentile is the fraction of 500
% null distances, from that non-synchronous pool, that fall below the
% synchronous distance. percentile_time is the fraction of those null
% means that lie strictly below the synchronous mean at each channel and
% time. A high percentile_time means the synchronous mean is larger than
% the nulls; a low value means it is smaller (a deeper trough).
%
% neuron_rows and pair_records are left in the workspace.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

data_root = sync_analysis_options().data_root;
saving_folder = fullfile(data_root, 'Analysis results', 'synchronous spike waveforms');
session_files = dir(fullfile(saving_folder, '*_sync_spike_waveforms.mat'));
if isempty(session_files)
    error('Analyze_sync_spike_waveforms:files', ...
        'No per-session waveform files in %s.', saving_folder);
end

neuron_rows = struct([]);
pair_records = struct([]);
for file_index = 1:numel(session_files)
    session_path = fullfile(session_files(file_index).folder, session_files(file_index).name);
    loaded = load(session_path, 'this_waveforms');
    if ~isfield(loaded, 'this_waveforms') || isempty(loaded.this_waveforms.rows)
        continue
    end
    if isempty(neuron_rows)
        neuron_rows = loaded.this_waveforms.rows;
    else
        neuron_rows = [neuron_rows, orderfields(loaded.this_waveforms.rows, neuron_rows)]; %#ok<AGROW>
    end
    if isfield(loaded.this_waveforms, 'pairs') && ~isempty(loaded.this_waveforms.pairs)
        if isempty(pair_records)
            pair_records = loaded.this_waveforms.pairs;
        else
            pair_records = [pair_records, orderfields(loaded.this_waveforms.pairs, pair_records)]; %#ok<AGROW>
        end
    end
end
if isempty(neuron_rows)
    error('Analyze_sync_spike_waveforms:rows', ...
        'The waveform files in %s have no neuron rows.', saving_folder);
end

[~, order] = sort([neuron_rows.percentile], 'descend');
neuron_rows = neuron_rows(order);
n_neurons = numel(neuron_rows);
n_high = sum([neuron_rows.percentile] >= 0.95);
n_high_channel = sum([neuron_rows.percentile_channel] >= 0.95);
fprintf(['%d neurons. Median distance percentile %.2f. ' ...
    '%d at or above 0.95 on the channel neighborhood (about %.1f if synchronous spikes match other spikes). ' ...
    '%d at or above 0.95 on the neuron channel alone.\n'], ...
    n_neurons, median([neuron_rows.percentile]), n_high, 0.05 * n_neurons, n_high_channel);
for row_index = 1:n_neurons
    row = neuron_rows(row_index);
    fprintf(['%s  pair %d  cluster %g  partner %g  %s  ch %g  ' ...
        'spikes %d  interval %d  sync %d  nonsync %d  ' ...
        'percentile %.3f  channel percentile %.3f  ' ...
        'extreme sample %.0f  abs percentile %.3f\n'], ...
        row.session_id, row.pair_index, row.cluster_id, row.partner_id, ...
        row.area, row.channel, row.n_spikes, row.n_interval, row.n_sync, row.n_nonsync, ...
        row.percentile, row.percentile_channel, ...
        row.extreme_time_samples, row.extreme_abs_percentile);
end

figure('Color', 'w', 'Name', 'Synchronous waveform percentiles');
histogram([neuron_rows.percentile], 0:0.05:1, ...
    'FaceColor', [0.25 0.25 0.25], 'EdgeColor', 'none');
hold on
plot([0.95 0.95], ylim, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5);
xlabel('Percentile of synchronous-spike waveform distance')
ylabel('Neurons')
title('Distance from the non-synchronous in-window mean')
box off

if isempty(pair_records)
    pair_records = pairs_from_rows(neuron_rows);
end
scores = pair_scores(pair_records);
[~, pair_order] = sort(scores, 'descend');
pair_records = pair_records(pair_order);
n_plot = min(4, numel(pair_records));
figure('Color', 'w', 'Name', 'Synchronous pairs, both neurons', ...
    'Position', [40 40 1280 240 * n_plot]);
for plot_index = 1:n_plot
    pair = pair_records(plot_index);
    pair.neuron_1.session_id = pair.session_id;
    pair.neuron_2.session_id = pair.session_id;
    subplot(n_plot, 4, (plot_index - 1) * 4 + 1);
    draw_means(pair.neuron_1, plot_index == 1);
    subplot(n_plot, 4, (plot_index - 1) * 4 + 2);
    draw_percentile_time(pair.neuron_1);
    subplot(n_plot, 4, (plot_index - 1) * 4 + 3);
    draw_means(pair.neuron_2, false);
    subplot(n_plot, 4, (plot_index - 1) * 4 + 4);
    draw_percentile_time(pair.neuron_2);
end

function pairs = pairs_from_rows(rows)
pairs = struct('session_id', {}, 'pair_index', {}, 'neuron_1', {}, 'neuron_2', {});
keys = strcat({rows.session_id}, '|', cellstr(num2str([rows.pair_index]')));
unique_keys = unique(keys, 'stable');
for key_index = 1:numel(unique_keys)
    members = rows(strcmp(keys, unique_keys{key_index}));
    record.session_id = members(1).session_id;
    record.pair_index = members(1).pair_index;
    record.neuron_1 = side_from_row(members(1));
    if numel(members) >= 2
        record.neuron_2 = side_from_row(members(2));
    else
        record.neuron_2 = blank_side();
    end
    if isempty(pairs)
        pairs = record;
    else
        pairs(end + 1) = orderfields(record, pairs); %#ok<AGROW>
    end
end
end

function side = side_from_row(row)
side = blank_side();
side.ok = true;
side.cluster_id = row.cluster_id;
side.partner_id = row.partner_id;
side.area = row.area;
side.channel = row.channel;
side.n_spikes = row.n_spikes;
side.n_interval = row.n_interval;
side.n_sync = row.n_sync;
side.n_nonsync = row.n_nonsync;
side.channels_kept = row.channels_kept;
side.time_samples = row.time_samples;
side.nonsync_mean = row.nonsync_mean;
side.sync_mean = row.sync_mean;
side.percentile_time = row.percentile_time;
side.percentile_abs_time = row.percentile_abs_time;
side.null_p025 = row.null_p025;
side.null_p975 = row.null_p975;
side.percentile = row.percentile;
side.percentile_channel = row.percentile_channel;
side.extreme_time_samples = row.extreme_time_samples;
side.extreme_abs_percentile = row.extreme_abs_percentile;
side.session_id = row.session_id;
end

function side = blank_side()
side.ok = false;
side.reason = 'missing neuron';
side.cluster_id = NaN;
side.area = '';
side.channel = NaN;
side.channels_kept = [];
side.time_samples = [];
side.nonsync_mean = [];
side.sync_mean = [];
side.percentile_time = [];
side.null_p025 = [];
side.null_p975 = [];
side.percentile = NaN;
side.percentile_channel = NaN;
side.extreme_time_samples = NaN;
side.extreme_abs_percentile = NaN;
side.session_id = '';
end

function scores = pair_scores(pairs)
scores = -inf(1, numel(pairs));
for pair_index = 1:numel(pairs)
    scores(pair_index) = max([ ...
        side_percentile(pairs(pair_index).neuron_1), ...
        side_percentile(pairs(pair_index).neuron_2)]);
end
end

function value = side_percentile(side)
value = -inf;
if isfield(side, 'ok') && side.ok && isfield(side, 'percentile') && isfinite(side.percentile)
    value = side.percentile;
end
end

function draw_means(side, show_legend)
if ~isfield(side, 'ok') || ~side.ok || isempty(side.sync_mean) || isempty(side.nonsync_mean)
    title(missing_title(side), 'Interpreter', 'none');
    axis off
    return
end
nonsync = side.nonsync_mean;
sync_mean = side.sync_mean;
n_channels = size(nonsync, 1);
time_ms = side.time_samples(:)' / 30000 * 1000;
scale = max(abs([nonsync(:); sync_mean(:)]));
if isempty(scale) || ~isfinite(scale) || scale == 0
    scale = 1;
end
gap = 2.6 * scale;
hold on
neuron_nonsync = [];
neuron_sync = [];
for channel_index = 1:n_channels
    y_offset = (channel_index - 1) * gap;
    line_width = 1.2;
    is_neuron = side.channels_kept(channel_index) == side.channel;
    if is_neuron
        line_width = 2.2;
        if ~isempty(side.null_p025) && ~isempty(side.null_p975)
            low = side.null_p025(channel_index, :) + y_offset;
            high = side.null_p975(channel_index, :) + y_offset;
            fill([time_ms, fliplr(time_ms)], [low, fliplr(high)], ...
                [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        end
    end
    nonsync_line = plot(time_ms, nonsync(channel_index, :) + y_offset, ...
        'Color', [0 0 0], 'LineWidth', line_width);
    sync_line = plot(time_ms, sync_mean(channel_index, :) + y_offset, ...
        'Color', [0.85 0.1 0.1], 'LineWidth', line_width);
    if is_neuron
        neuron_nonsync = nonsync_line;
        neuron_sync = sync_line;
    end
end
padding = 1.25 * scale;
ylim([-padding, (n_channels - 1) * gap + padding]);
xlim([time_ms(1) time_ms(end)]);
y_ticks = (0:n_channels - 1) * gap;
tick_labels = cell(1, n_channels);
for channel_index = 1:n_channels
    tick_labels{channel_index} = sprintf('%g', side.channels_kept(channel_index));
end
set(gca, 'YTick', y_ticks, 'YTickLabel', tick_labels);
xlabel('Time (ms)')
ylabel('Channel')
title(sprintf('%s   cluster %g   %s   ch %g   p %.2f', ...
    side_session(side), side.cluster_id, side.area, side.channel, side.percentile), ...
    'Interpreter', 'none');
if show_legend && ~isempty(neuron_nonsync)
    legend([neuron_nonsync, neuron_sync], ...
        {'Non-synchronous mean', 'Synchronous mean'}, 'Location', 'northeast');
end
box off
end

function draw_percentile_time(side)
if ~isfield(side, 'ok') || ~side.ok || isempty(side.percentile_time)
    title(missing_title(side), 'Interpreter', 'none');
    axis off
    return
end
time_ms = side.time_samples(:)' / 30000 * 1000;
n_channels = size(side.percentile_time, 1);
hold on
plot([time_ms(1) time_ms(end)], [0.5 0.5], 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
for channel_index = 1:n_channels
    is_neuron = side.channels_kept(channel_index) == side.channel;
    if is_neuron
        plot(time_ms, side.percentile_time(channel_index, :), ...
            'Color', [0.85 0.1 0.1], 'LineWidth', 2.2);
    else
        plot(time_ms, side.percentile_time(channel_index, :), ...
            'Color', [0.55 0.55 0.55], 'LineWidth', 1);
    end
end
ylim([-0.02 1.02]);
xlim([time_ms(1) time_ms(end)]);
xlabel('Time (ms)')
ylabel('Percentile')
title(sprintf('cluster %g   time percentile   extreme %.0f samples', ...
    side.cluster_id, side.extreme_time_samples), 'Interpreter', 'none');
box off
end

function label = missing_title(side)
if isfield(side, 'reason') && ~isempty(side.reason)
    label = side.reason;
else
    label = 'No waveform';
end
end

function label = side_session(side)
if isfield(side, 'session_id') && ~isempty(side.session_id)
    label = side.session_id;
else
    label = '';
end
end
