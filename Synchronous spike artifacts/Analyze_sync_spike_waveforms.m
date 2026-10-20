%% Synchronous-spike waveform distances
% Run Estimate_sync_spike_waveforms first. This loads the per-session mats,
% prints one row per neuron, and draws the percentile histogram plus the
% template and synchronous mean for the neurons with the highest percentiles.
%
% neuron_rows is one struct per neuron in a pair:
%   percentile            fraction of 500 null distances below the
%                         synchronous distance, all kept channels
%   percentile_channel    the same comparison on the neuron channel alone
% A high percentile means the synchronous mean is more unlike the
% all-spike mean than a typical draw of the same number of spikes.

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
    fprintf(['%s  cluster %g  partner %g  %s  ch %g  spikes %d  sync %d  ' ...
        'percentile %.3f  channel percentile %.3f\n'], ...
        row.session_id, row.cluster_id, row.partner_id, row.area, row.channel, ...
        row.n_spikes, row.n_sync, row.percentile, row.percentile_channel);
end

figure('Color', 'w', 'Name', 'Synchronous waveform percentiles');
histogram([neuron_rows.percentile], 0:0.05:1, ...
    'FaceColor', [0.25 0.25 0.25], 'EdgeColor', 'none');
hold on
plot([0.95 0.95], ylim, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5);
xlabel('Percentile of synchronous-spike waveform distance')
ylabel('Neurons')
title('Distance from the all-spike mean, all kept channels')
box off

n_plot = min(8, n_neurons);
figure('Color', 'w', 'Name', 'Highest synchronous-waveform percentiles', ...
    'Position', [80 40 760 180 * n_plot]);
for plot_index = 1:n_plot
    subplot(n_plot, 1, plot_index);
    draw_template_and_sync(neuron_rows(plot_index), plot_index == 1);
end

function draw_template_and_sync(row, show_legend)
template = row.template_mean;
sync_mean = row.sync_mean;
n_channels = size(template, 1);
time_ms = row.time_samples(:)' / 30000 * 1000;
scale = max(abs([template(:); sync_mean(:)]));
if isempty(scale) || ~isfinite(scale) || scale == 0
    scale = 1;
end
gap = 2.6 * scale;
hold on
neuron_template = [];
neuron_sync = [];
for channel_index = 1:n_channels
    y_offset = (channel_index - 1) * gap;
    line_width = 1.2;
    if row.channels_kept(channel_index) == row.channel
        line_width = 2.2;
    end
    template_line = plot(time_ms, template(channel_index, :) + y_offset, ...
        'Color', [0 0 0], 'LineWidth', line_width);
    sync_line = plot(time_ms, sync_mean(channel_index, :) + y_offset, ...
        'Color', [0.85 0.1 0.1], 'LineWidth', line_width);
    if row.channels_kept(channel_index) == row.channel
        neuron_template = template_line;
        neuron_sync = sync_line;
    end
end
padding = 1.25 * scale;
ylim([-padding, (n_channels - 1) * gap + padding]);
xlim([time_ms(1) time_ms(end)]);
y_ticks = (0:n_channels - 1) * gap;
tick_labels = cell(1, n_channels);
for channel_index = 1:n_channels
    tick_labels{channel_index} = sprintf('%g', row.channels_kept(channel_index));
end
set(gca, 'YTick', y_ticks, 'YTickLabel', tick_labels);
xlabel('Time (ms)')
ylabel('Channel')
title(sprintf('%s   cluster %g   %s   ch %g   percentile %.2f   channel %.2f', ...
    row.session_id, row.cluster_id, row.area, row.channel, ...
    row.percentile, row.percentile_channel));
if show_legend && ~isempty(neuron_template)
    legend([neuron_template, neuron_sync], ...
        {'All-spike mean', 'Synchronous mean'}, 'Location', 'northeast');
end
box off
end
