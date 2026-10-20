%% Raw-LFP power at synchronous play-bout spikes
% Run Estimate_sync_spike_lfp_power first. This loads the per-session mats,
% prints one row per neuron, and draws the percentile histogram plus the
% null distribution for the neurons with the highest percentiles.
%
% Power is the square of the Hilbert envelope of unfiltered LFP_PAG, one
% interpolated sample at each spike time, on the neuron's own channel.
% percentile is the fraction of 500 null means below the synchronous mean.
% A high percentile means those spikes sit on larger raw-LFP power than
% the same number of the neuron's other spikes.
%
% neuron_rows is left in the workspace.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

data_root = sync_analysis_options().data_root;
saving_folder = fullfile(data_root, 'Analysis results', 'synchronous spike waveforms');
session_files = dir(fullfile(saving_folder, '*_sync_spike_lfp_power.mat'));
if isempty(session_files)
    error('Analyze_sync_spike_lfp_power:files', ...
        'No per-session LFP-power files in %s.', saving_folder);
end

neuron_rows = struct([]);
for file_index = 1:numel(session_files)
    session_path = fullfile(session_files(file_index).folder, session_files(file_index).name);
    loaded = load(session_path, 'this_lfp_power');
    if ~isfield(loaded, 'this_lfp_power') || isempty(loaded.this_lfp_power.rows)
        continue
    end
    if isempty(neuron_rows)
        neuron_rows = loaded.this_lfp_power.rows;
    else
        neuron_rows = [neuron_rows, orderfields(loaded.this_lfp_power.rows, neuron_rows)]; %#ok<AGROW>
    end
end
if isempty(neuron_rows)
    error('Analyze_sync_spike_lfp_power:rows', ...
        'The LFP-power files in %s have no neuron rows.', saving_folder);
end

[~, order] = sort([neuron_rows.percentile], 'descend');
neuron_rows = neuron_rows(order);
n_neurons = numel(neuron_rows);
n_high = sum([neuron_rows.percentile] >= 0.95);
fprintf(['%d neurons. Median power percentile %.2f. ' ...
    '%d at or above 0.95 (about %.1f if synchronous spikes match other spikes).\n'], ...
    n_neurons, median([neuron_rows.percentile]), n_high, 0.05 * n_neurons);
for row_index = 1:n_neurons
    row = neuron_rows(row_index);
    fprintf(['%s  cluster %g  partner %g  %s  ch %g  spikes %d  sync %d  ' ...
        'sync power %.3g  all-spike power %.3g  percentile %.3f\n'], ...
        row.session_id, row.cluster_id, row.partner_id, row.area, row.channel, ...
        row.n_spikes, row.n_sync, row.mean_power_sync, row.mean_power_all, ...
        row.percentile);
end

figure('Color', 'w', 'Name', 'Synchronous LFP-power percentiles');
histogram([neuron_rows.percentile], 0:0.05:1, ...
    'FaceColor', [0.25 0.25 0.25], 'EdgeColor', 'none');
hold on
plot([0.95 0.95], ylim, 'Color', [0.85 0.1 0.1], 'LineWidth', 1.5);
xlabel('Percentile of mean raw-LFP power at synchronous spikes')
ylabel('Neurons')
title('Unfiltered LFP envelope power versus 500 matched draws')
box off

n_plot = min(8, n_neurons);
n_columns = min(4, n_plot);
n_figure_rows = ceil(n_plot / n_columns);
figure('Color', 'w', 'Name', 'Highest synchronous LFP-power percentiles', ...
    'Position', [80 80 320 * n_columns 260 * n_figure_rows]);
for plot_index = 1:n_plot
    subplot(n_figure_rows, n_columns, plot_index);
    row = neuron_rows(plot_index);
    histogram(row.null_mean_power, 20, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
    hold on
    y_limits = ylim;
    plot([row.mean_power_sync row.mean_power_sync], y_limits, ...
        'Color', [0.85 0.1 0.1], 'LineWidth', 2);
    xlabel('Mean power (count^2)')
    ylabel('Null draws')
    title(sprintf('%s  ch %g  p %.2f', row.session_id, row.channel, row.percentile), ...
        'Interpreter', 'none');
    box off
end
