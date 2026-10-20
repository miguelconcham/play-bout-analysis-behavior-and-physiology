function result = plot_cluster_waveform_and_acg( ...
    continuous_file, spike_samples, neuron_channel, varargin)
%PLOT_CLUSTER_WAVEFORM_AND_ACG  Multi-channel waveform and autocorrelogram.
%
%   result = plot_cluster_waveform_and_acg(continuous_file, spike_samples, ...
%       neuron_channel)
%
% Loads snippets with load_waveforms_from_continuous and plots every kept
% channel, not the neuron channel alone. ChannelRange is the number of
% channels on each side (default 5, so up to 11 channels, fewer at the ends
% of the probe). The neuron channel is drawn in red.
%
% The autocorrelogram uses the same lags as autocorrelogram.m.
%
% Name-value pairs: ChannelRange, NWaveforms, NChannels, SampleRate,
% AcgWindow, AcgBin, SessionName, ClusterID, Visible, Plot.

options = parse_plot_options(varargin);
spike_samples = double(spike_samples(:));
if isempty(spike_samples)
    error('plot_cluster_waveform_and_acg:spikes', 'No spike samples to plot.');
end

n_draw = min(options.n_waveforms, numel(spike_samples));
draw_index = randperm(numel(spike_samples), n_draw);
drawn_samples = spike_samples(draw_index);
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    continuous_file, drawn_samples, neuron_channel, options.channel_range, ...
    'NChannels', options.n_channels);

lags = autocorrelogram_lags(spike_samples / options.sample_rate, options.acg_window);

result.waveforms = waveforms;
result.channels_kept = channels_kept;
result.time_samples = time_samples;
result.time_ms = time_samples / options.sample_rate * 1000;
result.spike_samples = drawn_samples;
result.lags = lags;
result.neuron_channel = neuron_channel;
result.continuous_file = continuous_file;
result.figure = [];
if options.plot
    result.figure = draw_cluster_figure(result, options);
end
end

function options = parse_plot_options(name_value)
options.channel_range = 5;
options.n_waveforms = 40;
options.n_channels = 384;
options.sample_rate = 30000;
options.acg_window = [-0.05 0.05];
options.acg_bin = 0.001;
options.session_name = '';
options.cluster_id = [];
options.visible = 'on';
options.plot = true;

if mod(numel(name_value), 2) ~= 0
    error('plot_cluster_waveform_and_acg:nameValue', ...
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
        case 'nchannels'
            options.n_channels = value;
        case 'samplerate'
            options.sample_rate = value;
        case 'acgwindow'
            options.acg_window = value;
        case 'acgbin'
            options.acg_bin = value;
        case 'sessionname'
            options.session_name = char(value);
        case 'clusterid'
            options.cluster_id = value;
        case 'visible'
            options.visible = char(value);
        case 'plot'
            options.plot = logical(value);
        otherwise
            error('plot_cluster_waveform_and_acg:nameValue', ...
                'Unknown parameter "%s".', name);
    end
end
end

function fig = draw_cluster_figure(result, options)
fig = figure('Color', 'w', 'Visible', options.visible);
subplot(1, 2, 1);
draw_waveforms(result);
subplot(1, 2, 2);
draw_autocorrelogram(result.lags, options);
set(fig, 'Name', figure_name(result, options));
end

function draw_waveforms(result)
waveforms = result.waveforms;
n_channels = size(waveforms, 2);
n_times = size(waveforms, 3);
mean_wf = nan(n_channels, n_times);
for channel_index = 1:n_channels
    mean_wf(channel_index, :) = mean_omitnan(reshape( ...
        waveforms(:, channel_index, :), size(waveforms, 1), n_times), 1);
end
peak = max(abs(mean_wf(:)));
if isempty(peak) || ~isfinite(peak) || peak == 0
    peak = 1;
end
gap = 2.2 * peak;
time_ms = result.time_ms(:)';
hold on
y_ticks = (0:n_channels - 1) * gap;
for channel_index = 1:n_channels
    offset = y_ticks(channel_index);
    traces = reshape(waveforms(:, channel_index, :), size(waveforms, 1), n_times);
    plot(time_ms, traces' + offset, 'Color', [0.75 0.75 0.75]);
end
for channel_index = 1:n_channels
    offset = y_ticks(channel_index);
    mean_color = [0 0 0];
    mean_width = 1.5;
    if result.channels_kept(channel_index) == result.neuron_channel
        mean_color = [0.85 0.1 0.1];
        mean_width = 2.5;
    end
    plot(time_ms, mean_wf(channel_index, :) + offset, ...
        'Color', mean_color, 'LineWidth', mean_width);
end
tick_labels = arrayfun(@(channel) sprintf('%g', channel), ...
    result.channels_kept, 'UniformOutput', false);
set(gca, 'YTick', y_ticks, 'YTickLabel', tick_labels);
xlabel('Time (ms)')
ylabel('Channel')
title('Waveform')
box off
end

function draw_autocorrelogram(lags, options)
edges = options.acg_window(1):options.acg_bin:options.acg_window(2);
if numel(edges) < 2
    edges = [options.acg_window(1) options.acg_window(2)];
end
[counts, edges] = histcounts(lags, edges);
centers = (edges(1:end - 1) + edges(2:end)) / 2;
bar(centers, counts, 1, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none');
xlim(options.acg_window);
xlabel('Lag (s)')
ylabel('Count')
title('Autocorrelogram')
box off
end

function name = figure_name(result, options)
[parent_dir, file_name, file_ext] = fileparts(result.continuous_file);
[~, parent_name] = fileparts(parent_dir);
file_label = [parent_name filesep file_name file_ext];
if isempty(options.session_name)
    session_label = 'cluster';
else
    session_label = options.session_name;
end
if isempty(options.cluster_id)
    cluster_label = 'cluster';
else
    cluster_label = sprintf('cluster %g', options.cluster_id);
end
name = sprintf('%s  %s  ch %g  %s', session_label, cluster_label, ...
    result.neuron_channel, file_label);
end

function lags = autocorrelogram_lags(spike_seconds, window)
% Same pairs as autocorrelogram.m: lag of every other spike inside the window.
spike_seconds = sort(double(spike_seconds(:)));
n_spikes = numel(spike_seconds);
lags = zeros(0, 1);
stored = 0;
cursor = 1;
for spike_index = 1:n_spikes
    spike_time = spike_seconds(spike_index);
    while cursor < n_spikes && spike_seconds(cursor) < spike_time + window(1)
        cursor = cursor + 1;
    end
    other_index = cursor;
    while other_index <= n_spikes && spike_seconds(other_index) <= spike_time + window(2)
        lag = spike_seconds(other_index) - spike_time;
        if other_index ~= spike_index && lag ~= 0
            stored = stored + 1;
            if stored > numel(lags)
                lags(max(stored * 2, 16), 1) = 0; %#ok<AGROW>
            end
            lags(stored) = lag;
        end
        other_index = other_index + 1;
    end
end
lags = lags(1:stored);
end

function values = mean_omitnan(values, dimension)
missing = isnan(values);
values(missing) = 0;
counts = sum(~missing, dimension);
values = sum(values, dimension) ./ counts;
values(counts == 0) = NaN;
end
