function comparison = compare_sync_waveforms(waveforms, time_samples, sync_index, ...
    neuron_channel, channels_kept, n_repeats)
%COMPARE_SYNC_WAVEFORMS  Distance of the synchronous mean from the all-spike mean.
%
% waveforms is n_spikes x n_channels x n_times, in the same spike order as
% sync_index. Every snippet is baseline-subtracted (samples at or before
% -10). The template is the mean of every finite spike. The synchronous
% mean uses every finite synchronous spike. Each null mean uses that same
% count of the neuron's other finite spikes, without replacement inside
% the draw.
%
% Distance is the mean of abs(candidate mean - template) over channels and
% time, and the same quantity on the neuron channel alone. percentile is
% the fraction of the null distances that are smaller than the synchronous
% distance. A high percentile means the synchronous waveform is more
% unlike the all-spike mean than a typical draw of the same size.

if nargin < 6 || isempty(n_repeats)
    n_repeats = 500;
end
n_repeats = round(double(n_repeats));
channels_kept = double(channels_kept(:)');
time_samples = double(time_samples(:)');

comparison.ok = false;
comparison.reason = '';
comparison.n_sync = 0;
comparison.channels_kept = channels_kept;
comparison.time_samples = time_samples;
comparison.template_mean = [];
comparison.sync_mean = [];
comparison.distance_sync = NaN;
comparison.null_distances = nan(max(n_repeats, 0), 1);
comparison.percentile = NaN;
comparison.distance_sync_channel = NaN;
comparison.null_distances_channel = nan(max(n_repeats, 0), 1);
comparison.percentile_channel = NaN;

if ~isscalar(n_repeats) || ~isfinite(n_repeats) || n_repeats < 1
    comparison.reason = 'n_repeats must be a positive integer';
    return
end
neuron_row = find(channels_kept == neuron_channel, 1);
if isempty(neuron_row)
    comparison.reason = 'neuron channel is outside the loaded channel window';
    return
end
if isempty(waveforms)
    comparison.reason = 'no waveforms';
    return
end

centered = center_waveform_snippets(waveforms, time_samples);
n_spikes = size(centered, 1);
usable = neuron_channel_is_finite(centered, neuron_row);
is_sync = false(n_spikes, 1);
sync_index = unique(round(double(sync_index(:)')));
sync_index = sync_index(sync_index >= 1 & sync_index <= n_spikes);
if ~isempty(sync_index)
    is_sync(sync_index) = true;
end

sync_keep = find(is_sync & usable);
other_keep = find(~is_sync & usable);
comparison.n_sync = numel(sync_keep);
comparison.template_mean = mean_snippets(centered(usable, :, :));
if isempty(sync_keep)
    comparison.reason = 'no finite synchronous waveforms';
    return
end
if numel(other_keep) < numel(sync_keep)
    comparison.reason = 'fewer other spikes than synchronous spikes';
    return
end

comparison.sync_mean = mean_snippets(centered(sync_keep, :, :));
comparison.distance_sync = mean_abs_distance( ...
    comparison.sync_mean, comparison.template_mean);
comparison.distance_sync_channel = mean_abs_distance( ...
    comparison.sync_mean(neuron_row, :), comparison.template_mean(neuron_row, :));

n_sync = numel(sync_keep);
for repeat_index = 1:n_repeats
    draw = other_keep(randperm(numel(other_keep), n_sync));
    draw_mean = mean_snippets(centered(draw, :, :));
    comparison.null_distances(repeat_index) = mean_abs_distance( ...
        draw_mean, comparison.template_mean);
    comparison.null_distances_channel(repeat_index) = mean_abs_distance( ...
        draw_mean(neuron_row, :), comparison.template_mean(neuron_row, :));
end
if ~isfinite(comparison.distance_sync) || ~isfinite(comparison.distance_sync_channel)
    comparison.reason = 'synchronous waveform distance is not finite';
    return
end
comparison.percentile = mean(comparison.null_distances < comparison.distance_sync);
comparison.percentile_channel = mean( ...
    comparison.null_distances_channel < comparison.distance_sync_channel);
comparison.ok = true;
end

function usable = neuron_channel_is_finite(snippets, neuron_row)
n_spikes = size(snippets, 1);
n_times = size(snippets, 3);
traces = reshape(snippets(:, neuron_row, :), n_spikes, n_times);
usable = any(isfinite(traces), 2);
end

function mean_waveform = mean_snippets(snippets)
n_channels = size(snippets, 2);
n_times = size(snippets, 3);
if isempty(snippets)
    mean_waveform = nan(n_channels, n_times);
    return
end
finite_samples = isfinite(snippets);
snippets(~finite_samples) = 0;
counts = sum(finite_samples, 1);
sums = sum(snippets, 1);
mean_waveform = reshape(sums ./ max(counts, 1), n_channels, n_times);
missing = reshape(counts, n_channels, n_times) == 0;
mean_waveform(missing) = NaN;
end

function distance = mean_abs_distance(candidate, template)
difference = abs(candidate - template);
difference = difference(isfinite(difference));
if isempty(difference)
    distance = NaN;
else
    distance = mean(difference);
end
end
