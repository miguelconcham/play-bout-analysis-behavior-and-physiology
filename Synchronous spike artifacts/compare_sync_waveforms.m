function comparison = compare_sync_waveforms(waveforms, time_samples, sync_index, ...
    neuron_channel, channels_kept, n_repeats)
%COMPARE_SYNC_WAVEFORMS  Synchronous mean versus non-synchronous spikes.
%
% waveforms is n_spikes x n_channels x n_times for the comparison pool
% (spikes inside the play-bout windows), in the same order as sync_index.
% Every snippet is baseline-subtracted (samples at or before -10).
% nonsync_mean is the mean of the finite spikes that are not synchronous
% for this pair. sync_mean uses every finite synchronous spike.
%
% Each of n_repeats null means draws that same count from the
% non-synchronous spikes, without replacement inside the draw. The same
% non-synchronous spike may appear in more than one draw. Spikes that are
% synchronous for this pair are excluded. The 500 waveform tensors are not
% returned.
%
% Distance is the mean of abs(candidate mean - nonsync_mean) over channels
% and time, and the same quantity on the neuron channel alone. percentile
% is the fraction of the null distances strictly smaller than the
% synchronous distance.
%
% percentile_time is channels x time: the fraction of null means strictly
% below sync_mean at that sample. A negative trough that is below every
% null has a low percentile. percentile_abs_time is the fraction of null
% means whose absolute departure from nonsync_mean is strictly smaller, so
% a high value means the synchronous mean is farther from the
% non-synchronous mean regardless of sign. null_p025 and null_p975 are the
% pointwise 2.5% and 97.5% waveforms of the null means.

if nargin < 6 || isempty(n_repeats)
    n_repeats = 500;
end
n_repeats = round(double(n_repeats));
channels_kept = double(channels_kept(:)');
time_samples = double(time_samples(:)');

comparison = empty_comparison(channels_kept, time_samples, n_repeats);
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
n_channels = size(centered, 2);
n_times = size(centered, 3);
comparison.n_pool = n_spikes;
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
comparison.n_nonsync = numel(other_keep);
if isempty(sync_keep)
    comparison.reason = 'no finite synchronous waveforms';
    return
end
if numel(other_keep) < numel(sync_keep)
    comparison.reason = 'fewer other spikes than synchronous spikes';
    return
end

comparison.sync_mean = mean_snippets(centered(sync_keep, :, :));
comparison.nonsync_mean = mean_snippets(centered(other_keep, :, :));
comparison.distance_sync = mean_abs_distance( ...
    comparison.sync_mean, comparison.nonsync_mean);
comparison.distance_sync_channel = mean_abs_distance( ...
    comparison.sync_mean(neuron_row, :), comparison.nonsync_mean(neuron_row, :));

n_sync = numel(sync_keep);
null_means = nan(n_repeats, n_channels, n_times);
for repeat_index = 1:n_repeats
    draw = other_keep(randperm(numel(other_keep), n_sync));
    draw_mean = mean_snippets(centered(draw, :, :));
    null_means(repeat_index, :, :) = reshape(draw_mean, 1, n_channels, n_times);
    comparison.null_distances(repeat_index) = mean_abs_distance( ...
        draw_mean, comparison.nonsync_mean);
    comparison.null_distances_channel(repeat_index) = mean_abs_distance( ...
        draw_mean(neuron_row, :), comparison.nonsync_mean(neuron_row, :));
end
if ~isfinite(comparison.distance_sync) || ~isfinite(comparison.distance_sync_channel)
    comparison.reason = 'synchronous waveform distance is not finite';
    return
end

sync_page = reshape(comparison.sync_mean, 1, n_channels, n_times);
nonsync_page = reshape(comparison.nonsync_mean, 1, n_channels, n_times);
comparison.percentile_time = squeeze(mean(null_means < sync_page, 1));
comparison.percentile_abs_time = squeeze(mean( ...
    abs(null_means - nonsync_page) < abs(sync_page - nonsync_page), 1));
comparison.null_p025 = pointwise_quantile(null_means, 0.025);
comparison.null_p975 = pointwise_quantile(null_means, 0.975);
if n_channels == 1
    comparison.percentile_time = reshape(comparison.percentile_time, 1, n_times);
    comparison.percentile_abs_time = reshape(comparison.percentile_abs_time, 1, n_times);
    comparison.null_p025 = reshape(comparison.null_p025, 1, n_times);
    comparison.null_p975 = reshape(comparison.null_p975, 1, n_times);
end
clear null_means

comparison.percentile = mean(comparison.null_distances < comparison.distance_sync);
comparison.percentile_channel = mean( ...
    comparison.null_distances_channel < comparison.distance_sync_channel);
channel_abs = comparison.percentile_abs_time(neuron_row, :);
[comparison.extreme_abs_percentile, extreme_at] = max(channel_abs);
comparison.extreme_time_samples = time_samples(extreme_at);
comparison.extreme_percentile = comparison.percentile_time(neuron_row, extreme_at);
finite_abs = channel_abs(isfinite(channel_abs));
if isempty(finite_abs)
    comparison.neuron_channel_tail_fraction = NaN;
else
    comparison.neuron_channel_tail_fraction = mean(finite_abs >= 0.975);
end
comparison.ok = true;
end

function comparison = empty_comparison(channels_kept, time_samples, n_repeats)
comparison.ok = false;
comparison.reason = '';
comparison.n_pool = 0;
comparison.n_sync = 0;
comparison.n_nonsync = 0;
comparison.channels_kept = channels_kept;
comparison.time_samples = time_samples;
comparison.nonsync_mean = [];
comparison.sync_mean = [];
comparison.percentile_time = [];
comparison.percentile_abs_time = [];
comparison.null_p025 = [];
comparison.null_p975 = [];
comparison.distance_sync = NaN;
comparison.null_distances = nan(max(n_repeats, 0), 1);
comparison.percentile = NaN;
comparison.distance_sync_channel = NaN;
comparison.null_distances_channel = nan(max(n_repeats, 0), 1);
comparison.percentile_channel = NaN;
comparison.extreme_time_samples = NaN;
comparison.extreme_percentile = NaN;
comparison.extreme_abs_percentile = NaN;
comparison.neuron_channel_tail_fraction = NaN;
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

function quantile_waveform = pointwise_quantile(null_means, probability)
sorted = sort(null_means, 1);
n_repeats = size(sorted, 1);
rank = max(1, min(n_repeats, ceil(probability * n_repeats)));
quantile_waveform = squeeze(sorted(rank, :, :));
end
