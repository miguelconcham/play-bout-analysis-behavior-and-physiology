function [pairs, kept_neurons] = select_synchronous_pairs(neurons, bout_starts, options)
%SELECT_SYNCHRONOUS_PAIRS  Rank and cap play-bout coincident pairs.
%
% neurons have cluster_id, channel, area, spike_times (seconds), and
% spike_samples aligned with those times. Both neurons need at least
% min_spikes spikes in the session and at least min_sync coincident spikes
% inside the bout windows. The null pool is the other spikes inside those
% same windows, [bout_start + hist_range(1), bout_start + hist_range(2)]
% for every play bout. A pair is kept only when each neuron has at least
% as many non-synchronous spikes in the windows as synchronous spikes, so
% a null of the same size can be drawn. Every coincident spike is kept.
% Pairs are taken from the larger session spike counts first, with at most
% max_pairs_per_neuron, max_pairs_per_session, and max_pairs pairs.
%
% sync_index_1 and sync_index_2 index the sorted spike list stored on the
% returned neurons. interval_index on each neuron is the same list's spikes
% that fall inside any bout window.

if nargin < 3 || isempty(options)
    options = struct();
end
options = sync_analysis_options(options);
pairs = empty_pairs();
kept_neurons = empty_neurons();
if numel(neurons) < 2 || isempty(bout_starts) || options.max_pairs < 1
    return
end

neurons = sort_neuron_spikes(neurons);
neurons = attach_interval_index(neurons, bout_starts, options.hist_range);
eligible = empty_pairs();
for first = 1:numel(neurons) - 1
    for second = first + 1:numel(neurons)
        sync_index_1 = find(sync_spike_mask(neurons(first).spike_times, ...
            neurons(second).spike_times, bout_starts, ...
            options.hist_range, options.time_precision));
        sync_index_2 = find(sync_spike_mask(neurons(second).spike_times, ...
            neurons(first).spike_times, bout_starts, ...
            options.hist_range, options.time_precision));
        n1 = numel(neurons(first).spike_times);
        n2 = numel(neurons(second).spike_times);
        n_sync_1 = numel(sync_index_1);
        n_sync_2 = numel(sync_index_2);
        n_interval_1 = numel(neurons(first).interval_index);
        n_interval_2 = numel(neurons(second).interval_index);
        n_nonsync_1 = n_interval_1 - numel(intersect(sync_index_1, neurons(first).interval_index));
        n_nonsync_2 = n_interval_2 - numel(intersect(sync_index_2, neurons(second).interval_index));
        if n1 < options.min_spikes || n2 < options.min_spikes
            continue
        end
        if n_sync_1 < options.min_sync || n_sync_2 < options.min_sync
            continue
        end
        if n_nonsync_1 < n_sync_1 || n_nonsync_2 < n_sync_2
            continue
        end
        eligible(end + 1).neuron_index_1 = first; %#ok<AGROW>
        eligible(end).neuron_index_2 = second;
        eligible(end).cluster_id_1 = neurons(first).cluster_id;
        eligible(end).cluster_id_2 = neurons(second).cluster_id;
        eligible(end).channel_1 = neurons(first).channel;
        eligible(end).channel_2 = neurons(second).channel;
        eligible(end).area_1 = neurons(first).area;
        eligible(end).area_2 = neurons(second).area;
        eligible(end).n_spikes_1 = n1;
        eligible(end).n_spikes_2 = n2;
        eligible(end).n_interval_1 = n_interval_1;
        eligible(end).n_interval_2 = n_interval_2;
        eligible(end).n_sync_1 = n_sync_1;
        eligible(end).n_sync_2 = n_sync_2;
        eligible(end).sync_index_1 = sync_index_1(:);
        eligible(end).sync_index_2 = sync_index_2(:);
        eligible(end).rank_spikes = min(n1, n2);
    end
end
if isempty(eligible)
    return
end

[~, order] = sort([eligible.rank_spikes], 'descend');
eligible = eligible(order);
use_count = zeros(numel(neurons), 1);
taken = 0;
limit = min(options.max_pairs, options.max_pairs_per_session);
for pair_index = 1:numel(eligible)
    if taken >= limit
        break
    end
    pair = eligible(pair_index);
    if use_count(pair.neuron_index_1) >= options.max_pairs_per_neuron ...
            || use_count(pair.neuron_index_2) >= options.max_pairs_per_neuron
        continue
    end
    pairs(end + 1) = pair; %#ok<AGROW>
    use_count(pair.neuron_index_1) = use_count(pair.neuron_index_1) + 1;
    use_count(pair.neuron_index_2) = use_count(pair.neuron_index_2) + 1;
    taken = taken + 1;
end

kept_index = find(use_count > 0);
for index = 1:numel(kept_index)
    kept_neurons(end + 1) = neurons(kept_index(index)); %#ok<AGROW>
end
end

function neurons = attach_interval_index(neurons, bout_starts, hist_range)
for index = 1:numel(neurons)
    neurons(index).interval_index = bout_window_index( ...
        neurons(index).spike_times, bout_starts, hist_range);
end
end

function index = bout_window_index(spike_times, bout_starts, hist_range)
spike_times = double(spike_times(:));
bout_starts = double(bout_starts(:));
in_window = false(numel(spike_times), 1);
for bout_index = 1:numel(bout_starts)
    t0 = bout_starts(bout_index);
    in_window = in_window | (spike_times >= t0 + hist_range(1) ...
        & spike_times <= t0 + hist_range(2));
end
index = find(in_window);
end

function neurons = sort_neuron_spikes(neurons)
for index = 1:numel(neurons)
    times = double(neurons(index).spike_times(:));
    samples = double(neurons(index).spike_samples(:));
    if numel(samples) ~= numel(times)
        error('select_synchronous_pairs:samples', ...
            'spike_samples must have one entry for every spike time.');
    end
    keep = isfinite(times);
    times = times(keep);
    samples = samples(keep);
    [times, order] = sort(times);
    neurons(index).spike_times = times;
    neurons(index).spike_samples = round(samples(order));
end
end

function pairs = empty_pairs()
pairs = struct('neuron_index_1', {}, 'neuron_index_2', {}, ...
    'cluster_id_1', {}, 'cluster_id_2', {}, ...
    'channel_1', {}, 'channel_2', {}, ...
    'area_1', {}, 'area_2', {}, ...
    'n_spikes_1', {}, 'n_spikes_2', {}, ...
    'n_interval_1', {}, 'n_interval_2', {}, ...
    'n_sync_1', {}, 'n_sync_2', {}, ...
    'sync_index_1', {}, 'sync_index_2', {}, ...
    'rank_spikes', {});
end

function neurons = empty_neurons()
neurons = struct('cluster_id', {}, 'channel', {}, 'area', {}, ...
    'spike_times', {}, 'spike_samples', {}, 'interval_index', {});
end
