function [pairs, kept_neurons, coverage] = select_synchronous_pairs(neurons, bout_starts, options)
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
% Pairs are chosen to observe neurons, not to multiply combinations. A
% pair is added only when at least one of its neurons is not yet in the
% sample and that neuron's area is still below min_neuron_fraction. Among
% those pairs, one that covers two new neurons comes first, then the area
% furthest from the target, then the larger spike count. A neuron already
% in the sample is reused only as the partner that lets a new neuron in.
% max_pairs_per_neuron, max_pairs_per_session, and max_pairs still cap the
% sample when they are finite. The default is no cap, so the limit is the
% coverage target.
%
% coverage counts good units in the areas (n_area), how many have at least
% min_spikes, how many have a partner that can support a null, and how many
% were kept. The target is applied per area to the neurons with at least
% min_spikes.
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
coverage = blank_coverage(options);
if numel(neurons) < 2 || isempty(bout_starts) || options.max_pairs < 1
    coverage = coverage_report(neurons, false(numel(neurons), 1), ...
        false(numel(neurons), 1), options);
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
has_partner = false(numel(neurons), 1);
for pair_index = 1:numel(eligible)
    has_partner(eligible(pair_index).neuron_index_1) = true;
    has_partner(eligible(pair_index).neuron_index_2) = true;
end
covered = false(numel(neurons), 1);
if ~isempty(eligible)
    use_count = zeros(numel(neurons), 1);
    taken = 0;
    limit = min(options.max_pairs, options.max_pairs_per_session);
    while taken < limit && ~coverage_met(covered, neurons, options)
        best = 0;
        best_new = -1;
        best_deficit = -inf;
        best_rank = -inf;
        for pair_index = 1:numel(eligible)
            pair = eligible(pair_index);
            if use_count(pair.neuron_index_1) >= options.max_pairs_per_neuron ...
                    || use_count(pair.neuron_index_2) >= options.max_pairs_per_neuron
                continue
            end
            [n_new, deficit, useful] = cover_value( ...
                pair, covered, neurons, options.min_neuron_fraction, options.min_spikes);
            if ~useful
                continue
            end
            if n_new > best_new || (n_new == best_new && deficit > best_deficit + 1e-12) ...
                    || (n_new == best_new && abs(deficit - best_deficit) <= 1e-12 ...
                    && pair.rank_spikes > best_rank)
                best = pair_index;
                best_new = n_new;
                best_deficit = deficit;
                best_rank = pair.rank_spikes;
            end
        end
        if best == 0
            break
        end
        pair = eligible(best);
        pairs(end + 1) = pair; %#ok<AGROW>
        use_count(pair.neuron_index_1) = use_count(pair.neuron_index_1) + 1;
        use_count(pair.neuron_index_2) = use_count(pair.neuron_index_2) + 1;
        covered(pair.neuron_index_1) = true;
        covered(pair.neuron_index_2) = true;
        taken = taken + 1;
    end
end

kept_index = find(covered);
for index = 1:numel(kept_index)
    kept_neurons(end + 1) = neurons(kept_index(index)); %#ok<AGROW>
end
coverage = coverage_report(neurons, covered, has_partner, options);
end

function [n_new, deficit, useful] = cover_value(pair, covered, neurons, target, min_spikes)
ends = [pair.neuron_index_1, pair.neuron_index_2];
n_new = 0;
deficit = 0;
useful = false;
for end_index = 1:numel(ends)
    neuron_index = ends(end_index);
    if covered(neuron_index)
        continue
    end
    n_new = n_new + 1;
    area_gap = area_deficit(neurons(neuron_index).area, covered, neurons, target, min_spikes);
    deficit = deficit + area_gap;
    if area_gap > 0
        useful = true;
    end
end
end

function gap = area_deficit(area_name, covered, neurons, target, min_spikes)
n_spiking = 0;
n_covered = 0;
for index = 1:numel(neurons)
    if ~area_is(neurons(index).area, area_name)
        continue
    end
    if numel(neurons(index).spike_times) < min_spikes
        continue
    end
    n_spiking = n_spiking + 1;
    if covered(index)
        n_covered = n_covered + 1;
    end
end
if n_spiking < 1 || n_covered >= target * n_spiking - 1e-9
    gap = 0;
else
    gap = target - n_covered / n_spiking;
end
end

function met = coverage_met(covered, neurons, options)
names = area_names_of(neurons);
if isempty(names)
    met = true;
    return
end
met = true;
for index = 1:numel(names)
    if area_deficit(names{index}, covered, neurons, ...
            options.min_neuron_fraction, options.min_spikes) > 0
        met = false;
        return
    end
end
end

function names = area_names_of(neurons)
names = {};
for index = 1:numel(neurons)
    area_name = char(neurons(index).area);
    if isempty(area_name)
        continue
    end
    if ~any(strcmp(names, area_name))
        names{end + 1} = area_name; %#ok<AGROW>
    end
end
end

function same = area_is(area_name, target_name)
same = strcmp(char(area_name), char(target_name));
end

function coverage = blank_coverage(options)
coverage = coverage_report(empty_neurons(), false(0, 1), false(0, 1), options);
end

function coverage = coverage_report(neurons, covered, has_partner, options)
if numel(covered) ~= numel(neurons)
    covered = false(numel(neurons), 1);
end
if numel(has_partner) ~= numel(neurons)
    has_partner = false(numel(neurons), 1);
end
is_spiking = false(numel(neurons), 1);
for index = 1:numel(neurons)
    is_spiking(index) = numel(neurons(index).spike_times) >= options.min_spikes;
end
coverage.n_area = numel(neurons);
coverage.n_with_spikes = nnz(is_spiking);
coverage.n_with_partner = nnz(has_partner & is_spiking);
coverage.n_covered = nnz(covered & is_spiking);
coverage.n_below_min_spikes = coverage.n_area - coverage.n_with_spikes;
coverage.n_spiking_without_partner = coverage.n_with_spikes - coverage.n_with_partner;
coverage.fraction_of_area = fraction_or_nan(coverage.n_covered, coverage.n_area);
coverage.fraction_of_spiking = fraction_or_nan(coverage.n_covered, coverage.n_with_spikes);
coverage.target = options.min_neuron_fraction;
coverage.target_met = coverage_met(covered, neurons, options);
coverage.areas = area_rows(neurons, covered, is_spiking);
end

function rows = area_rows(neurons, covered, is_spiking)
rows = struct('area', {}, 'n_area', {}, 'n_with_spikes', {}, ...
    'n_covered', {}, 'fraction_of_spiking', {});
names = area_names_of(neurons);
for name_index = 1:numel(names)
    area_name = names{name_index};
    in_area = false(numel(neurons), 1);
    for index = 1:numel(neurons)
        in_area(index) = area_is(neurons(index).area, area_name);
    end
    n_spiking = nnz(in_area & is_spiking);
    n_covered = nnz(in_area & is_spiking & covered);
    rows(end + 1).area = area_name; %#ok<AGROW>
    rows(end).n_area = nnz(in_area);
    rows(end).n_with_spikes = n_spiking;
    rows(end).n_covered = n_covered;
    rows(end).fraction_of_spiking = fraction_or_nan(n_covered, n_spiking);
end
end

function value = fraction_or_nan(numerator, denominator)
if denominator < 1
    value = NaN;
else
    value = numerator / denominator;
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
