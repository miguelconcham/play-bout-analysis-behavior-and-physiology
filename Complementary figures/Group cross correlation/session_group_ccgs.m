function pairs = session_group_ccgs(group_spikes, n_neurons, windows, bin_edges, min_reference_spikes)
%SESSION_GROUP_CCGS  One directed correlogram for every pair of distinct groups.
%
% group_spikes{g} holds the pooled spike times of group g. n_neurons(g) is
% how many neurons entered that pool. A pair is usable when both groups
% have a neuron and the interior reference count reaches
% min_reference_spikes. An unusable pair keeps its reference and target
% names, and its rate curves are NaN.

names = group_ccg_names();
if numel(group_spikes) ~= numel(names.group_id) || numel(n_neurons) ~= numel(names.group_id)
    error('group_ccg:pairs', 'Provide one spike pool and one neuron count per group.');
end
if ~(isscalar(min_reference_spikes) && min_reference_spikes >= 1)
    error('group_ccg:pairs', 'min_reference_spikes must be at least 1.');
end
n_bins = max(numel(bin_edges) - 1, 0);
n_pairs = numel(names.pair_reference);
pairs = repmat(blank_pair(names, 1, n_bins), n_pairs, 1);
for k = 1:n_pairs
    pair = blank_pair(names, k, n_bins);
    raw = group_pair_ccg( ...
        group_spikes{pair.reference_index}, group_spikes{pair.target_index}, ...
        windows, bin_edges);
    pair.counts = raw.counts;
    pair.n_reference = raw.n_reference;
    pair.n_target = raw.n_target;
    pair.duration = raw.duration;
    pair.n_windows = raw.n_windows;
    pair.n_neurons_reference = n_neurons(pair.reference_index);
    pair.n_neurons_target = n_neurons(pair.target_index);
    pair.usable = n_neurons(pair.reference_index) >= 1 ...
        && n_neurons(pair.target_index) >= 1 ...
        && raw.n_reference >= min_reference_spikes ...
        && any(isfinite(raw.excess_rate));
    if pair.usable
        pair.conditional_rate = raw.conditional_rate;
        pair.target_rate = raw.target_rate;
        pair.excess_rate = raw.excess_rate;
        pair.rate_ratio = raw.rate_ratio;
    end
    pairs(k) = pair;
end
end


function pair = blank_pair(names, index, n_bins)
reference = names.pair_reference(index);
target = names.pair_target(index);
pair.reference = names.group_id{reference};
pair.target = names.group_id{target};
pair.reference_index = reference;
pair.target_index = target;
pair.counts = nan(1, n_bins);
pair.n_reference = 0;
pair.n_target = 0;
pair.duration = 0;
pair.n_windows = 0;
pair.conditional_rate = nan(1, n_bins);
pair.target_rate = nan;
pair.excess_rate = nan(1, n_bins);
pair.rate_ratio = nan(1, n_bins);
pair.n_neurons_reference = 0;
pair.n_neurons_target = 0;
pair.usable = false;
end
