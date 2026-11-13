function combinations = session_neuron_pair_ccgs(group_trains, group_ids, windows, bin_edges, min_reference_spikes)
%SESSION_NEURON_PAIR_CCGS  Every neuron pair in each directed group combination.
%
% group_trains{g} is a cell of spike trains. group_ids{g} holds the neuron
% id of each train. A group is not paired with itself. combination(k).counts
% has one row per neuron pair and one column per lag bin. A row whose
% reference count is below min_reference_spikes stays in the matrix and
% usable is false, so the mean later leaves it out.

names = group_ccg_names();
n_groups = numel(names.group_id);
if numel(group_trains) ~= n_groups || numel(group_ids) ~= n_groups
    error('group_ccg:pairs', 'Provide one train list and one id list per group.');
end
if ~(isscalar(min_reference_spikes) && min_reference_spikes >= 1)
    error('group_ccg:pairs', 'min_reference_spikes must be at least 1.');
end
n_bins = max(numel(bin_edges) - 1, 0);
n_combinations = numel(names.pair_reference);
combinations = repmat(blank_combination(names, 1, n_bins), n_combinations, 1);
for k = 1:n_combinations
    combination = blank_combination(names, k, n_bins);
    reference_trains = group_trains{combination.reference_index};
    target_trains = group_trains{combination.target_index};
    reference_ids = double(group_ids{combination.reference_index}(:));
    target_ids = double(group_ids{combination.target_index}(:));
    n_reference_neurons = numel(reference_trains);
    n_target_neurons = numel(target_trains);
    n_pairs = n_reference_neurons * n_target_neurons;
    combination.counts = nan(n_pairs, n_bins);
    combination.n_reference = zeros(n_pairs, 1);
    combination.n_target = zeros(n_pairs, 1);
    combination.duration = zeros(n_pairs, 1);
    combination.target_rate = nan(n_pairs, 1);
    combination.reference_id = zeros(n_pairs, 1);
    combination.target_id = zeros(n_pairs, 1);
    combination.usable = false(n_pairs, 1);
    row = 0;
    for a = 1:n_reference_neurons
        for b = 1:n_target_neurons
            row = row + 1;
            raw = neuron_pair_ccg(reference_trains{a}, target_trains{b}, windows, bin_edges);
            combination.reference_id(row) = reference_ids(a);
            combination.target_id(row) = target_ids(b);
            combination.n_reference(row) = raw.n_reference;
            combination.n_target(row) = raw.n_target;
            combination.duration(row) = raw.duration;
            combination.target_rate(row) = raw.target_rate;
            if any(isfinite(raw.counts))
                combination.counts(row, :) = raw.counts;
            end
            combination.usable(row) = raw.n_reference >= min_reference_spikes ...
                && isfinite(raw.target_rate) && raw.target_rate > 0 ...
                && any(isfinite(raw.conditional_rate));
        end
    end
    combination.n_pairs = n_pairs;
    combination.n_usable = sum(combination.usable);
    combinations(k) = combination;
end
end


function combination = blank_combination(names, index, n_bins)
reference = names.pair_reference(index);
target = names.pair_target(index);
combination.reference = names.group_id{reference};
combination.target = names.group_id{target};
combination.reference_index = reference;
combination.target_index = target;
combination.counts = zeros(0, n_bins);
combination.n_reference = zeros(0, 1);
combination.n_target = zeros(0, 1);
combination.duration = zeros(0, 1);
combination.target_rate = zeros(0, 1);
combination.reference_id = zeros(0, 1);
combination.target_id = zeros(0, 1);
combination.usable = false(0, 1);
combination.n_pairs = 0;
combination.n_usable = 0;
end
