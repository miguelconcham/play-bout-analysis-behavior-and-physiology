function result = neuron_pair_ccg(reference_spikes, target_spikes, windows, bin_edges)
%NEURON_PAIR_CCG  One neuron pair, counted with crosscorrelogram.
%
% Each Partners-session interval is passed to crosscorrelogram_bin_counts.
% Counts are added across intervals. conditional_rate is target spikes per
% second around one interior reference spike. target_rate is the target
% neuron's own rate in the intervals that supplied those references.

bin_edges = double(bin_edges(:)');
n_bins = numel(bin_edges) - 1;
result = empty_pair_result(n_bins);
if n_bins < 1
    return
end
windows = merge_time_intervals(windows);
if isempty(windows)
    return
end
gaps = diff(bin_edges);
bin_w = gaps(1);
margin = max(abs(bin_edges(1)), abs(bin_edges(end)));
reference_spikes = double(reference_spikes(:));
target_spikes = double(target_spikes(:));
reference_spikes = sort(reference_spikes(isfinite(reference_spikes)));
target_spikes = sort(target_spikes(isfinite(target_spikes)));
counts = zeros(1, n_bins);
n_reference = 0;
n_target = 0;
duration = 0;
expected_counts = 0;
n_windows = 0;
for w = 1:size(windows, 1)
    t0 = windows(w, 1);
    t1 = windows(w, 2);
    window_duration = t1 - t0;
    if ~(window_duration > 2 * margin)
        continue
    end
    reference = reference_spikes( ...
        reference_spikes >= (t0 + margin) & reference_spikes <= (t1 - margin));
    if isempty(reference)
        continue
    end
    target = target_spikes(target_spikes >= t0 & target_spikes <= t1);
    [these_counts, these_reference] = crosscorrelogram_bin_counts( ...
        reference, target, bin_edges);
    if ~any(isfinite(these_counts))
        continue
    end
    these_counts(~isfinite(these_counts)) = 0;
    this_rate = numel(target) / window_duration;
    counts = counts + these_counts;
    n_reference = n_reference + these_reference;
    n_target = n_target + numel(target);
    duration = duration + window_duration;
    expected_counts = expected_counts + these_reference * this_rate * bin_w;
    n_windows = n_windows + 1;
end
result.counts = counts;
result.n_reference = n_reference;
result.n_target = n_target;
result.duration = duration;
result.n_windows = n_windows;
if n_reference < 1 || ~(duration > 0)
    result.counts = nan(1, n_bins);
    return
end
result.conditional_rate = counts / (n_reference * bin_w);
result.target_rate = expected_counts / (n_reference * bin_w);
end


function result = empty_pair_result(n_bins)
result.counts = nan(1, n_bins);
result.n_reference = 0;
result.n_target = 0;
result.duration = 0;
result.n_windows = 0;
result.conditional_rate = nan(1, n_bins);
result.target_rate = nan;
end
