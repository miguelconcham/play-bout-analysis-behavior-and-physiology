function rate = interval_spike_acg(spikes, t0, t1, bin_edges, closed_end)
%INTERVAL_SPIKE_ACG  Edge-corrected autocorrelogram on one time interval.
%
% rate is pair counts divided by (T - |lag|). closed_end true keeps spikes
% at t1. Closed_end false keeps [t0, t1), which is how adjoining windows
% avoid counting a boundary spike twice.

if nargin < 5 || isempty(closed_end)
    closed_end = true;
end
bin_edges = double(bin_edges(:))';
n_bins = max(numel(bin_edges) - 1, 0);
rate = nan(1, n_bins);
duration = double(t1) - double(t0);
if n_bins < 1 || ~isfinite(duration) || duration <= 0
    return
end
spikes = double(spikes(:));
if closed_end
    spikes = spikes(spikes >= t0 & spikes <= t1);
else
    spikes = spikes(spikes >= t0 & spikes < t1);
end
counts = pair_acg_counts(spikes, bin_edges);
rate = edge_correct_acg(counts, bin_edges, duration);
end
