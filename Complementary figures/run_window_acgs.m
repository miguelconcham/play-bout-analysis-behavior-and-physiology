function acg = run_window_acgs(spike_times, run_starts, run_ends, bin_edges)
%RUN_WINDOW_ACGS  Edge-corrected autocorrelogram of each firing run.
%
% spike_times are in seconds. Lags are limited to the bin edges, and each
% row is pair counts divided by the overlap (T - |lag|). Bins past the end
% of a run are NaN.

spike_times = double(spike_times(:));
run_starts = double(run_starts(:));
run_ends = double(run_ends(:));
bin_edges = bin_edges(:)';
n_runs = numel(run_starts);
n_bins = max(numel(bin_edges) - 1, 0);
acg = nan(n_runs, n_bins);
if n_runs < 1 || n_bins < 1 || numel(run_ends) ~= n_runs
    return
end
range_limits = [bin_edges(1), bin_edges(end)];
for j = 1:n_runs
    if ~isfinite(run_starts(j)) || ~isfinite(run_ends(j)) || run_ends(j) <= run_starts(j)
        continue
    end
    in_run = spike_times >= run_starts(j) & spike_times <= run_ends(j);
    lags = autocorrelogram(spike_times(in_run), range_limits);
    counts = histcounts(lags, bin_edges);
    acg(j, :) = edge_correct_acg(counts, bin_edges, run_ends(j) - run_starts(j));
end
end
