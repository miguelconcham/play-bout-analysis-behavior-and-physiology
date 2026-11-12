function [acg, n_eligible, already_corrected] = firing_run_acgs(spike_times, starts, ends, bin_edges, method, saved_acg, saved_corrected)
%FIRING_RUN_ACGS  One autocorrelogram row per closed firing run.
%
% 'eligible' always counts partners of spikes that can see each bin.
% Saved histograms cannot be converted to that count. 'overlap' uses
% saved_acg when it has one row per run and the same number of bins.
% Otherwise it rebuilds pair counts / (T-|lag|) from the spikes inside
% each run.

method = canonical_acg_method(method);
starts = double(starts(:));
ends = double(ends(:));
n_runs = numel(starts);
n_bins = max(numel(bin_edges) - 1, 0);
n_eligible = [];
if strcmp(method, 'eligible')
    [acg, n_eligible] = eligible_window_counts(spike_times, starts, ends, bin_edges, true);
    already_corrected = false;
    return
end
if nargin >= 6 && ~isempty(saved_acg) && size(saved_acg, 1) == n_runs && size(saved_acg, 2) == n_bins
    acg = saved_acg;
    if nargin < 7 || isempty(saved_corrected)
        already_corrected = false;
    else
        already_corrected = logical(saved_corrected);
    end
    return
end
acg = run_window_acgs(spike_times, starts, ends, bin_edges);
already_corrected = true;
end
