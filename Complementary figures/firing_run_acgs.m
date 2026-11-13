function [acg, n_eligible, already_corrected] = firing_run_acgs(spike_times, starts, ends, bin_edges, method, saved_acg, saved_corrected, count_method)
%FIRING_RUN_ACGS  One autocorrelogram row per closed firing run.
%
% 'eligible' always counts partners of spikes that can see each bin.
% Saved histograms cannot be converted to that count. 'overlap' with
% count 'pairs' uses saved_acg when it has one row per run and the same
% number of bins. Otherwise it rebuilds the requested count from the
% spikes inside each run.
%
% count_method is 'pairs', 'fft', or 'crosscorr'. It applies to overlap.
% 'pairs' is the default: pair counts / (T-|lag|). 'fft' uses that same
% correction on a binned FFT. 'crosscorr' insets the reference spikes by
% the lag limit. Eligible estimates do not take a second count. Saved
% parsing histograms are pair counts, so fft and crosscorr always rebuild.

method = canonical_acg_method(method);
if nargin < 8 || isempty(count_method)
    count_method = 'pairs';
else
    count_method = canonical_count_method(count_method);
end
starts = double(starts(:));
ends = double(ends(:));
n_runs = numel(starts);
n_bins = max(numel(bin_edges) - 1, 0);
n_eligible = [];
if strcmp(method, 'eligible')
    if ~strcmp(count_method, 'pairs')
        error('firing_run_acgs:count', ...
            ['Eligible-spike autocorrelograms keep their own count. ' ...
            'acg_count applies to overlap estimates.']);
    end
    [acg, n_eligible] = eligible_window_counts(spike_times, starts, ends, bin_edges, true);
    already_corrected = false;
    return
end
if strcmp(count_method, 'pairs') && nargin >= 6 && ~isempty(saved_acg) ...
        && size(saved_acg, 1) == n_runs && size(saved_acg, 2) == n_bins
    acg = saved_acg;
    if nargin < 7 || isempty(saved_corrected)
        already_corrected = false;
    else
        already_corrected = logical(saved_corrected);
    end
    return
end
acg = run_window_acgs(spike_times, starts, ends, bin_edges, count_method, true);
already_corrected = true;
end
