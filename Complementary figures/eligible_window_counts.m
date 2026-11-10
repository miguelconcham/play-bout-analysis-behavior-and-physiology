function [partner_counts, n_eligible] = eligible_window_counts(spike_times, window_starts, window_ends, bin_edges, closed_end)
%ELIGIBLE_WINDOW_COUNTS  Partner counts from spikes that can see each bin.
%
% A spike at time t contributes to bin [L, R) only when [t+L, t+R] lies
% inside the window. Spikes that cannot see the bin are left out of that
% bin. An eligible spike with no partner in the bin contributes zero.
% partner_counts(w, b) is the sum of those partner counts in window w.
% n_eligible(w, b) is how many spikes in that window could see the bin.
% A window with no eligible spike keeps a zero, and the pool turns a bin
% whose total is zero into missing rather than a rate of zero.
%
% closed_end false uses [start, end), which is the play-epoch rule.
% closed_end true also keeps a spike at the window end, which is the
% firing-run rule. This count is not divided by T-|lag|.

if nargin < 5 || isempty(closed_end)
    closed_end = false;
end
spike_times = double(spike_times(:));
spike_times = spike_times(isfinite(spike_times));
spike_times = sort(spike_times);
window_starts = double(window_starts(:));
window_ends = double(window_ends(:));
bin_edges = double(bin_edges(:)');
n_bins = max(numel(bin_edges) - 1, 0);
n_win = numel(window_starts);
partner_counts = zeros(n_win, n_bins);
n_eligible = zeros(n_win, n_bins);
if n_win < 1 || n_bins < 1 || numel(window_ends) ~= n_win
    return
end
left = bin_edges(1:end - 1);
right = bin_edges(2:end);
tol = 1e-9;
for w = 1:n_win
    w0 = window_starts(w);
    w1 = window_ends(w);
    if ~isfinite(w0) || ~isfinite(w1) || w1 <= w0
        partner_counts(w, :) = NaN;
        n_eligible(w, :) = NaN;
        continue
    end
    if closed_end
        in_window = spike_times >= w0 & spike_times <= w1;
    else
        in_window = spike_times >= w0 & spike_times < w1;
    end
    spikes = spike_times(in_window);
    n_spikes = numel(spikes);
    if n_spikes < 1
        continue
    end
    eligible = false(n_spikes, n_bins);
    for i = 1:n_spikes
        eligible(i, :) = (spikes(i) + left) >= (w0 - tol) ...
            & (spikes(i) + right) <= (w1 + tol);
    end
    n_eligible(w, :) = sum(eligible, 1);
    if n_spikes < 2 || ~any(eligible(:))
        continue
    end
    for i = 1:n_spikes
        if ~any(eligible(i, :))
            continue
        end
        lags = spikes - spikes(i);
        lags(i) = [];
        counts = lag_bin_counts(lags, bin_edges);
        counts = counts(:)';
        partner_counts(w, :) = partner_counts(w, :) + counts .* eligible(i, :);
    end
end
end


function counts = lag_bin_counts(lags, edges)
if exist('histcounts', 'builtin') == 5 || exist('histcounts', 'file') == 2
    counts = histcounts(lags, edges);
    return
end
counts = zeros(1, numel(edges) - 1);
lags = double(lags(:));
lags = lags(isfinite(lags));
for b = 1:numel(counts)
    if b < numel(counts)
        in_bin = lags >= edges(b) & lags < edges(b + 1);
    else
        in_bin = lags >= edges(b) & lags <= edges(b + 1);
    end
    counts(b) = sum(in_bin);
end
end
