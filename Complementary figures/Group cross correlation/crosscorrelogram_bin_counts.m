function [counts, n_reference] = crosscorrelogram_bin_counts(reference_spikes, target_spikes, bin_edges)
%CROSSCORRELOGRAM_BIN_COUNTS  1 ms histogram of repository crosscorrelogram lags.
%
% Lags come from Custom functions/crosscorrelogram, not from xcorr.
% That function skips the first reference spike and stops one sample short
% of the last partner. A time is placed before the first reference and
% after the last partner, as crosscorr_acg_counts does, so every real
% spike stays in the correlogram. The lags are then counted on bin_edges.
% A lag of zero is kept: the two trains belong to different neurons.

bin_edges = double(bin_edges(:)');
n_bins = max(numel(bin_edges) - 1, 0);
reference_spikes = finite_column(reference_spikes);
target_spikes = finite_column(target_spikes);
n_reference = numel(reference_spikes);
counts = zeros(1, n_bins);
if n_bins < 1 || n_reference < 1 || isempty(target_spikes)
    if n_reference < 1 || n_bins < 1
        counts = nan(1, n_bins);
    end
    return
end
gaps = diff(bin_edges);
if any(abs(gaps - gaps(1)) > 1e-9) || ~(gaps(1) > 0)
    error('group_ccg:bins', 'CCG bin edges must be evenly spaced and increasing.');
end
if exist('crosscorrelogram', 'file') ~= 2
    error('group_ccg:crosscorrelogram', ...
        'crosscorrelogram.m is not on the path. Run add_repo_paths.');
end
margin = max(abs(bin_edges(1)), abs(bin_edges(end)));
reference_spikes = [reference_spikes(1) - (margin + 1); reference_spikes];
target_spikes = [target_spikes; target_spikes(end) + (margin + 1)];
offsets = crosscorrelogram(reference_spikes, target_spikes, [bin_edges(1), bin_edges(end)]);
offsets = double(offsets(:));
offsets = offsets(isfinite(offsets));
if exist('histcounts', 'builtin') == 5 || exist('histcounts', 'file') == 2
    counts = histcounts(offsets, bin_edges);
else
    counts = bin_lags(offsets, bin_edges);
end
end


function counts = bin_lags(offsets, bin_edges)
% Half-open bins, and the last bin includes its right edge, as histcounts does.
n_bins = numel(bin_edges) - 1;
counts = zeros(1, n_bins);
if isempty(offsets)
    return
end
bin_w = bin_edges(2) - bin_edges(1);
scaled = (offsets - bin_edges(1)) / bin_w;
bin = floor(scaled + 1e-8) + 1;
bin(bin == n_bins + 1) = n_bins;
keep = bin >= 1 & bin <= n_bins;
if any(keep)
    counts = counts + accumarray(bin(keep), 1, [n_bins, 1])';
end
end


function spikes = finite_column(spikes)
spikes = double(spikes(:));
spikes = sort(spikes(isfinite(spikes)));
end
