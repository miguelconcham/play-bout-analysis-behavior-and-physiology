function counts = pool_crosscorr_acgs(spikes, starts, ends, bin_edges)
%POOL_CROSSCORR_ACGS  Sum interior-reference autocorrelograms across windows.
%
% Each window uses crosscorr_acg_counts. The windows are [start, end) and
% do not overlap. Counts are added, not divided by (T-|lag|), because each
% reference spike already has the full lag range inside its window.

bin_edges = double(bin_edges(:))';
n_bins = max(numel(bin_edges) - 1, 0);
counts = nan(1, n_bins);
if n_bins < 1
    return
end
starts = double(starts(:));
ends = double(ends(:));
if numel(starts) ~= numel(ends)
    error('pool_crosscorr_acgs:windows', ...
        'Window starts and ends must have the same length.');
end
if isempty(starts)
    return
end
total = zeros(1, n_bins);
used = false;
for w = 1:numel(starts)
    these = crosscorr_acg_counts(spikes, starts(w), ends(w), bin_edges, false);
    if ~any(isfinite(these))
        continue
    end
    these(~isfinite(these)) = 0;
    total = total + these;
    used = true;
end
if used
    counts = total;
end
end
