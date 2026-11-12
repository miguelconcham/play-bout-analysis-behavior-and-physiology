function pooled = pool_selected_window_acgs(acg, bin_edges, window_length, already_corrected)
%POOL_SELECTED_WINDOW_ACGS  Overlap-weighted mean of window autocorrelograms.
%
% The value in a lag bin is the sum of pair counts divided by the sum of
% overlap (T - |lag|). A short window then does not count as much as a
% long one. already_corrected true means each row is already counts /
% overlap, and the counts are recovered before the sum. Bins whose
% overlap is shorter than one bin are ignored. This is the same pool used
% for the all-runs curve of the locking-period autocorrelogram.

bin_edges = bin_edges(:)';
n_bins = max(numel(bin_edges) - 1, 0);
pooled = nan(1, n_bins);
if n_bins < 1 || isempty(acg)
    return
end
acg = double(acg);
window_length = double(window_length(:));
if size(acg, 2) ~= n_bins
    error('pool_selected_window_acgs:width', ...
        'ACG columns (%d) do not match the bin edges.', size(acg, 2));
end
if isscalar(window_length) && size(acg, 1) > 1
    window_length = repmat(window_length, size(acg, 1), 1);
end
if numel(window_length) ~= size(acg, 1)
    error('pool_selected_window_acgs:rows', ...
        'Window length must have one entry per ACG row.');
end
bin_w = mean(diff(bin_edges));
centers = bin_edges(1:end - 1) + 0.5 * bin_w;
overlap = window_length - abs(centers);
overlap(overlap < bin_w | ~isfinite(overlap)) = NaN;
counts = acg;
if already_corrected
    counts = counts .* overlap;
end
counts(~isfinite(overlap)) = NaN;
den = sum(nan_to_zero(overlap), 1);
num = sum(nan_to_zero(counts), 1);
n_finite = sum(isfinite(counts) & isfinite(overlap), 1);
ok = den > 0 & n_finite > 0;
pooled(ok) = num(ok) ./ den(ok);
end


function filled = nan_to_zero(x)
filled = x;
filled(isnan(filled)) = 0;
end
