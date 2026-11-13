function mu = pool_crosscorr_rows(rate, window_length, bin_edges)
%POOL_CROSSCORR_ROWS  Pool crosscorr rates by inner-window duration.
%
% Each row is pair counts divided by the inner duration, T - 2*margin,
% where margin is the lag limit. A short window contributes in proportion
% to that inner duration. A window shorter than two margins stays out.

rate = double(rate);
bin_edges = double(bin_edges(:)');
n_bins = size(rate, 2);
mu = nan(1, n_bins);
if n_bins < 1 || isempty(rate)
    return
end
window_length = double(window_length(:));
if isscalar(window_length) && size(rate, 1) > 1
    window_length = repmat(window_length, size(rate, 1), 1);
end
if numel(window_length) ~= size(rate, 1)
    error('pool_crosscorr_rows:rows', ...
        'Window length must have one entry per autocorrelogram row.');
end
margin = max(abs(bin_edges(1)), abs(bin_edges(end)));
weight = window_length - 2 * margin;
finite_row = isfinite(weight) & weight > 0 & any(isfinite(rate), 2);
weight(~finite_row) = 0;
covered = isfinite(rate) & weight > 0;
counts = rate;
counts(~covered) = 0;
den = sum(weight);
if ~(isfinite(den) && den > 0)
    return
end
num = sum(counts .* weight, 1);
mu = num ./ den;
mu(~any(covered, 1)) = NaN;
end
