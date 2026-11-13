function band = session_acg_null_band(null_rows)
%SESSION_ACG_NULL_BAND  2.5%, median, and 97.5% of saved shuffle rows.
%
% null_rows is shuffles by lag bins. NaNs are left out of that bin. A bin
% with no finite shuffle stays NaN. low, median, and high are one row.

null_rows = double(null_rows);
if isempty(null_rows)
    band = empty_band(0);
    return
end
if ~ismatrix(null_rows)
    error('session_acg:band', 'Shuffle rows must be a matrix.');
end
n_bins = size(null_rows, 2);
band = empty_band(n_bins);
for bin = 1:n_bins
    column = null_rows(:, bin);
    column = column(isfinite(column));
    if isempty(column)
        continue
    end
    band.low(bin) = prctile(column, 2.5);
    band.median(bin) = median(column);
    band.high(bin) = prctile(column, 97.5);
end
end


function band = empty_band(n_bins)
band.low = nan(1, n_bins);
band.median = nan(1, n_bins);
band.high = nan(1, n_bins);
end
