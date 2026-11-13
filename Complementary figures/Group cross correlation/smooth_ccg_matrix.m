function smoothed = smooth_ccg_matrix(curves, centers, smooth_ms)
%SMOOTH_CCG_MATRIX  Gaussian-smooth every correlogram row.
%
% smooth_ms is the Gaussian standard deviation in milliseconds. A width
% that is not larger than the lag bin, in milliseconds, leaves every row
% as it was. Missing bins stay missing.

smoothed = curves;
if isempty(curves)
    return
end
bin_ms = median(diff(centers)) * 1000;
if ~isfinite(smooth_ms) || smooth_ms <= bin_ms
    return
end
sigma_bins = smooth_ms / bin_ms;
for row = 1:size(curves, 1)
    if any(isfinite(curves(row, :)))
        smoothed(row, :) = gaussian_smooth_acg(curves(row, :), sigma_bins);
    end
end
end
