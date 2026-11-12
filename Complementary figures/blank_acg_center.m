function Y = blank_acg_center(X, centers, flank_s)
%BLANK_ACG_CENTER  Set autocorrelogram bins that cover lag 0 to NaN.
%
% flank_s is the half-width, in seconds, of the refractory gap. A bin is
% blanked when its interval overlaps [-flank_s, flank_s]. The same gap is
% removed before the locking-period curves are normalized.

Y = X;
if isempty(X) || isempty(centers)
    return
end
t = double(centers(:)');
dt = median(diff(t));
if ~isfinite(dt) || dt <= 0
    dt = flank_s;
end
half = abs(dt) / 2;
lo = t - half;
hi = t + half;
Y(:, lo < flank_s & hi > -flank_s) = NaN;
end
