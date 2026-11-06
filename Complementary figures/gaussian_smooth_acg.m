function y = gaussian_smooth_acg(y, sigma_bins)
%GAUSSIAN_SMOOTH_ACG  Gaussian smooth that leaves NaN bins unchanged.
%
% sigma_bins is the standard deviation in samples. The kernel is zero past
% three standard deviations. Smoothing stops at a NaN, so a missing lag is
% not filled from the other side and is not averaged into a finite bin.
% The center of a peak is weighted more than the shoulders, unlike a
% moving average of the same width.

y = double(y(:))';
if isempty(y) || ~isfinite(sigma_bins) || sigma_bins <= 0
    return
end
half = max(1, ceil(3 * sigma_bins));
offsets = -half:half;
kernel = exp(-0.5 * (offsets / sigma_bins) .^ 2);
out = y;
n = numel(y);
for i = 1:n
    if ~isfinite(y(i))
        continue
    end
    left = i;
    while left > 1 && (i - left) < half && isfinite(y(left - 1))
        left = left - 1;
    end
    right = i;
    while right < n && (right - i) < half && isfinite(y(right + 1))
        right = right + 1;
    end
    idx = left:right;
    w = kernel(idx - i + half + 1);
    out(i) = sum(w .* y(idx)) / sum(w);
end
y = out;
end
