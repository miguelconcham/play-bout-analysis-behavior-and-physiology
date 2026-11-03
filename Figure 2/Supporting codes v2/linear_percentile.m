function value = linear_percentile(x, percentile)
%LINEAR_PERCENTILE  Inclusive linear percentile.
%
% Sorted values are indexed at 1 + (n - 1) * (p/100). This is the same
% percentile used by linear_iqr. An empty input is NaN.

x = sort(double(x(:)));
n = numel(x);
if n == 0
    value = NaN;
    return
end
if n == 1
    value = x(1);
    return
end
rank = 1 + (percentile / 100) * (n - 1);
low = floor(rank);
high = ceil(rank);
low = min(max(low, 1), n);
high = min(max(high, 1), n);
if low == high
    value = x(low);
else
    value = x(low) + (rank - low) * (x(high) - x(low));
end
end
