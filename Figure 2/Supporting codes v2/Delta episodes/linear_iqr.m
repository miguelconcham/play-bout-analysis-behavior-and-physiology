function value = linear_iqr(x)
%LINEAR_IQR  Interquartile range with linear percentile interpolation.
%
% Sorted values are indexed at 1 + (n - 1) * (p/100). This is the
% linear (numpy/MATLAB-style inclusive) percentile. One value has IQR 0.
% An empty input is NaN.

x = double(x(:));
if isempty(x)
    value = NaN;
    return
end
value = linear_percentile(x, 75) - linear_percentile(x, 25);
end
