function [high_percent, low_percent] = percentile_tail_fractions(stack, high_cut, low_cut)
%PERCENTILE_TAIL_FRACTIONS  Percent of neurons past the two percentile cuts.
%
% stack is neurons by lags, on a 0-100 percentile scale. high_percent is
% the percent of finite neurons whose percentile is at least high_cut
% (default 95). low_percent is the percent at or below low_cut (default 5).
% A NaN percentile is left out of that lag's denominator.

if nargin < 2 || isempty(high_cut)
    high_cut = 95;
end
if nargin < 3 || isempty(low_cut)
    low_cut = 5;
end
stack = double(stack);
finite = isfinite(stack);
n = sum(finite, 1);
high_percent = 100 * sum(finite & stack >= high_cut, 1) ./ n;
low_percent = 100 * sum(finite & stack <= low_cut, 1) ./ n;
high_percent(n == 0) = NaN;
low_percent(n == 0) = NaN;
end
