function y_lim = ccg_mean_ylim(mean_curve, ci_low, ci_high, centers, x_lim)
%CCG_MEAN_YLIM  Y limits that sit just outside the mean and its interval.
%
% Only lags inside x_lim are used. A small pad keeps the line off the
% frame. The pad is a fraction of that range, so a flat mean stays in a
% narrow window.

shown = centers >= x_lim(1) & centers <= x_lim(2);
values = [mean_curve(shown), ci_low(shown), ci_high(shown)];
values = values(isfinite(values));
if isempty(values)
    y_lim = [-1 1];
    return
end
low = min(values);
high = max(values);
span = high - low;
pad = 0.08 * span;
if ~(pad > 0)
    pad = 0.05 * max(abs([low, high]));
end
if ~(pad > 0)
    pad = 0.05;
end
y_lim = [low - pad, high + pad];
end
