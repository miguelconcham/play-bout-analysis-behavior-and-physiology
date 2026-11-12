function tails = play_power_tails(power, low_percent, high_percent)
%PLAY_POWER_TAILS  Low and high tails of one neuron or channel's play power.
%
% The cut is the same one used for firing runs in the locking-period
% autocorrelogram: at least four finite values, low at or below the
% low_percent percentile, high at or above high_percent. The default is
% 25 and 75. The middle is left out, and the tails are left empty when
% the high percentile is not strictly above the low one.

if nargin < 2 || isempty(low_percent)
    low_percent = 25;
end
if nargin < 3 || isempty(high_percent)
    high_percent = 75;
end
power = double(power(:));
n = numel(power);
tails.low = false(n, 1);
tails.high = false(n, 1);
tails.threshold_low = NaN;
tails.threshold_high = NaN;
finite_power = power(isfinite(power));
if numel(finite_power) < 4
    return
end
tails.threshold_low = linear_percentile(finite_power, low_percent);
tails.threshold_high = linear_percentile(finite_power, high_percent);
if isfinite(tails.threshold_low) && isfinite(tails.threshold_high) ...
        && tails.threshold_high > tails.threshold_low
    tails.low = isfinite(power) & power <= tails.threshold_low;
    tails.high = isfinite(power) & power >= tails.threshold_high;
end
end
