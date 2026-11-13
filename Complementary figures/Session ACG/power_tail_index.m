function [low_index, high_index, low_threshold, high_threshold] = power_tail_index(score, tail_percent)
%POWER_TAIL_INDEX  Bottom and top tails of one score per window.
%
% tail_percent is the low tail. The high tail is 100 minus that value, so
% 25 keeps the bottom 25% and the top 25%. Percentiles use
% linear_percentile, on the finite scores only. The tails are left empty
% unless at least four windows have a finite score and the high threshold
% is strictly above the low one.

score = double(score(:));
n_windows = numel(score);
low_index = false(n_windows, 1);
high_index = false(n_windows, 1);
low_threshold = NaN;
high_threshold = NaN;
if ~isfinite(tail_percent) || tail_percent <= 0 || tail_percent >= 50
    error('power_tail_index:tail', 'tail_percent must be between 0 and 50.');
end
finite = isfinite(score);
if sum(finite) < 4
    return
end
low_threshold = linear_percentile(score(finite), tail_percent);
high_threshold = linear_percentile(score(finite), 100 - tail_percent);
if ~(isfinite(low_threshold) && isfinite(high_threshold) && high_threshold > low_threshold)
    return
end
low_index = finite & score <= low_threshold;
high_index = finite & score >= high_threshold;
end
