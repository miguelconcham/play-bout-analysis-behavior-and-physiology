function epochs = fixed_duration_epochs(intervals, duration)
%FIXED_DURATION_EPOCHS  Non-overlapping epochs that sit fully inside intervals.
%
% Each epoch has length duration. A leftover shorter than duration is
% dropped, so an epoch never crosses from a play bout into the time
% outside it. Intervals that touch are joined first.

epochs = zeros(0, 2);
if isempty(intervals) || ~isscalar(duration) || ~isfinite(duration) || duration <= 0
    return
end
intervals = merge_time_intervals(intervals);
for i = 1:size(intervals, 1)
    n_epochs = floor((intervals(i, 2) - intervals(i, 1)) / duration);
    if n_epochs < 1
        continue
    end
    starts = intervals(i, 1) + (0:n_epochs - 1)' * duration;
    epochs = [epochs; [starts, starts + duration]]; %#ok<AGROW>
end
end
