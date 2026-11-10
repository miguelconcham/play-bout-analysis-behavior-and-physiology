function merged = merge_time_intervals(intervals)
%MERGE_TIME_INTERVALS  Sort intervals and join any that touch or overlap.
%
% intervals is N x 2, [start, end] in seconds. A row with a nonfinite
% endpoint, or with end <= start, is dropped. Touching intervals are
% joined so a boundary spike is not given two homes and a zero-length
% gap is not left between them.

merged = zeros(0, 2);
if isempty(intervals)
    return
end
intervals = double(intervals);
if size(intervals, 2) < 2
    error('merge_time_intervals:size', 'Intervals must have start and end columns.');
end
intervals = intervals(:, 1:2);
keep = all(isfinite(intervals), 2) & intervals(:, 2) > intervals(:, 1);
intervals = sortrows(intervals(keep, :), 1);
if isempty(intervals)
    return
end
merged = intervals(1, :);
for i = 2:size(intervals, 1)
    if intervals(i, 1) <= merged(end, 2)
        merged(end, 2) = max(merged(end, 2), intervals(i, 2));
    else
        merged(end + 1, :) = intervals(i, :); %#ok<AGROW>
    end
end
end
