function pieces = intersect_time_intervals(a, b)
%INTERSECT_TIME_INTERVALS  Overlap of two sets of time intervals.
%
% a and b are N x 2 and M x 2, [start, end] in seconds. The result is
% sorted and merged. A 1 x 2 row is one window, so a time limit can be
% passed directly.

if isempty(a) || isempty(b)
    pieces = zeros(0, 2);
    return
end
a = double(a);
b = double(b);
if isvector(a) && numel(a) == 2
    a = a(:)';
end
if isvector(b) && numel(b) == 2
    b = b(:)';
end
pieces = zeros(0, 2);
for i = 1:size(a, 1)
    for j = 1:size(b, 1)
        start_t = max(a(i, 1), b(j, 1));
        end_t = min(a(i, 2), b(j, 2));
        if end_t > start_t && all(isfinite([start_t, end_t]))
            pieces(end + 1, :) = [start_t, end_t]; %#ok<AGROW>
        end
    end
end
pieces = merge_time_intervals(pieces);
end
