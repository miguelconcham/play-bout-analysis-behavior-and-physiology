function gaps = outside_play_intervals(windows, play_bouts)
%OUTSIDE_PLAY_INTERVALS  Time inside windows that is not inside a play bout.
%
% windows and play_bouts are N x 2, [start, end] in seconds. Play that
% falls outside a window is ignored. Play that crosses a window edge is
% clipped to the window before the gap is cut. Windows with no play are
% returned whole.

windows = merge_time_intervals(windows);
play_bouts = merge_time_intervals(play_bouts);
gaps = zeros(0, 2);
for i = 1:size(windows, 1)
    w0 = windows(i, 1);
    w1 = windows(i, 2);
    if isempty(play_bouts)
        inside = zeros(0, 2);
    else
        overlap = play_bouts(:, 2) > w0 & play_bouts(:, 1) < w1;
        inside = play_bouts(overlap, :);
    end
    if ~isempty(inside)
        inside(:, 1) = max(inside(:, 1), w0);
        inside(:, 2) = min(inside(:, 2), w1);
        inside = inside(inside(:, 2) > inside(:, 1), :);
    end
    cursor = w0;
    for k = 1:size(inside, 1)
        if inside(k, 1) > cursor
            gaps(end + 1, :) = [cursor, inside(k, 1)]; %#ok<AGROW>
        end
        cursor = max(cursor, inside(k, 2));
    end
    if cursor < w1
        gaps(end + 1, :) = [cursor, w1]; %#ok<AGROW>
    end
end
end
