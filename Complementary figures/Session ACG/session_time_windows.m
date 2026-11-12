function [starts, ends] = session_time_windows(t0, t1, window_s)
%SESSION_TIME_WINDOWS  Non-overlapping windows that cover a session.
%
% Windows are [t0, t0+T), [t0+T, t0+2T), and so on. A leftover shorter
% than T is dropped, so every kept window has the same length.

t0 = double(t0);
t1 = double(t1);
window_s = double(window_s);
if ~(isfinite(t0) && isfinite(t1) && isfinite(window_s) && window_s > 0 && t1 > t0)
    starts = zeros(0, 1);
    ends = zeros(0, 1);
    return
end
n_windows = floor((t1 - t0) / window_s);
starts = t0 + (0:n_windows - 1)' * window_s;
ends = starts + window_s;
end
