function mu = logit_window_mean(logit, lag, window_s)
%LOGIT_WINDOW_MEAN  Mean logit inside an open lag window.
%
% window_s is [start stop] in seconds. A lag is inside when it is greater
% than start and less than stop. mu is one value per neuron. A neuron
% whose window is all missing stays NaN.

logit = double(logit);
lag = double(lag(:).');
if numel(window_s) ~= 2 || ~all(isfinite(window_s))
    error('logit_window_mean:window', 'The lag window must be two finite times.');
end
inside = lag > window_s(1) & lag < window_s(2);
if ~any(inside)
    error('logit_window_mean:window', ...
        'No lag falls strictly inside %.3f to %.3f s.', window_s(1), window_s(2));
end
if isempty(logit)
    mu = zeros(0, 1);
    return
end
if size(logit, 2) ~= numel(lag)
    error('logit_window_mean:axis', ...
        'The logit matrix has %d lags and the axis has %d.', size(logit, 2), numel(lag));
end
mu = mean(logit(:, inside), 2, 'omitnan');
end
