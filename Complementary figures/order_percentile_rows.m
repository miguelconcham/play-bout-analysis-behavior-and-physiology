function [percentile, logit, score, order] = order_percentile_rows(percentile, lag, window_s, n_shuffles)
%ORDER_PERCENTILE_ROWS  Sort neurons by mean percentile inside a lag window.
%
% percentile is neurons by lags, on a 0-100 scale. score is the mean
% percentile at lags inside window_s (seconds, inclusive). The highest
% score becomes row 1, which is the top row of imagesc. logit is
% log10(p/(1-p)) for every neuron, with 0 and 100 pulled inside (0, 1).
% order indexes the input rows in that sorted order.

percentile = double(percentile);
lag = double(lag(:).');
if isempty(percentile)
    logit = percentile;
    score = zeros(0, 1);
    return
end
if size(percentile, 2) ~= numel(lag)
    error('order_percentile_rows:axis', ...
        'The percentile matrix has %d lags and the axis has %d.', ...
        size(percentile, 2), numel(lag));
end
in_window = lag >= window_s(1) & lag <= window_s(2);
if ~any(in_window)
    error('order_percentile_rows:window', ...
        'No lag falls inside %.3f to %.3f s.', window_s(1), window_s(2));
end
score = mean(percentile(:, in_window), 2, 'omitnan');
finite_score = isfinite(score);
[~, finite_order] = sort(score(finite_score), 'descend');
finite_idx = find(finite_score);
order = [finite_idx(finite_order); find(~finite_score)];
percentile = percentile(order, :);
score = score(order);
logit = correct_percentile_logit(percentile, n_shuffles);
end
