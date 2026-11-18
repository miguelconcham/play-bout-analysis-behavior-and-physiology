function [acg, centers] = reflect_positive_acg(positive, pos_lag)
%REFLECT_POSITIVE_ACG  Two-sided curve that folds back to the positive side.
%
% The displayed observed autocorrelogram is already the average of the
% positive and negative lags. Reflecting it gives folded_acg_spectrum the
% same positive side, so the spectrum is the spectrum of that curve.

pos_lag = double(pos_lag(:)).';
positive = double(positive);
if isempty(positive) || isempty(pos_lag)
    acg = zeros(0, 0);
    centers = zeros(1, 0);
    return
end
if isvector(positive)
    positive = positive(:).';
end
n = min(size(positive, 2), numel(pos_lag));
positive = positive(:, 1:n);
pos_lag = pos_lag(1:n);
keep = pos_lag > 0 & isfinite(pos_lag);
positive = positive(:, keep);
pos_lag = pos_lag(keep);
if isempty(pos_lag)
    acg = zeros(size(positive, 1), 0);
    centers = zeros(1, 0);
    return
end
[pos_lag, order] = sort(pos_lag);
positive = positive(:, order);
acg = [fliplr(positive), positive];
centers = [-fliplr(pos_lag), pos_lag];
end
