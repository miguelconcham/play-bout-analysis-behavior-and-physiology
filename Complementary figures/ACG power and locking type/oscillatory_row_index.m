function [rows, score] = oscillatory_row_index(logit, lag, peak_window_s, trough_window_s, threshold)
%OSCILLATORY_ROW_INDEX  Neurons whose peak window exceeds the trough window.
%
% score is peak-window mean logit minus trough-window mean logit, in the
% input row order. rows lists the input rows whose score is finite and
% greater than threshold, strongest score first. Those rows are the
% oscillatory neurons.

score = logit_window_mean(logit, lag, peak_window_s) ...
    - logit_window_mean(logit, lag, trough_window_s);
if nargin < 5 || isempty(threshold) || ~isfinite(threshold)
    error('oscillatory_row_index:threshold', 'The difference threshold must be finite.');
end
finite = find(isfinite(score));
[~, order] = sort(score(finite), 'descend');
finite = finite(order);
rows = finite(score(finite) > threshold);
end
