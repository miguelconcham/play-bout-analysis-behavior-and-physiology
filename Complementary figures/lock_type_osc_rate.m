function [percent, observed, total, label] = lock_type_osc_rate(neuron_type, osc_rows)
%LOCK_TYPE_OSC_RATE  Percent of each lock type that is oscillatory.
%
% percent is 100 * observed/total. observed counts oscillatory neurons of
% that type. total counts every neuron of that type. The types are trough,
% peak, unlocked (stored as nonentrained), and non-classified. A type with
% no neurons has a missing percent.

if nargin < 2
    osc_rows = zeros(0, 1);
end
osc_rows = osc_rows(:);
[total, ~, label] = lock_type_proportion(neuron_type);
if isempty(osc_rows)
    observed = zeros(1, 4);
else
    observed = lock_type_proportion(neuron_type(osc_rows));
end
percent = nan(1, 4);
has_cells = total > 0;
percent(has_cells) = 100 * observed(has_cells) ./ total(has_cells);
end
