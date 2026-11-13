function [lag, stack] = percentile_curve_stack(records, split_name, tail_name, field_name)
%PERCENTILE_CURVE_STACK  One row per neuron for a shuffle-percentile curve.
%
% split_name 'all' reads the all-runs curve. Otherwise the curve is
% records{i}.(split_name).(tail_name).(field_name). Rows with no finite
% percentile are left out. The remaining rows are sorted by the lag of
% their highest percentile, and a higher peak comes first when that lag
% is shared, so an imagesc reads from early peaks at the top.

lag = [];
stack = zeros(0, 0);
if isempty(records)
    return
end
rows = {};
for i = 1:numel(records)
    tail = percentile_tail(records{i}, split_name, tail_name);
    if isempty(tail) || ~isfield(tail, field_name) || isempty(tail.(field_name))
        continue
    end
    y = double(tail.(field_name)(:).');
    if ~any(isfinite(y)) || ~isfield(tail, 'pctl_lag') || numel(tail.pctl_lag) ~= numel(y)
        continue
    end
    this_lag = double(tail.pctl_lag(:).');
    if isempty(lag)
        lag = this_lag;
    elseif numel(this_lag) ~= numel(lag) || max(abs(this_lag - lag)) > 1e-9
        continue
    end
    rows{end + 1, 1} = y; %#ok<AGROW>
end
if isempty(rows)
    lag = [];
    return
end
stack = cat(1, rows{:});
[peak_val, peak_bin] = row_peak(stack);
[~, order] = sortrows([peak_bin, -peak_val]);
stack = stack(order, :);
end


function [peak_val, peak_bin] = row_peak(stack)
n = size(stack, 1);
peak_val = nan(n, 1);
peak_bin = ones(n, 1);
for i = 1:n
    finite = find(isfinite(stack(i, :)));
    if isempty(finite)
        continue
    end
    [peak_val(i), local] = max(stack(i, finite));
    peak_bin(i) = finite(local);
end
end


function tail = percentile_tail(record, split_name, tail_name)
tail = [];
if strcmp(split_name, 'all')
    if isfield(record, 'all')
        tail = record.all;
    end
    return
end
if ~isfield(record, split_name)
    return
end
block = record.(split_name);
if isfield(block, tail_name)
    tail = block.(tail_name);
end
end
