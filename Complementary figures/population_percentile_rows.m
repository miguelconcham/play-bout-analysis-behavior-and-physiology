function [lag, stack, kept] = population_percentile_rows(records, split_name, tail_name, field_name)
%POPULATION_PERCENTILE_ROWS  Percentile curves in record order.
%
% split_name 'all' reads the all-runs curve. Otherwise the curve is
% records{i}.(split_name).(tail_name).(field_name). kept is the record
% index of each row. Rows stay in record order. A row with no finite
% percentile, or a lag axis that differs from the first kept row, is
% left out. This does not sort, so several calls can be aligned by kept.

lag = [];
stack = zeros(0, 0);
kept = zeros(0, 1);
if isempty(records)
    return
end
rows = {};
for i = 1:numel(records)
    tail = population_percentile_tail(records{i}, split_name, tail_name);
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
    kept(end + 1, 1) = i; %#ok<AGROW>
end
if isempty(rows)
    lag = [];
    return
end
stack = cat(1, rows{:});
end


function tail = population_percentile_tail(record, split_name, tail_name)
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
