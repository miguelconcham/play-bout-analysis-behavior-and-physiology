function [lag, percentile, neuron_type, record_index] = align_population_percentiles(records)
%ALIGN_POPULATION_PERCENTILES  One row per neuron across the population curves.
%
% percentile.sum and percentile.ratio each contain all, power_high, and
% power_low. A neuron is kept only when all six curves share one lag axis.
% Rows are in record order. neuron_type is the checkpoint lock group:
% 'trough', 'peak', 'nonentrained', or '' when the cell was not classified.
% record_index points back into records.

column_id = {'all', 'power_high', 'power_low'};
split_name = {'all', 'power', 'power'};
tail_name = {'', 'high', 'low'};
scale_id = {'sum', 'ratio'};
field_name = {'sum_pctl', 'ratio_pctl'};

stack = struct();
kept = struct();
lag = [];
for s = 1:numel(scale_id)
    for c = 1:numel(column_id)
        [this_lag, this_stack, this_kept] = population_percentile_rows( ...
            records, split_name{c}, tail_name{c}, field_name{s});
        if ~isempty(this_lag)
            if isempty(lag)
                lag = this_lag;
            elseif numel(this_lag) ~= numel(lag) || max(abs(this_lag - lag)) > 1e-9
                error('align_population_percentiles:axis', ...
                    'The %s %s curve is not on the same lag axis.', ...
                    scale_id{s}, column_id{c});
            end
        end
        stack.(scale_id{s}).(column_id{c}) = this_stack;
        kept.(scale_id{s}).(column_id{c}) = this_kept;
    end
end

record_index = kept.sum.all;
for s = 1:numel(scale_id)
    for c = 1:numel(column_id)
        record_index = intersect(record_index, kept.(scale_id{s}).(column_id{c}), 'stable');
    end
end

percentile = struct();
n_lag = numel(lag);
if isempty(record_index)
    lag = [];
    n_lag = 0;
end
for s = 1:numel(scale_id)
    for c = 1:numel(column_id)
        [found, loc] = ismember(record_index, kept.(scale_id{s}).(column_id{c}));
        if ~all(found)
            error('align_population_percentiles:index', ...
                'A shared neuron was missing from %s %s.', scale_id{s}, column_id{c});
        end
        if isempty(loc)
            percentile.(scale_id{s}).(column_id{c}) = zeros(0, n_lag);
        else
            percentile.(scale_id{s}).(column_id{c}) = stack.(scale_id{s}).(column_id{c})(loc, :);
        end
    end
end

neuron_type = cell(numel(record_index), 1);
for i = 1:numel(record_index)
    neuron_type{i} = lock_group_label(records{record_index(i)});
end
end


function label = lock_group_label(record)
label = '';
if ~isfield(record, 'lock_group') || isempty(record.lock_group)
    return
end
value = record.lock_group;
if iscell(value)
    if isempty(value)
        return
    end
    value = value{1};
end
label = strtrim(char(value(:).'));
end
