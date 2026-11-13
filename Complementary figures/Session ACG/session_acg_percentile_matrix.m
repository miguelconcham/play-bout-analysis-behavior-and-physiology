function packed = session_acg_percentile_matrix(records, group_name, curve_name, n_shuffles, sort_window_s)
%SESSION_ACG_PERCENTILE_MATRIX  Neurons by lags for one lock group and curve.
%
% group_name is trough, peak, delta, nonentrained, or all. Delta-locked
% is trough-locked and peak-locked together. Rows are sorted by the mean
% percentile inside sort_window_s, highest first, which is the top row of
% imagesc. Curves on a different lag axis are counted in n_dropped and
% left out. The axis with the most curves is the one that is drawn.

if nargin < 5 || isempty(sort_window_s)
    sort_window_s = [0.2 0.5];
end
group_name = session_acg_char(group_name);
curve_name = session_acg_char(curve_name);
known = {'trough', 'peak', 'delta', 'nonentrained', 'all'};
if ~any(strcmp(group_name, known))
    error('session_acg:percentile', ...
        'Lock group must be trough, peak, delta, nonentrained, or all.');
end
sort_window_s = double(sort_window_s(:))';
if numel(sort_window_s) ~= 2 || ~all(isfinite(sort_window_s)) || sort_window_s(2) < sort_window_s(1)
    error('session_acg:percentile', ...
        'sort_window_s must be [start end] in seconds.');
end
packed = empty_matrix();
if isempty(records)
    return
end
if ~isfield(records, 'curves') || ~isfield(records(1).curves, curve_name)
    error('session_acg:percentile', ...
        'Saved neurons have no %s percentile.', curve_name);
end

n = numel(records);
chosen_y = cell(n, 1);
chosen_lag = cell(n, 1);
n_chosen = 0;
for i = 1:n
    if ~neuron_in_group(records(i).lock_group, group_name)
        continue
    end
    curve = records(i).curves.(curve_name);
    y = curve.pctl;
    lag = curve.pctl_lag;
    if isempty(y) || numel(y) ~= numel(lag) || ~any(isfinite(y))
        continue
    end
    n_chosen = n_chosen + 1;
    chosen_y{n_chosen} = double(y(:))';
    chosen_lag{n_chosen} = double(lag(:))';
end
if n_chosen < 1
    return
end
chosen_y = chosen_y(1:n_chosen);
chosen_lag = chosen_lag(1:n_chosen);
[lag, keep] = majority_axis(chosen_lag);
stack = zeros(sum(keep), numel(lag));
row = 0;
for i = 1:n_chosen
    if ~keep(i)
        continue
    end
    row = row + 1;
    stack(row, :) = chosen_y{i};
end
packed.n_dropped = n_chosen - row;
[packed.percentile, packed.logit, packed.score, ~] = order_percentile_rows( ...
    stack, lag, sort_window_s, n_shuffles);
packed.lag = lag;
packed.n = size(packed.percentile, 1);
end


function packed = empty_matrix()
packed.lag = [];
packed.percentile = zeros(0, 0);
packed.logit = zeros(0, 0);
packed.score = zeros(0, 1);
packed.n = 0;
packed.n_dropped = 0;
end


function keep = neuron_in_group(lock_group, group_name)
lock_group = session_acg_char(lock_group);
if strcmp(group_name, 'all')
    keep = true;
elseif strcmp(group_name, 'delta')
    keep = any(strcmp(lock_group, {'trough', 'peak'}));
else
    keep = strcmp(lock_group, group_name);
end
end


function [lag, keep] = majority_axis(chosen_lag)
n = numel(chosen_lag);
keys = cell(n, 1);
for i = 1:n
    keys{i} = sprintf('%.8g ', chosen_lag{i});
end
[unique_keys, ~] = unique(keys, 'stable');
n_keys = numel(unique_keys);
counts = zeros(n_keys, 1);
index = zeros(n, 1);
for i = 1:n
    index(i) = find(strcmp(unique_keys, keys{i}), 1, 'first');
    counts(index(i)) = counts(index(i)) + 1;
end
[~, winner] = max(counts);
lag = chosen_lag{find(index == winner, 1, 'first')};
keep = index == winner;
end
