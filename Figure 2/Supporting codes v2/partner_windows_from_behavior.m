function windows = partner_windows_from_behavior(behavior)
%PARTNER_WINDOWS_FROM_BEHAVIOR  Partner-session intervals, earliest first.
%
%   Partner k is the k-th Session_structure row of type Partners session.
%   When that type is absent, non-tickling Session_structure rows are used.

windows = zeros(0, 2);
[animal, t0, t1, type] = behavior_columns(behavior);
if isempty(animal) || isempty(t0)
    return
end

animal = as_labels(animal);
is_block = strcmp(animal, 'Session_structure');
if ~isempty(type)
    type = as_labels(type);
    is_partner = is_block & strcmp(type, 'Partners session');
    if ~any(is_partner)
        is_partner = is_block & ~strcmp(type, 'Tickling');
    end
else
    is_partner = is_block;
end
if ~any(is_partner)
    return
end

windows = [t0(is_partner), t1(is_partner)];
windows = sortrows(windows, 1);
end

function [animal, t0, t1, type] = behavior_columns(behavior)
animal = {};
t0 = [];
t1 = [];
type = {};
if isempty(behavior)
    return
end

if is_table(behavior)
    animal = behavior.Animal;
    t0 = behavior.Start;
    t1 = behavior.End;
    names = behavior.Properties.VariableNames;
    if any(strcmp(names, 'Type'))
        type = behavior.Type;
    end
elseif isstruct(behavior) && isfield(behavior, 'Animal')
    if isscalar(behavior)
        animal = behavior.Animal;
        t0 = behavior.Start;
        t1 = behavior.End;
        if isfield(behavior, 'Type')
            type = behavior.Type;
        end
    else
        animal = {behavior.Animal};
        animal = animal(:);
        t0 = [behavior.Start]';
        t1 = [behavior.End]';
        if isfield(behavior, 'Type')
            type = {behavior.Type};
            type = type(:);
        end
    end
end

t0 = as_double_column(t0);
t1 = as_double_column(t1);
end

function tf = is_table(behavior)
tf = false;
if exist('istable', 'builtin') == 5 || exist('istable', 'file') == 2
    tf = istable(behavior);
end
end

function labels = as_labels(value)
if isempty(value)
    labels = {};
    return
end
if ischar(value)
    labels = cellstr(value);
elseif iscell(value)
    labels = cell(numel(value), 1);
    for k = 1:numel(value)
        item = value{k};
        if ischar(item)
            labels{k} = item;
        else
            labels{k} = char(item);
        end
    end
elseif exist('isstring', 'builtin') == 5 || exist('isstring', 'file') == 2
    if isstring(value)
        labels = cellstr(value);
    else
        labels = cellstr(num2str(value));
    end
else
    labels = cellstr(num2str(value));
end
labels = labels(:);
end

function x = as_double_column(x)
if isempty(x)
    x = [];
    return
end
if iscell(x)
    x = cell2mat(x(:));
end
x = double(x(:));
end
