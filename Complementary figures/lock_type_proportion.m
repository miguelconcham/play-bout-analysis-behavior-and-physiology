function [count, proportion, label] = lock_type_proportion(neuron_type)
%LOCK_TYPE_PROPORTION  Trough, peak, unlocked, and non-classified shares.
%
% neuron_type holds the checkpoint labels. 'nonentrained' is unlocked.
% An empty label, or any other label, is non-classified. proportion sums
% to 1 when neuron_type is not empty.

label = {'Trough', 'Peak', 'Unlocked', 'Non-classified'};
stored = {'trough', 'peak', 'nonentrained'};
n = numel(neuron_type);
count = zeros(1, 4);
for i = 1:n
    value = neuron_type{i};
    if iscell(value)
        if isempty(value)
            value = '';
        else
            value = value{1};
        end
    end
    name = lower(strtrim(char(value(:).')));
    hit = find(strcmp(stored, name), 1);
    if isempty(hit)
        count(4) = count(4) + 1;
    else
        count(hit) = count(hit) + 1;
    end
end
if n < 1
    proportion = nan(1, 4);
else
    proportion = count / n;
end
end
