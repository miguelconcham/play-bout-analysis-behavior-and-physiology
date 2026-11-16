function panels = significant_group_bars(values, group_index, titles, families, holm_p)
%SIGNIFICANT_GROUP_BARS  Significant unpaired features, z-scored for one axis.
%
% values has one column per feature and one row per neuron that entered
% the split. group_index is 1 for the oscillatory cluster and 2 for the
% other cluster. A feature is kept when holm_p is below 0.05. Panel 1
% holds the oscillatory features and panel 2 the non-oscillatory features.
% Each kept column is z-scored across the finite values of both groups.
% The stored means are the means of those z-scores.

titles = titles(:).';
families = families(:).';
holm_p = holm_p(:);
n_features = numel(titles);
if size(values, 2) ~= n_features || numel(families) ~= n_features ...
        || numel(holm_p) ~= n_features
    error('significant_group_bars:size', ...
        'Each feature needs one column, one title, one family, and one Holm-Sidak p.');
end
if numel(group_index) ~= size(values, 1)
    error('significant_group_bars:size', ...
        'group_index needs one entry per neuron.');
end

family_names = {'oscillatory', 'non-oscillatory'};
panels = struct([]);
for g = 1:numel(family_names)
    keep = [];
    for j = 1:n_features
        if strcmp(families{j}, family_names{g}) && isfinite(holm_p(j)) && holm_p(j) < 0.05
            keep(end + 1) = j; %#ok<AGROW>
        end
    end
    panels(g).family = family_names{g};
    panels(g).titles = titles(keep);
    panels(g).holm = holm_p(keep);
    panels(g).stars = cell(1, numel(keep));
    panels(g).z1 = cell(1, numel(keep));
    panels(g).z2 = cell(1, numel(keep));
    panels(g).mean1 = nan(1, numel(keep));
    panels(g).mean2 = nan(1, numel(keep));
    for k = 1:numel(keep)
        j = keep(k);
        y1 = values(group_index(:) == 1, j);
        y2 = values(group_index(:) == 2, j);
        y1 = y1(isfinite(y1));
        y2 = y2(isfinite(y2));
        [z1, z2, mean1, mean2] = zscore_two_groups(y1, y2);
        panels(g).z1{k} = z1;
        panels(g).z2{k} = z2;
        panels(g).mean1(k) = mean1;
        panels(g).mean2(k) = mean2;
        panels(g).stars{k} = holm_star_label(holm_p(j));
    end
end
end


function [z1, z2, mean1, mean2] = zscore_two_groups(y1, y2)
pooled = [y1(:); y2(:)];
if isempty(pooled)
    z1 = zeros(0, 1);
    z2 = zeros(0, 1);
    mean1 = nan;
    mean2 = nan;
    return
end
center = mean(pooled);
scale = std(pooled, 0);
if ~isfinite(scale) || scale == 0
    z1 = zeros(numel(y1), 1);
    z2 = zeros(numel(y2), 1);
else
    z1 = (y1(:) - center) / scale;
    z2 = (y2(:) - center) / scale;
end
mean1 = group_mean(z1);
mean2 = group_mean(z2);
end


function value = group_mean(z)
if isempty(z)
    value = nan;
else
    value = mean(z);
end
end


function label = holm_star_label(p)
if ~isfinite(p) || p >= 0.05
    label = '';
elseif p < 0.001
    label = '***';
elseif p < 0.01
    label = '**';
else
    label = '*';
end
end
