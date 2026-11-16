function panels = significant_paired_bars(high_values, low_values, titles, families, holm_p)
%SIGNIFICANT_PAIRED_BARS  Significant high-versus-low features, z-scored.
%
% high_values and low_values have one column per feature and one row per
% neuron. A feature is kept when holm_p is below 0.05. Only neurons with
% both power values enter that feature. Panel 1 holds the oscillatory
% features and panel 2 the non-oscillatory features. Each kept column is
% z-scored across those complete pairs. The stored means are the means of
% those z-scores.

titles = titles(:).';
families = families(:).';
holm_p = holm_p(:);
n_features = numel(titles);
if size(high_values, 2) ~= n_features || size(low_values, 2) ~= n_features ...
        || size(high_values, 1) ~= size(low_values, 1) ...
        || numel(families) ~= n_features || numel(holm_p) ~= n_features
    error('significant_paired_bars:size', ...
        'Each feature needs one high column, one low column, one title, one family, and one Holm-Sidak p.');
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
        pair = isfinite(high_values(:, j)) & isfinite(low_values(:, j));
        y_high = high_values(pair, j);
        y_low = low_values(pair, j);
        [z_high, z_low, mean_high, mean_low] = zscore_paired_groups(y_high, y_low);
        panels(g).z1{k} = z_high;
        panels(g).z2{k} = z_low;
        panels(g).mean1(k) = mean_high;
        panels(g).mean2(k) = mean_low;
        panels(g).stars{k} = paired_star_label(holm_p(j));
    end
end
end


function [z1, z2, mean1, mean2] = zscore_paired_groups(y1, y2)
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
mean1 = paired_group_mean(z1);
mean2 = paired_group_mean(z2);
end


function value = paired_group_mean(z)
if isempty(z)
    value = nan;
else
    value = mean(z);
end
end


function label = paired_star_label(p)
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
