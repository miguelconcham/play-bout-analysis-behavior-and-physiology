function [assignment, max_posterior] = assign_bump_centers_to_mixture( ...
        centers, component_table, minimum_posterior)
%ASSIGN_BUMP_CENTERS_TO_MIXTURE  Map each bump center to a mixture component.
%
% component_table is the selected mixture, with MeanHz, SDHz, and Weight.
% The posterior of a center is the normalized Gaussian likelihood of those
% components. assignment is the component row when that posterior is at
% least minimum_posterior and higher than every other component. It is 0
% when the center is left unassigned. max_posterior is the largest
% posterior, including for unassigned centers.

if nargin < 3 || isempty(minimum_posterior)
    minimum_posterior = 0.5;
end
if ~isscalar(minimum_posterior) || ~isfinite(minimum_posterior) || ...
        minimum_posterior <= 0 || minimum_posterior > 1
    error('regionGaussianCenters:posterior', ...
        'minimum_assignment_posterior must lie in (0, 1].');
end

centers = double(centers(:));
n_centers = numel(centers);
assignment = zeros(n_centers, 1);
max_posterior = zeros(n_centers, 1);
if n_centers == 0 || ~istable(component_table) || ...
        height(component_table) == 0
    return
end
required = {'MeanHz', 'SDHz', 'Weight'};
if ~all(ismember(required, component_table.Properties.VariableNames))
    error('regionGaussianCenters:posterior', ...
        'The mixture table needs MeanHz, SDHz, and Weight.');
end

n_components = height(component_table);
log_density = -inf(n_centers, n_components);
for component_n = 1:n_components
    component_sd = component_table.SDHz(component_n);
    component_weight = component_table.Weight(component_n);
    if ~isfinite(component_sd) || component_sd <= 0 || ...
            ~isfinite(component_weight) || component_weight <= 0
        continue
    end
    z = (centers - component_table.MeanHz(component_n)) ./ component_sd;
    log_density(:, component_n) = log(component_weight) - ...
        0.5 * z.^2 - log(component_sd) - 0.5 * log(2 * pi);
end

posterior = zeros(n_centers, n_components);
finite_row = any(isfinite(log_density), 2) & isfinite(centers);
if any(finite_row)
    shifted = log_density(finite_row, :);
    shifted(~isfinite(shifted)) = -inf;
    shifted = shifted - max(shifted, [], 2);
    density = exp(shifted);
    posterior(finite_row, :) = density ./ sum(density, 2);
end

ordered = sort(posterior, 2, 'descend');
max_posterior = ordered(:, 1);
if n_components >= 2
    second_posterior = ordered(:, 2);
else
    second_posterior = zeros(n_centers, 1);
end
[~, assignment] = max(posterior, [], 2);
unassigned = ~isfinite(centers) | max_posterior < minimum_posterior | ...
    (max_posterior - second_posterior) <= 1e-10;
assignment(unassigned) = 0;
end
