function [component_table, density_frequency, component_density, ...
        mixture_density] = gaussian_mixture_density_curves( ...
        means, sds, weights, frequency_limits, effective_n, hard_count)
%GAUSSIAN_MIXTURE_DENSITY_CURVES  Mixture curves from stored parameters.
%
% means, sds, and weights describe one component count. The curves use
% the same 1,000-point grid as the estimate step.

means = double(means(:));
sds = double(sds(:));
weights = double(weights(:));
if numel(means) ~= numel(sds) || numel(means) ~= numel(weights)
    error('regionGaussianCenters:curves', ...
        'Means, standard deviations, and weights must match.');
end
if nargin < 5 || isempty(effective_n)
    effective_n = nan(size(means));
else
    effective_n = double(effective_n(:));
end
if nargin < 6 || isempty(hard_count)
    hard_count = nan(size(means));
else
    hard_count = double(hard_count(:));
end
[means, order] = sort(means);
sds = sds(order);
weights = weights(order);
if numel(effective_n) == numel(means)
    effective_n = effective_n(order);
end
if numel(hard_count) == numel(means)
    hard_count = hard_count(order);
end
n_components = numel(means);
component_table = table((1:n_components).', means, sds, weights, ...
    hard_count, effective_n, ...
    'VariableNames', {'Component', 'MeanHz', 'SDHz', 'Weight', ...
    'NAssigned', 'EffectiveN'});

density_frequency = linspace( ...
    frequency_limits(1), frequency_limits(2), 1000).';
component_density = zeros(numel(density_frequency), n_components);
for component_n = 1:n_components
    if ~isfinite(sds(component_n)) || sds(component_n) <= 0
        continue
    end
    component_density(:, component_n) = weights(component_n) * ...
        exp(-0.5 * ((density_frequency - means(component_n)) ./ ...
        sds(component_n)).^2) ./ ...
        (sds(component_n) * sqrt(2 * pi));
end
mixture_density = sum(component_density, 2);
end
