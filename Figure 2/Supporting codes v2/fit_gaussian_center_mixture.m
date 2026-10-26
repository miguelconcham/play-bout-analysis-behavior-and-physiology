function result = fit_gaussian_center_mixture( ...
        centers, frequency_limits, options)
%FIT_GAUSSIAN_CENTER_MIXTURE  BIC or AIC mixture on bump centers.
%
% The fit matches the mid-LPAG Fourier-center mixture: fitgmdist on the
% finite center frequencies, a small covariance ridge, and several random
% restarts. A component stays eligible when it converges and its soft
% membership is at least minimum_effective_membership. selection_criterion
% is 'BIC' unless options asks for 'AIC'. Both criteria are stored.

if nargin < 3 || isempty(options)
    options = struct();
end
options = mixture_options(options);
centers = centers(:);
centers = centers(isfinite(centers));

if nargin < 2 || isempty(frequency_limits)
    if isempty(centers)
        frequency_limits = [0 1];
    else
        frequency_limits = [min(centers) max(centers)];
    end
end
frequency_limits = double(frequency_limits(:).');
if numel(frequency_limits) < 2 || ~all(isfinite(frequency_limits(1:2)))
    error('regionGaussianCenters:frequency', ...
        'frequency_limits must be [low high].');
end

result = empty_mixture_result(options.selection_criterion);
result.center_frequencies = centers;
result.frequency_limits = frequency_limits(1:2);

if numel(centers) < options.minimum_observations_per_component
    result.status = sprintf( ...
        'fewer than %d finite center frequencies', ...
        options.minimum_observations_per_component);
    return
end

max_components = min([options.max_components, numel(centers), ...
    numel(unique(centers)), floor(numel(centers) / ...
    options.minimum_observations_per_component)]);
if max_components < 1
    result.status = 'not enough distinct center frequencies';
    return
end

component_counts = (1:max_components).';
gmm_models = cell(max_components, 1);
gmm_aic = nan(max_components, 1);
gmm_bic = nan(max_components, 1);
gmm_negative_log_likelihood = nan(max_components, 1);
gmm_converged = false(max_components, 1);
gmm_minimum_effective_n = nan(max_components, 1);
gmm_eligible = false(max_components, 1);

previous_rng = rng;
rng_cleaner = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(options.random_seed, 'twister');
fit_options = statset('MaxIter', options.max_iter, 'Display', 'off');
for component_n = component_counts.'
    try
        this_gmm = fitgmdist(centers, component_n, ...
            'RegularizationValue', options.regularization, ...
            'Replicates', options.replicates, ...
            'Options', fit_options);
        gmm_models{component_n} = this_gmm;
        gmm_aic(component_n) = this_gmm.AIC;
        gmm_bic(component_n) = this_gmm.BIC;
        gmm_negative_log_likelihood(component_n) = ...
            this_gmm.NegativeLogLikelihood;
        gmm_converged(component_n) = this_gmm.Converged;
        this_posterior = posterior(this_gmm, centers);
        gmm_minimum_effective_n(component_n) = min(sum(this_posterior, 1));
        gmm_eligible(component_n) = this_gmm.Converged && ...
            gmm_minimum_effective_n(component_n) >= ...
            options.minimum_effective_membership && ...
            isfinite(this_gmm.AIC) && isfinite(this_gmm.BIC);
    catch ME
        warning('regionGaussianCenters:componentFitFailed', ...
            'The %d-component Gaussian-center mixture failed: %s', ...
            component_n, ME.message);
    end
end

result.model_selection = table(component_counts, gmm_aic, gmm_bic, ...
    gmm_negative_log_likelihood, gmm_converged, ...
    gmm_minimum_effective_n, gmm_eligible, ...
    'VariableNames', {'NComponents', 'AIC', 'BIC', ...
    'NegativeLogLikelihood', 'Converged', ...
    'MinimumEffectiveN', 'Eligible'});

if strcmp(options.selection_criterion, 'AIC')
    criterion_value = gmm_aic;
else
    criterion_value = gmm_bic;
end
best_component_count = select_gaussian_mixture_count( ...
    component_counts, criterion_value, gmm_eligible);
result.best_n_components = best_component_count;
if ~isfinite(best_component_count)
    result.status = 'no eligible model fits';
    return
end

best_gmm = gmm_models{best_component_count};
[component_table, density_frequency, component_density, mixture_density] = ...
    selected_mixture_curves(best_gmm, centers, frequency_limits(1:2));
result.status = 'fitted';
result.best_component_table = component_table;
result.density_frequency = density_frequency;
result.component_density = component_density;
result.mixture_density = mixture_density;
end

function options = mixture_options(options)
options = fill_mixture_option(options, 'selection_criterion', 'BIC');
options = fill_mixture_option(options, 'max_components', 10);
options = fill_mixture_option(options, 'replicates', 50);
options = fill_mixture_option(options, 'max_iter', 2000);
options = fill_mixture_option(options, 'regularization', 1e-4);
options = fill_mixture_option(options, ...
    'minimum_observations_per_component', 5);
options = fill_mixture_option(options, 'minimum_effective_membership', 3);
options = fill_mixture_option(options, 'random_seed', 1);

criterion = upper(strtrim(char(string(options.selection_criterion))));
if ~ismember(criterion, {'BIC', 'AIC'})
    error('regionGaussianCenters:criterion', ...
        'selection_criterion must be BIC or AIC.');
end
options.selection_criterion = criterion;
if options.max_components < 1 || options.replicates < 1 || ...
        options.max_iter < 1 || ...
        options.minimum_observations_per_component < 1 || ...
        options.minimum_effective_membership < 1
    error('regionGaussianCenters:options', ...
        ['max_components, replicates, max_iter, and the membership ' ...
        'minima must be at least 1.']);
end
end

function options = fill_mixture_option(options, name, default_value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = default_value;
end
end

function result = empty_mixture_result(criterion)
result = struct( ...
    'status', 'not fitted', ...
    'selection_criterion', criterion, ...
    'center_frequencies', zeros(0, 1), ...
    'frequency_limits', [NaN NaN], ...
    'model_selection', empty_model_selection(), ...
    'best_n_components', NaN, ...
    'best_component_table', empty_component_table(), ...
    'density_frequency', zeros(0, 1), ...
    'component_density', zeros(0, 0), ...
    'mixture_density', zeros(0, 1));
end

function model_selection = empty_model_selection()
model_selection = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), false(0, 1), zeros(0, 1), false(0, 1), ...
    'VariableNames', {'NComponents', 'AIC', 'BIC', ...
    'NegativeLogLikelihood', 'Converged', ...
    'MinimumEffectiveN', 'Eligible'});
end

function component_table = empty_component_table()
component_table = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'Component', 'MeanHz', 'SDHz', 'Weight', ...
    'NAssigned', 'EffectiveN'});
end

function [component_table, density_frequency, component_density, ...
        mixture_density] = selected_mixture_curves( ...
        this_gmm, centers, frequency_limits)
component_mean = this_gmm.mu(:);
component_sd = sqrt(squeeze(this_gmm.Sigma));
component_sd = component_sd(:);
component_weight = this_gmm.ComponentProportion(:);
[component_mean, component_order] = sort(component_mean);
component_sd = component_sd(component_order);
component_weight = component_weight(component_order);
n_components = numel(component_mean);

raw_assignment = cluster(this_gmm, centers);
sorted_index = zeros(n_components, 1);
sorted_index(component_order) = 1:n_components;
assignment = sorted_index(raw_assignment);
posterior_weight = posterior(this_gmm, centers);
posterior_weight = posterior_weight(:, component_order);
n_assigned = accumarray(assignment, 1, [n_components, 1]);
effective_n = sum(posterior_weight, 1).';

component_table = table((1:n_components).', component_mean, ...
    component_sd, component_weight, n_assigned, effective_n, ...
    'VariableNames', {'Component', 'MeanHz', 'SDHz', 'Weight', ...
    'NAssigned', 'EffectiveN'});

density_frequency = linspace( ...
    frequency_limits(1), frequency_limits(2), 1000).';
component_density = zeros(numel(density_frequency), n_components);
for component_n = 1:n_components
    if ~isfinite(component_sd(component_n)) || ...
            component_sd(component_n) <= 0
        continue
    end
    component_density(:, component_n) = ...
        component_weight(component_n) * exp(-0.5 * ( ...
        (density_frequency - component_mean(component_n)) ./ ...
        component_sd(component_n)).^2) ./ ...
        (component_sd(component_n) * sqrt(2 * pi));
end
mixture_density = sum(component_density, 2);
end
