function [scalars, components] = gaussian_mixture_measurements( ...
        means, sds, weights, centers, posterior_weight, fit_info)
%GAUSSIAN_MIXTURE_MEASUREMENTS  Summaries of one fitted center mixture.
%
% The extra summaries are computed from the fitted parameters and the
% posterior that eligibility already requires. No additional EM search is
% run. Goodness-of-fit uses every center for the 1 Hz histogram and at
% most 20,000 sorted centers for the Kolmogorov-Smirnov and
% Anderson-Darling statistics. The full posterior matrix is not stored.
%
% No inputs return the empty tables. A fit_info.failed row records the
% component count without measurements.

if nargin == 0
    scalars = empty_scalar_table();
    components = empty_component_parameter_table();
    return
end
if nargin < 6 || isempty(fit_info)
    fit_info = struct();
end
if isfield(fit_info, 'failed') && fit_info.failed
    scalars = failed_scalar_table(fit_info);
    components = empty_component_parameter_table();
    return
end

means = double(means(:));
sds = double(sds(:));
weights = double(weights(:));
centers = double(centers(:));
if numel(means) ~= numel(sds) || numel(means) ~= numel(weights)
    error('regionGaussianCenters:measurements', ...
        'Means, standard deviations, and weights must match.');
end
n_components = numel(means);
n_observations = numel(centers);
if size(posterior_weight, 1) ~= n_observations || ...
        size(posterior_weight, 2) ~= n_components
    error('regionGaussianCenters:measurements', ...
        'The posterior must be observations by components.');
end
posterior_weight = double(posterior_weight);

aic = info_number(fit_info, 'aic', NaN);
bic = info_number(fit_info, 'bic', NaN);
negative_log_likelihood = info_number( ...
    fit_info, 'negative_log_likelihood', NaN);
converged = logical(info_number(fit_info, 'converged', false));
num_iterations = info_number(fit_info, 'num_iterations', NaN);
replicates = info_number(fit_info, 'replicates', NaN);
regularization = info_number(fit_info, 'regularization', NaN);
minimum_membership = info_number( ...
    fit_info, 'minimum_effective_membership', 3);
shared_covariance = info_number(fit_info, 'shared_covariance', NaN);

effective_n = sum(posterior_weight, 1).';
[max_posterior, map_index] = max(posterior_weight, [], 2);
hard_count = accumarray(map_index, 1, [n_components, 1]);
positive = posterior_weight > 0;
entropy = -sum(posterior_weight(positive) .* log(posterior_weight(positive)));
if n_observations > 0
    mean_entropy = entropy / n_observations;
else
    mean_entropy = NaN;
end

log_component = component_log_density(centers, means, sds, weights);
log_density = logsumexp_rows(log_component);
map_linear = sub2ind([n_observations, n_components], ...
    (1:n_observations).', map_index);
hard_completed = sum(log_component(map_linear));
if ~isfinite(hard_completed)
    hard_completed = NaN;
end
if isfinite(negative_log_likelihood)
    completed_log_likelihood = -negative_log_likelihood - entropy;
else
    finite_density = log_density(isfinite(log_density));
    completed_log_likelihood = sum(finite_density) - entropy;
    negative_log_likelihood = -sum(finite_density);
end

n_parameters = parameter_count(aic, negative_log_likelihood, n_components);
if isfinite(aic) && isfinite(bic)
    eligible = converged && min(effective_n) >= minimum_membership;
else
    eligible = false;
end

finite_density = log_density(isfinite(log_density));
if numel(finite_density) >= 2
    pointwise_sd = std(finite_density);
    pointwise_p05 = sample_percentile(finite_density, 0.05);
    pointwise_p50 = sample_percentile(finite_density, 0.50);
    pointwise_p95 = sample_percentile(finite_density, 0.95);
else
    pointwise_sd = NaN;
    pointwise_p05 = NaN;
    pointwise_p50 = NaN;
    pointwise_p95 = NaN;
end

[order, geometry] = component_geometry(means, sds, weights);
effective_n = effective_n(order);
hard_count = hard_count(order);
means = means(order);
sds = sds(order);
weights = weights(order);

try
    gof = goodness_of_fit(centers, means, sds, weights);
catch
    gof = empty_goodness();
end

row = measurement_template();
row.NComponents = n_components;
row.NObservations = n_observations;
row.NParameters = n_parameters;
row.NegativeLogLikelihood = negative_log_likelihood;
row.LogLikelihood = -negative_log_likelihood;
row.AIC = aic;
row.BIC = bic;
row.Converged = converged;
row.NumIterations = num_iterations;
row.Replicates = replicates;
row.RegularizationValue = regularization;
row.SharedCovariance = shared_covariance;
row.Entropy = entropy;
row.MeanEntropy = mean_entropy;
row.CompletedLogLikelihood = completed_log_likelihood;
row.HardCompletedLogLikelihood = hard_completed;
row.MeanMaxPosterior = mean(max_posterior);
row.MedianMaxPosterior = median(max_posterior);
row.FractionMaxPosteriorAbove08 = mean(max_posterior >= 0.8);
row.FractionMaxPosteriorAbove095 = mean(max_posterior >= 0.95);
row.MinimumEffectiveN = min(effective_n);
row.MinimumHardCount = min(hard_count);
row.NComponentsBelowMinimum = sum(effective_n < minimum_membership);
row.MinimumWeight = min(weights);
row.MaximumWeight = max(weights);
row.WeightEntropy = weight_entropy(weights);
row.MinimumSDHz = min(sds);
row.MaximumSDHz = max(sds);
row.MeanSDHz = mean(sds);
row.MaximumComponentHeight = geometry.maximum_height;
row.MinimumNeighborSeparationHz = geometry.minimum_separation;
row.MinimumNeighborAshmanD = geometry.minimum_ashman_d;
row.MaximumPairBhattacharyya = geometry.maximum_bhattacharyya;
row.MeanNeighborBhattacharyya = geometry.mean_neighbor_bhattacharyya;
row.KSDistance = gof.ks_distance;
row.KSPoints = gof.ks_points;
row.AndersonDarling = gof.anderson_darling;
row.ChiSquare = gof.chi_square;
row.ChiSquareBins = gof.chi_square_bins;
row.HistogramL1 = gof.histogram_l1;
row.PointwiseLogLikSD = pointwise_sd;
row.PointwiseLogLikP05 = pointwise_p05;
row.PointwiseLogLikP50 = pointwise_p50;
row.PointwiseLogLikP95 = pointwise_p95;
row.Eligible = eligible;
scalars = struct2table(row, 'AsArray', true);

component_rows = repmat(component_template(), n_components, 1);
for component_n = 1:n_components
    component_rows(component_n).NComponents = n_components;
    component_rows(component_n).Component = component_n;
    component_rows(component_n).MeanHz = means(component_n);
    component_rows(component_n).SDHz = sds(component_n);
    component_rows(component_n).Weight = weights(component_n);
    component_rows(component_n).EffectiveN = effective_n(component_n);
    component_rows(component_n).HardCount = hard_count(component_n);
end
components = struct2table(component_rows);
end

function value = info_number(fit_info, name, default_value)
value = default_value;
if isstruct(fit_info) && isfield(fit_info, name) && ...
        ~isempty(fit_info.(name)) && isfinite(fit_info.(name)(1))
    value = double(fit_info.(name)(1));
end
end

function count = parameter_count(aic, negative_log_likelihood, n_components)
if isfinite(aic) && isfinite(negative_log_likelihood)
    count = (aic - 2 * negative_log_likelihood) / 2;
    if abs(count - round(count)) < 1e-4
        count = round(count);
    end
    return
end
count = 3 * n_components - 1;
end

function entropy = weight_entropy(weights)
weights = weights(weights > 0);
if isempty(weights)
    entropy = NaN;
    return
end
entropy = -sum(weights .* log(weights));
end

function [order, geometry] = component_geometry(means, sds, weights)
[~, order] = sort(means);
means = means(order);
sds = sds(order);
weights = weights(order);
height = weights ./ (sds * sqrt(2 * pi));
height = height(isfinite(height) & sds > 0);
if isempty(height)
    geometry.maximum_height = NaN;
else
    geometry.maximum_height = max(height);
end
n_components = numel(means);
if n_components < 2
    geometry.minimum_separation = NaN;
    geometry.minimum_ashman_d = NaN;
    geometry.maximum_bhattacharyya = NaN;
    geometry.mean_neighbor_bhattacharyya = NaN;
    return
end
separation = diff(means);
geometry.minimum_separation = min(separation);
pair_variance = sds(1:end - 1).^2 + sds(2:end).^2;
ashman = abs(separation) ./ sqrt(pair_variance / 2);
ashman = ashman(isfinite(ashman));
if isempty(ashman)
    geometry.minimum_ashman_d = NaN;
else
    geometry.minimum_ashman_d = min(ashman);
end
neighbor_bc = bhattacharyya_coeff( ...
    means(1:end - 1), sds(1:end - 1), means(2:end), sds(2:end));
neighbor_bc = neighbor_bc(isfinite(neighbor_bc));
if isempty(neighbor_bc)
    geometry.mean_neighbor_bhattacharyya = NaN;
else
    geometry.mean_neighbor_bhattacharyya = mean(neighbor_bc);
end
pair_bc = [];
for first_n = 1:n_components - 1
    for second_n = first_n + 1:n_components
        pair_bc(end + 1) = bhattacharyya_coeff( ... %#ok<AGROW>
            means(first_n), sds(first_n), ...
            means(second_n), sds(second_n));
    end
end
pair_bc = pair_bc(isfinite(pair_bc));
if isempty(pair_bc)
    geometry.maximum_bhattacharyya = NaN;
else
    geometry.maximum_bhattacharyya = max(pair_bc);
end
end

function coeff = bhattacharyya_coeff(mean_one, sd_one, mean_two, sd_two)
pair_variance = sd_one.^2 + sd_two.^2;
coeff = nan(size(pair_variance));
valid = pair_variance > 0 & sd_one > 0 & sd_two > 0;
coeff(valid) = sqrt(2 .* sd_one(valid) .* sd_two(valid) ./ ...
    pair_variance(valid)) .* exp(-0.25 * ...
    (mean_one(valid) - mean_two(valid)).^2 ./ pair_variance(valid));
end

function gof = goodness_of_fit(centers, means, sds, weights)
gof = empty_goodness();
n_observations = numel(centers);
if n_observations < 2
    return
end
bin_width = 1;
lower_edge = floor(min(centers));
upper_edge = ceil(max(centers));
if upper_edge <= lower_edge
    upper_edge = lower_edge + bin_width;
end
edges = lower_edge:bin_width:upper_edge;
observed = histcounts(centers, edges).';
expected = zeros(numel(observed), 1);
for bin_n = 1:numel(observed)
    expected(bin_n) = n_observations * ( ...
        univariate_mixture_cdf(edges(bin_n + 1), means, sds, weights) - ...
        univariate_mixture_cdf(edges(bin_n), means, sds, weights));
end
expected = max(expected, 0);
usable = expected >= 5;
gof.chi_square_bins = sum(usable);
if any(usable)
    gof.chi_square = sum((observed(usable) - expected(usable)).^2 ./ ...
        expected(usable));
end
gof.histogram_l1 = sum(abs(observed / n_observations - ...
    expected / n_observations));

sorted_centers = sort(centers);
point_cap = 20000;
stride = max(1, ceil(n_observations / point_cap));
sample = sorted_centers(1:stride:n_observations);
cdf_value = univariate_mixture_cdf(sample, means, sds, weights);
cdf_value = min(max(cdf_value, 1e-12), 1 - 1e-12);
n_points = numel(sample);
empirical_high = (1:n_points).' / n_points;
empirical_low = (0:n_points - 1).' / n_points;
gof.ks_distance = max([ ...
    max(abs(empirical_high - cdf_value)), ...
    max(abs(empirical_low - cdf_value))]);
gof.ks_points = n_points;
index = (1:n_points).';
gof.anderson_darling = -n_points - sum((2 * index - 1) ./ n_points .* ( ...
    log(cdf_value) + log(1 - cdf_value(n_points + 1 - index))));
end

function gof = empty_goodness()
gof = struct( ...
    'ks_distance', NaN, ...
    'ks_points', 0, ...
    'anderson_darling', NaN, ...
    'chi_square', NaN, ...
    'chi_square_bins', 0, ...
    'histogram_l1', NaN);
end

function probability = univariate_mixture_cdf(x, means, sds, weights)
probability = zeros(size(x));
for component_n = 1:numel(means)
    if ~(weights(component_n) > 0) || ~(sds(component_n) > 0) || ...
            ~isfinite(means(component_n))
        continue
    end
    probability = probability + weights(component_n) * normcdf( ...
        x, means(component_n), sds(component_n));
end
end

function log_component = component_log_density(x, means, sds, weights)
n_observations = numel(x);
n_components = numel(means);
log_component = -inf(n_observations, n_components);
for component_n = 1:n_components
    if ~(weights(component_n) > 0) || ~(sds(component_n) > 0) || ...
            ~isfinite(means(component_n))
        continue
    end
    z = (x - means(component_n)) ./ sds(component_n);
    log_component(:, component_n) = log(weights(component_n)) - ...
        log(sds(component_n)) - 0.5 * log(2 * pi) - 0.5 * z.^2;
end
end

function total = logsumexp_rows(log_component)
peak = max(log_component, [], 2);
safe_peak = peak;
safe_peak(~isfinite(safe_peak)) = 0;
total = safe_peak + log(sum(exp(log_component - safe_peak), 2));
total(all(~isfinite(log_component), 2)) = -inf;
end

function value = sample_percentile(x, probability)
x = sort(x(:));
if isempty(x)
    value = NaN;
    return
end
rank = 1 + (numel(x) - 1) * probability;
low = floor(rank);
high = ceil(rank);
if low < 1
    low = 1;
end
if high > numel(x)
    high = numel(x);
end
if low == high
    value = x(low);
else
    value = x(low) * (high - rank) + x(high) * (rank - low);
end
end

function row = measurement_template()
row = struct( ...
    'NComponents', NaN, ...
    'NObservations', NaN, ...
    'NParameters', NaN, ...
    'NegativeLogLikelihood', NaN, ...
    'LogLikelihood', NaN, ...
    'AIC', NaN, ...
    'BIC', NaN, ...
    'Converged', false, ...
    'NumIterations', NaN, ...
    'Replicates', NaN, ...
    'RegularizationValue', NaN, ...
    'SharedCovariance', NaN, ...
    'Entropy', NaN, ...
    'MeanEntropy', NaN, ...
    'CompletedLogLikelihood', NaN, ...
    'HardCompletedLogLikelihood', NaN, ...
    'MeanMaxPosterior', NaN, ...
    'MedianMaxPosterior', NaN, ...
    'FractionMaxPosteriorAbove08', NaN, ...
    'FractionMaxPosteriorAbove095', NaN, ...
    'MinimumEffectiveN', NaN, ...
    'MinimumHardCount', NaN, ...
    'NComponentsBelowMinimum', NaN, ...
    'MinimumWeight', NaN, ...
    'MaximumWeight', NaN, ...
    'WeightEntropy', NaN, ...
    'MinimumSDHz', NaN, ...
    'MaximumSDHz', NaN, ...
    'MeanSDHz', NaN, ...
    'MaximumComponentHeight', NaN, ...
    'MinimumNeighborSeparationHz', NaN, ...
    'MinimumNeighborAshmanD', NaN, ...
    'MaximumPairBhattacharyya', NaN, ...
    'MeanNeighborBhattacharyya', NaN, ...
    'KSDistance', NaN, ...
    'KSPoints', 0, ...
    'AndersonDarling', NaN, ...
    'ChiSquare', NaN, ...
    'ChiSquareBins', 0, ...
    'HistogramL1', NaN, ...
    'PointwiseLogLikSD', NaN, ...
    'PointwiseLogLikP05', NaN, ...
    'PointwiseLogLikP50', NaN, ...
    'PointwiseLogLikP95', NaN, ...
    'Eligible', false);
end

function row = component_template()
row = struct( ...
    'NComponents', NaN, ...
    'Component', NaN, ...
    'MeanHz', NaN, ...
    'SDHz', NaN, ...
    'Weight', NaN, ...
    'EffectiveN', NaN, ...
    'HardCount', NaN);
end

function scalars = empty_scalar_table()
row = measurement_template();
scalars = struct2table(row, 'AsArray', true);
scalars(1, :) = [];
end

function components = empty_component_parameter_table()
row = component_template();
components = struct2table(row, 'AsArray', true);
components(1, :) = [];
end

function scalars = failed_scalar_table(fit_info)
row = measurement_template();
row.NComponents = info_number(fit_info, 'n_components', NaN);
row.NObservations = info_number(fit_info, 'n_observations', NaN);
row.Replicates = info_number(fit_info, 'replicates', NaN);
row.RegularizationValue = info_number(fit_info, 'regularization', NaN);
scalars = struct2table(row, 'AsArray', true);
end
