function result = compare_play_power_acgs(side_a, side_b, centers, rate_a, rate_b, scale, options)
%COMPARE_PLAY_POWER_ACGS  Paired time-wise and property tests of two ACGs.
%
% side_a and side_b are neurons by lag bins. Each row is one neuron in
% one condition. For acg_method 'overlap', or when no method is set, a
% row is pair counts / overlap and the rate ratio divides by firing rate
% squared times the bin width. For 'eligible' a row is the mean partner
% count of spikes that can see the bin, and the rate ratio divides by
% firing rate times the bin width. Bins covering the refractory gap are
% blanked. scale 'sum' divides each neuron by its own finite mass. scale
% 'ratio' restores that mass and applies the rate ratio, so a flat
% Poisson autocorrelogram sits at 1.
%
% The cluster test is the sign-flip test on the Gaussian-smoothed curves.
% Difference is side B minus side A. Properties are the locking-period
% burst, prominences, latencies, and folded delta power.

if nargin < 7 || isempty(options)
    options = struct();
end
options = fill_compare_options(options);
side_a = double(side_a);
side_b = double(side_b);
centers = double(centers(:)');
if isvector(side_a)
    side_a = side_a(:)';
end
if isvector(side_b)
    side_b = side_b(:)';
end
n = size(side_a, 1);
if size(side_b, 1) ~= n || size(side_a, 2) ~= numel(centers) || size(side_b, 2) ~= numel(centers)
    error('compare_play_power_acgs:size', ...
        'The two autocorrelograms must share neurons and lag bins.');
end
rate_a = double(rate_a(:));
rate_b = double(rate_b(:));
if numel(rate_a) ~= n || numel(rate_b) ~= n
    error('compare_play_power_acgs:rate', 'Each neuron needs a firing rate on both sides.');
end

blank_a = blank_acg_center(side_a, centers, options.refractory_s);
blank_b = blank_acg_center(side_b, centers, options.refractory_s);
mass_a = sum_omitnan_rows(blank_a);
mass_b = sum_omitnan_rows(blank_b);
curve_a = normalize_acg_rows(blank_a);
curve_b = normalize_acg_rows(blank_b);
if strcmp(scale, 'ratio')
    bin_w = median(diff(centers));
    rate_exponent = acg_rate_exponent(options.acg_method);
    curve_a = scale_acg_by_rate(curve_a, mass_a, rate_a, bin_w, rate_exponent);
    curve_b = scale_acg_by_rate(curve_b, mass_b, rate_b, bin_w, rate_exponent);
elseif ~strcmp(scale, 'sum')
    error('compare_play_power_acgs:scale', 'scale must be ''sum'' or ''ratio''.');
end
result.keep = any(isfinite(curve_a), 2) & any(isfinite(curve_b), 2);
result.curve_a = curve_a(result.keep, :);
result.curve_b = curve_b(result.keep, :);
result.n = sum(result.keep);
smooth_n = max(1, round(options.smooth_half / median(diff(centers))));
result.smooth_a = smooth_acg_away_from_zero(result.curve_a, centers, options.skip_abs_lag, smooth_n);
result.smooth_b = smooth_acg_away_from_zero(result.curve_b, centers, options.skip_abs_lag, smooth_n);
result.scores_a = score_acg_curves(result.curve_a, centers, smooth_n, options.delta_lag);
result.scores_b = score_acg_curves(result.curve_b, centers, smooth_n, options.delta_lag);
result.clusters = paired_acg_clusters(result.smooth_a, result.smooth_b, centers, ...
    options.cluster_alpha, options.cluster_perm);
result.properties = score_property_pairs(result.scores_a, result.scores_b);
result.scale = scale;
end


function options = fill_compare_options(options)
if ~isfield(options, 'refractory_s') || isempty(options.refractory_s)
    options.refractory_s = 0.005;
end
if ~isfield(options, 'smooth_half') || isempty(options.smooth_half)
    options.smooth_half = 0.06;
end
if ~isfield(options, 'skip_abs_lag') || isempty(options.skip_abs_lag)
    options.skip_abs_lag = 0.02;
end
if ~isfield(options, 'delta_lag') || isempty(options.delta_lag)
    options.delta_lag = [0.2 0.5];
end
if ~isfield(options, 'cluster_alpha') || isempty(options.cluster_alpha)
    options.cluster_alpha = 0.05;
end
if ~isfield(options, 'cluster_perm') || isempty(options.cluster_perm)
    options.cluster_perm = 1000;
end
if ~isfield(options, 'acg_method') || isempty(options.acg_method)
    options.acg_method = 'overlap';
end
end


function Y = scale_acg_by_rate(Y, mass, rate, bin_w, rate_exponent)
if nargin < 5 || isempty(rate_exponent)
    rate_exponent = 2;
end
factor = mass(:) ./ ((rate(:) .^ rate_exponent) * bin_w);
factor(~isfinite(factor) | rate(:) <= 0 | mass(:) <= 0) = NaN;
Y = Y .* factor;
end


function total = sum_omitnan_rows(X)
filled = X;
filled(isnan(filled)) = 0;
total = sum(filled, 2);
end


function Y = smooth_acg_away_from_zero(X, centers, skip_abs_lag, smooth_n)
Y = X;
if isempty(X) || smooth_n <= 1
    return
end
pos = centers > skip_abs_lag;
neg = centers < -skip_abs_lag;
for j = 1:size(X, 1)
    if any(pos)
        Y(j, pos) = gaussian_smooth_acg(X(j, pos), smooth_n / 4);
    end
    if any(neg)
        Y(j, neg) = gaussian_smooth_acg(X(j, neg), smooth_n / 4);
    end
end
end


function properties = score_property_pairs(scores_a, scores_b)
specs = acg_property_specs();
properties = struct();
for m = 1:numel(specs)
    name = specs(m).name;
    ok = acg_measure_mask(scores_a, scores_b, name);
    ya = scores_a.(name);
    yb = scores_b.(name);
    stats.n = sum(ok);
    if stats.n < 1
        stats.median_a = NaN;
        stats.median_b = NaN;
        stats.median_delta = NaN;
    else
        stats.median_a = median(ya(ok));
        stats.median_b = median(yb(ok));
        stats.median_delta = median(yb(ok) - ya(ok));
    end
    stats.p_two = paired_signrank_p(yb(ok), ya(ok), 'both');
    if strcmp(specs(m).kind, 'level')
        stats.p_increase = paired_signrank_p(yb(ok), ya(ok), 'right');
    else
        stats.p_increase = NaN;
    end
    stats.label = specs(m).label;
    stats.kind = specs(m).kind;
    properties.(name) = stats;
end
end
