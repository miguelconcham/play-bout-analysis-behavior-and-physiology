function fit = speed_peak_ols(x, y, min_n)
%SPEED_PEAK_OLS  Ordinary least squares of y on one predictor.
%
% fit is y ~ 1 + x. Rows with a non-finite value are omitted. Fewer than
% min_n remaining rows, or a constant predictor, returns ok = false and
% NaN coefficients. Standard errors use the residual variance with n - 2
% degrees of freedom. AIC matches fitlm: the Gaussian MLE counts the
% intercept, the slope, and the error variance.

x = double(x(:));
y = double(y(:));
if numel(x) ~= numel(y)
    error('speed_peak_ols:size', 'x and y must have one value per row.');
end
if nargin < 3 || isempty(min_n)
    min_n = 3;
end
if ~isscalar(min_n) || min_n < 3 || min_n ~= round(min_n)
    error('speed_peak_ols:min', 'min_n must be an integer of at least 3.');
end

fit = empty_ols_fit();
keep = isfinite(x) & isfinite(y);
x = x(keep);
y = y(keep);
fit.n = numel(x);
if fit.n < min_n || std(x, 0) == 0
    return
end

x_mean = mean(x);
y_mean = mean(y);
x_centered = x - x_mean;
sum_xx = x_centered.' * x_centered;
slope = (x_centered.' * (y - y_mean)) / sum_xx;
intercept = y_mean - slope * x_mean;
residual = y - (intercept + slope * x);
sse = residual.' * residual;
if sse < 0
    sse = 0;
end
df = fit.n - 2;
sigma2 = sse / df;
sigma2_mle = sse / fit.n;
fit.residual_df = df;
fit.intercept = intercept;
fit.slope = slope;
fit.sigma = sqrt(sigma2);
fit.r_squared = r_squared_of(y, y_mean, sse);

sum_x = sum(x);
sum_x2 = sum(x.^2);
denominator = fit.n * sum_x2 - sum_x * sum_x;
fit.intercept_se = sqrt(sigma2 * sum_x2 / denominator);
fit.slope_se = sqrt(sigma2 * fit.n / denominator);
[fit.intercept_t, fit.intercept_p] = t_and_p(intercept, fit.intercept_se, df);
[fit.slope_t, fit.slope_p] = t_and_p(slope, fit.slope_se, df);

if sigma2_mle > 0
    fit.log_likelihood = -fit.n / 2 * log(2 * pi) ...
        - fit.n / 2 * log(sigma2_mle) - fit.n / 2;
    fit.aic = -2 * fit.log_likelihood + 2 * 3;
else
    fit.log_likelihood = Inf;
    fit.aic = -Inf;
end
fit.ok = true;
end

function fit = empty_ols_fit()
fit = struct( ...
    'ok', false, ...
    'n', 0, ...
    'intercept', NaN, ...
    'slope', NaN, ...
    'intercept_se', NaN, ...
    'slope_se', NaN, ...
    'intercept_t', NaN, ...
    'slope_t', NaN, ...
    'intercept_p', NaN, ...
    'slope_p', NaN, ...
    'r_squared', NaN, ...
    'sigma', NaN, ...
    'log_likelihood', NaN, ...
    'aic', NaN, ...
    'residual_df', NaN);
end

function r_squared = r_squared_of(y, y_mean, sse)
sst = (y - y_mean).' * (y - y_mean);
if sst <= 0
    r_squared = NaN;
else
    r_squared = 1 - sse / sst;
end
end

function [t_stat, p_value] = t_and_p(estimate, se, df)
if ~isfinite(se) || se <= 0
    if ~isfinite(estimate)
        t_stat = NaN;
        p_value = NaN;
    elseif estimate == 0
        t_stat = 0;
        p_value = 1;
    else
        t_stat = sign(estimate) * Inf;
        p_value = 0;
    end
    return
end
t_stat = estimate / se;
p_value = betainc(df / (df + t_stat * t_stat), df / 2, 0.5);
end
