function features = acg_nonoscillatory_components(acg_matrix, lags)
%ACG_NONOSCILLATORY_COMPONENTS  Gamma and difference-of-exponentials fits.
%
% Each row of acg_matrix is one autocorrelogram. lags are in seconds and
% match the columns. Only positive lags are used. The gamma kernel is
% amplitude * lag^(shape-1) * exp(-lag/scale), divided by its own peak,
% plus a baseline. Amplitude is then the height of that bump. The
% difference of exponentials is amplitude * (exp(-lag/decay) -
% exp(-lag/rise)), plus a baseline, with decay longer than rise. The table
% holds each fit's parameters and its goodness of fit. R2 is 1 minus the
% residual sum of squares over the total sum of squares. RMSE is the root
% mean squared residual. PeakAmplitude is the maximum of the positive-lag
% curve, and AUC is its integral over lag in seconds. Those two do not use
% the fit. A neuron whose curve cannot be fit is NaN in the fit columns.

acg_matrix = double(acg_matrix);
lags = double(lags(:)).';
if size(acg_matrix, 2) ~= numel(lags)
    error('acg_nonoscillatory_components:lag', ...
        'The autocorrelogram width must match the lag axis.');
end
use = lags > 0 & isfinite(lags);
lags = lags(use);
acg_matrix = acg_matrix(:, use);
names = { ...
    'GammaAmplitude', 'GammaShape', 'GammaScale', 'GammaBaseline', ...
    'GammaR2', 'GammaRMSE', ...
    'DiffExpAmplitude', 'DiffExpRise', 'DiffExpDecay', 'DiffExpBaseline', ...
    'DiffExpR2', 'DiffExpRMSE', 'PeakAmplitude', 'AUC'};
n = size(acg_matrix, 1);
values = nan(n, numel(names));
options = optimset('Display', 'off', 'MaxIter', 160, 'MaxFunEvals', 500);
for i = 1:n
    y = acg_matrix(i, :);
    keep = isfinite(y) & isfinite(lags);
    if nnz(keep) < 2
        continue
    end
    t = lags(keep);
    y = y(keep);
    values(i, 13) = max(y);
    values(i, 14) = trapz(t(:), y(:));
    if nnz(keep) < 8
        continue
    end
    values(i, 1:6) = fit_gamma(t, y, options);
    values(i, 7:12) = fit_difference_of_exponentials(t, y, options);
end
features = array2table(values, 'VariableNames', names);
end


function row = fit_gamma(t, y, options)
% y = baseline + amplitude * t^(shape-1) * exp(-t/scale)
row = nan(1, 6);
baseline0 = max(median(y(max(1, end - 4):end)), 0);
peak = max(y) - baseline0;
if ~isfinite(peak) || peak <= 0
    peak = max(max(y), 0) + 1e-6;
end
[~, peak_at] = max(y);
scale0 = max(t(peak_at), median(diff(t)));
starts = [ ...
    gamma_start(peak, 2, scale0, baseline0); ...
    gamma_start(peak, 1.2, max(scale0 / 2, 0.01), baseline0); ...
    gamma_start(peak, 5, max(scale0 * 2, 0.05), baseline0)];
best = inf;
best_x = starts(1, :);
for s = 1:size(starts, 1)
    try
        x = fminsearch(@(z) gamma_loss(z, t, y), starts(s, :), options);
    catch
        continue
    end
    loss = gamma_loss(x, t, y);
    if loss < best
        best = loss;
        best_x = x;
    end
end
if ~isfinite(best) || best >= 1e5
    return
end
amplitude = exp(best_x(1));
shape = exp(best_x(2));
scale = exp(best_x(3));
baseline = exp(best_x(4));
yhat = baseline + amplitude * gamma_kernel(t, shape, scale);
[r2, rmse] = goodness_of_fit(y, yhat);
row = [amplitude, shape, scale, baseline, r2, rmse];
end


function start = gamma_start(peak, shape, scale, baseline)
mode_t = max((shape - 1) * scale, scale);
kernel = gamma_kernel(mode_t, shape, scale);
amplitude = peak / max(kernel, 1e-8);
start = log([max(amplitude, 1e-8), shape, scale, max(baseline, 1e-6)]);
end


function loss = gamma_loss(x, t, y)
amplitude = exp(x(1));
shape = exp(x(2));
scale = exp(x(3));
baseline = exp(x(4));
if shape < 0.3 || shape > 20 || scale < 0.002 || scale > 1.5
    loss = 1e6;
    return
end
yhat = baseline + amplitude * gamma_kernel(t, shape, scale);
if any(~isfinite(yhat))
    loss = 1e6;
    return
end
loss = mean((y - yhat) .^ 2);
end


function g = gamma_kernel(t, shape, scale)
log_g = (shape - 1) * log(t) - t ./ scale;
log_g(t <= 0) = -inf;
log_g(log_g > 40) = 40;
g = exp(log_g);
g(~isfinite(g)) = 0;
peak = max(g);
if peak > 0
    g = g / peak;
end
end


function row = fit_difference_of_exponentials(t, y, options)
% y = baseline + amplitude * (exp(-t/decay) - exp(-t/rise)), decay > rise
row = nan(1, 6);
baseline0 = max(median(y(max(1, end - 4):end)), 0);
peak = max(y) - baseline0;
if ~isfinite(peak) || peak <= 0
    peak = max(max(y), 0) + 1e-6;
end
rise0 = max(min(t(t > 0)), 0.005);
decay0 = max(0.05, rise0 * 4);
height0 = exp_difference_height(rise0, decay0);
amplitude0 = peak / max(height0, 1e-6);
start = log([max(amplitude0, 1e-6), rise0, decay0 - rise0, max(baseline0, 1e-6)]);
try
    x = fminsearch(@(z) diffexp_loss(z, t, y), start, options);
catch
    return
end
if diffexp_loss(x, t, y) >= 1e5
    return
end
amplitude = exp(x(1));
rise = exp(x(2));
decay = rise + exp(x(3));
baseline = exp(x(4));
yhat = baseline + amplitude * (exp(-t ./ decay) - exp(-t ./ rise));
[r2, rmse] = goodness_of_fit(y, yhat);
row = [amplitude, rise, decay, baseline, r2, rmse];
end


function loss = diffexp_loss(x, t, y)
amplitude = exp(x(1));
rise = exp(x(2));
decay = rise + exp(x(3));
baseline = exp(x(4));
if rise < 0.001 || decay > 3
    loss = 1e6;
    return
end
yhat = baseline + amplitude * (exp(-t ./ decay) - exp(-t ./ rise));
if any(~isfinite(yhat))
    loss = 1e6;
    return
end
loss = mean((y - yhat) .^ 2);
end


function height = exp_difference_height(rise, decay)
% Maximum of exp(-t/decay) - exp(-t/rise) for decay > rise.
if ~(decay > rise)
    height = 0;
    return
end
t_peak = log(decay / rise) / (1 / rise - 1 / decay);
height = exp(-t_peak / decay) - exp(-t_peak / rise);
end


function [r2, rmse] = goodness_of_fit(y, yhat)
residual = y - yhat;
ss_res = sum(residual .^ 2);
ss_tot = sum((y - mean(y)) .^ 2);
if ss_tot <= eps(max(abs(y)))
    if ss_res <= eps(max(abs(y)))
        r2 = 1;
    else
        r2 = 0;
    end
else
    r2 = 1 - ss_res / ss_tot;
end
rmse = sqrt(mean(residual .^ 2));
end
