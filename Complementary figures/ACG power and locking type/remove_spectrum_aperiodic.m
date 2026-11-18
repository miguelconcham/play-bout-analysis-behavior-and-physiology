function flat = remove_spectrum_aperiodic(stack, f_grid)
%REMOVE_SPECTRUM_APERIODIC  Subtract the fixed aperiodic fit from each row.
%
% The fit is the FOOOF fixed model, offset minus exponent times
% log10(frequency). A first fit uses every positive bin. A second fit
% keeps the bins that sit on that background, so a peak is left in the
% residual. The residual is log10 power minus the second fit. A row with
% fewer than four positive bins stays NaN.

stack = double(stack);
f_grid = double(f_grid(:)).';
if isempty(stack)
    flat = stack;
    return
end
if isvector(stack)
    stack = stack(:).';
end
flat = nan(size(stack));
if numel(f_grid) ~= size(stack, 2)
    return
end
for row = 1:size(stack, 1)
    flat(row, :) = flatten_one_spectrum(stack(row, :), f_grid);
end
end


function residual = flatten_one_spectrum(power, f_grid)
residual = nan(size(power));
use = isfinite(power) & isfinite(f_grid) & power > 0 & f_grid > 0;
if nnz(use) < 4
    return
end
frequency = f_grid(use);
log_power = log10(power(use));
params = fit_fixed_aperiodic(frequency, log_power);
if isempty(params)
    return
end
params = refit_background(frequency, log_power, params);
if isempty(params)
    return
end
fit = params(1) - params(2) .* log10(frequency);
residual(use) = log_power - fit;
end


function params = refit_background(frequency, log_power, params)
%REFIT_BACKGROUND  Repeat the fit on bins that do not rise above it.
%
% Negative residuals are set to zero before the 2.5th percentile, matching
% the FOOOF robust aperiodic step. The peak then stays out of the refit.

fit = params(1) - params(2) .* log10(frequency);
flattened = log_power - fit;
flattened(flattened < 0) = 0;
threshold = low_percentile(flattened, 2.5);
keep = flattened <= threshold;
if nnz(keep) < 4
    return
end
refit = fit_fixed_aperiodic(frequency(keep), log_power(keep));
if ~isempty(refit)
    params = refit;
end
end


function params = fit_fixed_aperiodic(frequency, log_power)
%FIT_FIXED_APERIODIC  offset - exponent * log10(frequency).

frequency = double(frequency(:));
log_power = double(log_power(:));
log_frequency = log10(frequency);
if numel(log_power) < 2 || max(log_frequency) <= min(log_frequency)
    params = [];
    return
end
predictors = [ones(numel(log_power), 1), -log_frequency];
params = predictors \ log_power;
params = params(:).';
if numel(params) ~= 2 || any(~isfinite(params))
    params = [];
end
end


function value = low_percentile(samples, percent)
samples = sort(double(samples(:)));
rank = 1 + (percent / 100) * (numel(samples) - 1);
low = floor(rank);
high = ceil(rank);
low = min(max(low, 1), numel(samples));
high = min(max(high, 1), numel(samples));
if low == high
    value = samples(low);
    return
end
value = samples(low) + (rank - low) * (samples(high) - samples(low));
end
