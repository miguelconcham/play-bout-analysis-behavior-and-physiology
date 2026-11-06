function fit = irasa_power_law_fit(spec, frequency_range)
%IRASA_POWER_LAW_FIT  Log-log line on an IRASA fractal spectrum.
%
% Wen & Liu (2016), section 2.4, after the separation. The fractal
% spectrum is interpolated onto equally spaced log-frequency points so
% the high-frequency bins do not dominate the least-squares line. The
% slope of log10(power) versus log10(frequency) is negative for a
% descending spectrum. exponent is minus that slope, the positive beta
% in P(f) proportional to f^(-beta). log_intercept is the intercept in
% log10 power. broadband_log_power is the mean of the log-resampled
% log10 power, the broadband power in Miller et al. (2009) as used in
% the paper.
%
% FREQUENCY_RANGE is [fmin fmax] in Hz. DC cannot enter the log-log fit.

if nargin < 2 || ~isnumeric(frequency_range) || numel(frequency_range) ~= 2
    error('irasa_power_law_fit:range', ...
        'frequency_range must be [fmin fmax] in Hz.');
end
if frequency_range(1) <= 0 || frequency_range(2) <= frequency_range(1)
    error('irasa_power_law_fit:range', ...
        'The fit range must satisfy 0 < fmin < fmax.');
end
if ~isstruct(spec) || ~isfield(spec, 'frequency') || ~isfield(spec, 'fractal')
    error('irasa_power_law_fit:spec', ...
        'spec must come from irasa_separate and contain frequency and fractal.');
end

keep = spec.frequency >= frequency_range(1) & ...
    spec.frequency <= frequency_range(2);
frequency = spec.frequency(keep);
power = spec.fractal(keep, :);
if numel(frequency) < 3
    error('irasa_power_law_fit:range', ...
        'The fit range contains fewer than 3 frequency bins.');
end
if any(~isfinite(power(:)) | power(:) <= 0)
    error('irasa_power_law_fit:power', ...
        'The fractal spectrum in the fit range must be positive.');
end

log_frequency = log10(frequency);
log_power = log10(power);
log_grid = linspace(log_frequency(1), log_frequency(end), ...
    numel(log_frequency)).';
log_resampled = interp1(log_frequency, log_power, log_grid, 'linear');

n_channels = size(power, 2);
exponent = zeros(1, n_channels);
log_intercept = zeros(1, n_channels);
fitted_power = zeros(size(power));
for i_channel = 1:n_channels
    coefficients = polyfit(log_grid, log_resampled(:, i_channel), 1);
    exponent(i_channel) = -coefficients(1);
    log_intercept(i_channel) = coefficients(2);
    fitted_power(:, i_channel) = 10 .^ polyval(coefficients, log_frequency);
end

fit = struct();
fit.frequency = frequency;
fit.exponent = exponent;
fit.log_intercept = log_intercept;
fit.broadband_log_power = mean(log_resampled, 1);
fit.fitted_power = fitted_power;
fit.log_frequency_grid = log_grid;
fit.log_resampled_power = log_resampled;
fit.frequency_range = frequency_range(:).';
end
