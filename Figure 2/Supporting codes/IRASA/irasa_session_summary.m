function summary = irasa_session_summary(window_specs, fit_range)
%IRASA_SESSION_SUMMARY  Average IRASA windows, then fit the fractal line.
%
% WINDOW_SPECS is a struct array of irasa_separate outputs, one per
% sliding window of one Partners session. Spectra are averaged
% in linear power, the same domain used for the Welch play-session
% spectra. The oscillatory spectrum is then mixed minus fractal, which
% matches the mean of the window residuals. The power-law fit is applied
% to that session-mean fractal spectrum.

if nargin < 2
    error('irasa_session_summary:input', ...
        'window_specs and fit_range are required.');
end
if isempty(window_specs)
    error('irasa_session_summary:empty', ...
        'At least one IRASA window is required.');
end

n_windows = numel(window_specs);
frequency = window_specs(1).frequency(:);
n_frequency = numel(frequency);
window_mixed = zeros(n_windows, n_frequency);
window_fractal = zeros(n_windows, n_frequency);
window_oscillatory = zeros(n_windows, n_frequency);
window_exponent = zeros(n_windows, 1);
window_log_intercept = zeros(n_windows, 1);
window_broadband_log_power = zeros(n_windows, 1);

for window_n = 1:n_windows
    this_window = window_specs(window_n);
    if numel(this_window.frequency) ~= n_frequency || ...
            max(abs(this_window.frequency(:) - frequency)) > 1e-9
        error('irasa_session_summary:frequency', ...
            'IRASA windows in one session do not share a frequency grid.');
    end
    window_mixed(window_n, :) = this_window.mixed(:).';
    window_fractal(window_n, :) = this_window.fractal(:).';
    window_oscillatory(window_n, :) = this_window.oscillatory(:).';
    this_fit = irasa_power_law_fit(this_window, fit_range);
    window_exponent(window_n) = this_fit.exponent;
    window_log_intercept(window_n) = this_fit.log_intercept;
    window_broadband_log_power(window_n) = this_fit.broadband_log_power;
end

mixed = mean(window_mixed, 1).';
fractal = mean(window_fractal, 1).';
mean_spec = struct('frequency', frequency, 'fractal', fractal);
session_fit = irasa_power_law_fit(mean_spec, fit_range);

summary = struct();
summary.frequency = frequency;
summary.mixed = mixed;
summary.fractal = fractal;
summary.oscillatory = mixed - fractal;
summary.window_mixed = window_mixed;
summary.window_fractal = window_fractal;
summary.window_oscillatory = window_oscillatory;
summary.window_exponent = window_exponent;
summary.window_log_intercept = window_log_intercept;
summary.window_broadband_log_power = window_broadband_log_power;
summary.exponent = session_fit.exponent;
summary.log_intercept = session_fit.log_intercept;
summary.broadband_log_power = session_fit.broadband_log_power;
summary.fit_frequency = session_fit.frequency;
summary.fitted_power = session_fit.fitted_power;
summary.fit_range = fit_range(:).';
summary.n_windows = n_windows;
summary.averaging_domain = 'linear power';
end
