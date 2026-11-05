function components = irasa_session_components(window_specs, raw_frequency, welch_windows, fourier_windows, fit_range)
%IRASA_SESSION_COMPONENTS  Window-by-window raw, fractal, and periodic spectra.
%
% WINDOW_SPECS is one irasa_separate output per window. RAW_FREQUENCY,
% WELCH_WINDOWS, and FOURIER_WINDOWS come from irasa_raw_window_spectrum,
% one row per window. FIT_RANGE is [fmin fmax] in Hz for the power-law
% line of each window's fractal spectrum.
%
% The returned spectra are the windows themselves. This function does not
% average across windows. components.raw.welch and components.raw.fourier
% are the power estimates. components.fractal and components.periodic are
% the IRASA decomposition of the same windows.

if nargin < 5
    error('irasa_session_components:input', ...
        'window_specs, raw spectra, and fit_range are required.');
end
if isempty(window_specs)
    error('irasa_session_components:empty', ...
        'At least one IRASA window is required.');
end

raw_frequency = raw_frequency(:);
n_windows = numel(window_specs);
if size(welch_windows, 1) ~= n_windows || size(fourier_windows, 1) ~= n_windows || ...
        size(welch_windows, 2) ~= numel(raw_frequency) || ...
        size(fourier_windows, 2) ~= numel(raw_frequency)
    error('irasa_session_components:raw', ...
        'Raw Welch and Fourier windows must share one row per IRASA window.');
end

frequency = window_specs(1).frequency(:);
n_frequency = numel(frequency);
fractal = zeros(n_windows, n_frequency);
periodic = zeros(n_windows, n_frequency);
window_exponent = zeros(n_windows, 1);
window_log_intercept = zeros(n_windows, 1);
window_broadband_log_power = zeros(n_windows, 1);

for window_n = 1:n_windows
    this_window = window_specs(window_n);
    if numel(this_window.frequency) ~= n_frequency || ...
            max(abs(this_window.frequency(:) - frequency)) > 1e-9
        error('irasa_session_components:frequency', ...
            'IRASA windows in one session do not share a frequency grid.');
    end
    fractal(window_n, :) = this_window.fractal(:).';
    periodic(window_n, :) = this_window.oscillatory(:).';
    this_fit = irasa_power_law_fit(this_window, fit_range);
    window_exponent(window_n) = this_fit.exponent;
    window_log_intercept(window_n) = this_fit.log_intercept;
    window_broadband_log_power(window_n) = this_fit.broadband_log_power;
end

components = struct();
components.frequency = frequency;
components.n_windows = n_windows;
components.raw.frequency = raw_frequency;
components.raw.welch = welch_windows;
components.raw.fourier = fourier_windows;
components.fractal = fractal;
components.periodic = periodic;
components.window_exponent = window_exponent;
components.window_log_intercept = window_log_intercept;
components.window_broadband_log_power = window_broadband_log_power;
components.fit_range = fit_range(:).';
end
