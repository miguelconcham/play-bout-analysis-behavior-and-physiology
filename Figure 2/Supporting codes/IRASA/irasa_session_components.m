function components = irasa_session_components(summary, raw_frequency, welch_windows, fourier_windows)
%IRASA_SESSION_COMPONENTS  Raw spectrogram plus IRASA decomposition.
%
% RAW_FREQUENCY and the window matrices come from
% irasa_raw_window_spectrum, one row per 10 s window. They are the raw
% spectra estimated before decomposition, on the same grid as the FOOOF
% play-session pipeline. SUMMARY comes from irasa_session_summary.
%
% components.raw holds the Welch and Fourier spectrograms and their
% linear-power means. components.aperiodic is the IRASA fractal spectrum.
% components.periodic is the IRASA oscillatory spectrum. Means follow
% the FOOOF script: average the windows in linear power, then floor the
% positive spectra at realmin before the log10 copy.

if nargin < 4
    error('irasa_session_components:input', ...
        'summary, raw_frequency, welch_windows, and fourier_windows are required.');
end
raw_frequency = raw_frequency(:);
n_windows = size(welch_windows, 1);
if size(fourier_windows, 1) ~= n_windows || ...
        size(welch_windows, 2) ~= numel(raw_frequency) || ...
        size(fourier_windows, 2) ~= numel(raw_frequency)
    error('irasa_session_components:raw', ...
        'Raw Welch and Fourier windows must share the raw frequency grid.');
end
if summary.n_windows ~= n_windows
    error('irasa_session_components:windows', ...
        'The raw spectrogram and the IRASA summary have different window counts.');
end

welch_mean = max(mean(welch_windows, 1, 'omitnan').', realmin('double'));
fourier_mean = max(mean(fourier_windows, 1, 'omitnan').', realmin('double'));
aperiodic_mean = summary.fractal(:);
periodic_mean = summary.oscillatory(:);

components = struct();
components.raw.frequency = raw_frequency;
components.raw.welch.window_power = welch_windows;
components.raw.welch.mean_power = welch_mean;
components.raw.welch.log10_mean_power = log10(welch_mean);
components.raw.fourier.window_power = fourier_windows;
components.raw.fourier.mean_power = fourier_mean;
components.raw.fourier.log10_mean_power = log10(fourier_mean);

components.aperiodic.frequency = summary.frequency(:);
components.aperiodic.window_power = summary.window_fractal;
components.aperiodic.mean_power = aperiodic_mean;
components.aperiodic.log10_mean_power = log10(max(aperiodic_mean, realmin('double')));
components.aperiodic.exponent = summary.exponent;
components.aperiodic.log_intercept = summary.log_intercept;
components.aperiodic.broadband_log_power = summary.broadband_log_power;
components.aperiodic.fit_frequency = summary.fit_frequency;
components.aperiodic.fitted_power = summary.fitted_power;
components.aperiodic.fit_range = summary.fit_range;
components.aperiodic.window_exponent = summary.window_exponent;
components.aperiodic.window_log_intercept = summary.window_log_intercept;
components.aperiodic.window_broadband_log_power = ...
    summary.window_broadband_log_power;

components.periodic.frequency = summary.frequency(:);
components.periodic.window_power = summary.window_oscillatory;
components.periodic.mean_power = periodic_mean;
end
