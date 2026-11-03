function options = speed_peak_options(options)
%SPEED_PEAK_OPTIONS  Settings for speed-binned delta and theta peak frequency.
%
% Pass a struct to override any field. Empty fields keep the default.
% The same struct is copied into the saved recordings so a result file
% can be regenerated from the parameters it carries.

if nargin < 1 || isempty(options)
    options = struct();
end

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];

options = fill_option(options, 'data_root', data_root);
options = fill_option(options, 'sampling_rate', 2500);
options = fill_option(options, 'n_lfp_channels', 384);
options = fill_option(options, 'area_name', 'LPAG');
options = fill_option(options, 'video_fps', 30);
options = fill_option(options, 'tracking_smooth', 'loess');
options = fill_option(options, 'tracking_smooth_samples', 5);
options = fill_option(options, 'position_interpolation', 'cubic');

% One non-overlapping window, matching the partner-session FOOOF spectra.
% 10 s gives 0.1 Hz bins, the grid used by fooof_matlab in this repository.
options = fill_option(options, 'window_seconds', 10);
options = fill_option(options, 'window_overlap_seconds', 0);
options = fill_option(options, 'frequency_max', 40);
options = fill_option(options, 'welch_window', 'periodic Hann');
options = fill_option(options, 'fourier_window', 'rectangular');
options = fill_option(options, 'detrend', 'constant');
options = fill_option(options, 'averaging_domain', 'linear power');

options = fill_option(options, 'fooof_aperiodic_mode', 'fixed');
options = fill_option(options, 'fooof_peak_width_limits', [0.5 12]);
options = fill_option(options, 'fooof_max_n_peaks', Inf);
options = fill_option(options, 'fooof_min_peak_height', 0);
options = fill_option(options, 'fooof_peak_threshold', 2);

% Each speed bin holds at most the upper fraction of that session's finite
% speed samples, and at least the lower fraction when the sample count
% allows it. 10% is preferred because it is the largest share inside 5-10%.
options = fill_option(options, 'speed_bin_fraction_limits', [0.05 0.10]);
options = fill_option(options, 'min_window_speed_coverage', 0.8);

options = fill_option(options, 'delta_band', [1 5]);
options = fill_option(options, 'theta_band', [6 12]);
options = fill_option(options, 'min_bins_for_correlation', 3);
options = fill_option(options, 'speed_units', 'pixels per second');
options = fill_option(options, 'peak_selection', ...
    ['highest FOOOF peak power inside the band; ' ...
    'ties take the lower center frequency']);

validate_options(options);
options = derive_speed_peak_parameters(options);
end

function options = derive_speed_peak_parameters(options)
% Derived quantities are overwritten on every call so a changed window
% length cannot leave a stale frequency grid in the struct.
window_samples = options.window_seconds * options.sampling_rate;
if window_samples ~= round(window_samples)
    error('speed_peak_options:window', ...
        'window_seconds times sampling_rate must be an integer sample count.');
end
options.window_samples = window_samples;
options.frequency_resolution = options.sampling_rate / window_samples;
if abs(options.frequency_resolution - 1 / options.window_seconds) > 1e-9
    error('speed_peak_options:resolution', ...
        'The window does not give an integer number of frequency bins per Hz.');
end
options.n_frequency = round(options.frequency_max / options.frequency_resolution);
options.frequency = (1:options.n_frequency).' * options.frequency_resolution;
options.frequency_range = options.frequency([1 end]).';
options.min_speed_samples_per_window = ceil( ...
    options.min_window_speed_coverage * options.window_seconds * options.video_fps);
options.speed_bin_rule = ['Equal-count bins of the median speed of each spectral ' ...
    'window. Each bin holds at most the upper entry of speed_bin_fraction_limits ' ...
    'of the windows that have at least min_speed_samples_per_window finite speed ' ...
    'samples, and at least the lower entry when the window count allows it.'];
options.spectrum_speed_definition = ['median of the per-window median speeds, ' ...
    'in pixels per second, of the windows in the bin'];
options.speed_definition = ['frame-to-frame pixel displacement of the ' ...
    'loess-smoothed, cubic-interpolated focal animal, divided by the frame time step'];
options.session_definition = ['Session_structure / Partners session row overlapped ' ...
    'with one partner tracking file. Rows are not pooled.'];
options.spectral_methods = {'Welch', 'Fourier'};
options.correlation_method = ['Pearson and Spearman across speed bins with both ' ...
    'values finite. p-values use the two-sided Student-t transform of r ' ...
    'with n-2 degrees of freedom. A constant series or fewer than ' ...
    'min_bins_for_correlation pairs returns NaN.'];
end

function validate_options(options)
if options.sampling_rate <= 0 || options.video_fps <= 0 || options.window_seconds <= 0
    error('speed_peak_options:time', ...
        'sampling_rate, video_fps, and window_seconds must be positive.');
end
if options.window_overlap_seconds ~= 0
    error('speed_peak_options:overlap', ...
        'Windows are non-overlapping, matching the other FOOOF spectra.');
end
if options.frequency_max <= 0 || options.frequency_max >= options.sampling_rate / 2
    error('speed_peak_options:frequency', ...
        'frequency_max must lie between 0 and the Nyquist frequency.');
end
limits = options.speed_bin_fraction_limits;
if numel(limits) ~= 2 || limits(1) <= 0 || limits(2) < limits(1) || limits(2) >= 1
    error('speed_peak_options:bins', ...
        'speed_bin_fraction_limits must be [min max] with 0 < min <= max < 1.');
end
if ~isscalar(options.min_window_speed_coverage) || ...
        options.min_window_speed_coverage <= 0 || options.min_window_speed_coverage > 1
    error('speed_peak_options:coverage', ...
        'min_window_speed_coverage must lie in (0, 1].');
end
if numel(options.delta_band) ~= 2 || numel(options.theta_band) ~= 2 || ...
        options.delta_band(2) < options.delta_band(1) || ...
        options.theta_band(2) < options.theta_band(1)
    error('speed_peak_options:band', ...
        'delta_band and theta_band must be [low high] in Hz.');
end
if ~isscalar(options.min_bins_for_correlation) || ...
        options.min_bins_for_correlation < 3 || ...
        options.min_bins_for_correlation ~= round(options.min_bins_for_correlation)
    error('speed_peak_options:correlation', ...
        'min_bins_for_correlation must be an integer of at least 3.');
end
if ~strcmp(options.fooof_aperiodic_mode, 'fixed') && ...
        ~strcmp(options.fooof_aperiodic_mode, 'knee')
    error('speed_peak_options:fooof', ...
        'fooof_aperiodic_mode must be fixed or knee.');
end
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
