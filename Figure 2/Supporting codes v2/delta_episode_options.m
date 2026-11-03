function options = delta_episode_options(options)
%DELTA_EPISODE_OPTIONS  Settings for continuous mid-LPAG delta episodes.
%
% Pass a struct to override any field. Empty fields keep the default.
% A numeric zero is kept; it is not treated as empty.
%
% The detector looks for candidate runs of elevated 1-5 Hz amplitude in
% the continuous LFP. It does not decide that a run is a neural
% oscillation: a bandpass filter can ring around an isolated transient.

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
options = fill_option(options, 'delta_freq', [1 5]);
options = fill_option(options, 'delta_filter_order', 2000);
options = fill_option(options, 'filter_coefficients', []);

% 'std':      mean(amplitude) + amplitude_threshold * std(amplitude)
% 'median':   amplitude_threshold * median(amplitude)
% 'absolute': amplitude_threshold, in the same units as the LFP
options = fill_option(options, 'amplitude_threshold_mode', 'std');
options = fill_option(options, 'amplitude_threshold', 1);
options = fill_option(options, 'amplitude_smooth_seconds', 0);

% Runs separated by at most this gap are one episode. Zero keeps every
% contiguous above-threshold run separate.
options = fill_option(options, 'max_gap_seconds', 0);
options = fill_option(options, 'min_cycles', 3);
options = fill_option(options, 'min_duration_seconds', 0);

% Short enough that cycles faster than 5 Hz stay visible in the histogram.
% The delta-peak CSD code uses 1/5 Hz instead, which would hide them.
options = fill_option(options, 'min_peak_distance_seconds', 0.05);

% Samples within this margin of each filtered interval are not searched.
% Empty uses the FIR order divided by the sampling rate (filtfilt edges).
options = fill_option(options, 'edge_margin_seconds', []);

% 'partners_session': each Session_structure / Partners session row.
% 'full_recording': the whole LFP file as one interval.
options = fill_option(options, 'interval_mode', 'partners_session');

options = fill_option(options, 'example_pad_seconds', 2);
options = fill_option(options, 'example_max_seconds', 12);
options = fill_option(options, 'histogram_ms', [50 1000]);
options = fill_option(options, 'histogram_bin_ms', 25);

% Cole & Voytek (2019), J Neurophysiol 122:849-861. Extrema are read on
% a low-pass at about four times the oscillation frequency. Their theta
% example (4-10 Hz) uses 40 Hz, so the default here is four times the
% upper edge of delta_freq. The simulation low-pass is described as a
% 100 ms filter; cycle_lowpass_order is that duration in samples.
if ~isfield(options, 'cycle_lowpass_hz') || isempty(options.cycle_lowpass_hz)
    options.cycle_lowpass_hz = 4 * options.delta_freq(2);
end
options = fill_option(options, 'cycle_lowpass_order', ...
    round(0.1 * options.sampling_rate));
options = fill_option(options, 'cycle_lowpass_coefficients', []);

% Published sequence thresholds. Not the later bycycle defaults
% (amplitude fraction 0, consistency 0.5, period consistency 0.5,
% monotonicity 0.8). 'motor_beta' is the empirical LFP set: amplitude
% consistency 0.3, period consistency 0.5, monotonicity 0.6, and
% amplitude at or above the 20th percentile of cycles in the interval.
% 'simulation' is 0.6, 0.6, and 0.9, with no amplitude percentile.
% Cole and Voytek chose both sets by eye. The beta set was for 13-30 Hz
% motor cortex, not for delta; it is the published real-LFP suggestion.
% Either symmetry, period CV, and amplitude CV stay continuous and are
% not thresholded. A missing numeric field takes the preset; a value
% passed in the struct replaces the preset.
options = fill_option(options, 'cycle_threshold_preset', 'motor_beta');
preset = cycle_threshold_preset(options.cycle_threshold_preset);
options = fill_option(options, 'amp_consistency_threshold', preset.amp_consistency);
options = fill_option(options, 'period_consistency_threshold', preset.period_consistency);
options = fill_option(options, 'monotonicity_threshold', preset.monotonicity);
options = fill_option(options, 'cycle_amplitude_percentile', preset.amplitude_percentile);
options = fill_option(options, 'min_consecutive_cycles', 3);
options.cycle_threshold_note = [ ...
    'Cole & Voytek 2019 J Neurophysiol 122:849-861. ' ...
    'motor_beta: amplitude consistency 0.3, period consistency 0.5, ' ...
    'monotonicity 0.6, amplitude at or above the 20th percentile. ' ...
    'simulation: 0.6, 0.6, 0.9, no amplitude percentile. ' ...
    'At least 3 consecutive cycles. Chosen by those authors by eye; ' ...
    'the beta numbers are not a delta-specific recommendation. ' ...
    'Bycycle library defaults (0.5, 0.5, 0.8) are not used. ' ...
    'Period CV, amplitude CV, and symmetry are not gated.'];

options = fill_option(options, 'interpretation', ...
    ['Candidate episodes of elevated 1-5 Hz amplitude in the continuous ' ...
     'LFP. Filtering can produce apparent oscillations around an isolated ' ...
     'transient, so compare each episode with the minimally processed trace. ' ...
     'Cycle-by-cycle measures describe trough-to-trough structure inside ' ...
     'those episodes. They do not prove a neural origin and do not exclude ' ...
     'movement or electrode artifacts.']);

if numel(options.delta_freq) ~= 2 || any(options.delta_freq(:) <= 0) || ...
        options.delta_freq(1) >= options.delta_freq(2)
    error('delta_episode_options:delta_freq', ...
        'delta_freq must be [low high] in Hz.');
end
if ~ismember(options.amplitude_threshold_mode, {'std', 'median', 'absolute'})
    error('delta_episode_options:threshold', ...
        'amplitude_threshold_mode must be std, median, or absolute.');
end
if ~ismember(options.interval_mode, {'partners_session', 'full_recording'})
    error('delta_episode_options:interval', ...
        'interval_mode must be partners_session or full_recording.');
end
if options.min_cycles < 0 || options.min_duration_seconds < 0 || ...
        options.max_gap_seconds < 0 || options.min_peak_distance_seconds < 0
    error('delta_episode_options:criterion', ...
        'Cycle, duration, gap, and peak-distance criteria must be >= 0.');
end
if options.cycle_lowpass_hz <= 0 || ...
        options.cycle_lowpass_hz >= options.sampling_rate / 2
    error('delta_episode_options:lowpass', ...
        'cycle_lowpass_hz must lie between 0 and the Nyquist frequency.');
end
if options.cycle_lowpass_order < 1 || ...
        options.cycle_lowpass_order ~= round(options.cycle_lowpass_order)
    error('delta_episode_options:lowpass', ...
        'cycle_lowpass_order must be a positive integer.');
end
if options.min_consecutive_cycles < 1
    error('delta_episode_options:consecutive', ...
        'min_consecutive_cycles must be at least 1.');
end
check_unit_interval(options.amp_consistency_threshold, 'amp_consistency_threshold');
check_unit_interval(options.period_consistency_threshold, 'period_consistency_threshold');
check_unit_interval(options.monotonicity_threshold, 'monotonicity_threshold');
percentile = options.cycle_amplitude_percentile;
if ~(isscalar(percentile) && (isnan(percentile) || ...
        (percentile >= 0 && percentile <= 100)))
    error('delta_episode_options:percentile', ...
        'cycle_amplitude_percentile must be NaN or lie in [0, 100].');
end
end

function preset = cycle_threshold_preset(name)
switch name
    case 'motor_beta'
        preset = struct( ...
            'amp_consistency', 0.3, ...
            'period_consistency', 0.5, ...
            'monotonicity', 0.6, ...
            'amplitude_percentile', 20);
    case 'simulation'
        preset = struct( ...
            'amp_consistency', 0.6, ...
            'period_consistency', 0.6, ...
            'monotonicity', 0.9, ...
            'amplitude_percentile', NaN);
    otherwise
        error('delta_episode_options:preset', ...
            'cycle_threshold_preset must be motor_beta or simulation.');
end
end

function check_unit_interval(value, name)
if ~isscalar(value) || ~isfinite(value) || value < 0 || value > 1
    error('delta_episode_options:threshold', ...
        '%s must lie in [0, 1].', name);
end
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
