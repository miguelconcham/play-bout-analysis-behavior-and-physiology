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

options = fill_option(options, 'interpretation', ...
    ['Candidate episodes of elevated 1-5 Hz amplitude in the continuous ' ...
     'LFP. Filtering can produce apparent oscillations around an isolated ' ...
     'transient, so compare each episode with the minimally processed trace. ' ...
     'These episodes are not proof of a neural oscillation.']);

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
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
