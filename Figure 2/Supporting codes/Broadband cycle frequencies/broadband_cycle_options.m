function options = broadband_cycle_options(options)
%BROADBAND_CYCLE_OPTIONS  Settings for low-pass cycle frequencies.
%
% Cycles are trough-to-trough intervals on a broad low-pass. No delta or
% theta bandpass decides where a cycle starts. Tiny ripples are removed
% by the consistency and amplitude-percentile gates.

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
options = fill_option(options, 'lowpass_hz', 40);
options = fill_option(options, 'lowpass_order', round(0.1 * options.sampling_rate));
options = fill_option(options, 'lowpass_coefficients', []);
options = fill_option(options, 'frequency_limits_hz', [0.5 20]);
options = fill_option(options, 'histogram_edges_hz', 0.5:0.5:20);
options = fill_option(options, 'amplitude_consistency_threshold', 0.3);
options = fill_option(options, 'period_consistency_threshold', 0.5);
options = fill_option(options, 'monotonicity_threshold', 0.6);
options = fill_option(options, 'amplitude_percentile', 20);
options = fill_option(options, 'interval_mode', 'partners_session');
options = fill_option(options, 'example_half_seconds', 1);
options = fill_option(options, 'output_folder', fullfile(data_root, ...
    'Analysis results', 'broadband cycle frequencies'));
options = fill_option(options, 'output_filename', 'broadband_cycle_frequencies.mat');

if options.lowpass_hz <= 0 || options.lowpass_hz >= options.sampling_rate / 2
    error('broadband_cycle_options:lowpass', ...
        'lowpass_hz must lie between 0 and the Nyquist frequency.');
end
if numel(options.frequency_limits_hz) ~= 2 || ...
        options.frequency_limits_hz(1) <= 0 || ...
        options.frequency_limits_hz(1) >= options.frequency_limits_hz(2)
    error('broadband_cycle_options:frequency', ...
        'frequency_limits_hz must be [low high] with low > 0.');
end
if ~ismember(options.interval_mode, {'partners_session', 'full_recording'})
    error('broadband_cycle_options:interval', ...
        'interval_mode must be partners_session or full_recording.');
end
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
