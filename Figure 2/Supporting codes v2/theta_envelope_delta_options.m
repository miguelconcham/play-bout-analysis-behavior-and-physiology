function options = theta_envelope_delta_options(options)
%THETA_ENVELOPE_DELTA_OPTIONS  Settings for theta-envelope delta modulation.
%
% The pipeline asks two questions. The session question is whether the
% Hilbert amplitude of 6-12 Hz LPAG activity oscillates between 0.1 and 6 Hz,
% and whether that oscillation is coherent with the LFP. The onset question
% is whether a 2-3 Hz rise and fall of that amplitude, visible in the average
% aligned to play onset, is present in each animal and in single behaviors.
%
% Frequency grid: centers run from modulation_band(1) to modulation_band(2).
% Each center is an average over bandwidth_hz. test_step_hz is the grid used
% for surrogate tests. fine_step_hz localizes the peak. A 1-2 s pre-play
% window cannot resolve 0.2 Hz, so the onset test uses test_step_hz.

if nargin < 1 || isempty(options)
    options = struct();
end

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];

options = fill_option(options, 'data_root', data_root);
options = fill_option(options, 'sampling_rate', 2500);
options = fill_option(options, 'envelope_rate', 50);
options = fill_option(options, 'n_lfp_channels', 384);
options = fill_option(options, 'area_name', 'LPAG');
options = fill_option(options, 'theta_freq', [6 12]);
options = fill_option(options, 'theta_transition_hz', 1.5);
options = fill_option(options, 'theta_filter_order', []);
options = fill_option(options, 'modulation_band', [0.1 6]);
options = fill_option(options, 'bandwidth_hz', 0.5);
options = fill_option(options, 'test_step_hz', 0.5);
options = fill_option(options, 'fine_step_hz', 0.2);
options = fill_option(options, 'confirmatory_band', [2 3]);
options = fill_option(options, 'welch_window_seconds', 10);
options = fill_option(options, 'welch_overlap', 0.5);
options = fill_option(options, 'welch_resolution_hz', 0.1);
options = fill_option(options, 'n_perm', 500);
options = fill_option(options, 'n_group_flip', 5000);
options = fill_option(options, 'min_shift_seconds', 5);
options = fill_option(options, 'min_events', 15);
options = fill_option(options, 'clean_pre_seconds', 2);
options = fill_option(options, 'pre_windows', [-2 0; -1 0]);
options = fill_option(options, 'epoch_window', [-3 1]);
options = fill_option(options, 'play_behaviors', {'Pounce', 'CC', 'Boxing', ...
    'Evasion', 'Pin', 'Escape', 'CB', 'CD'});
options = fill_option(options, 'interval_mode', 'partners_session');
options = fill_option(options, 'random_seed', NaN);
options = fill_option(options, 'output_folder', fullfile(data_root, ...
    'Analysis results', 'theta envelope delta'));
options = fill_option(options, 'output_filename', 'theta_envelope_delta.mat');
options = fill_option(options, 'interpretation', [ ...
    'Session peaks are the 0.2 Hz center with the largest envelope power or ' ...
    'coherence. Z-scores ask whether that value exceeds the surrogate null. ' ...
    'The 0.5 Hz grid is the surrogate test. Onset spectra of 1-2 s windows ' ...
    'are tests on the 0.5 Hz grid; the 0.2 Hz curve of the average waveform ' ...
    'is descriptive. Coherence nulls circularly shift each Welch segment of ' ...
    'the LFP, because one lag of the whole recording leaves magnitude-squared ' ...
    'coherence unchanged. The envelope-spectrum null phase-randomizes the LFP ' ...
    'and recomputes the theta envelope. Onset nulls shift play times inside ' ...
    'each partner session.']);

if numel(options.theta_freq) ~= 2 || options.theta_freq(1) >= options.theta_freq(2)
    error('theta_envelope_delta_options:theta', ...
        'theta_freq must be [low high] in Hz.');
end
if numel(options.modulation_band) ~= 2 || ...
        options.modulation_band(1) >= options.modulation_band(2) || ...
        options.modulation_band(1) <= 0
    error('theta_envelope_delta_options:band', ...
        'modulation_band must be [low high] in Hz with low > 0.');
end
if options.bandwidth_hz <= 0 || options.test_step_hz <= 0 || options.fine_step_hz <= 0
    error('theta_envelope_delta_options:step', ...
        'Bandwidth and frequency steps must be positive.');
end
if options.envelope_rate <= 2 * options.modulation_band(2) || ...
        options.envelope_rate <= 2 * options.theta_freq(2)
    error('theta_envelope_delta_options:rate', ...
        'envelope_rate must sit above twice the theta and modulation bands.');
end
if options.sampling_rate < options.envelope_rate || ...
        mod(options.sampling_rate, options.envelope_rate) ~= 0
    error('theta_envelope_delta_options:rate', ...
        'sampling_rate must be an integer multiple of envelope_rate.');
end
if ~ismember(options.interval_mode, {'partners_session', 'full_recording'})
    error('theta_envelope_delta_options:interval', ...
        'interval_mode must be partners_session or full_recording.');
end
if size(options.pre_windows, 2) ~= 2 || any(options.pre_windows(:, 1) >= options.pre_windows(:, 2))
    error('theta_envelope_delta_options:window', ...
        'pre_windows must be [start end] rows, in seconds relative to onset.');
end
if numel(options.epoch_window) ~= 2 || options.epoch_window(1) >= options.epoch_window(2)
    error('theta_envelope_delta_options:window', ...
        'epoch_window must be [start end] in seconds.');
end
if options.n_perm < 0 || options.n_perm ~= round(options.n_perm)
    error('theta_envelope_delta_options:perm', 'n_perm must be a non-negative integer.');
end
if options.welch_window_seconds <= 0 || options.welch_overlap < 0 || ...
        options.welch_overlap >= 1 || options.welch_resolution_hz <= 0
    error('theta_envelope_delta_options:welch', ...
        'Welch window, overlap, and resolution are not usable.');
end
if isempty(options.theta_filter_order)
    options.theta_filter_order = 2 * ceil( ...
        (3.3 * options.envelope_rate / options.theta_transition_hz) / 2);
end
if options.theta_filter_order < 8 || mod(options.theta_filter_order, 2) ~= 0
    error('theta_envelope_delta_options:filter', ...
        'theta_filter_order must be an even integer of at least 8.');
end
end

function options = fill_option(options, name, value)
if ~isfield(options, name) || isempty(options.(name))
    options.(name) = value;
end
end
