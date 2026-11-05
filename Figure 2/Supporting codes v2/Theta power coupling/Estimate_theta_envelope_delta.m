%% Estimate_theta_envelope_delta
% Theta-envelope modulation from 0.1 to 6 Hz, and its alignment to each behavior.
%
% Run add_repo_paths, then run this script. It calls
% GENERATE_THETA_ENVELOPE_DELTA once per Neuropixels folder and saves
% Data/Analysis results/theta envelope delta/theta_envelope_delta.mat.
%
% Then run Analyze_theta_envelope_delta.m from the top, or open
% Analyze_theta_envelope_session.m and Analyze_theta_envelope_onset.m
% and run one %% section at a time. Each plotting section opens one figure.
%
% Knobs live in theta_envelope_delta_options.m.
%   theta_freq, theta_transition_hz
%   modulation_band, bandwidth_hz, test_step_hz, fine_step_hz
%   confirmatory_band, welch_window_seconds, n_perm
%   pre_windows, epoch_window, min_shift_seconds, min_events
%
% The session peak is the 0.2 Hz center with the largest envelope power or
% coherence. Surrogate z-scores use the 0.5 Hz grid. The envelope-spectrum null
% phase-randomizes the LFP and rebuilds the theta envelope. The coherence
% null shifts each Welch segment of the LFP by its own lag. A single lag of
% the whole recording would leave magnitude-squared coherence unchanged.
%
% Each behavior is aligned on its own, for the implanted animal and for the
% partner. Play types are not pooled into a bout. Non-play behaviors are
% included. Theta amplitude is the Hilbert envelope and the 6-12 Hz
% spectrogram power. The primary window is -2 to 0 s. A shift of that
% behavior's onsets is the null. The -1 to 0 s window is saved beside it.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Theta power coupling\Estimate_theta_envelope_delta.m'];
end
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
run(fullfile(repo_root, 'add_repo_paths.m'));

options = theta_envelope_delta_options();
npx_raw_data = fullfile(options.data_root, 'NPX data', 'NPX raw data');
saving_folder = options.output_folder;

% Leave animals_to_run empty to analyze every recording folder.
animals_to_run = {};
animals_to_exclude = {};

%% Recordings
animal_list = dir(npx_raw_data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));
if ~isempty(animals_to_run)
    animal_list = animal_list(ismember({animal_list.name}, animals_to_run));
end
if ~isempty(animals_to_exclude)
    animal_list = animal_list(~ismember({animal_list.name}, animals_to_exclude));
end

fprintf(['Theta-envelope delta scan. %d recordings. ' ...
    'Grid %.1f-%.1f Hz, bandwidth %.1f Hz, test step %.1f Hz, fine step %.1f Hz. ' ...
    '%d surrogates.\n%s\n'], ...
    numel(animal_list), options.modulation_band(1), options.modulation_band(2), ...
    options.bandwidth_hz, options.test_step_hz, options.fine_step_hz, ...
    options.n_perm, options.interpretation);

recordings = [];
animal_names = cell(0, 3);
failures = cell(0, 2);

for recording_n = 1:numel(animal_list)
    animal_name = animal_list(recording_n).name;
    animal_path = fullfile(npx_raw_data, animal_name);
    fprintf('Recording %d/%d: %s\n', recording_n, numel(animal_list), animal_name);
    try
        new_recordings = GENERATE_THETA_ENVELOPE_DELTA(animal_path, options);
    catch exception
        warning('Estimate_theta_envelope_delta:Fail', ...
            '%s failed: %s', animal_name, exception.message);
        failures(end + 1, :) = {animal_name, exception.message}; %#ok<AGROW>
        continue
    end
    if isempty(new_recordings)
        failures(end + 1, :) = {animal_name, 'skipped'}; %#ok<AGROW>
        continue
    end
    if isempty(recordings)
        recordings = new_recordings(:);
    else
        recordings(end + (1:numel(new_recordings)), 1) = new_recordings(:);
    end
    for channel_n = 1:numel(new_recordings)
        animal_names(end + 1, :) = {animal_name, ...
            new_recordings(channel_n).lpag_probe_number, ...
            new_recordings(channel_n).lpag_channel}; %#ok<AGROW>
        n_self = 0;
        n_partner = 0;
        catalog = new_recordings(channel_n).behavior;
        for behavior_n = 1:numel(catalog)
            n_self = n_self + catalog(behavior_n).self.envelope.pre(1).n_valid;
            n_partner = n_partner + catalog(behavior_n).partner.envelope.pre(1).n_valid;
        end
        fprintf(['  probe %g channel %g: envelope peak %.2f Hz, ' ...
            'coherence peak %.2f Hz, %d implanted-animal onsets, ' ...
            '%d partner onsets\n'], ...
            new_recordings(channel_n).lpag_probe_number, ...
            new_recordings(channel_n).lpag_channel, ...
            new_recordings(channel_n).session.peak_power_hz, ...
            new_recordings(channel_n).session.peak_coherence_hz, ...
            n_self, n_partner);
    end
end

%% Save
if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
parameters = options;
output_file = fullfile(saving_folder, options.output_filename);
fprintf('Saving %s\n', output_file);
save(output_file, 'recordings', 'animal_names', 'failures', 'parameters', '-v7.3');
fprintf('Theta-envelope estimation finished. %d channels, %d skipped or failed.\n', ...
    size(animal_names, 1), size(failures, 1));
