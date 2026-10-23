%% Estimate_delta_oscillatory_episodes
% Continuous 1-5 Hz episodes on the mid-LPAG channel of every recording.
%
% Run add_repo_paths, then run this script. It calls
% GENERATE_DELTA_OSCILLATORY_EPISODES once per Neuropixels folder and
% saves Data/Analysis results/delta oscillatory episodes/
% delta_oscillatory_episodes.mat.
%
% Then run Analyze_delta_oscillatory_episodes.m for traces, cycle-period
% histograms, and the session summary.
%
% Knobs live in delta_episode_options.m. The ones that change detection:
%   amplitude_threshold_mode, amplitude_threshold
%   amplitude_smooth_seconds, max_gap_seconds
%   min_cycles, min_duration_seconds, min_peak_distance_seconds
%   edge_margin_seconds, interval_mode, delta_freq, delta_filter_order
%
% interval_mode 'partners_session' (default) uses each
% Session_structure / Partners session row. 'full_recording' uses the
% whole LFP file. This is not a play-bout average.
%
% Detected episodes are candidates. The 1-5 Hz filter can ring around a
% transient. Keep the minimally processed trace in view when judging them.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Estimate_delta_oscillatory_episodes.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

options = delta_episode_options();
options.filter_coefficients = delta_bandpass_coefficients(options);
npx_raw_data = fullfile(options.data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(options.data_root, 'Analysis results', ...
    'delta oscillatory episodes');

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

fprintf(['Delta-episode detection. %d recordings. Mode %s. ' ...
    'Threshold %s = %g. Minimum %g cycles.\n%s\n'], ...
    numel(animal_list), options.interval_mode, ...
    options.amplitude_threshold_mode, options.amplitude_threshold, ...
    options.min_cycles, options.interpretation);

delta_episodes = struct([]);
animal_names = cell(0, 3);
failures = cell(0, 2);

for recording_n = 1:numel(animal_list)
    animal_name = animal_list(recording_n).name;
    animal_path = fullfile(npx_raw_data, animal_name);
    fprintf('Recording %d/%d: %s\n', recording_n, numel(animal_list), animal_name);
    try
        new_recordings = GENERATE_DELTA_OSCILLATORY_EPISODES(animal_path, options);
    catch exception
        warning('Estimate_delta_oscillatory_episodes:Fail', ...
            '%s failed: %s', animal_name, exception.message);
        failures(end + 1, :) = {animal_name, exception.message}; %#ok<AGROW>
        continue
    end
    if isempty(new_recordings)
        failures(end + 1, :) = {animal_name, 'skipped'}; %#ok<AGROW>
        continue
    end

    if isempty(delta_episodes)
        delta_episodes = new_recordings(:);
    else
        delta_episodes(end + (1:numel(new_recordings)), 1) = new_recordings(:);
    end
    for channel_n = 1:numel(new_recordings)
        animal_names(end + 1, :) = {animal_name, ...
            new_recordings(channel_n).lpag_probe_number, ...
            new_recordings(channel_n).lpag_channel}; %#ok<AGROW>
        channel_summary = new_recordings(channel_n).summary;
        fprintf(['  probe %g channel %g: %d episodes, ' ...
            'median duration %.2f s, median period %.0f ms, ' ...
            'median period IQR %.0f ms\n'], ...
            channel_summary.probe_number, channel_summary.lpag_channel, ...
            channel_summary.n_episodes, channel_summary.median_duration_s, ...
            1000 * channel_summary.median_cycle_period_s, ...
            1000 * channel_summary.median_cycle_period_iqr_s);
    end
end

%% Save
if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
if isempty(delta_episodes)
    session_summary = table();
else
    session_summary = delta_episode_summary_table(delta_episodes);
end
detection_parameters = options;
output_file = fullfile(saving_folder, 'delta_oscillatory_episodes.mat');
fprintf('Saving %s\n', output_file);
save(output_file, 'delta_episodes', 'animal_names', 'session_summary', ...
    'failures', 'detection_parameters', '-v7.3');
if ~isempty(session_summary)
    writetable(session_summary, fullfile(saving_folder, 'session_summary.csv'));
    disp(session_summary);
end
fprintf('Delta-episode estimation finished. %d channels, %d skipped or failed.\n', ...
    size(animal_names, 1), size(failures, 1));
