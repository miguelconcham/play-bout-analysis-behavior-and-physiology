%% Estimate_speed_peak_frequency
% Delta and theta peak frequency against locomotion speed, one session at
% a time.
%
% Run add_repo_paths, then run this script. It calls
% GENERATE_SPEED_PEAK_FREQUENCY once per Neuropixels folder. Each partner
% tracking file is split by Session_structure / Partners session rows.
% Each 10 s window is labeled by its median speed. Those speeds are divided
% into equal-count bins of 5-10% of the windows in that session. Welch and
% Fourier spectra in a bin are averaged in linear power and parameterized
% with FOOOF.
%
% Saved to Data/Analysis results/speed peak frequency/:
%   speed_peak_frequency.mat
%       recordings, animal_names, peak_table, correlation_table,
%       analysis_parameters, failures
%   speed_peak_by_session.csv
%   speed_peak_correlations.csv
%
% analysis_parameters is the full settings struct, including the frequency
% grid, FOOOF settings, speed-bin fractions, and band limits.
%
% Then run Analyze_speed_peak_frequency.m. Sessions are not pooled before
% the peak-versus-speed correlation.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Estimate_speed_peak_frequency.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

analysis_parameters = speed_peak_options();
npx_raw_data = fullfile(analysis_parameters.data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(analysis_parameters.data_root, 'Analysis results', ...
    'speed peak frequency');

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

fprintf(['Speed-binned peak frequency. %d recordings. ' ...
    'Bins %.0f-%.0f%% of each session. Windows %.1f s. ' ...
    'Delta %.1f-%.1f Hz, theta %.1f-%.1f Hz.\n'], ...
    numel(animal_list), ...
    100 * analysis_parameters.speed_bin_fraction_limits(1), ...
    100 * analysis_parameters.speed_bin_fraction_limits(2), ...
    analysis_parameters.window_seconds, ...
    analysis_parameters.delta_band(1), analysis_parameters.delta_band(2), ...
    analysis_parameters.theta_band(1), analysis_parameters.theta_band(2));

recordings = struct([]);
animal_names = cell(0, 4);
failures = cell(0, 2);

for recording_n = 1:numel(animal_list)
    animal_name = animal_list(recording_n).name;
    animal_path = fullfile(npx_raw_data, animal_name);
    fprintf('Recording %d/%d: %s\n', recording_n, numel(animal_list), animal_name);
    try
        [new_recordings, analysis_parameters] = ...
            GENERATE_SPEED_PEAK_FREQUENCY(animal_path, analysis_parameters);
    catch exception
        warning('Estimate_speed_peak_frequency:Fail', ...
            '%s failed: %s', animal_name, exception.message);
        failures(end + 1, :) = {animal_name, exception.message}; %#ok<AGROW>
        continue
    end
    if isempty(new_recordings) || all(arrayfun(@(x) isempty(x.sessions), new_recordings))
        failures(end + 1, :) = {animal_name, 'no speed-binned session'}; %#ok<AGROW>
        if ~isempty(new_recordings) && ~isempty(new_recordings(1).notes)
            fprintf('  %s\n', new_recordings(1).notes);
        end
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
            new_recordings(channel_n).lpag_channel, ...
            numel(new_recordings(channel_n).sessions)}; %#ok<AGROW>
        fprintf('  probe %g channel %g: %d sessions\n', ...
            new_recordings(channel_n).lpag_probe_number, ...
            new_recordings(channel_n).lpag_channel, ...
            numel(new_recordings(channel_n).sessions));
    end
end

%% Save
if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
[peak_table, correlation_table] = speed_peak_session_table(recordings);
output_file = fullfile(saving_folder, 'speed_peak_frequency.mat');
fprintf('Saving %s\n', output_file);
save(output_file, 'recordings', 'animal_names', 'peak_table', ...
    'correlation_table', 'analysis_parameters', 'failures', '-v7.3');
writetable(peak_table, fullfile(saving_folder, 'speed_peak_by_session.csv'));
writetable(correlation_table, fullfile(saving_folder, 'speed_peak_correlations.csv'));
if ~isempty(correlation_table)
    disp(correlation_table);
end
fprintf(['Speed peak-frequency estimation finished. %d channels, ' ...
    '%d skipped or failed.\n'], size(animal_names, 1), size(failures, 1));
