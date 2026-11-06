%% Estimate_broadband_cycle_frequencies
% Trough-to-trough frequencies on a 40 Hz low-pass of every mid-LPAG channel.
%
% Run add_repo_paths, then run this script. It saves
% Data/Analysis results/broadband cycle frequencies/
% broadband_cycle_frequencies.mat.
%
% Then run Analyze_broadband_cycle_frequencies.m. Set figure_folder at the
% top of that script before running it. The analysis script writes SVG and
% JPEG figures and result_summary.txt.
%
% This pipeline does not apply a delta or theta bandpass. A cycle is one
% trough-to-trough interval on the low-pass waveform.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes\Broadband cycle frequencies\' ...
        'Estimate_broadband_cycle_frequencies.m'];
end
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
run(fullfile(repo_root, 'add_repo_paths.m'));

options = broadband_cycle_options();
options.lowpass_coefficients = broadband_lowpass_coefficients(options);
npx_raw_data = fullfile(options.data_root, 'NPX data', 'NPX raw data');
saving_folder = options.output_folder;
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

fprintf(['Broadband cycle frequencies. %d recordings. Low-pass %.4g Hz. ' ...
    'Accepted cycles %.4g-%.4g Hz.\n'], numel(animal_list), ...
    options.lowpass_hz, options.frequency_limits_hz(1), ...
    options.frequency_limits_hz(2));

recordings = struct([]);
failures = cell(0, 2);
for recording_n = 1:numel(animal_list)
    animal_name = animal_list(recording_n).name;
    fprintf('Recording %d/%d: %s\n', recording_n, numel(animal_list), animal_name);
    try
        new_recordings = GENERATE_BROADBAND_CYCLE_FREQUENCIES( ...
            fullfile(npx_raw_data, animal_name), options);
    catch exception
        warning('Estimate_broadband_cycle_frequencies:Fail', ...
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
        fprintf('  probe %g channel %g: %d accepted of %d candidate cycles\n', ...
            new_recordings(channel_n).lpag_probe_number, ...
            new_recordings(channel_n).lpag_channel, ...
            new_recordings(channel_n).n_accepted_cycles, ...
            new_recordings(channel_n).n_candidate_cycles);
    end
end

%% Save
if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
detection_parameters = options;
output_file = fullfile(saving_folder, options.output_filename);
fprintf('Saving %s\n', output_file);
save(output_file, 'recordings', 'failures', 'detection_parameters', '-v7.3');
if ~isempty(recordings)
    write_broadband_cycle_summary_text(recordings, ...
        fullfile(saving_folder, 'result_summary.txt'));
end
fprintf('Broadband cycle-frequency estimation finished. %d channels, %d failures.\n', ...
    numel(recordings), size(failures, 1));
