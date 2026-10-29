%% Estimate_behavior_probability
% Build binary behavior-presence matrices around selected behavior offsets.
%
% Edit target_behaviors, time_window, and bin_size below. Every
% Session_structure / Partners session is retained separately. Within each
% session:
%   self_behavior    = target events x behavior types x time bins
%   partner_behavior = target events x behavior types x time bins
%
% "Self" is the animal that performed that row's target behavior, even
% when the target performer is the partner rather than the focal animal.

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 1\' ...
        'Behavior probability\Estimate_behavior_probability.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

%% Parameters to select
% Include the two targets that will be compared in Analyze.
target_behaviors = {'Pounce', 'PounceI'};
% Another possible comparison:
% target_behaviors = {'PounceI', 'Bite'};
time_window = [-5 5];
bin_size = 0.01;
% Gray analysis area used by both the window-average and bin-wise tests.
% Time zero is displayed but excluded from inference.
comparison_window = [0 2];

analysis_parameters = behavior_probability_options(struct( ...
    'target_behaviors', {target_behaviors}, ...
    'time_window', time_window, ...
    'bin_size', bin_size, ...
    'comparison_window', comparison_window));

% Leave animals_to_run empty to process every behavior file.
animals_to_run = {};
animals_to_exclude = {};

%% Generate one recording structure, with play sessions kept separate
behavior_files = dir(fullfile(analysis_parameters.behavior_folder, '*.txt'));
if ~isempty(animals_to_run)
    [~, names] = cellfun(@fileparts, ...
        fullfile({behavior_files.folder}, {behavior_files.name}), ...
        'UniformOutput', false);
    behavior_files = behavior_files(ismember(names, animals_to_run));
end
if ~isempty(animals_to_exclude)
    [~, names] = cellfun(@fileparts, ...
        fullfile({behavior_files.folder}, {behavior_files.name}), ...
        'UniformOutput', false);
    behavior_files = behavior_files(~ismember(names, animals_to_exclude));
end

recordings = struct([]);
failures = cell(0, 2);
fprintf(['Behavior probability around %s offsets. %d files, ' ...
    'window [%g %g] s, bin %.4f s.\n'], ...
    strjoin(analysis_parameters.target_behaviors, ' + '), ...
    numel(behavior_files), analysis_parameters.time_window, ...
    analysis_parameters.bin_size);

for file_n = 1:numel(behavior_files)
    behavior_file = fullfile(behavior_files(file_n).folder, ...
        behavior_files(file_n).name);
    fprintf('Recording %d/%d: %s\n', file_n, numel(behavior_files), ...
        behavior_files(file_n).name);
    try
        new_recording = GENERATE_BEHAVIOR_PROBABILITY( ...
            behavior_file, '', analysis_parameters);
    catch exception
        warning('Estimate_behavior_probability:Failed', ...
            '%s failed: %s', behavior_files(file_n).name, ...
            exception.message);
        failures(end + 1, :) = ...
            {behavior_files(file_n).name, exception.message}; %#ok<AGROW>
        continue
    end
    if isempty(recordings)
        recordings = new_recording;
    else
        recordings(end + 1, 1) = new_recording; %#ok<SAGROW>
    end
    n_targets = sum(arrayfun(@(x) size(x.self_behavior, 1), ...
        new_recording.sessions));
    fprintf('  %d play sessions, %d target events\n', ...
        numel(new_recording.sessions), n_targets);
end

%% Save
if ~isfolder(analysis_parameters.output_folder)
    mkdir(analysis_parameters.output_folder);
end
output_file = fullfile(analysis_parameters.output_folder, ...
    analysis_parameters.output_filename);
save(output_file, 'recordings', 'analysis_parameters', 'failures', '-v7.3');
fprintf('Saved %d recordings to %s\n', numel(recordings), output_file);
fprintf('Skipped or failed files: %d\n', size(failures, 1));
