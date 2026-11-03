%% Estimate_actor_behavior_probability
% Build an independent actor-relative behavior-probability data set.
%
% Existing files under Figure 1/Behavior probability are not read or
% overwritten. Each play session is split by the animal that performed the
% target. Every target-actor block stores separate self and other matrices.

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 1\' ...
        'Behavior probability by actor\Estimate_actor_behavior_probability.m'];
end
this_dir = fileparts(this_file);
addpath(this_dir);

%% Parameters to select
target_behaviors = {'Pounce', 'PounceI'};
time_window = [-5 5];
bin_size = 0.01;
comparison_window = [0 2];
flank_seconds = 2;
wrapped_bins = 20;

analysis_parameters = actor_behavior_probability_options(struct( ...
    'target_behaviors', {target_behaviors}, ...
    'time_window', time_window, ...
    'bin_size', bin_size, ...
    'comparison_window', comparison_window, ...
    'flank_seconds', flank_seconds, ...
    'wrapped_bins', wrapped_bins));

% Leave empty to process every behavior file.
animals_to_run = {};
animals_to_exclude = {};

%% Generate one recording structure, retaining sessions and target actors
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
fprintf(['Actor-relative behavior probability around %s offsets. %d files, ' ...
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
        new_recording = GENERATE_ACTOR_BEHAVIOR_PROBABILITY( ...
            behavior_file, analysis_parameters);
    catch exception
        warning('Estimate_actor_behavior_probability:Failed', ...
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
    n_actor_blocks = sum(arrayfun( ...
        @(session) numel(session.actor_data), new_recording.sessions));
    n_targets = 0;
    for session_n = 1:numel(new_recording.sessions)
        n_targets = n_targets + sum(arrayfun( ...
            @(actor) size(actor.self_behavior, 1), ...
            new_recording.sessions(session_n).actor_data));
    end
    fprintf('  %d sessions, %d actor blocks, %d target events\n', ...
        numel(new_recording.sessions), n_actor_blocks, n_targets);
end

%% Save to the actor-pipeline output only
if ~isfolder(analysis_parameters.output_folder)
    mkdir(analysis_parameters.output_folder);
end
output_file = fullfile(analysis_parameters.output_folder, ...
    analysis_parameters.output_filename);
save(output_file, 'recordings', 'analysis_parameters', 'failures', '-v7.3');
fprintf('Saved %d recordings to %s\n', numel(recordings), output_file);
fprintf('Skipped or failed files: %d\n', size(failures, 1));
