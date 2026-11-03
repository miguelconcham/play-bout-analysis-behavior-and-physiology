%% Analyze_theta_envelope_onset
% Theta amplitude aligned to each behavior, not to a play bout.
%
% Run Estimate_theta_envelope_delta.m first. Then run this file from the top,
% or run the sections in order. Each plotting section opens one figure, with
% one panel per behavior.
%
% Play and non-play behaviors are separate panels. The implanted animal and
% the playing partner are separate figures. Hilbert envelope and spectrogram
% power are separate figures.
%
% Gray line: one animal. Black line: mean of animals with at least 15 onsets
% of that behavior. Blue band: -2 to 0 s. Red line: behavior onset.
% The p-value in each panel is the 2-3 Hz onset-shift test for that panel.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Analyze_theta_envelope_onset.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

options = theta_envelope_delta_options();
result_file = fullfile(options.output_folder, options.output_filename);

%% Load
fprintf('Loading %s\n', result_file);
loaded = load(result_file, 'recordings', 'parameters');
recordings = loaded.recordings;
if isempty(recordings)
    error('No recordings in %s.', result_file);
end
if isfield(loaded, 'parameters') && isfield(loaded.parameters, 'min_events')
    min_events = loaded.parameters.min_events;
else
    min_events = 15;
end
grouped = theta_envelope_group_behaviors(recordings, min_events);
if isempty(grouped)
    error('The saved recordings have no per-behavior alignment.');
end

%% Print one test per behavior, actor, and amplitude
fprintf(['Behavior-aligned theta amplitude. %d animals. ' ...
    'Minimum %d onsets.\n'], numel(grouped(1).animals), min_events);
fprintf('%-12s %-10s %-10s %8s %10s %10s\n', ...
    'behavior', 'actor', 'measure', 'onsets', 'locked_p', 'single_p');
for item_n = 1:numel(grouped)
    item = grouped(item_n);
    if strcmp(item.actor, 'self')
        actor_label = 'implanted';
    else
        actor_label = 'partner';
    end
    fprintf('%-12s %-10s %-10s %8d %10.4f %10.4f\n', ...
        item.behavior, actor_label, item.measure, sum(item.n_trials), ...
        item.locked_p, item.loo_p);
end

%% Hilbert envelope, implanted animal
theta_envelope_plot_aligned(grouped, 'self', 'envelope')

%% Hilbert envelope, playing partner
theta_envelope_plot_aligned(grouped, 'partner', 'envelope')

%% Spectrogram power, implanted animal
theta_envelope_plot_aligned(grouped, 'self', 'power')

%% Spectrogram power, playing partner
theta_envelope_plot_aligned(grouped, 'partner', 'power')

%% Save animal table
rows = cell(0, 8);
for item_n = 1:numel(grouped)
    item = grouped(item_n);
    for animal_n = 1:numel(item.animals)
        rows(end + 1, :) = {item.animals{animal_n}, item.behavior, item.actor, ...
            item.measure, item.is_play, item.n_trials(animal_n), ...
            item.included(animal_n), item.animal_locked_p(animal_n)}; %#ok<AGROW>
    end
end
summary = cell2table(rows, 'VariableNames', {'animal', 'behavior', 'actor', ...
    'measure', 'is_play', 'n_onsets', 'included', 'locked_p'});
[folder, ~, ~] = fileparts(result_file);
writetable(summary, fullfile(folder, 'theta_envelope_onset_animals.csv'));
disp(summary)
