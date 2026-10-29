%% Analyze_behavior_probability
% Convert binary event matrices into per-session probability estimates.
%
% Run Estimate_behavior_probability.m first. The saved analysis.sessions
% array remains separated by play session. mean_session_* averages those
% session estimates with equal session weight; event_weighted_* instead
% pools all target events.

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 1\' ...
        'Behavior probability\Analyze_behavior_probability.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

default_parameters = behavior_probability_options();
input_file = fullfile(default_parameters.output_folder, ...
    default_parameters.output_filename);
loaded = load(input_file, 'recordings', 'analysis_parameters', 'failures');
if ~isfield(loaded, 'recordings') || isempty(loaded.recordings)
    error('Analyze_behavior_probability:Input', ...
        'No recordings were found in %s.', input_file);
end

analysis = summarize_behavior_probability(loaded.recordings);
analysis_file = fullfile(default_parameters.output_folder, ...
    'behavior_probability_analysis.mat');
save(analysis_file, 'analysis', '-v7.3');

%% Equal-session-weight probability maps
figure_handle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1200 700]);
maps = {analysis.mean_session_self, analysis.mean_session_partner};
map_titles = {'Target performer (self)', 'Other animal (partner)'};
for map_n = 1:2
    subplot(1, 2, map_n);
    imagesc(analysis.time, 1:numel(analysis.behavior_types), maps{map_n});
    axis xy
    hold on
    plot([0 0], [0.5 numel(analysis.behavior_types) + 0.5], ...
        'w--', 'LineWidth', 1);
    yticks(1:numel(analysis.behavior_types));
    yticklabels(analysis.behavior_types);
    xlabel('Time from target behavior offset (s)');
    title(map_titles{map_n});
    caxis([0 1]);
    colorbar
end
sgtitle(sprintf('%s: %d events in %d/%d play sessions', ...
    strjoin(analysis.target_behaviors, ' + '), ...
    analysis.total_target_events, analysis.n_sessions_with_targets, ...
    analysis.n_sessions_total));

figure_base = fullfile(default_parameters.output_folder, ...
    'behavior_probability_equal_session_weight');
savefig(figure_handle, [figure_base, '.fig']);
saveas(figure_handle, [figure_base, '.png']);
close(figure_handle);

fprintf(['Saved per-session behavior probabilities to %s\n' ...
    '%d target events; %d of %d play sessions contained a target.\n'], ...
    analysis_file, analysis.total_target_events, ...
    analysis.n_sessions_with_targets, analysis.n_sessions_total);
