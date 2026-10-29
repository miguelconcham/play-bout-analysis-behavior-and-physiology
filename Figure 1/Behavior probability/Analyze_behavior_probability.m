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

%% Compare two target behaviors after pooling all sessions
% Both names must have been included in target_behaviors when running
% Estimate_behavior_probability.m.
behavior_pair = {'Pounce', 'PounceI'};
comparison_options = struct();
% Choose 'self', 'partner', or 'either'. 'either' is present when the
% behavior was performed by the target animal or its partner.
comparison_options.actor_relation = 'either';
% The permutation statistic is the difference in mean occupancy across
% this offset-relative window. Full unsmoothed time courses are retained.
comparison_options.test_window = [0 2];
comparison_options.n_permutations = 10000;
comparison_options.random_seed = 1;
plot_smoothing_seconds = 0.10;

comparison = compare_behavior_pair_probability( ...
    loaded.recordings, behavior_pair, comparison_options);
comparison_table = table(comparison.group_names(:), ...
    comparison.observed_difference(:), comparison.p_value(:), ...
    comparison.p_value_holm(:), ...
    'VariableNames', {'BehaviorGroup', 'ProbabilityDifference', ...
    'PermutationP', 'HolmP'});
disp(comparison_table);

figure_handle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 900 1000]);
colors = lines(2);
smooth_bins = max(1, round(plot_smoothing_seconds / analysis.bin_size));
for group_n = 1:numel(comparison.group_names)
    subplot(3, 1, group_n);
    hold on
    patch(comparison.test_window([1 2 2 1]), [0 0 1 1], ...
        [0.85 0.85 0.85], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.35, 'HandleVisibility', 'off');
    for target_n = 1:2
        curve = reshape(comparison.probability(group_n, target_n, :), ...
            1, []);
        curve = movmean(curve, smooth_bins);
        plot(comparison.time, curve, 'Color', colors(target_n, :), ...
            'LineWidth', 1.8, 'DisplayName', sprintf('%s (N=%d)', ...
            comparison.behavior_pair{target_n}, ...
            comparison.n_events_plot(target_n)));
    end
    plot([0 0], [0 1], 'k--', 'HandleVisibility', 'off');
    xlim(comparison.time([1 end]));
    ylim([0 1]);
    ylabel('Probability');
    title(sprintf('%s: permutation p = %.4g, Holm p = %.4g', ...
        comparison.group_names{group_n}, comparison.p_value(group_n), ...
        comparison.p_value_holm(group_n)));
    legend('Location', 'best');
end
xlabel('Time from target behavior offset (s)');
sgtitle(sprintf(['Grouped behavior probability (%s); test [%g, %g] s ' ...
    '(strictly post-offset); ' ...
    '%d shared sessions; %.2f s display smoothing'], ...
    comparison.actor_relation, comparison.test_time([1 end]), ...
    comparison.n_eligible_sessions, plot_smoothing_seconds));

pair_tag = sprintf('%s_vs_%s', comparison.behavior_pair{:});
figure_base = fullfile(default_parameters.output_folder, ...
    ['behavior_group_probability_', pair_tag]);
savefig(figure_handle, [figure_base, '.fig']);
saveas(figure_handle, [figure_base, '.png']);
close(figure_handle);

save(analysis_file, 'comparison', 'comparison_table', '-append');
fprintf(['Saved per-session behavior probabilities to %s\n' ...
    '%d target events; %d of %d play sessions contained a target.\n'], ...
    analysis_file, analysis.total_target_events, ...
    analysis.n_sessions_with_targets, analysis.n_sessions_total);
