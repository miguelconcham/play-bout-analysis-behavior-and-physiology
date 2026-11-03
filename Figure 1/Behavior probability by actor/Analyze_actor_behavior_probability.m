%% Analyze_actor_behavior_probability
% Analyze the independent target-actor/self-other estimate.
%
% Run Estimate_actor_behavior_probability.m first. Outputs from the older
% behavior-probability pipeline are neither loaded nor overwritten.

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 1\' ...
        'Behavior probability by actor\Analyze_actor_behavior_probability.m'];
end
this_dir = fileparts(this_file);
addpath(this_dir);

parameters = actor_behavior_probability_options();
input_file = fullfile(parameters.output_folder, parameters.output_filename);
loaded = load(input_file, 'recordings', 'analysis_parameters', 'failures');
if ~isfield(loaded, 'recordings') || isempty(loaded.recordings)
    error('Analyze_actor_behavior_probability:Input', ...
        'No actor-relative recordings were found in %s.', input_file);
end
if isfield(loaded, 'analysis_parameters')
    parameters = loaded.analysis_parameters;
end

analysis = summarize_actor_behavior_probability(loaded.recordings);
analysis_file = fullfile(parameters.output_folder, ...
    'actor_behavior_probability_analysis.mat');
save(analysis_file, 'analysis', '-v7.3');

%% Equal actor-session-weight heatmaps
for target_n = 1:numel(analysis.target_behaviors)
    figure_handle = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 1200 700]);
    maps = { ...
        reshape(analysis.mean_stratum_self(:, target_n, :), ...
        numel(analysis.behavior_types), []), ...
        reshape(analysis.mean_stratum_other(:, target_n, :), ...
        numel(analysis.behavior_types), [])};
    map_titles = {'Target performer (self)', 'Any other animal'};
    for relation_n = 1:2
        subplot(1, 2, relation_n);
        imagesc(analysis.time, 1:numel(analysis.behavior_types), ...
            maps{relation_n});
        axis xy
        hold on
        plot([0 0], [0.5 numel(analysis.behavior_types) + 0.5], ...
            'w--', 'LineWidth', 1);
        yticks(1:numel(analysis.behavior_types));
        yticklabels(analysis.behavior_types);
        xlabel('Time from target behavior offset (s)');
        title(map_titles{relation_n});
        caxis([0 1]);
        colorbar
    end
    sgtitle(sprintf('%s: %d events in %d actor-session strata', ...
        analysis.target_behaviors{target_n}, analysis.n_events(target_n), ...
        analysis.n_strata_with_target(target_n)));
    target_tag = lower(regexprep(analysis.target_behaviors{target_n}, ...
        '[^A-Za-z0-9]+', '_'));
    figure_base = fullfile(parameters.output_folder, ...
        ['actor_behavior_probability_', target_tag]);
    savefig(figure_handle, [figure_base, '.fig']);
    saveas(figure_handle, [figure_base, '.png']);
    close(figure_handle);
end

%% Compare two target behaviors separately for self and other
if numel(parameters.target_behaviors) < 2
    error('Analyze_actor_behavior_probability:TargetPair', ...
        'Estimate at least two target behaviors for the pair comparison.');
end
behavior_pair = parameters.target_behaviors(1:2);
relations = {'self', 'other'};
pair_comparisons = struct();
pair_summary_tables = cell(1, numel(relations));
pair_bin_tables = cell(1, numel(relations));
for relation_n = 1:numel(relations)
    relation = relations{relation_n};
    settings = struct('relation', relation, ...
        'test_window', parameters.comparison_window, ...
        'n_permutations', 10000, 'random_seed', 1);
    comparison = compare_actor_behavior_pair_probability( ...
        loaded.recordings, behavior_pair, settings);
    pair_comparisons.(relation) = comparison;
    pair_summary_tables{relation_n} = pair_summary_table(comparison);
    pair_bin_tables{relation_n} = pair_bin_table(comparison);
    plot_pair_comparison(comparison, parameters.output_folder);
end
pair_summary_table_all = vertcat(pair_summary_tables{:});
pair_bin_table_all = vertcat(pair_bin_tables{:});
writetable(pair_summary_table_all, fullfile(parameters.output_folder, ...
    'actor_target_pair_comparison.csv'));
writetable(pair_bin_table_all, fullfile(parameters.output_folder, ...
    'actor_target_pair_comparison_binwise.csv'));

%% Compare post-offset bins with pre-onset baseline for each relation
baseline_comparisons = struct();
baseline_summary_tables = cell(1, numel(relations));
baseline_bin_tables = cell(1, numel(relations));
for relation_n = 1:numel(relations)
    relation = relations{relation_n};
    settings = struct('relation', relation, ...
        'pre_average_window', [], ...
        'n_permutations', 10000, 'random_seed', 1);
    baseline = compare_actor_behavior_baseline_probability( ...
        loaded.recordings, parameters.target_behaviors, settings);
    baseline_comparisons.(relation) = baseline;
    baseline_summary_tables{relation_n} = baseline_summary_table(baseline);
    baseline_bin_tables{relation_n} = baseline_bin_table(baseline);
    plot_baseline_comparison(baseline, parameters.output_folder);
end
baseline_summary_table_all = vertcat(baseline_summary_tables{:});
baseline_bin_table_all = vertcat(baseline_bin_tables{:});
writetable(baseline_summary_table_all, fullfile(parameters.output_folder, ...
    'actor_baseline_comparison.csv'));
writetable(baseline_bin_table_all, fullfile(parameters.output_folder, ...
    'actor_baseline_comparison_binwise.csv'));

save(analysis_file, 'pair_comparisons', 'pair_summary_table_all', ...
    'pair_bin_table_all', 'baseline_comparisons', ...
    'baseline_summary_table_all', 'baseline_bin_table_all', '-append');
fprintf(['Saved actor-relative probabilities and statistics to %s\n' ...
    '%d target events across %d actor-session strata.\n'], ...
    analysis_file, analysis.total_target_events, ...
    analysis.n_actor_session_strata);

function output = pair_summary_table(comparison)
n_groups = numel(comparison.group_names);
output = table(repmat({comparison.relation}, n_groups, 1), ...
    comparison.group_names(:), comparison.observed_difference(:), ...
    comparison.p_value(:), comparison.p_value_holm(:), ...
    repmat(comparison.n_eligible_strata, n_groups, 1), ...
    'VariableNames', {'Relation', 'BehaviorGroup', ...
    'ProbabilityDifference', 'PermutationP', 'HolmP', ...
    'NActorSessionStrata'});
end

function output = pair_bin_table(comparison)
n_groups = numel(comparison.group_names);
n_bins = numel(comparison.test_time);
output = table(repmat({comparison.relation}, n_groups * n_bins, 1), ...
    repmat(comparison.group_names(:), n_bins, 1), ...
    repelem(comparison.test_time(:), n_groups), ...
    comparison.observed_difference_by_bin(:), ...
    comparison.p_value_by_bin(:), ...
    comparison.p_value_by_bin_fwer_group(:), ...
    comparison.p_value_by_bin_fwer_all(:), ...
    'VariableNames', {'Relation', 'BehaviorGroup', 'Time', ...
    'ProbabilityDifference', 'PermutationP', ...
    'GroupMaxTP', 'GlobalMaxTP'});
end

function plot_pair_comparison(comparison, output_folder)
figure_handle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 900 1000]);
colors = lines(2);
for group_n = 1:numel(comparison.group_names)
    subplot(numel(comparison.group_names), 1, group_n);
    hold on
    patch(comparison.test_window([1 2 2 1]), [0 0 1 1], ...
        [0.85 0.85 0.85], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.35, 'HandleVisibility', 'off');
    for target_n = 1:2
        curve = reshape(comparison.probability(group_n, target_n, :), 1, []);
        plot(comparison.time, curve, 'Color', colors(target_n, :), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('%s (N=%d)', ...
            comparison.behavior_pair{target_n}, ...
            comparison.n_events_test(target_n)));
    end
    significant = comparison.p_value_by_bin_fwer_all(group_n, :) < 0.05;
    if any(significant)
        plot(comparison.test_time(significant), ...
            0.97 * ones(1, sum(significant)), 'k.', ...
            'MarkerSize', 9, 'HandleVisibility', 'off');
    end
    plot([0 0], [0 1], 'k--', 'HandleVisibility', 'off');
    xlim(comparison.time([1 end]));
    ylim([0 1]);
    ylabel('Probability');
    title(sprintf('%s: Holm p = %.4g', ...
        comparison.group_names{group_n}, comparison.p_value_holm(group_n)));
    legend('Location', 'best');
end
xlabel('Time from target behavior offset (s)');
sgtitle(sprintf('%s behavior; %d paired actor-session strata', ...
    comparison.relation, comparison.n_eligible_strata));
pair_tag = sprintf('%s_vs_%s', comparison.behavior_pair{:});
base = fullfile(output_folder, sprintf( ...
    'actor_target_pair_%s_%s', pair_tag, comparison.relation));
savefig(figure_handle, [base, '.fig']);
saveas(figure_handle, [base, '.png']);
close(figure_handle);
end

function output = baseline_summary_table(result)
n_groups = numel(result.group_names);
n_targets = numel(result.target_names);
relation = cell(n_groups * n_targets, 1);
target = cell(size(relation));
group = cell(size(relation));
row_n = 0;
for target_n = 1:n_targets
    for group_n = 1:n_groups
        row_n = row_n + 1;
        relation{row_n} = result.relation;
        target{row_n} = result.target_names{target_n};
        group{row_n} = result.group_names{group_n};
    end
end
output = table(relation, target, group, ...
    result.pre_probability(:), result.post_probability(:), ...
    result.observed_difference(:), result.p_value(:), ...
    result.p_value_holm(:), ...
    repelem(result.n_events_test(:), n_groups), ...
    repelem(result.n_actor_session_strata(:), n_groups), ...
    'VariableNames', {'Relation', 'Target', 'BehaviorGroup', ...
    'BaselineProbability', 'PostProbability', 'PostMinusBaseline', ...
    'PermutationP', 'HolmP', 'NEvents', 'NActorSessionStrata'});
end

function output = baseline_bin_table(result)
n_groups = numel(result.group_names);
n_targets = numel(result.target_names);
n_bins = numel(result.post_time);
n_rows = n_groups * n_targets * n_bins;
relation = cell(n_rows, 1);
target = cell(n_rows, 1);
group = cell(n_rows, 1);
time = zeros(n_rows, 1);
difference = zeros(n_rows, 1);
p = zeros(n_rows, 1);
group_p = zeros(n_rows, 1);
target_p = zeros(n_rows, 1);
row_n = 0;
for target_n = 1:n_targets
    for bin_n = 1:n_bins
        for group_n = 1:n_groups
            row_n = row_n + 1;
            relation{row_n} = result.relation;
            target{row_n} = result.target_names{target_n};
            group{row_n} = result.group_names{group_n};
            time(row_n) = result.post_time(bin_n);
            difference(row_n) = ...
                result.observed_difference_by_bin(group_n, bin_n, target_n);
            p(row_n) = result.p_value_by_bin(group_n, bin_n, target_n);
            group_p(row_n) = ...
                result.p_value_by_bin_fwer_group(group_n, bin_n, target_n);
            target_p(row_n) = ...
                result.p_value_by_bin_fwer_target(group_n, bin_n, target_n);
        end
    end
end
output = table(relation, target, group, time, difference, p, ...
    group_p, target_p, 'VariableNames', {'Relation', 'Target', ...
    'BehaviorGroup', 'TimeAfterOffset', 'PostMinusBaseline', ...
    'PermutationP', 'GroupMaxTP', 'TargetMaxTP'});
end

function plot_baseline_comparison(result, output_folder)
n_groups = numel(result.group_names);
n_targets = numel(result.target_names);
figure_handle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [80 80 1100 900]);
for group_n = 1:n_groups
    for target_n = 1:n_targets
        subplot(n_groups, n_targets, ...
            (group_n - 1) * n_targets + target_n);
        hold on
        curve = reshape(result.probability(group_n, target_n, :), 1, []);
        plot(result.wrapped_plot_time, curve, 'LineWidth', 1.5);
        yline(result.pre_probability(group_n, target_n), ':', ...
            'Baseline');
        offset_x = result.wrapped_plot_time(result.wrap_index(end));
        plot([0 0], [0 1], 'k--');
        plot([offset_x offset_x], [0 1], 'k--');
        significant = result.p_value_by_bin_fwer_target( ...
            group_n, :, target_n) < 0.05;
        if any(significant)
            plot(offset_x + result.post_time(significant), ...
                0.97 * ones(1, sum(significant)), 'k.', 'MarkerSize', 8);
        end
        ylim([0 1]);
        title(sprintf('%s / %s; Holm p=%.4g', ...
            result.target_names{target_n}, result.group_names{group_n}, ...
            result.p_value_holm(group_n, target_n)));
        if target_n == 1
            ylabel('Probability');
        end
    end
end
xlabel('Wrapped time (dashed lines mark target onset and offset)');
sgtitle(sprintf('%s behavior: post-offset versus pre-onset baseline', ...
    result.relation));
base = fullfile(output_folder, ...
    ['actor_baseline_probability_', result.relation]);
savefig(figure_handle, [base, '.fig']);
saveas(figure_handle, [base, '.png']);
close(figure_handle);
end
