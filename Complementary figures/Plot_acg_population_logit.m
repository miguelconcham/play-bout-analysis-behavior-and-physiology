% Plot_acg_population_logit
% All neurons together, on a logit scale, from a saved ACG shuffle.
%
% Overlap is the standard method. Section 1 loads
% locking_period_acg_shuffle_n500.mat and builds the matrices and
% feature_table. Section 2 runs PCA, t-SNE, and the two-group split.
% Section 3 keeps cluster 1 of cluster_method, the higher-modulation group.
% cluster_method is 'ward' or 'kmeans'. Ward is the default.
% Section 4 only draws. After section 3 has run, edit section 4 and run
% that section alone.
%
% output_folder is the first line of section 1. Change that path to choose
% where this script writes each figure as jpg and svg, and where it writes
% analysis_summary.txt. Later sections keep that folder when it is already
% in the workspace.
% n_example_neurons and n_example_subplots are set at the start of section 4.
% n_example_neurons chooses how many of the strongest oscillatory neurons
% are drawn. The examples also keep the strongest oscillatory neuron from
% each analysed area except VLPAG, so that count can be larger than
% n_example_neurons. n_example_subplots is how many panels go on each figure.
% Each neuron uses two panels, so the figures continue until every
% requested neuron is drawn.
%
% Every matrix shares one row order. Row i of imagesc_logit.sum.all is
% the same neuron as row i of imagesc_logit.sum.power_high, power_low,
% and the rate-ratio matrices. neuron_type{i} is that neuron's lock group:
% 'trough', 'peak', 'nonentrained' (unlocked), or '' when it was not
% classified. sort_score.* is the mean percentile from 200 to 400 ms.
% Row 1 is the highest all-runs sum-to-1 score.
%
% Section 2 builds ACG features from the high-power sum-to-1 logit and
% splits the neurons with Ward and k-means. cluster_method chooses which
% split is kept. Ward is the default. Cluster 1 of that split is the
% higher-modulation group.
% shown_logit is that group, strongest modulation first.
% imagesc_osc_sum_all is the same sum-to-1 image under a plain name.
% osc_neuron_type is its lock group. feature_table has one row per neuron.

%% 1. Load the checkpoint and build the imagesc matrices
output_folder = fullfile(fileparts(mfilename('fullpath')), 'acg_population_logit_figures');
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
addpath(this_dir);
addpath(repo_root);
if exist('align_population_percentiles', 'file') ~= 2
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

acg_method = 'overlap';
n_shuffles_default = 500;
% Mean percentile inside this lag window sets the stored row order.
% 0.2 s is 200 ms. Section 2 sorts the figure again.
sort_window_s = [0.200 0.400];

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\' ...
    'NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data'];
parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
checkpoint = fullfile(parsing_folder, shuffle_checkpoint_filename(n_shuffles_default, acg_method));
if exist(checkpoint, 'file') ~= 2
    error('Plot_acg_population_logit:checkpoint', ...
        'No shuffle checkpoint at %s.', checkpoint);
end

stored = load(checkpoint, 'shuffle_records', 'parameters');
if ~isfield(stored, 'shuffle_records') || ~isfield(stored, 'parameters') ...
        || ~isfield(stored.parameters, 'n_shuffles')
    error('Plot_acg_population_logit:checkpoint', ...
        '%s has no shuffle records.', checkpoint);
end
records = stored.shuffle_records;
if ~iscell(records)
    records = num2cell(records(:));
end
records = records(:);
n_shuffles = stored.parameters.n_shuffles;

column_id = {'all', 'power_high', 'power_low'};
column_title = {'All neurons', 'All neurons, high 25% power', 'All neurons, low 25% power'};
scale_id = {'sum', 'ratio'};
scale_title = {'sum-to-1', 'rate ratio (Poisson = 1)'};

[lag, percentile, neuron_type, neuron_record_index] = align_population_percentiles(records);
if isempty(neuron_type)
    error('Plot_acg_population_logit:neurons', ...
        'No neuron had all-runs, high-power, and low-power curves on one lag axis.');
end
[percentile.sum.all, imagesc_sum_all, sort_score_all, row_order] = order_percentile_rows( ...
    percentile.sum.all, lag, sort_window_s, n_shuffles);
neuron_type = neuron_type(row_order);
neuron_record_index = neuron_record_index(row_order);
[neuron_session, neuron_area, neuron_cluster_id] = ...
    neuron_identity_from_records(records, neuron_record_index);

imagesc_logit = struct();
sort_score = struct();
imagesc_logit.sum.all = imagesc_sum_all;
sort_score.sum.all = sort_score_all;
in_sort_window = lag >= sort_window_s(1) & lag <= sort_window_s(2);
for s = 1:numel(scale_id)
    for c = 1:numel(column_id)
        if s == 1 && c == 1
            continue
        end
        percentile.(scale_id{s}).(column_id{c}) = ...
            percentile.(scale_id{s}).(column_id{c})(row_order, :);
        imagesc_logit.(scale_id{s}).(column_id{c}) = correct_percentile_logit( ...
            percentile.(scale_id{s}).(column_id{c}), n_shuffles);
        sort_score.(scale_id{s}).(column_id{c}) = mean( ...
            percentile.(scale_id{s}).(column_id{c})(:, in_sort_window), 2, 'omitnan');
    end
end

% Plain names for the sum-to-1 imagesc matrices. Row order matches neuron_type.
imagesc_sum_all = imagesc_logit.sum.all;
imagesc_sum_power_high = imagesc_logit.sum.power_high;
imagesc_sum_power_low = imagesc_logit.sum.power_low;

prep_file = fullfile(parsing_folder, 'acg_population_logit_prep.mat');
% acg2features scores the lags it is given. The checkpoint axis is
% two-sided, and the rhythmicity function uses the positive side.
positive_lag = lag > 0;
feature_table = acg2features(imagesc_sum_power_high(:, positive_lag), lag(positive_lag));

save(prep_file, 'lag', 'imagesc_logit', 'percentile', 'sort_score', ...
    'imagesc_sum_all', 'imagesc_sum_power_high', 'imagesc_sum_power_low', ...
    'neuron_type', 'neuron_record_index', 'neuron_session', 'neuron_area', ...
    'neuron_cluster_id', 'feature_table', ...
    'column_id', 'column_title', 'scale_id', 'scale_title', ...
    'n_shuffles', 'sort_window_s', 'checkpoint', 'prep_file');
disp(['Section 1 saved the matrices to ', prep_file])
disp('imagesc_sum_all is neurons by lag. neuron_type is the lock group of each row.')
disp('feature_table is one row of ACG features per neuron. Run section 2 to classify.')
if exist(output_folder, 'dir') ~= 7
    mkdir(output_folder);
end
write_summary_section(output_folder, '1. Checkpoint and population matrices', strjoin({ ...
    ['Output folder: ', output_folder], ...
    ['Checkpoint: ', checkpoint], ...
    sprintf('Shuffles: %d', n_shuffles), ...
    sprintf('Neurons with all, high-power, and low-power curves: %d', numel(neuron_type)), ...
    sprintf('Lag axis: %.3f to %.3f s, %d bins', lag(1), lag(end), numel(lag)), ...
    sprintf('Sort window for the stored row order: %.3f to %.3f s', ...
        sort_window_s(1), sort_window_s(2)), ...
    ['Feature columns from the high-power sum-to-1 logit, positive lags: ', ...
        strjoin(feature_table.Properties.VariableNames, ', ')], ...
    ['Prep file: ', prep_file]}, newline))

%% 2. ACG features, PCA, t-SNE, and the two-group split
% Re-run this section after changing the feature list or cluster_method.
% Section 1 does not need to run again while feature_table is still in
% the workspace. cluster_method is 'ward' by default. Set it to 'kmeans'
% to keep the k-means split instead. Cluster 1 of that split is the
% higher-modulation group, and section 3 uses it. Both
% splits are still drawn. The violins after the split also show
% frequency, amplitude, and the peak-to-peak drop. Those extra columns
% do not enter the split.
cluster_method = 'ward';
if ~exist('output_folder', 'var') || isempty(output_folder)
    output_folder = fullfile(fileparts(mfilename('fullpath')), 'acg_population_logit_figures');
end
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
if ~exist('feature_table', 'var') || ~exist('neuron_type', 'var') ...
        || ~exist('imagesc_sum_power_high', 'var')
    data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\' ...
        'NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data'];
    parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
    prep_file = fullfile(parsing_folder, 'acg_population_logit_prep.mat');
    load(prep_file)
end
if ~exist('feature_table', 'var')
    positive_lag = lag > 0;
    feature_table = acg2features(imagesc_sum_power_high(:, positive_lag), lag(positive_lag));
end

rng(42);

candidate_features = { ...
    'Modulation1', 'Modulation2', 'PeakInterval', 'TroughInterval', ...
    'OscillationScore', 'Baseline'};
feature_names = candidate_features( ...
    ismember(candidate_features, feature_table.Properties.VariableNames));
if isempty(feature_names)
    error('Plot_acg_population_logit:features', ...
        'No candidate features were found in feature_table.');
end
X_raw = feature_table{:, feature_names};
if ~isnumeric(X_raw)
    error('Plot_acg_population_logit:features', 'Selected features must be numeric.');
end
fprintf('Selected features:\n');
disp(feature_names(:)')

% K-means below numbers only these complete rows. Section 3 maps them
% back onto the image rows.
valid_rows = all(isfinite(X_raw), 2);
X = X_raw(valid_rows, :);
fprintf('Total neurons: %d\n', size(X_raw, 1));
fprintf('Neurons retained: %d\n', size(X, 1));
fprintf('Neurons excluded: %d\n', sum(~valid_rows));
if size(X, 1) < 10
    error('Plot_acg_population_logit:features', ...
        'Too few complete neurons for a meaningful PCA or t-SNE analysis.');
end

feature_sd = std(X, 0, 1);
keep_features = feature_sd > eps(max(abs(X(:))));
X = X(:, keep_features);
feature_names = feature_names(keep_features);
if numel(feature_names) < 2
    error('Plot_acg_population_logit:features', ...
        'Fewer than two variable features remain.');
end
fprintf('Features retained after variance check:\n');
disp(feature_names(:)')

feature_correlation = corr(X, 'Rows', 'complete');
open_full_figure('ACG feature correlations')
heatmap(feature_names, feature_names, feature_correlation, 'ColorLimits', [-1 1])
title('Correlation between ACG features')
save_script_figure(output_folder)

feature_mean = mean(X, 1);
feature_sigma = std(X, 0, 1);
Xz = (X - feature_mean) ./ feature_sigma;

% pca_mu is the centering of this PCA. Section 4 uses it to place the
% low-power features in the same plane. t-SNE has no equivalent transform.
[pca_coeff, pca_score, pca_latent, ~, pca_explained, pca_mu] = pca(Xz);
open_full_figure('PCA variance explained')
subplot(1, 2, 1)
bar(pca_explained)
xlabel('Principal component')
ylabel('Variance explained (%)')
title('Variance explained')
subplot(1, 2, 2)
plot(cumsum(pca_explained), '-o', 'LineWidth', 1.5)
xlabel('Number of components')
ylabel('Cumulative variance explained (%)')
title('Cumulative variance')
ylim([0 100])
grid on
save_script_figure(output_folder)

open_full_figure('PCA population structure')
scatter(pca_score(:, 1), pca_score(:, 2), 35, 'filled', 'MarkerFaceAlpha', 0.7)
xlabel(sprintf('PC1 (%.1f%%)', pca_explained(1)))
ylabel(sprintf('PC2 (%.1f%%)', pca_explained(2)))
title('ACG features: PCA')
grid on
axis square
save_script_figure(output_folder)

open_full_figure('PCA feature loadings')
biplot(pca_coeff(:, 1:2), 'Scores', pca_score(:, 1:2), 'VarLabels', feature_names)
title('PCA loadings and neuron scores')
save_script_figure(output_folder)

n_complete = size(Xz, 1);
perplexity_values = [5 15 30 50];
perplexity_values = perplexity_values(perplexity_values < (n_complete - 1) / 3);
if isempty(perplexity_values)
    perplexity_values = max(2, floor((n_complete - 1) / 3));
end
Y_tsne = cell(numel(perplexity_values), 1);
for k = 1:numel(perplexity_values)
    rng(42);
    Y_tsne{k} = tsne(Xz, ...
        'NumDimensions', 2, ...
        'Perplexity', perplexity_values(k), ...
        'Standardize', false, ...
        'Distance', 'euclidean');
end

open_full_figure('ACG t-SNE')
for k = 1:numel(perplexity_values)
    subplot(1, numel(perplexity_values), k)
    scatter(Y_tsne{k}(:, 1), Y_tsne{k}(:, 2), 30, 'filled', 'MarkerFaceAlpha', 0.7)
    xlabel('t-SNE 1')
    ylabel('t-SNE 2')
    title(sprintf('Perplexity = %d', perplexity_values(k)))
    axis square
    grid on
end
save_script_figure(output_folder)

rng(42);
[idx_kmeans, kmeans_centers] = kmeans(Xz, 2, ...
    'Distance', 'sqeuclidean', ...
    'Replicates', 100, ...
    'MaxIter', 1000);
fprintf('\nK-means cluster sizes:\n');
disp(accumarray(idx_kmeans, 1))

ward_tree = linkage(Xz, 'ward');
idx_ward = cluster(ward_tree, 'maxclust', 2);
fprintf('\nWard cluster sizes:\n');
disp(accumarray(idx_ward, 1))
if numel(unique(idx_kmeans)) < 2 || numel(unique(idx_ward)) < 2
    error('Plot_acg_population_logit:cluster', ...
        'One of the two-group splits did not use both groups.');
end

open_full_figure('Comparison of exploratory partitions')
tiledlayout(2, 3, 'TileSpacing', 'compact')
nexttile
gscatter(pca_score(:, 1), pca_score(:, 2), idx_kmeans)
title('PCA: k-means')
xlabel('PC1')
ylabel('PC2')
grid on
nexttile
gscatter(pca_score(:, 1), pca_score(:, 2), idx_ward)
title('PCA: Ward')
xlabel('PC1')
ylabel('PC2')
grid on
nexttile
tsne_for_groups = Y_tsne{min(2, numel(Y_tsne))};
gscatter(tsne_for_groups(:, 1), tsne_for_groups(:, 2), idx_kmeans)
title('t-SNE: k-means assignments')
grid on
nexttile
gscatter(tsne_for_groups(:, 1), tsne_for_groups(:, 2), idx_ward)
title('t-SNE: Ward assignments')
grid on
nexttile
group_counts = crosstab(idx_kmeans, idx_ward);
imagesc(group_counts)
set(gca, 'XTick', 1:size(group_counts, 2), 'YTick', 1:size(group_counts, 1))
xlabel('Ward')
ylabel('K-means')
title('K-means vs Ward counts')
colorbar
disp(group_counts)
save_script_figure(output_folder)

open_full_figure('Silhouette comparison')
subplot(1, 2, 1)
sil_kmeans = silhouette(Xz, idx_kmeans, 'euclidean');
draw_silhouette(sil_kmeans, idx_kmeans)
title(sprintf('K-means, mean = %.3g', mean(sil_kmeans)))
subplot(1, 2, 2)
sil_ward = silhouette(Xz, idx_ward, 'euclidean');
draw_silhouette(sil_ward, idx_ward)
title(sprintf('Ward, mean = %.3g', mean(sil_ward)))
save_script_figure(output_folder)

modulation_names = {'Modulation1', 'Modulation2'};
modulation_column = find(ismember(feature_names, modulation_names));
if isempty(modulation_column)
    error('Plot_acg_population_logit:features', 'No modulation features remain.');
end
modulation_score = mean(Xz(:, modulation_column), 2, 'omitnan');
kmeans_modulation = [ ...
    mean(modulation_score(idx_kmeans == 1)); ...
    mean(modulation_score(idx_kmeans == 2))];
if kmeans_modulation(2) > kmeans_modulation(1)
    idx_kmeans = 3 - idx_kmeans;
end
ward_modulation = [ ...
    mean(modulation_score(idx_ward == 1)); ...
    mean(modulation_score(idx_ward == 2))];
if ward_modulation(2) > ward_modulation(1)
    idx_ward = 3 - idx_ward;
end
fprintf('\nAfter relabeling, cluster 1 is the higher-modulation group.\n')
if strcmpi(cluster_method, 'kmeans') || strcmpi(cluster_method, 'k-means')
    idx_cluster = idx_kmeans;
    cluster_method_label = 'k-means';
elseif strcmpi(cluster_method, 'ward')
    idx_cluster = idx_ward;
    cluster_method_label = 'Ward';
else
    error('Plot_acg_population_logit:cluster', ...
        'cluster_method must be ''kmeans'' or ''ward''.');
end
fprintf('The oscillatory group is %s cluster 1.\n', cluster_method_label)
counts_after = crosstab(idx_kmeans, idx_ward);
kmeans_sizes = accumarray(idx_kmeans, 1);
ward_sizes = accumarray(idx_ward, 1);
kmeans_modulation = [ ...
    mean(modulation_score(idx_kmeans == 1)); ...
    mean(modulation_score(idx_kmeans == 2))];
ward_modulation = [ ...
    mean(modulation_score(idx_ward == 1)); ...
    mean(modulation_score(idx_ward == 2))];
cluster_stats = empty_feature_stats();
cluster_stats(1) = [];

complete_rows = find(valid_rows);
% The first six names entered the split when they had variance. The rest
% are drawn so the two groups can be compared. Peak1 - Peak2 is positive
% when the first peak is larger than the second.
descriptor_names = { ...
    'Modulation1', 'Modulation2', 'PeakInterval', 'TroughInterval', ...
    'OscillationScore', 'Baseline', ...
    'DominantFrequency', 'PeakFrequency', 'TroughFrequency', ...
    'Persistence', 'Peak1', 'Peak2', 'Trough1', 'Trough2'};
descriptor_titles = { ...
    'Modulation1', 'Modulation2', 'PeakInterval', 'TroughInterval', ...
    'OscillationScore', 'Baseline', ...
    'DominantFrequency', 'PeakFrequency', 'TroughFrequency', ...
    'Persistence (2nd / 1st modulation)', 'Peak1', 'Peak2', 'Trough1', 'Trough2'};
table_names = feature_table.Properties.VariableNames;
descriptor_values = zeros(numel(complete_rows), 0);
shown_titles = {};
shown_roles = {};
for j = 1:numel(descriptor_names)
    this_name = descriptor_names{j};
    if ~ismember(this_name, table_names)
        continue
    end
    this_column = feature_table{complete_rows, this_name};
    if ~isnumeric(this_column)
        continue
    end
    descriptor_values(:, end + 1) = this_column; %#ok<AGROW>
    shown_titles{end + 1} = descriptor_titles{j}; %#ok<AGROW>
    if ismember(this_name, feature_names)
        shown_roles{end + 1} = 'used in the split'; %#ok<AGROW>
    else
        shown_roles{end + 1} = 'not used in the split'; %#ok<AGROW>
    end
end
if all(ismember({'Peak1', 'Peak2'}, table_names))
    descriptor_values(:, end + 1) = feature_table{complete_rows, 'Peak1'} ...
        - feature_table{complete_rows, 'Peak2'}; %#ok<AGROW>
    shown_titles{end + 1} = 'Peak1 - Peak2'; %#ok<AGROW>
    shown_roles{end + 1} = 'not used in the split'; %#ok<AGROW>
end
if ~isempty(shown_titles)
    n_shown = numel(shown_titles);
    cluster_stats = repmat(empty_feature_stats(), n_shown, 1);
    for j = 1:n_shown
        group_one = descriptor_values(idx_cluster == 1, j);
        group_two = descriptor_values(idx_cluster == 2, j);
        cluster_stats(j) = compare_two_samples(group_one, group_two, 'unpaired');
        cluster_stats(j).feature = shown_titles{j};
        cluster_stats(j).role = shown_roles{j};
    end
    cluster_raw_p = [cluster_stats.raw_p].';
    cluster_adjusted_p = holm_sidak_adjust(cluster_raw_p);
    for j = 1:n_shown
        cluster_stats(j).holm_p = cluster_adjusted_p(j);
    end
    open_full_figure(['Feature distributions by ', cluster_method_label, ' group'])
    n_descriptor_cols = 5;
    tiledlayout(ceil(numel(shown_titles) / n_descriptor_cols), n_descriptor_cols, ...
        'TileSpacing', 'compact')
    for j = 1:numel(shown_titles)
        nexttile
        hold on
        group_values = cell(1, 2);
        group_mean = nan(1, 2);
        for group = 1:2
            this_group = descriptor_values(idx_cluster == group, j);
            this_group = this_group(isfinite(this_group));
            group_values{group} = this_group;
            if ~isempty(this_group)
                group_mean(group) = mean(this_group);
            end
            draw_one_violin(group, this_group, [0.75 0.75 0.75])
        end
        if all(isfinite(group_mean))
            plot([1 2], group_mean, '-', 'Color', [0.1 0.1 0.1], ...
                'LineWidth', 2.4, 'HandleVisibility', 'off')
        end
        for group = 1:2
            draw_violin_mean_bar(group, group_values{group})
            if isempty(group_values{group})
                continue
            end
            point_x = group + violin_point_offset(numel(group_values{group}));
            scatter(point_x, group_values{group}, 12, [0 0 0], 'filled', ...
                'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none')
        end
        set(gca, 'XTick', [1 2])
        xlim([0.4 2.6])
        if ~isfinite(cluster_stats(j).holm_p)
            p_text = 'Too few neurons';
        else
            p_text = sprintf('Holm-Sidak p = %.3g', cluster_stats(j).holm_p);
        end
        title({shown_titles{j}, shown_roles{j}, p_text}, 'Interpreter', 'none')
        xlabel([cluster_method_label, ' group'])
        grid on
    end
    sgtitle({[cluster_method_label, ' cluster 1 is the higher modulation.'], ...
        'Black bars are means, joined to show the group difference.', ...
        'The title is the Holm-Sidak p across these panels.', ...
        'Peak1 - Peak2 is the drop from the first peak to the second.'})
    save_script_figure(output_folder)
end

ACG_PCA = struct( ...
    'feature_names', {feature_names}, ...
    'complete_rows', complete_rows, ...
    'X', X, ...
    'Xz', Xz, ...
    'coeff', pca_coeff, ...
    'score', pca_score, ...
    'explained', pca_explained, ...
    'latent', pca_latent, ...
    'Y_tsne', {Y_tsne}, ...
    'perplexity', perplexity_values, ...
    'idx_kmeans', idx_kmeans, ...
    'idx_ward', idx_ward, ...
    'idx_cluster', idx_cluster, ...
    'cluster_method', cluster_method, ...
    'kmeans_centers', kmeans_centers, ...
    'modulation_score', modulation_score, ...
    'valid_rows', valid_rows);
write_section_2_summary(output_folder, feature_names, size(X_raw, 1), size(X, 1), ...
    sum(~valid_rows), feature_correlation, pca_explained, perplexity_values, ...
    group_counts, counts_after, kmeans_sizes, ward_sizes, ...
    mean(sil_kmeans), mean(sil_ward), cluster_method_label, ...
    kmeans_modulation, ward_modulation, cluster_stats)
disp(['Run section 3 to keep ', cluster_method_label, ' cluster 1.'])

%% 3. Choose the oscillatory neurons
% Cluster 1 is the higher-modulation group of cluster_method. Its rows are
% the complete feature rows, not yet the image rows. The strongest
% modulation is placed first.
if ~exist('output_folder', 'var') || isempty(output_folder)
    output_folder = fullfile(fileparts(mfilename('fullpath')), 'acg_population_logit_figures');
end
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
if ~exist('idx_cluster', 'var') || ~exist('complete_rows', 'var') ...
        || ~exist('modulation_score', 'var') || ~exist('neuron_type', 'var')
    error('Plot_acg_population_logit:cluster', ...
        'The chosen split is not in the workspace. Run section 2 first.');
end

osc_score = nan(numel(neuron_type), 1);
osc_score(complete_rows) = modulation_score;
osc_member = false(numel(neuron_type), 1);
osc_member(complete_rows(idx_cluster == 1)) = true;
osc_rows = oscillatory_member_rows(osc_member, osc_score);
osc_neuron_type = neuron_type(osc_rows);
[type_percent, type_observed, type_total, lock_type_label] = ...
    lock_type_osc_rate(neuron_type, osc_rows);
n_osc = numel(osc_rows);
n_before_osc = numel(neuron_type);

shown_logit = struct();
shown_percentile = struct();
for s = 1:numel(scale_id)
    for c = 1:numel(column_id)
        shown_logit.(scale_id{s}).(column_id{c}) = ...
            imagesc_logit.(scale_id{s}).(column_id{c})(osc_rows, :);
        shown_percentile.(scale_id{s}).(column_id{c}) = ...
            percentile.(scale_id{s}).(column_id{c})(osc_rows, :);
    end
end
high_minus_low = shown_logit.sum.power_high - shown_logit.sum.power_low;
% Plain names for the oscillatory sum-to-1 images. Row 1 is the strongest
% modulation. The same rows are neuron_type(osc_rows).
imagesc_osc_sum_all = shown_logit.sum.all;
imagesc_osc_sum_power_high = shown_logit.sum.power_high;
imagesc_osc_sum_power_low = shown_logit.sum.power_low;

disp(sprintf('Oscillatory neurons: %d of %d. shown_logit uses that row order.', ...
    n_osc, n_before_osc))
disp('Run section 4 to draw.')
type_lines = cell(1, numel(lock_type_label));
for k = 1:numel(lock_type_label)
    type_lines{k} = sprintf('%-16s %d/%d   %.1f%%', lock_type_label{k}, ...
        type_observed(k), type_total(k), type_percent(k));
end
write_summary_section(output_folder, '3. Oscillatory neurons kept for the figures', join_summary_lines( ...
    ['Output folder: ', output_folder], ...
    sprintf('Oscillatory neurons: %d of %d', n_osc, n_before_osc), ...
    'Cluster 1 of the chosen split, strongest modulation first.', ...
    'Share of each lock type that is in that group:', ...
    type_lines))

%% 4. Figures
% Re-run this section after editing the figure. Section 3 does not need
% to run again while shown_logit is still in the workspace. The low-power
% PCA figure and the paired feature figure also need the section 2 PCA.
if ~exist('output_folder', 'var') || isempty(output_folder)
    output_folder = fullfile(fileparts(mfilename('fullpath')), 'acg_population_logit_figures');
end
n_example_neurons = 10;
n_example_subplots = 6;
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
if ~exist('shown_logit', 'var')
    error('Plot_acg_population_logit:shown', ...
        'shown_logit is not in the workspace. Run section 3 first.');
end

% Open intervals, in seconds, for the cycle scatter below.
trough_windows_s = [0.05 0.15; 0.35 0.45];
peak_windows_s = [0.25 0.35; 0.45 0.55];

% Symmetric color range for the logit image. 1 is about the 9th and 91st
% percentiles. 2 is about the 1st and 99th.
logit_color_limit = 2;
high_cut = 95;
low_cut = 5;
n_rows = 4;
n_cols = numel(column_id);
n_colors = 256;
half = n_colors / 2;
logit_map = [ ...
    linspace(0.15, 0.96, half)', linspace(0.28, 0.96, half)', linspace(0.70, 0.96, half)'; ...
    linspace(0.96, 0.70, half)', linspace(0.96, 0.18, half)', linspace(0.96, 0.18, half)'];
first_peak_lag = lag > peak_windows_s(1, 1) & lag < peak_windows_s(1, 2);
first_trough_lag = lag > trough_windows_s(1, 1) & lag < trough_windows_s(1, 2);
second_peak_lag = lag > peak_windows_s(2, 1) & lag < peak_windows_s(2, 2);
second_trough_lag = lag > trough_windows_s(2, 1) & lag < trough_windows_s(2, 2);
if ~any(first_peak_lag) || ~any(first_trough_lag) ...
        || ~any(second_peak_lag) || ~any(second_trough_lag)
    error('Plot_acg_population_logit:window', ...
        'A peak or trough window contains no lag bin.')
end
cycle_n = nan(1, n_cols);
cycle_r = nan(1, n_cols);
cycle_p = nan(1, n_cols);

for s = 1
    open_full_figure(['ACG population logit, ', scale_title{s}])
    for c = 1:n_cols
        this_logit = shown_logit.(scale_id{s}).(column_id{c});
        this_pctl = shown_percentile.(scale_id{s}).(column_id{c});
        n = size(this_logit, 1);
        if n < 1
            for row = 1:n_rows
                subplot(n_rows, n_cols, (row - 1) * n_cols + c)
                if row == 1
                    title({column_title{c}, sprintf('osc/n = 0/%d', n_before_osc)})
                end
                axis off
            end
            continue
        end

        subplot(n_rows, n_cols, c)
        if n == 1
            imagesc(lag, [0.5 1.5], this_logit)
            ylim([0.5 1.5])
        else
            imagesc(lag, 1:n, this_logit)
            set(gca, 'YTick', [1 n])
        end
        axis tight
        clim([-logit_color_limit logit_color_limit])
        colormap(gca, logit_map)
        ylabel('Neuron')
        title({column_title{c}, sprintf('osc/n = %d/%d', n_osc, n_before_osc)})
        box off
        if c == n_cols
            colorbar
        end

        first_diff = mean(this_logit(:, first_peak_lag), 2, 'omitnan') ...
            - mean(this_logit(:, first_trough_lag), 2, 'omitnan');
        second_diff = mean(this_logit(:, second_peak_lag), 2, 'omitnan') ...
            - mean(this_logit(:, second_trough_lag), 2, 'omitnan');
        subplot(n_rows, n_cols, n_cols + c)
        finite_pair = isfinite(first_diff) & isfinite(second_diff);
        cycle_n(c) = sum(finite_pair);
        if sum(finite_pair) < 3
            plot(first_diff, second_diff, '.k')
            title('Too few neurons')
        else
            cycle_fit = fitlm(first_diff(finite_pair), second_diff(finite_pair));
            plot(cycle_fit)
            legend off
            [r_mat, p_mat] = corrcoef(first_diff(finite_pair), second_diff(finite_pair));
            cycle_r(c) = r_mat(1, 2);
            cycle_p(c) = p_mat(1, 2);
            title(sprintf('r = %.2f, p = %.3g', r_mat(1, 2), p_mat(1, 2)))
        end
        xlabel('First cycle, peak - trough')
        ylabel('Second cycle, peak - trough')
        box off

        [high_percent, low_percent] = percentile_tail_fractions(this_pctl, high_cut, low_cut);
        subplot(n_rows, n_cols, 2 * n_cols + c)
        hold on
        plot(lag, high_percent, 'Color', [0.7 0.25 0.25], 'LineWidth', 1.8)
        plot(lag, low_percent, 'Color', [0.25 0.45 0.75], 'LineWidth', 1.8)
        ylim([0 100])
        ylabel('Neurons (%)')
        if c == 1
            legend({sprintf('Percentile >= %g', high_cut), ...
                sprintf('Percentile <= %g', low_cut)}, ...
                'Location', 'northeast', 'Box', 'off')
        end
        box off

        subplot(n_rows, n_cols, 3 * n_cols + c)
        hold on
        plot(lag([1 end]), [0 0], ':', 'Color', [0.45 0.45 0.45], ...
            'HandleVisibility', 'off')
        plot(lag, mean(this_logit, 1, 'omitnan'), 'k', 'LineWidth', 2.4)
        xlabel('Lag (s)')
        ylabel('Mean logit')
        ylim([-logit_color_limit logit_color_limit])
        box off
    end
    sgtitle({['Oscillatory neurons, ', scale_title{s}], ...
        [cluster_method_label, ' group with the higher modulation, strongest first.']})
    save_script_figure(output_folder)
end

open_full_figure('ACG population logit, high minus low power')
subplot(2, 1, 1)
if n_osc < 1
    title('No oscillatory neurons')
    axis off
else
    if n_osc == 1
        imagesc(lag, [0.5 1.5], high_minus_low)
        ylim([0.5 1.5])
    else
        imagesc(lag, 1:n_osc, high_minus_low)
        set(gca, 'YTick', [1 n_osc])
    end
    axis tight
    clim([-logit_color_limit logit_color_limit])
    colormap(gca, logit_map)
    ylabel('Neuron')
    title('High 25% power minus low 25% power')
    colorbar
    box off
end

subplot(2, 1, 2)
hold on
plot(lag([1 end]), [0 0], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
if n_osc >= 2
    [high_vs_low_h, high_vs_low_p] = ttest(shown_logit.sum.power_high, shown_logit.sum.power_low);
    mean_difference = mean(high_minus_low, 1, 'omitnan');
    plot(lag, mean_difference, 'k', 'LineWidth', 1.8)
    significant = high_vs_low_h == 1 & isfinite(mean_difference);
    plot(lag(significant), mean_difference(significant), '.', ...
        'Color', [0.7 0.25 0.25], 'MarkerSize', 14)
else
    high_vs_low_h = false(1, numel(lag));
    high_vs_low_p = nan(1, numel(lag));
    if n_osc == 1
        plot(lag, high_minus_low, 'k', 'LineWidth', 1.8)
    end
end
xlabel('Lag (s)')
ylabel('Mean logit difference')
title('Red marks a paired difference at p < 0.05')
box off
save_script_figure(output_folder)

% The PCA below keeps its original axes and only receives the low-power
% points. t-SNE cannot do that, so the t-SNE figure is a new embedding of
% the high-power neurons plus the oscillatory neurons at low power.
if ~exist('pca_coeff', 'var') || ~exist('pca_mu', 'var') ...
        || ~exist('pca_score', 'var') || ~exist('pca_explained', 'var') ...
        || ~exist('feature_mean', 'var') || ~exist('feature_sigma', 'var') ...
        || ~exist('feature_names', 'var') || ~exist('complete_rows', 'var') ...
        || ~exist('feature_table', 'var') || ~exist('imagesc_sum_power_low', 'var') ...
        || ~exist('osc_rows', 'var') || ~exist('Xz', 'var')
    error('Plot_acg_population_logit:pca', ...
        'The high-power PCA is not in the workspace. Run sections 2 and 3 first.');
end
positive_lag = lag > 0;
feature_table_low = acg2features( ...
    imagesc_sum_power_low(:, positive_lag), lag(positive_lag));
X_low = feature_table_low{complete_rows, feature_names};
Xz_low = (X_low - feature_mean) ./ feature_sigma;
score_low = (Xz_low - pca_mu) * pca_coeff;
[in_complete, osc_pos] = ismember(osc_rows, complete_rows);
osc_pos = osc_pos(in_complete);
high_xy = pca_score(osc_pos, 1:2);
low_xy = score_low(osc_pos, 1:2);
finite_projection = all(isfinite(high_xy), 2) & all(isfinite(low_xy), 2);
osc_in_complete = ismember(complete_rows(:), osc_rows(:));
group_color = [0.84 0.37 0.00; 0.00 0.62 0.45; 0.00 0.45 0.70];
group_label = { ...
    'Oscillatory, high power', ...
    'Non-oscillatory, high power', ...
    'Oscillatory, low power'};

pca_all = [pca_score(:, 1:2); low_xy(finite_projection, :)];
pca_groups = { ...
    pca_score(osc_in_complete, 1:2), ...
    pca_score(~osc_in_complete, 1:2), ...
    low_xy(finite_projection, :)};
open_full_figure('Low-power ACGs in the high-power PCA')
pca_axes = gobjects(1, 3);
for g = 1:3
    pca_axes(g) = subplot(1, 3, g);
    draw_iso_probability(pca_all, pca_groups{g}, group_color(g, :), ...
        low_xy(finite_projection, :), group_color(3, :))
    title(group_label{g}, 'Interpreter', 'none')
    xlabel(sprintf('PC1 (%.1f%%)', pca_explained(1)))
    ylabel(sprintf('PC2 (%.1f%%)', pca_explained(2)))
    grid on
    box off
end
pca_limits = embedding_axis_limits(pca_all);
for g = 1:3
    set(pca_axes(g), 'XLim', pca_limits(1:2), 'YLim', pca_limits(3:4))
end
linkaxes(pca_axes, 'xy')
sgtitle(['Black points are every neuron. Blue points are oscillatory neurons at low power. ', ...
    'Filled bands are that panel''s group.'])
save_script_figure(output_folder)

low_member = osc_in_complete & all(isfinite(Xz_low), 2);
X_embed = [Xz; Xz_low(low_member, :)];
group_id = ones(size(Xz, 1), 1);
group_id(~osc_in_complete) = 2;
group_id = [group_id; 3 * ones(sum(low_member), 1)];
n_embed = size(X_embed, 1);
candidate_perplexity = [5 15 30 50];
candidate_perplexity = candidate_perplexity( ...
    candidate_perplexity < (n_embed - 1) / 3);
if isempty(candidate_perplexity)
    candidate_perplexity = max(2, floor((n_embed - 1) / 3));
end
if numel(candidate_perplexity) >= 2
    embed_perplexity = candidate_perplexity(2);
else
    embed_perplexity = candidate_perplexity(1);
end
rng(42);
Y_power = tsne(X_embed, ...
    'NumDimensions', 2, ...
    'Perplexity', embed_perplexity, ...
    'Standardize', false, ...
    'Distance', 'euclidean');
fprintf(['Low-power oscillatory neurons in the new t-SNE: %d. ', ...
    'Omitted for missing low-power features: %d.\n'], ...
    sum(low_member), sum(osc_in_complete) - sum(low_member))

tsne_groups = { ...
    Y_power(group_id == 1, :), ...
    Y_power(group_id == 2, :), ...
    Y_power(group_id == 3, :)};
open_full_figure('t-SNE with low-power oscillatory neurons')
tsne_axes = gobjects(1, 3);
for g = 1:3
    tsne_axes(g) = subplot(1, 3, g);
    draw_iso_probability(Y_power, tsne_groups{g}, group_color(g, :), ...
        tsne_groups{3}, group_color(3, :))
    title({group_label{g}, sprintf('Perplexity = %g', embed_perplexity)}, ...
        'Interpreter', 'none')
    xlabel('t-SNE 1')
    ylabel('t-SNE 2')
    grid on
    box off
end
tsne_limits = embedding_axis_limits(Y_power);
for g = 1:3
    set(tsne_axes(g), 'XLim', tsne_limits(1:2), 'YLim', tsne_limits(3:4))
end
linkaxes(tsne_axes, 'xy')
sgtitle(['Black points are every neuron. Blue points are oscillatory neurons at low power. ', ...
    'Filled bands are that panel''s group.'])
save_script_figure(output_folder)

% Same columns as the cluster violins. High values come from the
% high-power feature table. Low values are the same neurons at low power.
comparison_names = { ...
    'Modulation1', 'Modulation2', 'PeakInterval', 'TroughInterval', ...
    'OscillationScore', 'Baseline', ...
    'DominantFrequency', 'PeakFrequency', 'TroughFrequency', ...
    'Persistence', 'Peak1', 'Peak2', 'Trough1', 'Trough2'};
comparison_titles = { ...
    'Modulation1', 'Modulation2', 'PeakInterval', 'TroughInterval', ...
    'OscillationScore', 'Baseline', ...
    'DominantFrequency', 'PeakFrequency', 'TroughFrequency', ...
    'Persistence (2nd / 1st modulation)', 'Peak1', 'Peak2', 'Trough1', 'Trough2'};
low_names = feature_table_low.Properties.VariableNames;
high_names = feature_table.Properties.VariableNames;
high_compare = zeros(numel(osc_rows), 0);
low_compare = zeros(numel(osc_rows), 0);
compare_titles = {};
for j = 1:numel(comparison_names)
    this_name = comparison_names{j};
    if ~ismember(this_name, high_names) || ~ismember(this_name, low_names)
        continue
    end
    high_column = feature_table{osc_rows, this_name};
    low_column = feature_table_low{osc_rows, this_name};
    if ~isnumeric(high_column) || ~isnumeric(low_column)
        continue
    end
    high_compare(:, end + 1) = high_column; %#ok<AGROW>
    low_compare(:, end + 1) = low_column; %#ok<AGROW>
    compare_titles{end + 1} = comparison_titles{j}; %#ok<AGROW>
end
if all(ismember({'Peak1', 'Peak2'}, high_names)) ...
        && all(ismember({'Peak1', 'Peak2'}, low_names))
    high_compare(:, end + 1) = feature_table{osc_rows, 'Peak1'} ...
        - feature_table{osc_rows, 'Peak2'}; %#ok<AGROW>
    low_compare(:, end + 1) = feature_table_low{osc_rows, 'Peak1'} ...
        - feature_table_low{osc_rows, 'Peak2'}; %#ok<AGROW>
    compare_titles{end + 1} = 'Peak1 - Peak2'; %#ok<AGROW>
end
if ~isempty(compare_titles)
    n_compare = numel(compare_titles);
    paired_stats = repmat(empty_feature_stats(), n_compare, 1);
    for j = 1:n_compare
        pair = isfinite(high_compare(:, j)) & isfinite(low_compare(:, j));
        paired_stats(j) = compare_two_samples( ...
            high_compare(pair, j), low_compare(pair, j), 'paired');
        paired_stats(j).feature = compare_titles{j};
        paired_stats(j).role = '';
    end
    paired_adjusted_p = holm_sidak_adjust([paired_stats.raw_p].');
    for j = 1:n_compare
        paired_stats(j).holm_p = paired_adjusted_p(j);
    end
    open_full_figure('High versus low power features, oscillatory neurons')
    n_compare_cols = 5;
    tiledlayout(ceil(numel(compare_titles) / n_compare_cols), n_compare_cols, ...
        'TileSpacing', 'compact')
    for j = 1:numel(compare_titles)
        nexttile
        high_vals = high_compare(:, j);
        low_vals = low_compare(:, j);
        pair = isfinite(high_vals) & isfinite(low_vals);
        hold on
        if any(pair)
            y_high = high_vals(pair);
            y_low = low_vals(pair);
            point_offset = violin_point_offset(numel(y_high));
            draw_one_violin(1, y_high, [0.84 0.37 0.00])
            draw_one_violin(2, y_low, [0.00 0.45 0.70])
            plot([1 + point_offset.'; 2 + point_offset.'], ...
                [y_high.'; y_low.'], '-', 'Color', [0.62 0.62 0.62], ...
                'LineWidth', 0.3)
            plot([1 2], [mean(y_high), mean(y_low)], '-', ...
                'Color', [0.1 0.1 0.1], 'LineWidth', 2.4, 'HandleVisibility', 'off')
            draw_violin_mean_bar(1, y_high)
            draw_violin_mean_bar(2, y_low)
            scatter(1 + point_offset, y_high, 12, [0 0 0], 'filled', ...
                'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none')
            scatter(2 + point_offset, y_low, 12, [0 0 0], 'filled', ...
                'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'none')
        end
        if ~isfinite(paired_stats(j).holm_p)
            title({compare_titles{j}, 'Too few neurons'}, 'Interpreter', 'none')
        else
            title({compare_titles{j}, sprintf('Holm-Sidak p = %.3g', ...
                paired_stats(j).holm_p)}, 'Interpreter', 'none')
        end
        xlim([0.4 2.6])
        set(gca, 'XTick', [1 2], 'XTickLabel', {'High power', 'Low power'})
        grid on
    end
    sgtitle({'Same oscillatory neurons at high and low power.', ...
        'Thin lines join one neuron. Thick bars are means and the line between them.', ...
        'The title is the Holm-Sidak p across these panels.'})
    save_script_figure(output_folder)
else
    paired_stats = empty_feature_stats();
    paired_stats(1) = [];
end

type_color = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.45 0.45 0.45; 0.75 0.75 0.75];
open_full_figure('Oscillatory cells within each lock type')
shown_percent = type_percent;
shown_percent(~isfinite(shown_percent)) = 0;
hold on
type_legend = cell(1, 4);
for k = 1:4
    bar(k, shown_percent(k), 'FaceColor', type_color(k, :), 'BarWidth', 0.8)
    type_legend{k} = sprintf('%s %d/%d', lock_type_label{k}, type_observed(k), type_total(k));
    text(k, min(shown_percent(k) + 4, 108), sprintf('%d/%d', type_observed(k), type_total(k)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end
legend(type_legend, 'Location', 'northeast', 'Box', 'off')
set(gca, 'XTick', 1:4, 'XTickLabel', lock_type_label)
ylabel('Oscillatory within type (%)')
ylim([0 115])
title(sprintf('Trough %d/%d, peak %d/%d, unlocked %d/%d, non-classified %d/%d', ...
    type_observed(1), type_total(1), type_observed(2), type_total(2), ...
    type_observed(3), type_total(3), type_observed(4), type_total(4)))
box off
save_script_figure(output_folder)
if ~exist('neuron_area', 'var') || numel(neuron_area) ~= numel(neuron_type)
    if exist('records', 'var') && exist('neuron_record_index', 'var')
        [neuron_session, neuron_area, neuron_cluster_id] = ...
            neuron_identity_from_records(records, neuron_record_index);
    elseif exist('checkpoint', 'var') && exist('neuron_record_index', 'var')
        stored_id = load(checkpoint, 'shuffle_records');
        id_records = stored_id.shuffle_records;
        if ~iscell(id_records)
            id_records = num2cell(id_records(:));
        end
        [neuron_session, neuron_area, neuron_cluster_id] = ...
            neuron_identity_from_records(id_records(:), neuron_record_index);
    else
        neuron_session = repmat({''}, numel(neuron_type), 1);
        neuron_area = repmat({''}, numel(neuron_type), 1);
        neuron_cluster_id = nan(numel(neuron_type), 1);
    end
end
[area_percent, area_observed, area_total, area_label] = ...
    area_osc_rate(neuron_area, osc_rows);
area_test = area_enrichment(area_label, area_observed, area_total);
open_full_figure('Oscillatory cells within each brain area')
shown_area = area_percent;
shown_area(~isfinite(shown_area)) = 0;
hold on
for k = 1:numel(area_label)
    bar(k, shown_area(k), 'FaceColor', [0.20 0.45 0.80], 'BarWidth', 0.8)
    if area_test.tested(k) && isfinite(area_test.holm_p(k))
        bar_text = sprintf('%d/%d\nHolm-Sidak p = %.3g', ...
            area_observed(k), area_total(k), area_test.holm_p(k));
    else
        bar_text = sprintf('%d/%d', area_observed(k), area_total(k));
    end
    text(k, min(shown_area(k) + 4, 100), bar_text, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Interpreter', 'none')
end
set(gca, 'XTick', 1:numel(area_label), 'XTickLabel', area_label)
try
    set(gca, 'TickLabelInterpreter', 'none')
catch
end
ylabel('Oscillatory within area (%)')
ylim([0 130])
title({'Share of analysed neurons in each area that are oscillatory', ...
    'Holm-Sidak p: more oscillatory than the other analysed areas'})
box off
save_script_figure(output_folder)
area_lines = cell(1, numel(area_label));
for k = 1:numel(area_label)
    if area_test.tested(k)
        area_lines{k} = sprintf('%-16s %4d/%-4d  expected %7.2f  raw p %8.3g  Holm-Sidak p %8.3g', ...
            area_label{k}, area_observed(k), area_total(k), area_test.expected(k), ...
            area_test.raw_p(k), area_test.holm_p(k));
    else
        area_lines{k} = sprintf('%-16s %4d/%-4d  not tested', ...
            area_label{k}, area_observed(k), area_total(k));
    end
end
missing_area = area_label(area_observed == 0);
if isempty(missing_area)
    missing_line = 'Every analysed area has at least one oscillatory neuron.';
else
    missing_line = ['Areas with no oscillatory neuron: ', strjoin(missing_area, ', '), '.'];
end
if all(strcmp(area_label, 'area unknown'))
    area_note = 'No brain area was stored on these neurons. Re-run section 1, or run section 4 while the shuffle checkpoint is available.';
else
    area_note = 'The denominator is every analysed neuron in that area. The numerator is how many of them are in the oscillatory group.';
end
if isfinite(area_test.chi2_p)
    chi_line = sprintf('Pearson chi-square across tested areas: chi2 = %.4g, df = %.0f, p = %.4g. Smallest expected cell = %.2f.', ...
        area_test.chi2, area_test.df, area_test.chi2_p, area_test.min_expected);
else
    chi_line = 'The chi-square test was not run. Fewer than two named areas had analysed neurons.';
end
if isfinite(area_test.min_expected) && area_test.min_expected < 5
    chi_caution = 'An expected cell is below 5, so the chi-square p is a rough check. Use the Fisher p-values for each area.';
else
    chi_caution = 'The chi-square p asks whether the areas share one oscillatory rate.';
end
write_summary_section(output_folder, '4. Oscillatory share by brain area', join_summary_lines( ...
    area_note, ...
    'Each tested area is compared with the other analysed neurons by a one-sided Fisher exact test.', ...
    'The alternative is a higher oscillatory count than expected from those other neurons.', ...
    'The expected count is the area size times the oscillatory rate of all tested neurons.', ...
    'Holm-Sidak adjusts the Fisher p-values across the tested areas. Area unknown is not tested.', ...
    chi_line, ...
    chi_caution, ...
    area_lines, ...
    missing_line))
area_feature_names = [comparison_names, {'Peak1 - Peak2'}];
area_feature_titles = [comparison_titles, {'Peak1 - Peak2'}];
n_area_features = numel(area_feature_names);
area_anova = repmat(empty_area_anova(), n_area_features, 1);
for j = 1:n_area_features
    [groups, names] = values_by_area(feature_table, neuron_area, ...
        area_feature_names{j}, {'VLPAG', 'area unknown'}, 2);
    area_anova(j).feature = area_feature_titles{j};
    area_anova(j).areas = names;
    area_anova(j).groups = groups;
    stats = one_way_anova(groups);
    area_anova(j).n = stats.n;
    area_anova(j).mean = stats.mean;
    area_anova(j).F = stats.F;
    area_anova(j).df1 = stats.df1;
    area_anova(j).df2 = stats.df2;
    area_anova(j).raw_p = stats.p;
    area_anova(j).eta = stats.eta;
end
area_anova_holm = holm_sidak_adjust([area_anova.raw_p]);
for j = 1:n_area_features
    area_anova(j).holm_p = area_anova_holm(j);
    names = area_anova(j).areas;
    groups = area_anova(j).groups;
    if ~isfinite(area_anova(j).holm_p) || area_anova(j).holm_p >= 0.05 || numel(names) < 2
        continue
    end
    pair_list = empty_feature_stats();
    pair_list(1) = [];
    pair_index = zeros(0, 2);
    for a = 1:numel(names)
        for b = a + 1:numel(names)
            pair = compare_two_samples(groups{a}, groups{b}, 'unpaired');
            pair.feature = sprintf('%s vs %s', names{a}, names{b});
            pair.role = area_anova(j).feature;
            pair_list(end + 1) = pair; %#ok<AGROW>
            pair_index(end + 1, :) = [a, b]; %#ok<AGROW>
        end
    end
    if isempty(pair_list)
        continue
    end
    pair_holm = holm_sidak_adjust([pair_list.raw_p]);
    for p = 1:numel(pair_list)
        pair_list(p).holm_p = pair_holm(p);
    end
    area_anova(j).pairs = pair_list;
    area_anova(j).pair_index = pair_index;
end
open_full_figure('Oscillatory features by brain area')
n_descriptor_cols = 5;
tiledlayout(ceil(n_area_features / n_descriptor_cols), n_descriptor_cols, ...
    'TileSpacing', 'compact')
for j = 1:n_area_features
    nexttile
    hold on
    names = area_anova(j).areas;
    groups = area_anova(j).groups;
    for a = 1:numel(names)
        y = groups{a}(:);
        y = y(isfinite(y));
        if isempty(y)
            continue
        end
        scatter(a + violin_point_offset(numel(y)), y, 12, [0 0 0], 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor', 'none')
    end
    for a = 1:numel(names)
        if ~isfinite(area_anova(j).mean(a))
            continue
        end
        bar(a, area_anova(j).mean(a), 0.55, 'FaceColor', area_bar_color(names{a}), ...
            'FaceAlpha', 0.45, 'EdgeColor', 'none')
        hold on
    end
    y_top = nan;
    if ~isempty(area_anova(j).pairs)
        show = find([area_anova(j).pairs.holm_p] < 0.05);
        if ~isempty(show)
            pooled = [];
            for a = 1:numel(groups)
                pooled = [pooled; groups{a}(:)]; %#ok<AGROW>
            end
            y_top = draw_sig_brackets(area_anova(j).pair_index(show, :), ...
                [area_anova(j).pairs(show).holm_p], pooled);
        end
    end
    if ~isempty(names)
        set(gca, 'XTick', 1:numel(names), 'XTickLabel', names)
        xlim([0.4, numel(names) + 0.6])
    end
    try
        set(gca, 'TickLabelInterpreter', 'none')
    catch
    end
    if isfinite(y_top)
        yl = ylim;
        ylim([yl(1), max(yl(2), y_top)])
    end
    if ~isfinite(area_anova(j).holm_p)
        p_text = 'Too few neurons';
    else
        p_text = sprintf('Holm-Sidak p = %.3g', area_anova(j).holm_p);
    end
    title({area_anova(j).feature, p_text}, 'Interpreter', 'none')
    grid on
    box off
end
sgtitle({'High-power oscillatory features by brain area. VLPAG is left out.', ...
    'Points are neurons. Bars are means. Stars mark pairwise differences after a significant ANOVA.'})
save_script_figure(output_folder)
write_summary_section(output_folder, '4. Oscillatory features by brain area', ...
    area_feature_report(area_anova))
lock_anova = repmat(empty_area_anova(), n_area_features, 1);
for j = 1:n_area_features
    [groups, names] = values_by_lock(feature_table, neuron_type, ...
        area_feature_names{j}, 2);
    lock_anova(j).feature = area_feature_titles{j};
    lock_anova(j).areas = names;
    lock_anova(j).groups = groups;
    stats = one_way_anova(groups);
    lock_anova(j).n = stats.n;
    lock_anova(j).mean = stats.mean;
    lock_anova(j).F = stats.F;
    lock_anova(j).df1 = stats.df1;
    lock_anova(j).df2 = stats.df2;
    lock_anova(j).raw_p = stats.p;
    lock_anova(j).eta = stats.eta;
end
lock_anova_holm = holm_sidak_adjust([lock_anova.raw_p]);
for j = 1:n_area_features
    lock_anova(j).holm_p = lock_anova_holm(j);
    names = lock_anova(j).areas;
    groups = lock_anova(j).groups;
    if ~isfinite(lock_anova(j).holm_p) || lock_anova(j).holm_p >= 0.05 || numel(names) < 2
        continue
    end
    pair_list = empty_feature_stats();
    pair_list(1) = [];
    pair_index = zeros(0, 2);
    for a = 1:numel(names)
        for b = a + 1:numel(names)
            pair = compare_two_samples(groups{a}, groups{b}, 'unpaired');
            pair.feature = sprintf('%s vs %s', names{a}, names{b});
            pair.role = lock_anova(j).feature;
            pair_list(end + 1) = pair; %#ok<AGROW>
            pair_index(end + 1, :) = [a, b]; %#ok<AGROW>
        end
    end
    if isempty(pair_list)
        continue
    end
    pair_holm = holm_sidak_adjust([pair_list.raw_p]);
    for p = 1:numel(pair_list)
        pair_list(p).holm_p = pair_holm(p);
    end
    lock_anova(j).pairs = pair_list;
    lock_anova(j).pair_index = pair_index;
end
open_full_figure('Oscillatory features by lock type')
tiledlayout(ceil(n_area_features / n_descriptor_cols), n_descriptor_cols, ...
    'TileSpacing', 'compact')
for j = 1:n_area_features
    nexttile
    hold on
    names = lock_anova(j).areas;
    groups = lock_anova(j).groups;
    for a = 1:numel(names)
        y = groups{a}(:);
        y = y(isfinite(y));
        if isempty(y)
            continue
        end
        scatter(a + violin_point_offset(numel(y)), y, 12, [0 0 0], 'filled', ...
            'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor', 'none')
    end
    for a = 1:numel(names)
        if ~isfinite(lock_anova(j).mean(a))
            continue
        end
        bar(a, lock_anova(j).mean(a), 0.55, 'FaceColor', lock_bar_color(names{a}), ...
            'FaceAlpha', 0.45, 'EdgeColor', 'none')
        hold on
    end
    y_top = nan;
    if ~isempty(lock_anova(j).pairs)
        show = find([lock_anova(j).pairs.holm_p] < 0.05);
        if ~isempty(show)
            pooled = [];
            for a = 1:numel(groups)
                pooled = [pooled; groups{a}(:)]; %#ok<AGROW>
            end
            y_top = draw_sig_brackets(lock_anova(j).pair_index(show, :), ...
                [lock_anova(j).pairs(show).holm_p], pooled);
        end
    end
    if ~isempty(names)
        set(gca, 'XTick', 1:numel(names), 'XTickLabel', names)
        xlim([0.4, numel(names) + 0.6])
    end
    try
        set(gca, 'TickLabelInterpreter', 'none')
    catch
    end
    if isfinite(y_top)
        yl = ylim;
        ylim([yl(1), max(yl(2), y_top)])
    end
    if ~isfinite(lock_anova(j).holm_p)
        p_text = 'Too few neurons';
    else
        p_text = sprintf('Holm-Sidak p = %.3g', lock_anova(j).holm_p);
    end
    title({lock_anova(j).feature, p_text}, 'Interpreter', 'none')
    grid on
    box off
end
sgtitle({'High-power oscillatory features by lock type. Non-classified cells are left out.', ...
    'Points are neurons. Bars are means. Stars mark pairwise differences after a significant ANOVA.'})
save_script_figure(output_folder)
write_summary_section(output_folder, '4. Oscillatory features by lock type', ...
    area_feature_report(lock_anova, { ...
    'One-way ANOVA of each high-power oscillatory feature, with lock type as the factor.', ...
    'The groups are trough, peak, and unlocked. Non-classified cells are excluded.', ...
    'A lock type is included in a feature only when it has at least two finite values.', ...
    'Holm-Sidak adjusts the ANOVA p-values across features.', ...
    'Pairwise equal-variance t-tests are run only when that ANOVA Holm-Sidak p is below 0.05.', ...
    'Those pairs are then Holm-Sidak adjusted within the feature. Stars use that adjusted p.'}))
if ~exist('acg_method', 'var') || isempty(acg_method)
    acg_method = 'overlap';
end
if exist('checkpoint', 'var')
    example_summary_path = fullfile(fileparts(checkpoint), ...
        locking_summary_filename(false, acg_method));
else
    example_summary_path = '';
end
example_bank = load_example_acg_bank(example_summary_path);
example_pos = example_positions(osc_rows, neuron_area, n_example_neurons, 'VLPAG');
n_example = numel(example_pos);
example_names = [comparison_names, {'Peak1 - Peak2'}];
neurons_per_figure = max(1, floor(n_example_subplots / 2));
n_example_cols = neurons_per_figure;
if n_example > 0
    example_lines = cell(n_example, 1);
    lag_limits = example_lag_limits(lag);
    n_example_figures = ceil(n_example / neurons_per_figure);
    for page = 1:n_example_figures
        i0 = (page - 1) * neurons_per_figure + 1;
        i1 = min(page * neurons_per_figure, n_example);
        open_full_figure(sprintf('Example oscillatory neurons %d', page))
        for k = 1:(i1 - i0 + 1)
            i = i0 + k - 1;
            rank = example_pos(i);
            pop_row = osc_rows(rank);
            this_area = area_display_name(neuron_area{pop_row});
            lock_label = display_lock_label(neuron_type{pop_row});
            [acg_lag, acg_y, acg_found] = example_acg_curve(example_bank, ...
                neuron_session{pop_row}, neuron_cluster_id(pop_row));
            subplot(2, n_example_cols, k)
            if acg_found
                plot(acg_lag, acg_y, 'k', 'LineWidth', 1.4)
                xlim(lag_limits)
                ylim(example_acg_ylim(acg_lag, acg_y, lag_limits))
            else
                text(0.5, 0.5, 'ACG not in summary', 'HorizontalAlignment', 'center')
                axis([0 1 0 1])
            end
            xlabel('Lag (s)')
            ylabel('Sum-to-1')
            title({sprintf('%s, %s', this_area, lock_label), ...
                sprintf('%s id %g', neuron_session{pop_row}, neuron_cluster_id(pop_row))}, ...
                'Interpreter', 'none')
            box off
            subplot(2, n_example_cols, n_example_cols + k)
            plot(lag, shown_percentile.sum.power_high(rank, :), 'Color', [0.15 0.35 0.7], ...
                'LineWidth', 1.4)
            xlim(lag_limits)
            ylim([0 100])
            xlabel('Lag (s)')
            ylabel('Percentile')
            box off
            example_lines{i} = example_neuron_summary(i, neuron_session{pop_row}, ...
                neuron_cluster_id(pop_row), this_area, lock_label, osc_score(pop_row), ...
                acg_found, feature_table, feature_table_low, pop_row, example_names);
        end
        sgtitle({sprintf('Neurons %d-%d of %d. VLPAG is left out.', i0, i1, n_example), ...
            'Top: autocorrelogram on its own scale. Bottom: high-power shuffle percentile.'}, ...
            'Interpreter', 'none')
        save_script_figure(output_folder)
    end
    write_summary_section(output_folder, '4. Example oscillatory neurons', join_summary_lines( ...
        sprintf('n_example_neurons = %d. Drawn: %d. Subplots per figure: %d.', ...
            n_example_neurons, n_example, 2 * neurons_per_figure), ...
        sprintf('Figures: %d. Each figure holds up to %d neurons, autocorrelogram above percentile.', ...
            n_example_figures, neurons_per_figure), ...
        'VLPAG is left out. Each other area that has an oscillatory neuron contributes its strongest one.', ...
        'Further neurons are the next highest modulation scores outside VLPAG, until n_example_neurons is reached.', ...
        'If those areas already exceed n_example_neurons, every area representative is still drawn.', ...
        'The drawn order is strongest modulation first.', ...
        example_bank.note, ...
        'Each autocorrelogram uses the same lag window as the percentile.', ...
        'Its y-axis is set from that neuron alone, using lags at and after 0.05 s, so the bin beside zero does not flatten the trace.', ...
        'The percentile is the high-power sum-to-1 shuffle percentile.', ...
        'Features are the high-power and low-power values for that neuron.', ...
        example_lines))
else
    write_summary_section(output_folder, '4. Example oscillatory neurons', join_summary_lines( ...
        sprintf('n_example_neurons = %d. Drawn: 0.', n_example_neurons), ...
        'No oscillatory neuron outside VLPAG was available.'))
end
write_section_4_summary(output_folder, column_title, cycle_n, cycle_r, cycle_p, ...
    lag, high_vs_low_h, high_vs_low_p, n_osc, ...
    sum(osc_in_complete), sum(~osc_in_complete), sum(finite_projection), ...
    sum(low_member), sum(osc_in_complete) - sum(low_member), embed_perplexity, ...
    paired_stats, lock_type_label, type_observed, type_total, type_percent)


function draw_iso_probability(all_xy, group_xy, face_color, highlight_xy, highlight_color)
%DRAW_ISO_PROBABILITY  Filled density bands for one group on a black cloud.
%
% The levels follow Figure 1 LDA: a kernel density on a shared grid, then
% the density values that enclose 10% to 90% of that density, in steps of
% 10%. Darker bands are the higher-density core. highlight_xy is drawn
% after the bands so those points stay visible on every panel.

if nargin < 4
    highlight_xy = zeros(0, 2);
    highlight_color = [0 0 0];
end
all_xy = all_xy(all(isfinite(all_xy), 2), :);
group_xy = group_xy(all(isfinite(group_xy), 2), :);
highlight_xy = highlight_xy(all(isfinite(highlight_xy), 2), :);
if isempty(all_xy) && isempty(highlight_xy)
    return
end
if ~isempty(all_xy)
    scatter(all_xy(:, 1), all_xy(:, 2), 12, [0 0 0], 'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none')
    hold on
end
if size(group_xy, 1) >= 3 && ~isempty(all_xy)
    span = max(all_xy, [], 1) - min(all_xy, [], 1);
    span(span <= 0) = 1;
    pad = 0.05 * span;
    x1 = linspace(min(all_xy(:, 1)) - pad(1), max(all_xy(:, 1)) + pad(1), 100);
    x2 = linspace(min(all_xy(:, 2)) - pad(2), max(all_xy(:, 2)) + pad(2), 100);
    [xg, yg] = meshgrid(x1, x2);
    density = [];
    try
        density = ksdensity(group_xy, [xg(:), yg(:)]);
    catch
        density = [];
    end
    if ~isempty(density)
        density_grid = reshape(density, size(xg));
        if any(density_grid(:) > 0)
            density_sorted = sort(density_grid(:), 'descend');
            cumulative = cumsum(density_sorted) / sum(density_sorted);
            probability_levels = 0.1:0.1:0.9;
            thresholds = nan(size(probability_levels));
            for level_n = 1:numel(probability_levels)
                level_index = find(cumulative >= probability_levels(level_n), 1, 'first');
                if ~isempty(level_index)
                    thresholds(level_n) = density_sorted(level_index);
                end
            end
            thresholds = unique(thresholds(isfinite(thresholds)));
            if numel(thresholds) >= 2
                sorted_thresholds = sort(thresholds(:));
                n_thresholds = numel(sorted_thresholds);
                for band_n = (n_thresholds - 1):-1:1
                    contourf(x1, x2, density_grid, ...
                        sorted_thresholds([band_n band_n + 1]), ...
                        'LineWidth', 0.5, 'FaceColor', face_color, ...
                        'EdgeColor', 'none', ...
                        'FaceAlpha', 0.28 + 0.42 * band_n / n_thresholds);
                end
            end
        end
    end
end
if ~isempty(highlight_xy)
    scatter(highlight_xy(:, 1), highlight_xy(:, 2), 18, highlight_color, ...
        'filled', 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor', 'none')
end
end


function limits = embedding_axis_limits(xy)
xy = xy(all(isfinite(xy), 2), :);
span = max(xy, [], 1) - min(xy, [], 1);
span(span <= 0) = 1;
pad = 0.08 * span;
limits = [ ...
    min(xy(:, 1)) - pad(1), max(xy(:, 1)) + pad(1), ...
    min(xy(:, 2)) - pad(2), max(xy(:, 2)) + pad(2)];
end


function draw_one_violin(x_center, y, face_color)
%DRAW_ONE_VIOLIN  Kernel-density violin centered on a chosen x.
%
% The shape does not place the observations. Callers scatter those points
% at x positions they keep, so a later line can join a pair exactly.
% A long black bar marks the mean.

y = y(:);
y = y(isfinite(y));
if isempty(y)
    return
end
if numel(y) >= 2 && max(y) > min(y)
    density = [];
    try
        [density, support] = ksdensity(y, 'NumPoints', 100);
    catch
        density = [];
    end
    if ~isempty(density)
        density = density(:).';
        support = support(:).';
        if max(density) > 0
            half_width = 0.32;
            density = density / max(density) * half_width;
            fill([x_center + density, x_center - density(end:-1:1)], ...
                [support, support(end:-1:1)], face_color, ...
                'FaceAlpha', 0.22, 'EdgeColor', face_color, 'LineWidth', 1, ...
                'HandleVisibility', 'off')
        end
    end
end
draw_violin_mean_bar(x_center, y)
end


function draw_violin_mean_bar(x_center, y)
%DRAW_VIOLIN_MEAN_BAR  Long thick bar at the mean of one violin.

y = y(:);
y = y(isfinite(y));
if isempty(y)
    return
end
mean_value = mean(y);
plot([x_center - 0.36, x_center + 0.36], [mean_value, mean_value], ...
    'k', 'LineWidth', 2.8, 'HandleVisibility', 'off')
end


function s = empty_area_anova()
s = struct('feature', '', 'areas', {{}}, 'groups', {{}}, ...
    'n', [], 'mean', [], 'F', nan, 'df1', nan, 'df2', nan, ...
    'raw_p', nan, 'holm_p', nan, 'eta', nan, ...
    'pairs', empty_feature_stats(), 'pair_index', zeros(0, 2));
s.pairs(1) = [];
end


function [groups, names] = values_by_area(feature_table, neuron_area, feature_name, skip_names, min_n)
values = feature_column(feature_table, feature_name);
labels = cell(numel(neuron_area), 1);
for i = 1:numel(neuron_area)
    labels{i} = area_display_name(neuron_area{i});
end
ordered = {'SC', 'LPAG', 'DLPAG', 'DR'};
for i = 1:numel(labels)
    if ~any(strcmp(ordered, labels{i}))
        ordered{end + 1} = labels{i}; %#ok<AGROW>
    end
end
groups = {};
names = {};
for k = 1:numel(ordered)
    if any(strcmpi(skip_names, ordered{k}))
        continue
    end
    y = values(strcmp(labels, ordered{k}));
    y = y(isfinite(y));
    if numel(y) < min_n
        continue
    end
    groups{end + 1} = y(:); %#ok<AGROW>
    names{end + 1} = ordered{k}; %#ok<AGROW>
end
end


function values = feature_column(feature_table, name)
if strcmp(name, 'Peak1 - Peak2')
    values = feature_column(feature_table, 'Peak1') - feature_column(feature_table, 'Peak2');
    return
end
n = size(feature_table, 1);
if ~ismember(name, feature_table.Properties.VariableNames)
    values = nan(n, 1);
    return
end
values = double(feature_table{:, name});
values = values(:);
end


function stats = one_way_anova(groups)
%ONE_WAY_ANOVA  Between-area F test for one feature.

n_groups = numel(groups);
ns = zeros(1, n_groups);
means = nan(1, n_groups);
ss_within = 0;
all_y = zeros(0, 1);
for g = 1:n_groups
    y = groups{g}(:);
    y = y(isfinite(y));
    ns(g) = numel(y);
    if ns(g) > 0
        means(g) = mean(y);
        all_y = [all_y; y]; %#ok<AGROW>
    end
    if ns(g) > 1
        ss_within = ss_within + (ns(g) - 1) * var(y, 0);
    end
end
n_total = numel(all_y);
stats.n = ns;
stats.mean = means;
stats.df1 = n_groups - 1;
stats.df2 = n_total - n_groups;
stats.F = nan;
stats.p = nan;
stats.eta = nan;
if n_groups < 2 || stats.df2 < 1 || n_total < 1
    return
end
grand = mean(all_y);
ss_between = sum(ns .* (means - grand) .^ 2);
ss_total = ss_between + ss_within;
if ss_total > 0
    stats.eta = ss_between / ss_total;
end
if ss_within <= 0
    if ss_between > 0
        stats.F = inf;
        stats.p = 0;
        stats.eta = 1;
    else
        stats.F = 0;
        stats.p = 1;
        stats.eta = 0;
    end
    return
end
stats.F = (ss_between / stats.df1) / (ss_within / stats.df2);
stats.p = f_survival(stats.F, stats.df1, stats.df2);
end


function p = f_survival(F, df1, df2)
if ~isfinite(df1) || ~isfinite(df2) || df1 <= 0 || df2 <= 0
    p = nan;
    return
end
if ~isfinite(F)
    if F > 0
        p = 0;
    else
        p = nan;
    end
    return
end
if F < 0
    p = nan;
    return
end
p = betainc(df2 / (df2 + df1 * F), df2 / 2, df1 / 2);
end


function [groups, names] = values_by_lock(feature_table, neuron_type, feature_name, min_n)
%VALUES_BY_LOCK  Feature values for trough, peak, and unlocked neurons.

values = feature_column(feature_table, feature_name);
kept = {'Trough', 'Peak', 'Unlocked'};
labels = cell(numel(neuron_type), 1);
for i = 1:numel(neuron_type)
    labels{i} = display_lock_label(neuron_type{i});
end
groups = {};
names = {};
for k = 1:numel(kept)
    y = values(strcmp(labels, kept{k}));
    y = y(isfinite(y));
    if numel(y) < min_n
        continue
    end
    groups{end + 1} = y(:); %#ok<AGROW>
    names{end + 1} = kept{k}; %#ok<AGROW>
end
end


function color = lock_bar_color(name)
switch name
    case 'Trough'
        color = [0.20 0.45 0.80];
    case 'Peak'
        color = [0.85 0.25 0.20];
    case 'Unlocked'
        color = [0.45 0.45 0.45];
    otherwise
        color = [0.75 0.75 0.75];
end
end


function color = area_bar_color(name)
switch name
    case 'SC'
        color = [0.20 0.45 0.80];
    case 'LPAG'
        color = [0.85 0.33 0.10];
    case 'DLPAG'
        color = [0.47 0.67 0.19];
    case 'DR'
        color = [0.49 0.18 0.56];
    otherwise
        color = [0.45 0.45 0.45];
end
end


function y_top = draw_sig_brackets(index_pairs, p_values, y_values)
%DRAW_SIG_BRACKETS  Stacked lines and stars for significant area pairs.

y_values = y_values(isfinite(y_values));
y_top = nan;
if isempty(y_values) || isempty(index_pairs)
    return
end
y_hi = max(y_values);
y_lo = min(y_values);
span = y_hi - y_lo;
if span <= 0
    span = max(abs(y_hi), 1);
end
step = 0.16 * span;
tick = 0.035 * span;
base = y_hi + 0.1 * span;
distance = index_pairs(:, 2) - index_pairs(:, 1);
[~, order] = sort(distance);
index_pairs = index_pairs(order, :);
p_values = p_values(order);
y_top = base;
for r = 1:size(index_pairs, 1)
    y = base + (r - 1) * step;
    x1 = index_pairs(r, 1);
    x2 = index_pairs(r, 2);
    plot([x1, x1, x2, x2], [y - tick, y, y, y - tick], 'k-', 'LineWidth', 1.1, ...
        'HandleVisibility', 'off')
    text((x1 + x2) / 2, y, star_text(p_values(r)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'Interpreter', 'none')
    y_top = y + 0.12 * span;
end
end


function label = star_text(p)
if ~isfinite(p) || p >= 0.05
    label = '';
elseif p < 0.001
    label = '***';
elseif p < 0.01
    label = '**';
else
    label = '*';
end
end


function text = area_feature_report(area_anova, intro)
%AREA_FEATURE_REPORT  ANOVA and, when it is significant, its pairwise tests.

if nargin < 2 || isempty(intro)
    intro = { ...
        'One-way ANOVA of each high-power oscillatory feature, with brain area as the factor.', ...
        'VLPAG is excluded. Neurons without a stored area are excluded.', ...
        'An area is included in a feature only when it has at least two finite values.', ...
        'Holm-Sidak adjusts the ANOVA p-values across features.', ...
        'Pairwise equal-variance t-tests are run only when that ANOVA Holm-Sidak p is below 0.05.', ...
        'Those pairs are then Holm-Sidak adjusted within the feature. Stars use that adjusted p.'};
end
lines = intro;
for j = 1:numel(area_anova)
    row = area_anova(j);
    if isempty(row.areas)
        n_text = 'none';
    else
        n_bits = cell(1, numel(row.areas));
        for a = 1:numel(row.areas)
            n_bits{a} = sprintf('%s %g', row.areas{a}, row.n(a));
        end
        n_text = strjoin(n_bits, ', ');
    end
    lines{end + 1} = sprintf('%-28s F(%s, %s) = %s   eta2 = %s   raw p = %s   Holm-Sidak p = %s   n: %s', ...
        row.feature, stat_num(row.df1, '%.0f'), stat_num(row.df2, '%.0f'), ...
        stat_num(row.F, '%.4g'), stat_num(row.eta, '%.4g'), ...
        stat_num(row.raw_p, '%.4g'), stat_num(row.holm_p, '%.4g'), n_text); %#ok<AGROW>
    if ~isfinite(row.holm_p) || row.holm_p >= 0.05
        lines{end + 1} = '  Pairwise tests were not run.'; %#ok<AGROW>
        continue
    end
    for p = 1:numel(row.pairs)
        pair = row.pairs(p);
        lines{end + 1} = sprintf('  %-20s t = %s   df = %s   d = %s   raw p = %s   Holm-Sidak p = %s', ...
            pair.feature, stat_num(pair.t, '%.4g'), stat_num(pair.df, '%.0f'), ...
            stat_num(pair.effect, '%.4g'), stat_num(pair.raw_p, '%.4g'), ...
            stat_num(pair.holm_p, '%.4g')); %#ok<AGROW>
    end
end
text = strjoin(lines, newline);
end


function draw_silhouette(score, cluster_id)
%DRAW_SILHOUETTE  Horizontal silhouette bars in the current axes.
%
% silhouette returns values and does not draw when those values are
% saved. This draws one bar per neuron, longest first within each cluster.

score = score(:);
cluster_id = cluster_id(:);
ok = isfinite(score) & isfinite(cluster_id);
score = score(ok);
cluster_id = cluster_id(ok);
hold on
if isempty(score)
    text(0.5, 0.5, 'No silhouette values', 'HorizontalAlignment', 'center')
    axis([0 1 0 1])
    return
end
groups = unique(cluster_id);
colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
y_cursor = 0;
tick_at = zeros(1, numel(groups));
tick_label = cell(1, numel(groups));
for g = 1:numel(groups)
    this_score = sort(score(cluster_id == groups(g)), 'descend');
    n_this = numel(this_score);
    y_pos = y_cursor + (1:n_this).';
    color = colors(mod(g - 1, size(colors, 1)) + 1, :);
    barh(y_pos, this_score, 1, 'FaceColor', color, 'EdgeColor', 'none')
    hold on
    tick_at(g) = y_cursor + (n_this + 1) / 2;
    tick_label{g} = sprintf('%g', groups(g));
    y_cursor = y_cursor + n_this;
end
plot([0 0], [0, y_cursor + 1], '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 0.8)
ylim([0, y_cursor + 1])
low = min(0, min(score));
xlim([low - 0.05 * max(1 + abs(low), eps), 1])
set(gca, 'YTick', tick_at, 'YTickLabel', tick_label, 'YDir', 'reverse')
xlabel('Silhouette value')
ylabel('Cluster')
box off
end


function result = area_enrichment(area_label, area_observed, area_total)
%AREA_ENRICHMENT  One-sided Fisher tests for extra oscillatory neurons.

n = numel(area_label);
result.expected = nan(1, n);
result.raw_p = nan(1, n);
result.holm_p = nan(1, n);
result.tested = false(1, n);
result.chi2 = nan;
result.df = nan;
result.chi2_p = nan;
result.min_expected = nan;
named = true(1, n);
for k = 1:n
    if strcmp(area_label{k}, 'area unknown') || area_total(k) < 1
        named(k) = false;
    end
end
osc_all = sum(area_observed(named));
n_all = sum(area_total(named));
non_all = n_all - osc_all;
for k = 1:n
    if ~named(k) || n_all < 1
        continue
    end
    result.tested(k) = true;
    a = area_observed(k);
    b = area_total(k) - a;
    c = osc_all - a;
    d = non_all - b;
    result.expected(k) = area_total(k) * osc_all / n_all;
    result.raw_p(k) = fisher_exact_greater(a, b, c, d);
end
result.holm_p = holm_sidak_adjust(result.raw_p);
if sum(named) >= 2 && osc_all > 0 && non_all > 0
    observed = [area_observed(named); area_total(named) - area_observed(named)];
    expected = [area_total(named) * osc_all / n_all; ...
        area_total(named) * non_all / n_all];
    result.min_expected = min(expected(:));
    result.chi2 = sum((observed(:) - expected(:)) .^ 2 ./ expected(:));
    result.df = sum(named) - 1;
    result.chi2_p = chi2_survival(result.chi2, result.df);
end
end


function p = fisher_exact_greater(a, b, c, d)
%FISHER_EXACT_GREATER  One-sided probability of at least a oscillatory neurons.

a = round(a);
b = round(b);
c = round(c);
d = round(d);
population = a + b + c + d;
oscillatory = a + c;
drawn = a + b;
if population < 1 || drawn < 0 || oscillatory < 0 || drawn > population
    p = nan;
    return
end
lo = max(0, drawn - (population - oscillatory));
hi = min(drawn, oscillatory);
logs = -inf(hi - lo + 1, 1);
for k = lo:hi
    logs(k - lo + 1) = log_choose(oscillatory, k) ...
        + log_choose(population - oscillatory, drawn - k) ...
        - log_choose(population, drawn);
end
largest = max(logs);
if ~isfinite(largest)
    p = nan;
    return
end
weights = exp(logs - largest);
weights = weights / sum(weights);
p = sum(weights((lo:hi).' >= a));
end


function value = log_choose(n, k)
if k < 0 || n < 0 || k > n
    value = -inf;
    return
end
value = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1);
end


function p = chi2_survival(x, df)
if ~isfinite(x) || ~isfinite(df) || df < 1 || x < 0
    p = nan;
    return
end
p = gammainc(x / 2, df / 2, 'upper');
end


function adjusted = holm_sidak_adjust(raw_p)
%HOLM_SIDAK_ADJUST  Step-down Holm-Sidak adjusted p-values.
%
% Non-finite raw p-values stay out of the family. For the rank-r test in
% a family of m, the Sidak value is 1 - (1 - p)^(m - r + 1). The adjusted
% values then cannot decrease as the raw p increases, and they are capped
% at 1. One finite p is returned unchanged.

adjusted = nan(size(raw_p));
usable = isfinite(raw_p);
p = raw_p(usable);
p = p(:);
n_tests = numel(p);
if n_tests == 0
    return
end
p = min(max(p, 0), 1);
[sorted_p, order] = sort(p);
sorted_adjusted = zeros(n_tests, 1);
running_max = 0;
for rank = 1:n_tests
    remaining = n_tests - rank + 1;
    this_adjusted = 1 - (1 - sorted_p(rank)) ^ remaining;
    if this_adjusted < running_max
        this_adjusted = running_max;
    end
    if this_adjusted > 1
        this_adjusted = 1;
    end
    running_max = this_adjusted;
    sorted_adjusted(rank) = this_adjusted;
end
unsorted = nan(n_tests, 1);
unsorted(order) = sorted_adjusted;
adjusted(usable) = unsorted;
end


function offset = violin_point_offset(n_points)
%VIOLIN_POINT_OFFSET  Fixed x offsets inside a violin.
%
% The same offset vector is used for both members of a pair, so the line
% and the two scatters share one x on each side.

if n_points <= 1
    offset = zeros(n_points, 1);
    return
end
offset = linspace(-0.12, 0.12, n_points).';
end


function [session, area, cluster_id] = neuron_identity_from_records(records, record_index)
%NEURON_IDENTITY_FROM_RECORDS  Session, area, and cluster id in population order.

n = numel(record_index);
session = cell(n, 1);
area = cell(n, 1);
cluster_id = nan(n, 1);
for i = 1:n
    record = records{record_index(i)};
    session{i} = record_text(record, 'session');
    area{i} = record_text(record, 'area');
    cluster_id(i) = record_number(record, 'cluster_id');
end
end


function text = record_text(record, field_name)
text = '';
if ~isfield(record, field_name) || isempty(record.(field_name))
    return
end
text = value_text(record.(field_name));
end


function number = record_number(record, field_name)
number = nan;
if ~isfield(record, field_name) || isempty(record.(field_name))
    return
end
value = record.(field_name);
if iscell(value)
    if isempty(value)
        return
    end
    value = value{1};
end
if isnumeric(value) && isfinite(value)
    number = double(value);
end
end


function text = value_text(value)
if iscell(value)
    if isempty(value)
        text = '';
        return
    end
    value = value{1};
end
text = strtrim(char(value(:).'));
end


function label = display_lock_label(raw)
name = lower(value_text(raw));
if strcmp(name, 'trough')
    label = 'Trough';
elseif strcmp(name, 'peak')
    label = 'Peak';
elseif strcmp(name, 'nonentrained')
    label = 'Unlocked';
elseif isempty(name)
    label = 'Non-classified';
else
    label = value_text(raw);
end
end


function bank = load_example_acg_bank(summary_path)
%LOAD_EXAMPLE_ACG_BANK  Mirrored sum-to-1 autocorrelograms keyed by neuron.

bank = empty_acg_bank('The locking summary file was not found, so the original autocorrelogram is not drawn.');
if isempty(summary_path) || exist(summary_path, 'file') ~= 2
    return
end
try
    stored = load(summary_path, 'power_acg');
catch
    bank.note = 'The locking summary file could not be read.';
    return
end
if isfield(stored, 'power_acg') && has_neuron_curve(stored.power_acg, 'high')
    bank = bank_from_curves(stored.power_acg.all_neurons, stored.power_acg.centers, 'high', ...
        'Top row: high-power sum-to-1 autocorrelogram, positive lags, both sides averaged.');
    return
end
try
    stored = load(summary_path, 'locking_acg');
catch
    bank.note = 'The locking summary file could not be read.';
    return
end
if isfield(stored, 'locking_acg') && has_neuron_curve(stored.locking_acg, 'all')
    bank = bank_from_curves(stored.locking_acg.all_neurons, stored.locking_acg.centers, 'all', ...
        'Top row: all-runs autocorrelogram, positive lags, both sides averaged. The high-power curve was not in the summary.');
    return
end
bank.note = 'The summary file has no per-neuron autocorrelogram.';
end


function ok = has_neuron_curve(summary, curve_name)
ok = isfield(summary, 'all_neurons') && isfield(summary.all_neurons, curve_name) ...
    && isfield(summary, 'centers') && ~isempty(summary.all_neurons.(curve_name));
end


function bank = empty_acg_bank(note)
bank.lag = zeros(1, 0);
bank.y = zeros(0, 0);
bank.session = {};
bank.cluster_id = zeros(0, 1);
bank.note = note;
end


function bank = bank_from_curves(block, centers, curve_name, note)
bank = empty_acg_bank(note);
centers = double(centers(:).');
curves = double(block.(curve_name));
if size(curves, 2) ~= numel(centers)
    bank.note = 'The autocorrelogram width does not match its lag axis.';
    return
end
[folded, pos_lag] = mirror_acg_matrix(curves, centers);
bank.lag = pos_lag;
bank.y = folded;
if iscell(block.session)
    bank.session = block.session(:);
else
    bank.session = cellstr(block.session);
end
bank.cluster_id = double(block.cluster_id(:));
end


function [folded, pos_lag] = mirror_acg_matrix(curves, centers)
%MIRROR_ACG_MATRIX  Average each positive lag with its negative-lag partner.

centers = double(centers(:).');
pos = find(centers > 0);
neg = find(centers < 0);
neg = neg(end:-1:1);
n = min(numel(pos), numel(neg));
pos = pos(1:n);
neg = neg(1:n);
if n < 1 || max(abs(centers(pos) + centers(neg))) > 1e-6
    folded = zeros(size(curves, 1), 0);
    pos_lag = zeros(1, 0);
    return
end
folded = (curves(:, pos) + curves(:, neg)) / 2;
pos_lag = centers(pos);
end


function [acg_lag, acg_y, found] = example_acg_curve(bank, session, cluster_id)
acg_lag = bank.lag;
acg_y = nan(1, numel(bank.lag));
found = false;
row = find_example_neuron(bank, session, cluster_id);
if isempty(row)
    return
end
acg_y = bank.y(row, :);
found = any(isfinite(acg_y));
end


function row = find_example_neuron(bank, session, cluster_id)
row = [];
if isempty(bank.cluster_id) || ~isfinite(cluster_id)
    return
end
session = value_text(session);
for i = 1:numel(bank.cluster_id)
    if abs(bank.cluster_id(i) - double(cluster_id)) < 1e-6 ...
            && strcmp(value_text(bank.session{i}), session)
        row = i;
        return
    end
end
end


function name = area_display_name(raw)
name = strtrim(value_text(raw));
if isempty(name)
    name = 'area unknown';
end
end


function tf = is_skipped_example_area(raw, skip_name)
tf = strcmpi(strtrim(value_text(raw)), strtrim(skip_name));
end


function [percent, observed, total, label] = area_osc_rate(neuron_area, osc_rows)
%AREA_OSC_RATE  Percent of analysed neurons in each area that are oscillatory.

preferred = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
names = cell(numel(neuron_area), 1);
for i = 1:numel(neuron_area)
    names{i} = area_display_name(neuron_area{i});
end
label = cell(1, 0);
for k = 1:numel(preferred)
    if any(strcmp(names, preferred{k}))
        label{end + 1} = preferred{k}; %#ok<AGROW>
    end
end
for i = 1:numel(names)
    if ~any(strcmp(label, names{i}))
        label{end + 1} = names{i}; %#ok<AGROW>
    end
end
n = numel(label);
observed = zeros(1, n);
total = zeros(1, n);
if isempty(osc_rows)
    osc_names = {};
else
    osc_names = names(osc_rows);
end
for k = 1:n
    total(k) = sum(strcmp(names, label{k}));
    observed(k) = sum(strcmp(osc_names, label{k}));
end
percent = nan(1, n);
has_cells = total > 0;
percent(has_cells) = 100 * observed(has_cells) ./ total(has_cells);
end


function position = example_positions(osc_rows, neuron_area, n_requested, skip_area)
%EXAMPLE_POSITIONS  Oscillatory rows to draw, skipping one area.
%
% position indexes osc_rows. The strongest oscillatory neuron in every
% other area is kept. More neurons are added in modulation order until
% n_requested, and the area representatives are kept even past that count.

if nargin < 3 || isempty(n_requested)
    n_requested = 0;
end
eligible = zeros(0, 1);
for i = 1:numel(osc_rows)
    if ~is_skipped_example_area(neuron_area{osc_rows(i)}, skip_area)
        eligible(end + 1, 1) = i; %#ok<AGROW>
    end
end
covered = cell(1, 0);
reps = zeros(0, 1);
for k = 1:numel(eligible)
    name = area_display_name(neuron_area{osc_rows(eligible(k))});
    if ~any(strcmp(covered, name))
        covered{end + 1} = name; %#ok<AGROW>
        reps(end + 1, 1) = eligible(k); %#ok<AGROW>
    end
end
position = reps(:);
target = max(n_requested, numel(reps));
for k = 1:numel(eligible)
    if numel(position) >= target
        break
    end
    if ~any(position == eligible(k))
        position(end + 1, 1) = eligible(k); %#ok<AGROW>
    end
end
position = sort(position(:));
end


function x_lim = example_lag_limits(lag)
lag = double(lag(isfinite(lag)));
if isempty(lag)
    x_lim = [0 1];
    return
end
span = max(max(lag) - min(lag), eps);
pad = 0.02 * span;
x_lim = [min(lag) - pad, max(lag) + pad];
if min(lag) >= 0
    x_lim(1) = 0;
end
end


function y_lim = example_acg_ylim(acg_lag, acg_y, lag_limits)
%EXAMPLE_ACG_YLIM  Scale one autocorrelogram without the near-zero peak.

acg_lag = double(acg_lag(:).');
acg_y = double(acg_y(:).');
shown = acg_lag >= lag_limits(1) & acg_lag <= lag_limits(2) & isfinite(acg_y);
tail = shown & acg_lag >= 0.05;
values = acg_y(tail);
if numel(values) < 5
    values = acg_y(shown);
end
if isempty(values)
    y_lim = [0 1];
    return
end
lo = min(values);
hi = max(values);
span = hi - lo;
if span <= 1e-6 * max(abs(hi), 1)
    span = max(abs(hi), eps);
end
pad = 0.08 * span;
y_lim = [lo - pad, hi + pad];
end


function text = example_neuron_summary(index, session, cluster_id, area_label, ...
    lock_label, modulation, acg_found, feature_table, feature_table_low, row, names)
header = sprintf('Neuron %d   session %s   id %g   area %s   lock %s   modulation %.4g', ...
    index, session, cluster_id, area_label, lock_label, modulation);
lines = {header};
if ~acg_found
    lines{end + 1} = 'Original autocorrelogram was not matched in the summary file.'; %#ok<AGROW>
end
lines{end + 1} = sprintf('%-42s %14s %14s', 'Feature', 'high power', 'low power'); %#ok<AGROW>
for j = 1:numel(names)
    high_value = table_feature_value(feature_table, row, names{j});
    low_value = table_feature_value(feature_table_low, row, names{j});
    lines{end + 1} = sprintf('%-42s %14s %14s', names{j}, ...
        stat_num(high_value, '%.4g'), stat_num(low_value, '%.4g')); %#ok<AGROW>
end
lines{end + 1} = ''; %#ok<AGROW>
text = strjoin(lines, newline);
end


function value = table_feature_value(feature_table, row, name)
value = nan;
if isempty(feature_table)
    return
end
if strcmp(name, 'Peak1 - Peak2')
    value = table_feature_value(feature_table, row, 'Peak1') ...
        - table_feature_value(feature_table, row, 'Peak2');
    return
end
names = feature_table.Properties.VariableNames;
if ~ismember(name, names)
    return
end
column = feature_table{row, name};
if isnumeric(column) && isscalar(column) && isfinite(column)
    value = double(column);
end
end


function fig = open_full_figure(name)
%OPEN_FULL_FIGURE  White figure that fills the screen.

fig = figure('Name', name, 'Color', 'w');
make_figure_fullscreen(fig)
end


function make_figure_fullscreen(fig)
set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
if isprop(fig, 'WindowState')
    set(fig, 'WindowState', 'maximized')
end
end


function square_plot_boxes(fig)
%SQUARE_PLOT_BOXES  Keep each panel square inside the full-screen figure.

ax = findall(fig, 'Type', 'axes');
for k = 1:numel(ax)
    this_class = class(ax(k));
    if ~isempty(strfind(this_class, 'ColorBar'))
        continue
    end
    tag = '';
    try
        tag = get(ax(k), 'Tag');
    catch
    end
    if strcmpi(tag, 'legend') || strcmpi(tag, 'Colorbar')
        continue
    end
    try
        axis(ax(k), 'square')
    catch
    end
end
charts = findall(fig, 'Type', 'heatmap');
for k = 1:numel(charts)
    square_chart_position(charts(k))
end
end


function square_chart_position(chart)
try
    pos = get(chart, 'Position');
catch
    return
end
if numel(pos) < 4 || min(pos(3:4)) <= 0
    return
end
side = min(pos(3), pos(4));
pos(1) = pos(1) + (pos(3) - side) / 2;
pos(2) = pos(2) + (pos(4) - side) / 2;
pos(3) = side;
pos(4) = side;
try
    set(chart, 'Position', pos)
catch
end
end


function text = join_summary_lines(varargin)
%JOIN_SUMMARY_LINES  One summary string from text and from a row or column of lines.

lines = cell(1, 0);
for k = 1:nargin
    piece = varargin{k};
    if iscell(piece)
        piece = piece(:).';
    else
        piece = {piece};
    end
    lines = [lines, piece]; %#ok<AGROW>
end
text = strjoin(lines, newline);
end


function save_script_figure(folder)
%SAVE_SCRIPT_FIGURE  Write the current figure as svg and jpg.

if nargin < 1 || isempty(folder)
    return
end
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
fig = gcf;
make_figure_fullscreen(fig)
square_plot_boxes(fig)
drawnow
name = regexprep(strtrim(get(fig, 'Name')), '[^\w]+', '_');
name = regexprep(name, '^_|_$', '');
if isempty(name)
    name = sprintf('figure_%d', double(get(fig, 'Number')));
end
set(fig, 'PaperPositionMode', 'auto')
svg_path = fullfile(folder, [name, '.svg']);
jpg_path = fullfile(folder, [name, '.jpg']);
print(fig, svg_path, '-dsvg')
print(fig, jpg_path, '-djpeg', '-r200')
fprintf('Saved %s\n', svg_path)
fprintf('Saved %s\n', jpg_path)
end


function write_summary_section(folder, title, body)
%WRITE_SUMMARY_SECTION  Replace one titled block in analysis_summary.txt.

if nargin < 1 || isempty(folder)
    return
end
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
if ~text_ends_with(body, newline)
    body = [body, newline];
end
block = ['-------', newline, title, newline, '--', newline, body, newline];
path = fullfile(folder, 'analysis_summary.txt');
if exist(path, 'file') == 2
    updated = replace_summary_block(fileread(path), title, block);
else
    updated = block;
end
fid = fopen(path, 'w');
if fid < 0
    error('Plot_acg_population_logit:summary', 'Could not write %s.', path);
end
fwrite(fid, updated);
fclose(fid);
fprintf('Wrote %s\n', path)
end


function updated = replace_summary_block(previous, title, block)
%REPLACE_SUMMARY_BLOCK  Swap the block whose title line matches.

if ~text_ends_with(previous, newline)
    previous = [previous, newline];
end
starts = regexp(previous, '(?m)^-------\r?\n');
if isempty(starts)
    updated = [previous, block];
    return
end
starts = [starts, numel(previous) + 1];
kept = cell(1, numel(starts) - 1);
replaced = false;
for k = 1:numel(starts) - 1
    chunk = previous(starts(k):starts(k + 1) - 1);
    lines = regexp(chunk, '\r?\n', 'split');
    this_title = '';
    if numel(lines) >= 2
        this_title = lines{2};
    end
    if strcmp(this_title, title)
        kept{k} = block;
        replaced = true;
    else
        kept{k} = chunk;
    end
end
if ~replaced
    kept{end + 1} = block;
end
updated = [kept{:}];
end


function tf = text_ends_with(text, suffix)
n_suffix = numel(suffix);
tf = numel(text) >= n_suffix && strcmp(text(end - n_suffix + 1:end), suffix);
end


function write_section_2_summary(folder, feature_names, n_neurons, n_complete, ...
    n_excluded, feature_correlation, pca_explained, perplexity_values, ...
    counts_before, counts_after, kmeans_sizes, ward_sizes, ...
    mean_sil_kmeans, mean_sil_ward, cluster_method_label, ...
    kmeans_modulation, ward_modulation, cluster_stats)

pca_lines = cell(numel(pca_explained), 1);
cumulative = cumsum(pca_explained(:));
for k = 1:numel(pca_explained)
    pca_lines{k} = sprintf('PC%-2d %7.2f%%    cumulative %7.2f%%', ...
        k, pca_explained(k), cumulative(k));
end
perplexity_text = sprintf('%g ', perplexity_values);
perplexity_text = strtrim(perplexity_text);
write_summary_section(folder, '2. Features kept for the split', join_summary_lines( ...
    ['Output folder: ', folder], ...
    'Figures from this section are saved there as jpg and svg.', ...
    sprintf('Neurons in the feature table: %d', n_neurons), ...
    sprintf('Complete rows used in PCA, t-SNE, and the split: %d', n_complete), ...
    sprintf('Rows excluded for a missing feature: %d', n_excluded), ...
    ['Features in the split: ', strjoin(feature_names, ', ')], ...
    'Correlation of those standardized features:', ...
    labeled_number_matrix('Feature', feature_names, feature_names, feature_correlation, '%.3f')))
write_summary_section(folder, '2. PCA', join_summary_lines( ...
    'PCA of the standardized features that entered the split.', ...
    'Variance explained:', pca_lines))
write_summary_section(folder, '2. t-SNE', join_summary_lines( ...
    't-SNE of the same standardized features. Distance: Euclidean. Seed: 42.', ...
    ['Perplexities: ', perplexity_text]))
write_summary_section(folder, '2. K-means and Ward', join_summary_lines( ...
    'Both splits use two groups. The count figure is drawn before relabeling.', ...
    'K-means by Ward, before relabeling:', ...
    labeled_number_matrix('K-means', {'cluster 1', 'cluster 2'}, ...
        {'Ward 1', 'Ward 2'}, counts_before, '%.0f'), ...
    'After relabeling, cluster 1 is the higher-modulation group in each split.', ...
    sprintf('K-means sizes: cluster 1 = %d, cluster 2 = %d', kmeans_sizes(1), kmeans_sizes(2)), ...
    sprintf('K-means mean standardized modulation: cluster 1 = %.4g, cluster 2 = %.4g', ...
        kmeans_modulation(1), kmeans_modulation(2)), ...
    sprintf('Ward sizes: cluster 1 = %d, cluster 2 = %d', ward_sizes(1), ward_sizes(2)), ...
    sprintf('Ward mean standardized modulation: cluster 1 = %.4g, cluster 2 = %.4g', ...
        ward_modulation(1), ward_modulation(2)), ...
    'K-means by Ward, after relabeling:', ...
    labeled_number_matrix('K-means', {'cluster 1', 'cluster 2'}, ...
        {'Ward 1', 'Ward 2'}, counts_after, '%.0f'), ...
    sprintf('Mean silhouette, Euclidean: k-means = %.4g, Ward = %.4g', ...
        mean_sil_kmeans, mean_sil_ward), ...
    ['The oscillatory group is ', cluster_method_label, ' cluster 1.']))
write_summary_section(folder, ...
    '2. Feature comparison, oscillatory versus non-oscillatory', ...
    feature_comparison_report(cluster_stats, 'unpaired'))
end


function write_section_4_summary(folder, column_title, cycle_n, cycle_r, cycle_p, ...
    lag, high_vs_low_h, high_vs_low_p, n_osc, ...
    n_osc_high, n_other_high, n_osc_low_pca, n_osc_low_tsne, n_omitted, ...
    embed_perplexity, paired_stats, lock_type_label, type_observed, type_total, type_percent)

cycle_lines = cell(1, numel(column_title));
for c = 1:numel(column_title)
    if ~isfinite(cycle_r(c))
        cycle_lines{c} = sprintf('%-32s n = %.0f    too few neurons for r', ...
            column_title{c}, cycle_n(c));
    else
        cycle_lines{c} = sprintf('%-32s n = %.0f    r = %.4g    uncorrected p = %.4g', ...
            column_title{c}, cycle_n(c), cycle_r(c), cycle_p(c));
    end
end
n_sig = sum(high_vs_low_h(:) == 1);
finite_p = isfinite(high_vs_low_p);
if any(finite_p)
    [min_p, local_index] = min(high_vs_low_p(finite_p));
    finite_index = find(finite_p);
    min_lag = lag(finite_index(local_index));
    lag_line = sprintf(['Paired t-tests across %d lag bins. %d bins have uncorrected p < 0.05. ', ...
        'Smallest uncorrected p = %.4g at lag %.3f s.'], numel(lag), n_sig, min_p, min_lag);
else
    lag_line = sprintf('No paired lag test was run. Oscillatory neurons: %d.', n_osc);
end
type_lines = cell(1, numel(lock_type_label));
for k = 1:numel(lock_type_label)
    type_lines{k} = sprintf('%-16s %d/%d   %.1f%%', lock_type_label{k}, ...
        type_observed(k), type_total(k), type_percent(k));
end
write_summary_section(folder, '4. Oscillatory logit figure', join_summary_lines( ...
    ['Output folder: ', folder], ...
    'Figures from this section are saved there as jpg and svg.', ...
    'Pearson correlation of the first cycle (peak minus trough) with the second cycle.', ...
    'The p on the figure is uncorrected. These three correlations are not in the Holm-Sidak family.', ...
    cycle_lines))
write_summary_section(folder, '4. High power minus low power across lag', join_summary_lines( ...
    'Each lag is a paired t-test of the same oscillatory neurons, high power minus low power.', ...
    'Red points on the figure mark uncorrected p < 0.05.', ...
    'These lag tests are not part of the Holm-Sidak correction used for the feature panels.', ...
    lag_line))
write_summary_section(folder, '4. Low-power points in the PCA and t-SNE', join_summary_lines( ...
    sprintf('PCA counts. Oscillatory high power: %d. Other high-power group: %d. Oscillatory low power with a finite projection: %d.', ...
        n_osc_high, n_other_high, n_osc_low_pca), ...
    sprintf('New t-SNE perplexity: %g.', embed_perplexity), ...
    sprintf('Low-power oscillatory neurons included: %d. Omitted for missing low-power features: %d.', ...
        n_osc_low_tsne, n_omitted)))
write_summary_section(folder, '4. Feature comparison, high versus low power', ...
    feature_comparison_report(paired_stats, 'paired'))
write_summary_section(folder, '4. Oscillatory share of each lock type', join_summary_lines( ...
    'Percent of each lock type that belongs to the oscillatory group.', ...
    'The label is observed/total.', ...
    type_lines))
end


function text = feature_comparison_report(rows, kind)
if strcmp(kind, 'paired')
    intro = { ...
        'Paired t-tests. The same oscillatory neurons are compared at high and low power.', ...
        't uses the mean of (high - low) and the standard deviation of those differences.', ...
        'df = n - 1.', ...
        'Cohen dz = mean(high - low) / standard deviation of (high - low).', ...
        'Power is the two-sided probability of p < 0.05 if the true effect equals the observed t.', ...
        'raw p is uncorrected. Holm-Sidak p is adjusted across every row with a finite raw p.', ...
        'The panel title shows only the Holm-Sidak p.', ...
        'A blank t means fewer than two complete pairs. That row is left out of the correction.'};
    header = sprintf('%-42s %6s %12s %12s %12s %10s %8s %10s %12s %14s %8s', ...
        'Feature', 'n', 'mean high', 'mean low', 'high-low', 't', 'df', ...
        'Cohen dz', 'raw p', 'Holm-Sidak p', 'power');
else
    intro = { ...
        'Unpaired equal-variance t-tests. Group 1 is the oscillatory cluster. Group 2 is the other cluster.', ...
        't = (mean 1 - mean 2) / pooled standard error.', ...
        'df = n1 + n2 - 2.', ...
        'Cohen d = (mean 1 - mean 2) / pooled standard deviation.', ...
        'Power is the two-sided probability of p < 0.05 if the true effect equals the observed t.', ...
        'raw p is uncorrected. Holm-Sidak p is adjusted across every row with a finite raw p.', ...
        'The panel title shows only the Holm-Sidak p.', ...
        'A blank t means one group has fewer than two finite values. That row is left out of the correction.'};
    header = sprintf('%-42s %-22s %6s %6s %12s %12s %12s %10s %8s %10s %12s %14s %8s', ...
        'Feature', 'role', 'n1', 'n2', 'mean 1', 'mean 2', 'mean 1-2', 't', 'df', ...
        'Cohen d', 'raw p', 'Holm-Sidak p', 'power');
end
lines = [intro, {header}];
if isempty(rows)
    lines{end + 1} = 'No feature was compared.';
else
    for j = 1:numel(rows)
        lines{end + 1} = format_feature_row(rows(j), kind); %#ok<AGROW>
    end
end
text = strjoin(lines, newline);
end


function line = format_feature_row(row, kind)
if strcmp(kind, 'paired')
    line = sprintf('%-42s %6s %12s %12s %12s %10s %8s %10s %12s %14s %8s', ...
        row.feature, stat_num(row.n1, '%.0f'), stat_num(row.mean1, '%.4g'), ...
        stat_num(row.mean2, '%.4g'), stat_num(row.difference, '%.4g'), ...
        stat_num(row.t, '%.4g'), stat_num(row.df, '%.0f'), ...
        stat_num(row.effect, '%.4g'), stat_num(row.raw_p, '%.4g'), ...
        stat_num(row.holm_p, '%.4g'), stat_num(row.power, '%.3f'));
else
    line = sprintf('%-42s %-22s %6s %6s %12s %12s %12s %10s %8s %10s %12s %14s %8s', ...
        row.feature, row.role, stat_num(row.n1, '%.0f'), stat_num(row.n2, '%.0f'), ...
        stat_num(row.mean1, '%.4g'), stat_num(row.mean2, '%.4g'), ...
        stat_num(row.difference, '%.4g'), stat_num(row.t, '%.4g'), ...
        stat_num(row.df, '%.0f'), stat_num(row.effect, '%.4g'), ...
        stat_num(row.raw_p, '%.4g'), stat_num(row.holm_p, '%.4g'), ...
        stat_num(row.power, '%.3f'));
end
end


function text = stat_num(value, fmt)
if isempty(value) || isnan(value)
    text = '';
    return
end
text = sprintf(fmt, value);
end


function text = labeled_number_matrix(corner, row_labels, col_labels, values, fmt)
header = sprintf('%-16s', corner);
for c = 1:numel(col_labels)
    header = [header, sprintf('%12s', col_labels{c})]; %#ok<AGROW>
end
lines = {header};
for r = 1:numel(row_labels)
    line = sprintf('%-16s', row_labels{r});
    for c = 1:numel(col_labels)
        line = [line, sprintf('%12s', sprintf(fmt, values(r, c)))]; %#ok<AGROW>
    end
    lines{end + 1} = line; %#ok<AGROW>
end
text = strjoin(lines, newline);
end


function s = empty_feature_stats()
s = struct('feature', '', 'role', '', 'n1', 0, 'n2', 0, ...
    'mean1', nan, 'mean2', nan, 'difference', nan, ...
    't', nan, 'df', nan, 'effect', nan, ...
    'raw_p', nan, 'holm_p', nan, 'power', nan);
end


function s = compare_two_samples(y1, y2, mode)
%COMPARE_TWO_SAMPLES  t, df, effect size, raw p, and observed power.

s = empty_feature_stats();
y1 = y1(:);
y2 = y2(:);
if strcmp(mode, 'paired')
    pair = isfinite(y1) & isfinite(y2);
    y1 = y1(pair);
    y2 = y2(pair);
else
    y1 = y1(isfinite(y1));
    y2 = y2(isfinite(y2));
end
s.n1 = numel(y1);
s.n2 = numel(y2);
if s.n1 > 0
    s.mean1 = mean(y1);
end
if s.n2 > 0
    s.mean2 = mean(y2);
end
if strcmp(mode, 'paired')
    if s.n1 < 2
        return
    end
    difference = y1 - y2;
    s.difference = mean(difference);
    s.df = s.n1 - 1;
    sd_difference = std(difference, 0);
    s = fill_t_effect(s, s.difference, sd_difference, s.df, 1 / s.n1);
    return
end
if isfinite(s.mean1) && isfinite(s.mean2)
    s.difference = s.mean1 - s.mean2;
end
if s.n1 < 2 || s.n2 < 2
    return
end
s.df = s.n1 + s.n2 - 2;
pooled_variance = ((s.n1 - 1) * var(y1, 0) + (s.n2 - 1) * var(y2, 0)) / s.df;
if pooled_variance < 0
    pooled_variance = 0;
end
s = fill_t_effect(s, s.difference, sqrt(pooled_variance), s.df, 1 / s.n1 + 1 / s.n2);
end


function s = fill_t_effect(s, difference, scale, df, variance_factor)
%FILL_T_EFFECT  Shared t, effect size, p, and power for one comparison.
%
% scale is the pooled standard deviation, or the standard deviation of
% the paired differences. variance_factor is 1/n1+1/n2, or 1/n for pairs.

if scale == 0
    if difference == 0
        s.t = 0;
        s.effect = 0;
        s.raw_p = 1;
        s.power = observed_t_power(0, df);
    else
        s.t = sign(difference) * inf;
        s.effect = s.t;
        s.raw_p = 0;
        s.power = 1;
    end
    return
end
s.effect = difference / scale;
s.t = difference / (scale * sqrt(variance_factor));
s.raw_p = two_sided_t_p(s.t, df);
s.power = observed_t_power(s.t, df);
end


function p = two_sided_t_p(tstat, df)
if ~isfinite(df) || df <= 0
    p = nan;
    return
end
if ~isfinite(tstat)
    p = 0;
    return
end
if tstat == 0
    p = 1;
    return
end
x = df / (df + tstat ^ 2);
p = betainc(x, df / 2, 0.5);
p = min(max(p, 0), 1);
end


function tc = t_critical(tail_p, df)
x = betaincinv(2 * tail_p, df / 2, 0.5);
tc = sqrt(df * (1 - x) / x);
end


function power = observed_t_power(tstat, df)
%OBSERVED_T_POWER  Two-sided power at alpha 0.05 for the observed t.

if ~isfinite(df) || df <= 0
    power = nan;
    return
end
if ~isfinite(tstat)
    power = 1;
    return
end
tcrit = t_critical(0.025, df);
power = noncentral_t_cdf(-tcrit, df, tstat) + (1 - noncentral_t_cdf(tcrit, df, tstat));
power = min(max(power, 0), 1);
end


function cdf = noncentral_t_cdf(t, df, ncp)
%NONCENTRAL_T_CDF  Probability integral of the noncentral t distribution.

if ~isfinite(df) || df <= 0 || ~isfinite(ncp)
    cdf = nan;
    return
end
if t == inf
    cdf = 1;
    return
end
if t == -inf
    cdf = 0;
    return
end
vmax = df + 20 * sqrt(2 * df) + 80;
n_points = 20000;
step = vmax / n_points;
v = ((1:n_points) - 0.5) * step;
log_pdf = (df / 2 - 1) * log(v) - v / 2 - (df / 2) * log(2) - gammaln(df / 2);
z = t * sqrt(v / df) - ncp;
phi = 0.5 * erfc(-z / sqrt(2));
cdf = sum(phi .* exp(log_pdf)) * step;
cdf = min(max(cdf, 0), 1);
end
