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
% Set output_folder in the workspace before running a section. That path is
% where each figure and analysis_summary.txt are written. This script does
% not choose it. The neuron table, the population matrices, and the
% lock-type model cache are written to
% Data\Analysis results\acg population logit. A cache or matrix saved by
% an earlier run is read from output_folder or from spike train parsing.
% n_example_neurons, n_neurons_per_figure, example_smooth_ms, and
% plot_example_neurons_only are set at the start of section 4.
% n_example_neurons chooses how many of the strongest oscillatory neurons
% are drawn. The examples also keep the strongest oscillatory neuron from
% each analysed area except VLPAG, so that count can be larger than
% n_example_neurons. n_neurons_per_figure is how many neurons go on each
% figure. Each neuron is one column: both autocorrelograms on the top row,
% both shuffle percentiles on the bottom row. example_smooth_ms is the
% Gaussian smoothing width in milliseconds. plot_example_neurons_only
% draws those neurons and leaves the other section 4 figures unrun.
% classifier_method is set just before the lock-type classifier.
% The default is naive Bayes, LDA, and a medium Gaussian SVM.
% 'all' adds the tree and multinomial logit. One name keeps that classifier.
% The features are those whose lock-type ANOVA has Holm-Sidak p below 0.05
% in all classified neurons or in the oscillatory neurons only.
% Each classifier is naive Bayes, LDA, or a medium Gaussian SVM. Each one is
% trained twice, with inverse-frequency weights and with equal weights, on
% every trough, peak, and unlocked neuron and again on the oscillatory
% neurons only. A joint model mixes the two class probabilities. The mix
% uses each model's training-fold accuracy and the ACG features only.
% n_accuracy_permutations label shuffles are drawn once per population and
% reused for all nine models. The default is 1000. The model accuracy is
% the mean of the three lock-type accuracies. A bar is filled when that
% permutation p is below 0.05. A finished fit is written to
% lock_type_model_cache.mat under Analysis results\acg population logit. The same code version,
% features, neurons, and permutation count load that file instead of
% fitting again. lock_model_cache_version is 1 for this fit.
%
% Backup. After an ACG run, paste the block below into the command window.
% It copies the ACG inputs and the saved models into
% Data\Analysis results\acg unin back up. Run it again to refresh that
% folder. To put a folder back, copy it from there into Analysis results.
% The large parts are spike train parsing and phase locking data.
% acg population logit holds lock_type_model_cache.mat,
% acg_population_logit_prep.mat, and oscillatory_neuron_labels.csv.
% The shuffle checkpoint and the power summary are in spike train parsing.
% A copy of those model files still sitting in output_folder is copied too.
%
% data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
% results_folder = fullfile(data_root, 'Analysis results');
% backup_folder = fullfile(results_folder, 'acg unin back up');
% if exist(backup_folder, 'dir') ~= 7
%     mkdir(backup_folder);
% end
% source_names = { ...
%     'spike train parsing', ...
%     'phase locking data', ...
%     'play power acg', ...
%     'session acg shuffle', ...
%     'acg population logit'};
% for i = 1:numel(source_names)
%     source = fullfile(results_folder, source_names{i});
%     if exist(source, 'dir') ~= 7
%         fprintf('Missing, left out: %s\n', source);
%         continue
%     end
%     destination = fullfile(backup_folder, source_names{i});
%     if exist(destination, 'dir') ~= 7
%         mkdir(destination);
%     end
%     copyfile(source, destination);
%     fprintf('Backed up %s\n', destination);
% end
% model_names = {'lock_type_model_cache.mat', 'acg_population_logit_prep.mat', 'oscillatory_neuron_labels.csv'};
% if exist('output_folder', 'var') && ~isempty(output_folder) && exist(output_folder, 'dir') == 7
%     model_destination = fullfile(backup_folder, 'acg population logit');
%     if exist(model_destination, 'dir') ~= 7
%         mkdir(model_destination);
%     end
%     for i = 1:numel(model_names)
%         model_file = fullfile(output_folder, model_names{i});
%         if exist(model_file, 'file') == 2
%             copyfile(model_file, fullfile(model_destination, model_names{i}));
%             fprintf('Backed up %s\n', model_file);
%         end
%     end
% end
% fprintf('Backup folder: %s\n', backup_folder);
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
% Section 5 stores the mean channel delta power of the runs that form
% each neuron's low-power autocorrelogram and each neuron's high-power
% autocorrelogram. For neurons with both tails, it regresses the
% high-minus-low power gap on the high-minus-low change in the six
% oscillatory features and on the high-minus-low change in firing rate.
% The intercept is the gap expected when those predictors do not change.
% The scale is the less skewed of the raw gap and the log gap. The same
% fit is repeated on the non-oscillatory neurons alone. Pearson
% correlations of the predictors are saved with each fit. Each predictor
% with p below 0.05 is drawn against the power gap with the other
% predictors fixed at their means. Each feature change is also compared
% with the firing-rate change, and each significant feature is scattered
% against that firing-rate change. Section 6 tests
% delta power across areas for the high-power tail, the low-power tail,
% and all runs. Neurons on one channel are averaged, then channels in one
% area and session are averaged. A session intercept is removed before
% area is tested, and the figure shows those session-adjusted values.
% When that area test is below 0.05, pairwise comparisons are Holm-Sidak
% adjusted and significant pairs are drawn with a line and stars.
% It also reads the saved
% power summary and compares, within each neuron, the firing rate of the
% low-power runs, the high-power runs, and the runs in the all-runs
% autocorrelogram. Each pair is tested with a paired t-test, a
% signed-rank test, a sign test, and a paired permutation
% Kolmogorov-Smirnov test.

%% 1. Load the checkpoint and build the imagesc matrices
% output_folder is an input. Set it in the workspace before running this
% section. Figures and analysis_summary.txt are written there.
if ~exist('output_folder', 'var') || isempty(output_folder)
    error('Plot_acg_population_logit:outputFolder', ...
        ['Set output_folder to the folder where figures should be saved, ' ...
        'then run this section.']);
end
[output_folder, estimate_folder] = logit_folders(output_folder);
this_dir = fileparts(mfilename('fullpath'));
figures_dir = fileparts(this_dir);
repo_root = fileparts(figures_dir);
addpath(repo_root);
addpath(figures_dir);
addpath(this_dir);
if exist('align_population_percentiles', 'file') ~= 2
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

acg_method = 'overlap';
acg_count = 'pairs';
n_shuffles_default = 500;
% Mean percentile inside this lag window sets the stored row order.
% 0.2 s is 200 ms. Section 2 sorts the figure again.
sort_window_s = [0.200 0.400];

data_root = experiment_data_root();
parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
checkpoint = fullfile(parsing_folder, shuffle_checkpoint_filename(n_shuffles_default, acg_method, acg_count));
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
[high_power_delta, low_power_delta, all_power_delta] = ...
    stratum_mean_delta(records, neuron_record_index);

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

prep_file = fullfile(estimate_folder, 'acg_population_logit_prep.mat');
% acg2features scores the lags it is given. The checkpoint axis is
% two-sided, and the rhythmicity function uses the positive side.
positive_lag = lag > 0;
feature_table = acg2features(imagesc_sum_power_high(:, positive_lag), lag(positive_lag));

save(prep_file, 'lag', 'imagesc_logit', 'percentile', 'sort_score', ...
    'imagesc_sum_all', 'imagesc_sum_power_high', 'imagesc_sum_power_low', ...
    'neuron_type', 'neuron_record_index', 'neuron_session', 'neuron_area', ...
    'neuron_cluster_id', 'feature_table', 'high_power_delta', 'low_power_delta', ...
    'all_power_delta', ...
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
    sprintf('Neurons with a finite high-power mean delta power: %d', sum(isfinite(high_power_delta))), ...
    sprintf('Neurons with a finite low-power mean delta power: %d', sum(isfinite(low_power_delta))), ...
    sprintf('Neurons with a finite all-runs mean delta power: %d', sum(isfinite(all_power_delta))), ...
    'high_power_delta is the mean, across that neuron''s high-power runs, of the squared median delta-band envelope on the neuron channel.', ...
    'low_power_delta is that mean for the runs at or below the neuron''s 25th percentile.', ...
    'all_power_delta is that mean over every run with a finite score. Older checkpoints leave it missing, and section 6 then reads score_all from the power summary.', ...
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
if ~exist('output_folder', 'var')
    output_folder = '';
end
[output_folder, estimate_folder] = logit_folders(output_folder);
this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
addpath(this_dir);
if ~exist('feature_table', 'var') || ~exist('neuron_type', 'var') ...
        || ~exist('imagesc_sum_power_high', 'var')
    prep_file = logit_prep_file(estimate_folder, output_folder);
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

% This PCA and the t-SNE below use high-power neurons only. Section 4 fits
% a separate PCA and t-SNE that also include the oscillatory neurons at
% low power, and the section 4 figures share that plane.
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
% modulation is placed first. This section also writes
% oscillatory_neuron_labels.csv, with the session and neuron id of every
% population neuron and whether that neuron was classified as oscillatory.
if ~exist('output_folder', 'var')
    output_folder = '';
end
[output_folder, estimate_folder] = logit_folders(output_folder);
this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
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
if ~exist('neuron_session', 'var') || ~exist('neuron_cluster_id', 'var') ...
        || numel(neuron_session) ~= numel(neuron_type) ...
        || numel(neuron_cluster_id) ~= numel(neuron_type)
    if ~exist('records', 'var') || ~exist('neuron_record_index', 'var')
        error('Plot_acg_population_logit:cluster', ...
            'Session and neuron id are not in the workspace. Run section 1 first.');
    end
    [neuron_session, ~, neuron_cluster_id] = ...
        neuron_identity_from_records(records, neuron_record_index);
end
oscillatory_label_file = save_oscillatory_label_table(estimate_folder, ...
    neuron_session, neuron_cluster_id, osc_member, complete_rows);
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
    ['Neuron table: ', oscillatory_label_file], ...
    'Each row is one population neuron. session and neuron_id match the original spike train.', ...
    'oscillatory is true only for cluster 1. classification is oscillatory, non-oscillatory, or not classified when the neuron did not enter the split.', ...
    'Share of each lock type that is in that group:', ...
    type_lines))

%% 4. Figures
% Re-run this section after editing the figure. Section 3 does not need
% to run again while shown_logit is still in the workspace. The low-power
% PCA figure and the paired feature figure also need the section 2 PCA.
if ~exist('output_folder', 'var')
    output_folder = '';
end
[output_folder, estimate_folder] = logit_folders(output_folder);
n_example_neurons = 10;
n_neurons_per_figure = 3;
example_smooth_ms = 0.05;
plot_example_neurons_only = false;
this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
addpath(this_dir);
if ~exist('shown_logit', 'var')
    error('Plot_acg_population_logit:shown', ...
        'shown_logit is not in the workspace. Run section 3 first.');
end
if ~plot_example_neurons_only
% Open intervals, in seconds, for the cycle scatter below.
trough_windows_s = [0.05 0.15; 0.35 0.45];
peak_windows_s = [0.25 0.35; 0.45 0.55];

% Symmetric color range for the logit image. 1 is about the 9th and 91st
% percentiles. 2 is about the 1st and 99th.
logit_color_limit = 2;
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
stratum_columns = {'power_high', 'power_low'};
stratum_titles = {'High 25% power', 'Low 25% power'};
if ~exist('osc_member', 'var') || numel(osc_member) ~= numel(neuron_type)
    osc_member = false(numel(neuron_type), 1);
    osc_member(osc_rows) = true;
end
if ~exist('osc_score', 'var') || numel(osc_score) ~= numel(neuron_type)
    osc_score = nan(numel(neuron_type), 1);
end
if ~exist('imagesc_logit', 'var')
    error('Plot_acg_population_logit:shown', ...
        'imagesc_logit is not in the workspace. Run section 1 first.');
end
non_osc_rows = oscillatory_member_rows(~osc_member, osc_score);
[cycle_n, cycle_r, cycle_p] = draw_power_logit_pair( ...
    imagesc_logit, osc_rows, stratum_columns, stratum_titles, ...
    'oscillatory', n_osc, n_before_osc, lag, logit_color_limit, logit_map, ...
    first_peak_lag, first_trough_lag, second_peak_lag, second_trough_lag, ...
    scale_title{1}, ...
    [cluster_method_label, ' group with the higher modulation, strongest first.'], ...
    output_folder);
[non_cycle_n, non_cycle_r, non_cycle_p] = draw_power_logit_pair( ...
    imagesc_logit, non_osc_rows, stratum_columns, stratum_titles, ...
    'non-oscillatory', numel(non_osc_rows), n_before_osc, lag, ...
    logit_color_limit, logit_map, ...
    first_peak_lag, first_trough_lag, second_peak_lag, second_trough_lag, ...
    scale_title{1}, ...
    'Neurons outside that group, strongest remaining modulation first.', ...
    output_folder);
non_cycle_lines = cycle_summary_lines(stratum_titles, non_cycle_n, non_cycle_r, non_cycle_p);
write_summary_section(output_folder, '4. Non-oscillatory logit figure', join_summary_lines( ...
    ['Output folder: ', output_folder], ...
    sprintf('Non-oscillatory neurons: %d of %d.', numel(non_osc_rows), n_before_osc), ...
    'Same two-by-two layout as the oscillatory logit figure: high power and low power, image above mean logit.', ...
    'Rows are ordered by the remaining modulation score, strongest first.', ...
    'Pearson correlation of the first cycle with the second cycle is reported here and is not drawn.', ...
    non_cycle_lines))

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

% One PCA and one t-SNE include the high-power neurons and the oscillatory
% neurons at low power. Every embedding figure below uses that plane.
% Figures that omit the low-power points keep the same coordinates.
if ~exist('feature_mean', 'var') || ~exist('feature_sigma', 'var') ...
        || ~exist('feature_names', 'var') || ~exist('complete_rows', 'var') ...
        || ~exist('feature_table', 'var') || ~exist('imagesc_sum_power_low', 'var') ...
        || ~exist('osc_rows', 'var') || ~exist('Xz', 'var')
    error('Plot_acg_population_logit:pca', ...
        'The high-power features are not in the workspace. Run sections 2 and 3 first.');
end
positive_lag = lag > 0;
feature_table_low = acg2features( ...
    imagesc_sum_power_low(:, positive_lag), lag(positive_lag));
X_low = feature_table_low{complete_rows, feature_names};
Xz_low = (X_low - feature_mean) ./ feature_sigma;
osc_in_complete = ismember(complete_rows(:), osc_rows(:));
low_member = osc_in_complete & all(isfinite(Xz_low), 2);
X_embed = [Xz; Xz_low(low_member, :)];
group_id = ones(size(Xz, 1), 1);
group_id(~osc_in_complete) = 2;
group_id = [group_id; 3 * ones(sum(low_member), 1)];
[~, power_score, ~, ~, power_explained] = pca(X_embed);
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
fprintf(['Shared embedding. Low-power oscillatory neurons: %d. ', ...
    'Omitted for missing low-power features: %d.\n'], ...
    sum(low_member), sum(osc_in_complete) - sum(low_member))
power_color = [0.84 0.37 0.00; 0.45 0.28 0.62; 0.00 0.45 0.70];
power_label = { ...
    'Oscillatory, high power', ...
    'Non-oscillatory, high power', ...
    'Oscillatory, low power'};
power_caption = { ...
    'Orange circles are oscillatory neurons at high power. Purple circles are the other high-power neurons.', ...
    'Bands use the same colors. Axis limits come from the shared plane, including the low-power points.'};
draw_power_embedding(power_score(:, 1:2), group_id, [1 2], power_color, power_label, ...
    output_folder, 'PCA oscillatory and non-oscillatory', ...
    sprintf('PC1 (%.1f%%)', power_explained(1)), ...
    sprintf('PC2 (%.1f%%)', power_explained(2)), ...
    [power_caption, {'Low-power oscillatory neurons are in this PCA and are not drawn.'}])
draw_power_embedding(power_score(:, 1:2), group_id, [1 2 3], power_color, power_label, ...
    output_folder, 'PCA oscillatory and non-oscillatory with low power', ...
    sprintf('PC1 (%.1f%%)', power_explained(1)), ...
    sprintf('PC2 (%.1f%%)', power_explained(2)), ...
    [power_caption, {'Blue triangles are the same oscillatory neurons at low power.'}])
draw_power_embedding(Y_power, group_id, [1 2], power_color, power_label, ...
    output_folder, 't-SNE oscillatory and non-oscillatory', ...
    't-SNE 1', 't-SNE 2', ...
    [power_caption, {sprintf('Perplexity = %g. Low-power oscillatory neurons are in this t-SNE and are not drawn.', ...
        embed_perplexity)}])
draw_power_embedding(Y_power, group_id, [1 2 3], power_color, power_label, ...
    output_folder, 't-SNE oscillatory and non-oscillatory with low power', ...
    't-SNE 1', 't-SNE 2', ...
    [power_caption, {sprintf('Perplexity = %g. Blue triangles are the same oscillatory neurons at low power.', ...
        embed_perplexity)}])

n_high = size(Xz, 1);
lock_embedding_labels = lock_labels_for_rows(neuron_type, complete_rows);
lock_embedding_keep = { ...
    true(n_high, 1), ...
    ismember(complete_rows(:), osc_rows(:))};
lock_embedding_names = {'all neurons', 'oscillatory neurons'};
lock_xy = {power_score(1:n_high, 1:2), Y_power(1:n_high, :)};
lock_x_label = { ...
    sprintf('PC1 (%.1f%%)', power_explained(1)), ...
    't-SNE 1'};
lock_y_label = { ...
    sprintf('PC2 (%.1f%%)', power_explained(2)), ...
    't-SNE 2'};
lock_plane = {'PCA', 't-SNE'};
for plane = 1:2
    for set_n = 1:2
        if ~any(lock_embedding_keep{set_n})
            fprintf('No %s were available for the lock-type %s.\n', ...
                lock_embedding_names{set_n}, lock_plane{plane})
            continue
        end
        draw_lock_type_embedding( ...
            lock_xy{plane}(lock_embedding_keep{set_n}, :), ...
            lock_embedding_labels(lock_embedding_keep{set_n}), output_folder, ...
            ['Lock type in the shared ', lock_plane{plane}, ', ', lock_embedding_names{set_n}], ...
            lock_x_label{plane}, lock_y_label{plane}, ...
            {['Black points are the ', lock_embedding_names{set_n}, ...
                ' in the shared ', lock_plane{plane}, '.'], ...
            'That plane also contains the oscillatory neurons at low power. Those points are not drawn here.', ...
            'Blue is trough, red is peak, and green is unlocked.'})
    end
end

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

open_full_figure('Oscillatory cells within each lock type')
shown_percent = type_percent;
shown_percent(~isfinite(shown_percent)) = 0;
hold on
type_legend = cell(1, 4);
for k = 1:4
    bar(k, shown_percent(k), 'FaceColor', lock_bar_color(lock_type_label{k}), 'BarWidth', 0.8)
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
lock_anova = lock_type_feature_anova(feature_table, neuron_type, ...
    area_feature_names, area_feature_titles);
draw_lock_feature_figure(lock_anova, output_folder, n_descriptor_cols, ...
    'Oscillatory features by lock type', ...
    {'High-power features by lock type, all classified neurons.', ...
    'Non-classified cells are left out.', ...
    'Points are neurons. Bars are means. Stars mark pairwise differences after a significant ANOVA.'})
write_summary_section(output_folder, '4. Oscillatory features by lock type', ...
    area_feature_report(lock_anova, { ...
    'One-way ANOVA of each high-power feature, with lock type as the factor.', ...
    'Every classified high-power neuron is included.', ...
    'The groups are trough, peak, and unlocked. Non-classified cells are excluded.', ...
    'A lock type is included in a feature only when it has at least two finite values.', ...
    'Holm-Sidak adjusts the ANOVA p-values across features.', ...
    'Pairwise equal-variance t-tests are run only when that ANOVA Holm-Sidak p is below 0.05.', ...
    'Those pairs are then Holm-Sidak adjusted within the feature. Stars use that adjusted p.'}))
osc_lock_anova = lock_type_feature_anova(feature_table, neuron_type, ...
    area_feature_names, area_feature_titles, osc_rows);
draw_lock_feature_figure(osc_lock_anova, output_folder, n_descriptor_cols, ...
    'Oscillatory features by lock type, oscillatory neurons', ...
    {'High-power features by lock type, oscillatory neurons only.', ...
    'Oscillatory neurons are cluster 1 from section 3. Non-classified cells are left out.', ...
    'Points are neurons. Bars are means. Stars mark pairwise differences after a significant ANOVA.'})
write_summary_section(output_folder, '4. Oscillatory features by lock type, oscillatory neurons', ...
    area_feature_report(osc_lock_anova, { ...
    'One-way ANOVA of each high-power feature, with lock type as the factor.', ...
    'Only oscillatory neurons are included. They are cluster 1 from section 3.', ...
    'The groups are trough, peak, and unlocked. Non-classified cells are excluded.', ...
    'A lock type is included in a feature only when it has at least two finite values.', ...
    'Holm-Sidak adjusts the ANOVA p-values across features.', ...
    'Pairwise equal-variance t-tests are run only when that ANOVA Holm-Sidak p is below 0.05.', ...
    'Those pairs are then Holm-Sidak adjusted within the feature. Stars use that adjusted p.'}))
% Subplot numbers on the feature panels, in reading order:
% 1 Modulation1, 2 Modulation2, 6 Baseline, 12 Peak2, 13 Trough1,
% 15 Peak1 - Peak2.
two_way_panels = [1 2 6 12 13 15];
n_two_way = numel(two_way_panels);
two_way_samples = cell(n_two_way, 1);
two_way_with = cell(n_two_way, 1);
two_way_without = cell(n_two_way, 1);
for j = 1:n_two_way
    feature_index = two_way_panels(j);
    two_way_samples{j} = two_way_sample(feature_table, neuron_area, neuron_type, ...
        area_feature_names{feature_index});
    two_way_samples{j}.feature = area_feature_titles{feature_index};
    [area_idx, type_idx] = two_way_indices(two_way_samples{j});
    n_area = numel(two_way_samples{j}.area_names);
    n_type = numel(two_way_samples{j}.type_names);
    two_way_with{j} = two_way_anova(two_way_samples{j}.y, area_idx, type_idx, ...
        n_area, n_type, true);
    two_way_without{j} = two_way_anova(two_way_samples{j}.y, area_idx, type_idx, ...
        n_area, n_type, false);
end
two_way_with = assign_two_way_holm(two_way_with);
two_way_without = assign_two_way_holm(two_way_without);
draw_two_way_panels(two_way_samples, two_way_with, output_folder, ...
    'Two-factor ANOVA with interaction', ...
    {'Area and cell type, with an interaction.', ...
    'Blue is trough, red is peak, green is unlocked.', ...
    'Each title line is the Holm-Sidak p of that term.'})
draw_two_way_panels(two_way_samples, two_way_without, output_folder, ...
    'Two-factor ANOVA without interaction', ...
    {'Area and cell type, without an interaction.', ...
    'Blue is trough, red is peak, green is unlocked.', ...
    'Each title line is the Holm-Sidak p of that term.'})
write_summary_section(output_folder, '4. Two-factor ANOVA with interaction', ...
    two_way_report(two_way_samples, two_way_with, { ...
    'Features are subplot numbers 1, 2, 6, 12, 13, and 15.', ...
    'They are Modulation1, Modulation2, Baseline, Peak2, Trough1, and Peak1 - Peak2.', ...
    'The model is area + cell type + area x cell type.', ...
    'Each term is tested after the other terms, with effect coding.', ...
    'VLPAG, neurons without a stored area, and non-classified cells are excluded.', ...
    'Holm-Sidak adjusts every term in this figure together.', ...
    'A term is listed with its statistic only when that adjusted p is below 0.05.'}))
write_summary_section(output_folder, '4. Two-factor ANOVA without interaction', ...
    two_way_report(two_way_samples, two_way_without, { ...
    'Features are subplot numbers 1, 2, 6, 12, 13, and 15.', ...
    'They are Modulation1, Modulation2, Baseline, Peak2, Trough1, and Peak1 - Peak2.', ...
    'The model is area + cell type, with no interaction.', ...
    'Each main effect is tested after the other main effect.', ...
    'VLPAG, neurons without a stored area, and non-classified cells are excluded.', ...
    'Holm-Sidak adjusts every term in this figure together.', ...
    'A term is listed with its statistic only when that adjusted p is below 0.05.'}))
classifier_method = {'svm', 'naivebayes', 'lda'};
classifier_method_name = classifier_method_text(classifier_method);
n_accuracy_permutations = 1000;
[classifier_names, classifier_titles, classifier_sources] = significant_lock_features( ...
    lock_anova, osc_lock_anova, area_feature_names, area_feature_titles);
lock_class_names = {'Trough', 'Peak', 'Unlocked'};
classifier_methods = {'svm', 'naivebayes', 'lda'};
fprintf('Classifier features (%d): %s\n', numel(classifier_titles), ...
    strjoin(classifier_titles, ', '))
if isempty(classifier_names)
    fprintf(['No lock-type feature had Holm-Sidak p below 0.05, ', ...
        'so the classifiers were not run.\n'])
    write_summary_section(output_folder, '4. Lock-type classifier', join_summary_lines( ...
        'No classifier was run.', ...
        'A feature is used when its lock-type ANOVA Holm-Sidak p is below 0.05 for all classified neurons or for oscillatory neurons only.', ...
        'Neither comparison marked a feature.'))
else
% Bump this when the fit itself changes. A file from another version is not loaded.
lock_model_cache_version = 1;
lock_model_cache_file = lock_cache_file(estimate_folder, output_folder);
classifier_populations = { ...
    struct('name', 'All classified neurons', 'rows', (1:size(feature_table, 1)).'), ...
    struct('name', 'Oscillatory neurons', 'rows', osc_rows(:))};
[classifier_populations, cache_loaded] = try_load_lock_model_cache( ...
    lock_model_cache_file, lock_model_cache_version, 42, n_accuracy_permutations, ...
    classifier_methods, lock_class_names, classifier_names, ...
    feature_table, neuron_type, classifier_populations);
if ~cache_loaded
    rng(42);
    fprintf(['Lock-type models: %d shuffles, shared by the weighted, unweighted, ', ...
        'and joint fits of each classifier.\n'], n_accuracy_permutations)
    for p = 1:numel(classifier_populations)
        [classifier_X, classifier_y] = lock_classifier_data(feature_table, neuron_type, ...
            classifier_names, classifier_populations{p}.rows);
        classifier_populations{p}.X = classifier_X;
        classifier_populations{p}.y = classifier_y;
        fprintf('%s: %d neurons, %d features.\n', ...
            classifier_populations{p}.name, numel(classifier_y), size(classifier_X, 2))
        shuffled_labels = shared_label_shuffles(classifier_y, n_accuracy_permutations);
        model_grid = cell(3, numel(classifier_methods));
        for m = 1:numel(classifier_methods)
            method_label = classifier_method_label(classifier_methods{m});
            fprintf('%s, %s: real labels.\n', classifier_populations{p}.name, method_label)
            observed = classifier_weight_trio(classifier_X, classifier_y, ...
                lock_class_names, classifier_methods{m});
            null_accuracy = nan(n_accuracy_permutations, 3);
            null_recall = nan(n_accuracy_permutations, 3, numel(lock_class_names));
            if observed{1}.ok || observed{2}.ok
                shuffle_clock = tic;
                for shuffle = 1:n_accuracy_permutations
                    shuffled = classifier_weight_trio(classifier_X, ...
                        shuffled_labels(shuffle, :).', lock_class_names, classifier_methods{m});
                    for row = 1:3
                        if shuffled{row}.ok
                            null_accuracy(shuffle, row) = shuffled{row}.accuracy;
                            null_recall(shuffle, row, :) = shuffled{row}.recall;
                        end
                    end
                    if shuffle == 1 || mod(shuffle, 25) == 0 || shuffle == n_accuracy_permutations
                        fprintf('%s, %s: shuffle %d of %d (%.0f s).\n', ...
                            classifier_populations{p}.name, method_label, shuffle, ...
                            n_accuracy_permutations, toc(shuffle_clock))
                    end
                end
            else
                fprintf('%s, %s: %s\n', classifier_populations{p}.name, ...
                    method_label, observed{1}.note)
            end
            for row = 1:3
                recall_null = reshape(null_recall(:, row, :), n_accuracy_permutations, []);
                model_grid{row, m} = apply_shared_null(observed{row}, ...
                    null_accuracy(:, row), recall_null);
            end
        end
        classifier_populations{p}.grid = model_grid;
        classifier_populations{p}.results = model_grid(:);
    end
    write_lock_model_cache(lock_model_cache_file, lock_model_cache_version, 42, ...
        n_accuracy_permutations, classifier_methods, lock_class_names, ...
        classifier_names, classifier_populations)
end
for p = 1:numel(classifier_populations)
    report_lock_model_population(classifier_populations{p}, lock_class_names)
    draw_lock_classifier(classifier_populations{p}.grid, output_folder, classifier_method_name, ...
        ['Lock type classifier, ', classifier_populations{p}.name], ...
        classifier_populations{p}.name)
end
draw_accuracy_null(classifier_populations, output_folder, classifier_method_name)
draw_class_accuracy_bars(classifier_populations, output_folder)
write_summary_section(output_folder, '4. Lock-type classifier', ...
    lock_classifier_report(classifier_titles, classifier_sources, classifier_method_name, ...
    n_accuracy_permutations, classifier_populations))
end
end

if ~exist('feature_table_low', 'var') || size(feature_table_low, 1) ~= numel(neuron_type)
    if exist('imagesc_sum_power_low', 'var') && exist('lag', 'var')
        positive_lag = lag > 0;
        feature_table_low = acg2features( ...
            imagesc_sum_power_low(:, positive_lag), lag(positive_lag));
    end
end
if ~exist('acg_method', 'var') || isempty(acg_method)
    acg_method = 'overlap';
end
if ~exist('acg_count', 'var') || isempty(acg_count)
    acg_count = 'pairs';
end
if exist('checkpoint', 'var')
    example_summary_path = fullfile(fileparts(checkpoint), ...
        locking_summary_filename(false, acg_method, acg_count));
else
    example_summary_path = '';
end
example_bank = load_example_acg_bank(example_summary_path);
example_pos = example_positions(osc_rows, neuron_area, n_example_neurons, 'VLPAG');
n_example = numel(example_pos);
example_names = { ...
    'Modulation1', 'Modulation2', 'PeakInterval', 'TroughInterval', ...
    'OscillationScore', 'Baseline', ...
    'DominantFrequency', 'PeakFrequency', 'TroughFrequency', ...
    'Persistence', 'Peak1', 'Peak2', 'Trough1', 'Trough2', 'Peak1 - Peak2'};
n_example_rows = 2;
neurons_per_figure = max(1, round(n_neurons_per_figure));
n_example_cols = neurons_per_figure;
high_color = [0.75 0.22 0.18];
low_color = [0.15 0.35 0.70];
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
            [acg_lag, acg_high, acg_low, found_high, found_low] = example_stratum_curves( ...
                example_bank, neuron_session{pop_row}, neuron_cluster_id(pop_row));
            acg_high = smooth_example_trace(acg_lag, acg_high, example_smooth_ms);
            acg_low = smooth_example_trace(acg_lag, acg_low, example_smooth_ms);
            pct_high = smooth_example_trace(lag, ...
                shown_percentile.sum.power_high(rank, :), example_smooth_ms);
            pct_low = smooth_example_trace(lag, ...
                shown_percentile.sum.power_low(rank, :), example_smooth_ms);
            shared_ylim = example_acg_ylim(acg_lag, [acg_high(:), acg_low(:)], lag_limits);
            subplot(n_example_rows, n_example_cols, k)
            hold on
            if found_high
                plot(acg_lag, acg_high, 'Color', high_color, 'LineWidth', 1.4)
            end
            if found_low
                plot(acg_lag, acg_low, 'Color', low_color, 'LineWidth', 1.4)
            end
            if found_high || found_low
                xlim(lag_limits)
                ylim(shared_ylim)
                if k == 1
                    legend_names = {};
                    if found_high
                        legend_names{end + 1} = 'High 25% power'; %#ok<AGROW>
                    end
                    if found_low
                        legend_names{end + 1} = 'Low 25% power'; %#ok<AGROW>
                    end
                    legend(legend_names, 'Location', 'northeast', 'Box', 'off')
                end
            else
                text(0.5, 0.5, 'ACG not in summary', 'HorizontalAlignment', 'center')
                axis([0 1 0 1])
            end
            ylabel('Sum-to-1')
            title({sprintf('%s, %s', this_area, lock_label), ...
                sprintf('%s id %g', neuron_session{pop_row}, neuron_cluster_id(pop_row))}, ...
                'Interpreter', 'none')
            box off
            subplot(n_example_rows, n_example_cols, n_example_cols + k)
            hold on
            plot(lag, pct_high, 'Color', high_color, 'LineWidth', 1.4)
            plot(lag, pct_low, 'Color', low_color, 'LineWidth', 1.4)
            xlim(lag_limits)
            ylim([0 100])
            xlabel('Lag (s)')
            ylabel('Percentile')
            box off
            example_lines{i} = example_neuron_summary(i, neuron_session{pop_row}, ...
                neuron_cluster_id(pop_row), this_area, lock_label, osc_score(pop_row), ...
                found_high, found_low, feature_table, feature_table_low, pop_row, example_names);
        end
        sgtitle({sprintf('Neurons %d-%d of %d. VLPAG is left out.', i0, i1, n_example), ...
            'Each column is one neuron. Top: high-power and low-power autocorrelograms.', ...
            'Bottom: high-power and low-power shuffle percentiles.', ...
            sprintf('Gaussian smooth, T = %g ms.', example_smooth_ms)}, ...
            'Interpreter', 'none')
        save_script_figure(output_folder)
    end
    write_summary_section(output_folder, '4. Example oscillatory neurons', join_summary_lines( ...
        sprintf('n_example_neurons = %d. Drawn: %d. Neurons per figure: %d.', ...
            n_example_neurons, n_example, neurons_per_figure), ...
        sprintf('Figures: %d. Each neuron is one column with two rows.', n_example_figures), ...
        'Top row: high-power and low-power autocorrelograms on one axes.', ...
        'Bottom row: high-power and low-power sum-to-1 shuffle percentiles on one axes.', ...
        sprintf('Gaussian smoothing width T = %g ms.', example_smooth_ms), ...
        'VLPAG is left out. Each other area that has an oscillatory neuron contributes its strongest one.', ...
        'Further neurons are the next highest modulation scores outside VLPAG, until n_example_neurons is reached.', ...
        'If those areas already exceed n_example_neurons, every area representative is still drawn.', ...
        'The drawn order is strongest modulation first.', ...
        example_bank.note, ...
        'The autocorrelogram y-axis is taken from that neuron, using lags at and after 0.05 s, so the bin beside zero does not flatten the trace.', ...
        'Features are the high-power and low-power values for that neuron.', ...
        example_lines))
else
    write_summary_section(output_folder, '4. Example oscillatory neurons', join_summary_lines( ...
        sprintf('n_example_neurons = %d. Drawn: 0.', n_example_neurons), ...
        'No oscillatory neuron outside VLPAG was available.'))
end
if ~plot_example_neurons_only
write_section_4_summary(output_folder, stratum_titles, cycle_n, cycle_r, cycle_p, ...
    lag, high_vs_low_h, high_vs_low_p, n_osc, ...
    sum(osc_in_complete), sum(~osc_in_complete), sum(low_member), ...
    sum(low_member), sum(osc_in_complete) - sum(low_member), embed_perplexity, ...
    paired_stats, lock_type_label, type_observed, type_total, type_percent)
end

%% 5. Within-neuron delta-power gap and the feature changes
% Re-run this section after section 3. One row per neuron that has both
% tails and a positive firing rate in each tail. The outcome is the
% high-power mean minus the low-power mean. Each predictor is that same
% high-minus-low change in an autocorrelogram feature, or in firing rate.
% The intercept is the gap expected when those predictors do not change.
% The same equation is fit again on the non-oscillatory neurons.
if ~exist('output_folder', 'var')
    output_folder = '';
end
[output_folder, estimate_folder] = logit_folders(output_folder);
this_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(this_dir));
addpath(this_dir);
if ~exist('feature_table', 'var') || ~exist('neuron_type', 'var') ...
        || ~exist('complete_rows', 'var') || ~exist('osc_rows', 'var')
    error('Plot_acg_population_logit:powerModel', ...
        'The oscillatory split is not in the workspace. Run sections 2 and 3 first.');
end
if ~exist('high_power_delta', 'var') || ~exist('low_power_delta', 'var') ...
        || numel(high_power_delta) ~= numel(neuron_type) ...
        || numel(low_power_delta) ~= numel(neuron_type)
    if ~exist('records', 'var') || ~exist('neuron_record_index', 'var')
        error('Plot_acg_population_logit:powerModel', ...
            'The low-power and high-power means are not in the workspace. Run section 1 first.');
    end
    [high_power_delta, low_power_delta, all_power_delta] = ...
        stratum_mean_delta(records, neuron_record_index);
end
if ~exist('feature_table_low', 'var') || size(feature_table_low, 1) ~= numel(neuron_type)
    if ~exist('imagesc_sum_power_low', 'var') || ~exist('lag', 'var')
        error('Plot_acg_population_logit:powerModel', ...
            'The low-power autocorrelogram is not in the workspace. Run section 1 first.');
    end
    if exist('acg2features', 'file') ~= 2
        run(fullfile(fileparts(fileparts(this_dir)), 'add_repo_paths.m'));
    end
    positive_lag = lag > 0;
    feature_table_low = acg2features( ...
        imagesc_sum_power_low(:, positive_lag), lag(positive_lag));
end
power_feature_names = { ...
    'Modulation1', 'Modulation2', 'PeakInterval', 'TroughInterval', ...
    'OscillationScore', 'Baseline'};
power_feature_names = power_feature_names( ...
    ismember(power_feature_names, feature_table.Properties.VariableNames) ...
    & ismember(power_feature_names, feature_table_low.Properties.VariableNames));
if isempty(power_feature_names)
    error('Plot_acg_population_logit:powerModel', ...
        'None of the oscillatory features are in both feature tables.');
end
power_rows = complete_rows(:);
power_rows = power_rows(power_rows >= 1 & power_rows <= numel(neuron_type));
x_high = feature_table{power_rows, power_feature_names};
x_low = feature_table_low{power_rows, power_feature_names};
y_high = high_power_delta(power_rows);
y_low = low_power_delta(power_rows);
osc = double(ismember(power_rows, osc_rows(:)));
high_ok = all(isfinite(x_high), 2) & isfinite(y_high) & y_high > 0;
low_ok = all(isfinite(x_low), 2) & isfinite(y_low) & y_low > 0;
if ~exist('acg_method', 'var') || isempty(acg_method)
    acg_method = 'overlap';
end
if ~exist('acg_count', 'var') || isempty(acg_count)
    acg_count = 'pairs';
end
if ~exist('checkpoint', 'var')
    checkpoint = '';
end
if ~exist('neuron_session', 'var') || ~exist('neuron_cluster_id', 'var') ...
        || numel(neuron_session) ~= numel(neuron_type) ...
        || numel(neuron_cluster_id) ~= numel(neuron_type)
    error('Plot_acg_population_logit:powerModel', ...
        'Session and cluster id are not in the workspace. Run section 1 first.');
end
rate_summary_path = stratum_summary_path(checkpoint, acg_method, acg_count);
[rate_low_pop, rate_high_pop] = population_stratum_rates( ...
    rate_summary_path, neuron_session, neuron_cluster_id);
rate_low_rows = rate_low_pop(power_rows);
rate_high_rows = rate_high_pop(power_rows);
rate_ok = isfinite(rate_low_rows) & isfinite(rate_high_rows);
both = high_ok & low_ok & rate_ok;
if sum(both) < 3
    error('Plot_acg_population_logit:powerModel', ...
        'Too few neurons with both power tails and a positive firing rate in each tail.');
end
gap_x = [x_high(both, :) - x_low(both, :), rate_high_rows(both) - rate_low_rows(both)];
gap_names = [power_feature_names, {'FiringRate'}];
gap_raw = y_high(both) - y_low(both);
gap_log = log(y_high(both)) - log(y_low(both));
gap_osc = osc(both);
n_rate_missing = sum(high_ok & low_ok & ~rate_ok);
power_gap = fit_within_neuron_gap(gap_raw, gap_log, gap_x, gap_names, ...
    'All neurons with both tails');
n_power_osc = sum(gap_osc == 1);
fprintf('Neurons with only a high-power tail: %d. Only a low-power tail: %d.\n', ...
    sum(high_ok & ~low_ok), sum(low_ok & ~high_ok))
fprintf('Both power tails but no positive firing rate in both tails: %d.\n', n_rate_missing)
fprintf('Oscillatory neurons: %d. Non-oscillatory neurons: %d.\n', ...
    n_power_osc, power_gap.n - n_power_osc)
print_gap_model(power_gap)
draw_power_gap_scale(gap_raw, gap_log, power_gap.scale, output_folder, ...
    'Delta power gap scale')
draw_power_model_coefficients(power_gap.fit, power_gap.outcome_label, output_folder, ...
    'Within-neuron delta power gap')
power_gap_partial = draw_ceteris_paribus(power_gap, output_folder, ...
    'Ceteris paribus, delta power gap')
draw_regressor_correlation(power_gap.correlation, power_gap.names, output_folder, ...
    'Delta power gap predictor correlations')
write_summary_section(output_folder, '5. Within-neuron delta-power gap', join_summary_lines( ...
    ['Output folder: ', output_folder], ...
    'One row per neuron that has both power tails and a positive firing rate in each tail.', ...
    'A run is high power when its delta power is at or above the neuron''s own 75th percentile.', ...
    'A run is low power when its delta power is at or below the neuron''s own 25th percentile.', ...
    'Delta power of a run is the square of the median delta-band envelope on the neuron channel.', ...
    'The outcome is the high-power mean minus the low-power mean.', ...
    'On the log scale that outcome is log(high) minus log(low).', ...
    'Each autocorrelogram predictor is the high-power value minus the low-power value of that feature.', ...
    'Features are Modulation1, Modulation2, PeakInterval, TroughInterval, OscillationScore, and Baseline.', ...
    'FiringRate is the high-power firing rate minus the low-power firing rate, in spikes per second.', ...
    'Each tail rate is the spike count in that tail divided by the summed duration of its runs.', ...
    'The intercept is the power gap expected when every predictor change, including the firing-rate change, is zero.', ...
    'A slope is the change in the power gap associated with a one-unit within-neuron change in that predictor.', ...
    'A predictor with p below 0.05 is plotted against the power gap while every other predictor is fixed at its mean.', ...
    power_gap_partial, ...
    'Each neuron has one row, so the standard errors are ordinary least squares.', ...
    'The model F tests the slopes together. The intercept is the reference gap and is left out of that test.', ...
    'The correlation table is the Pearson correlation of the predictors in this fit. The intercept is not included.', ...
    sprintf('Neurons in the model: %d. Oscillatory: %d. Non-oscillatory: %d.', ...
        power_gap.n, n_power_osc, power_gap.n - n_power_osc), ...
    sprintf('Neurons with only a high-power tail: %d. Only a low-power tail: %d.', ...
        sum(high_ok & ~low_ok), sum(low_ok & ~high_ok)), ...
    sprintf('Both power tails but no positive firing rate in both tails: %d.', n_rate_missing), ...
    sprintf('Mean raw gap, oscillatory: %.4g. Non-oscillatory: %.4g.', ...
        finite_mean(gap_raw(gap_osc == 1)), finite_mean(gap_raw(gap_osc == 0))), ...
    sprintf('Mean log gap, oscillatory: %.4g. Non-oscillatory: %.4g.', ...
        finite_mean(gap_log(gap_osc == 1)), finite_mean(gap_log(gap_osc == 0))), ...
    gap_fit_lines(power_gap)))
non_osc = gap_osc == 0;
non_gap = fit_within_neuron_gap(gap_raw(non_osc), gap_log(non_osc), ...
    gap_x(non_osc, :), gap_names, 'Non-oscillatory neurons with both tails');
fprintf('Non-oscillatory gap model. Oscillatory neurons left out: %d.\n', n_power_osc)
print_gap_model(non_gap)
draw_power_gap_scale(gap_raw(non_osc), gap_log(non_osc), non_gap.scale, output_folder, ...
    'Delta power gap scale, non-oscillatory')
draw_power_model_coefficients(non_gap.fit, non_gap.outcome_label, output_folder, ...
    'Within-neuron delta power gap, non-oscillatory')
non_gap_partial = draw_ceteris_paribus(non_gap, output_folder, ...
    'Ceteris paribus, delta power gap, non-oscillatory')
draw_regressor_correlation(non_gap.correlation, non_gap.names, output_folder, ...
    'Delta power gap predictor correlations, non-oscillatory')
write_summary_section(output_folder, ...
    '5c. Within-neuron delta-power gap, non-oscillatory neurons', join_summary_lines( ...
    ['Output folder: ', output_folder], ...
    'This is the same within-neuron equation as section 5, fit only on neurons outside the oscillatory group.', ...
    'The oscillatory group is cluster 1 from section 3.', ...
    'One row per non-oscillatory neuron that has both power tails and a positive firing rate in each tail.', ...
    'The outcome, the six feature changes, and FiringRate are defined as in section 5.', ...
    'The scale is chosen again from the skewness of this subset.', ...
    'A predictor with p below 0.05 is plotted against the power gap while every other predictor is fixed at its mean.', ...
    non_gap_partial, ...
    'The correlation table is the Pearson correlation of the predictors in this subset.', ...
    sprintf('Non-oscillatory neurons in this fit: %d. Oscillatory neurons left out: %d.', ...
        non_gap.n, n_power_osc), ...
    gap_fit_lines(non_gap)))

% Does the within-neuron firing-rate change track the feature changes?
rate_link = feature_change_vs_rate(gap_x, gap_names);
fprintf('Feature change versus firing-rate change, n = %d neurons.\n', rate_link.n)
for i = 1:numel(rate_link.lines)
    fprintf('%s\n', rate_link.lines{i})
end
draw_feature_rate_link(rate_link, output_folder)
rate_scatter_note = draw_feature_rate_scatter(rate_link, output_folder)
write_summary_section(output_folder, '5d. Feature changes and the firing-rate change', ...
    join_summary_lines( ...
    ['Output folder: ', output_folder], ...
    'Each row is one neuron with both power tails and a positive firing rate in each tail.', ...
    'A feature change is the high-power autocorrelogram value minus the low-power value.', ...
    'The firing-rate change is the high-power rate minus the low-power rate, in spikes per second.', ...
    'Pearson and Spearman correlations test whether those two changes move together.', ...
    'The slope is from the within-neuron regression of the feature change on the firing-rate change, with an intercept.', ...
    'Holm-Sidak adjusted p-values are for the Pearson tests across the features in this table.', ...
    'A feature whose Pearson p is below 0.05 is scattered against the firing-rate change.', ...
    rate_scatter_note, ...
    sprintf('Neurons: %d.', rate_link.n), ...
    rate_link.lines))

% Firing rates are already stored on the power-split summary. This block
% only reads them. It does not rebuild the autocorrelograms.
if ~exist('acg_method', 'var') || isempty(acg_method)
    acg_method = 'overlap';
end
if ~exist('acg_count', 'var') || isempty(acg_count)
    acg_count = 'pairs';
end
if ~exist('checkpoint', 'var')
    checkpoint = '';
end
if ~exist('rate_summary_path', 'var') || isempty(rate_summary_path)
    rate_summary_path = stratum_summary_path(checkpoint, acg_method, acg_count);
end
if ~exist('neuron_session', 'var') || ~exist('neuron_cluster_id', 'var') ...
        || numel(neuron_session) ~= numel(neuron_cluster_id)
    rate_session = {};
    rate_cluster_id = [];
else
    rate_session = neuron_session;
    rate_cluster_id = neuron_cluster_id;
end
rate_report = load_stratum_firing_rates(rate_summary_path, rate_session, rate_cluster_id);
for i = 1:numel(rate_report.lines)
    fprintf('%s\n', rate_report.lines{i})
end
if numel(rate_report.low) >= 2
    draw_stratum_rate_figure(rate_report.low, rate_report.high, output_folder)
end
write_summary_section(output_folder, '5b. Firing rates of the power strata', ...
    join_summary_lines(rate_report.lines))


%% 6. Nested comparison of delta power across areas
% Re-run this section after section 1. area_power_conditions chooses the
% set of runs. 'high' is the high-power mean, 'low' is the low-power mean,
% and 'all' is the mean over every run with a finite score.
area_power_conditions = {'high', 'low', 'all'};
if ~exist('output_folder', 'var')
    output_folder = '';
end
[output_folder, estimate_folder] = logit_folders(output_folder);
if ~exist('neuron_area', 'var') || ~exist('neuron_session', 'var') ...
        || ~exist('neuron_cluster_id', 'var') || ~exist('high_power_delta', 'var') ...
        || ~exist('low_power_delta', 'var')
    error('Plot_acg_population_logit:areaNest', ...
        'The population and its delta power are not in the workspace. Run section 1 first.');
end
if ~exist('parsing_folder', 'var') || isempty(parsing_folder)
    parsing_folder = default_parsing_folder();
end
if ~exist('acg_method', 'var') || isempty(acg_method)
    acg_method = 'overlap';
end
if ~exist('acg_count', 'var') || isempty(acg_count)
    acg_count = 'pairs';
end
if ~exist('checkpoint', 'var')
    checkpoint = '';
end
neuron_channel = channels_for_neurons(parsing_folder, neuron_session, neuron_cluster_id);
for condition_index = 1:numel(area_power_conditions)
    area_power_condition = area_power_conditions{condition_index};
    if strcmp(area_power_condition, 'low')
        area_power_values = low_power_delta;
    elseif strcmp(area_power_condition, 'high')
        area_power_values = high_power_delta;
    elseif strcmp(area_power_condition, 'all')
        if ~exist('all_power_delta', 'var') || numel(all_power_delta) ~= numel(neuron_session) ...
                || ~any(isfinite(all_power_delta))
            all_power_delta = population_summary_field( ...
                stratum_summary_path(checkpoint, acg_method, acg_count), ...
                neuron_session, neuron_cluster_id, 'score_all');
        end
        area_power_values = all_power_delta;
    else
        error('Plot_acg_population_logit:areaNest', ...
            'area_power_conditions accepts ''high'', ''low'', and ''all''.');
    end
    area_nest = nested_area_power(area_power_values, neuron_area, neuron_session, ...
        neuron_channel, area_power_condition);
    fprintf('Nested area comparison, %s power.\n', area_power_condition)
    for line_index = 1:numel(area_nest.lines)
        fprintf('%s\n', area_nest.lines{line_index})
    end
    if area_nest.ok
        draw_nested_area_power(area_nest, output_folder)
    end
    write_summary_section(output_folder, ...
        ['6. Nested area comparison of delta power, ', area_power_condition], ...
        join_summary_lines(area_nest.lines))
end


function draw_power_embedding(xy, group_id, shown, colors, labels, folder, figure_name, x_label, y_label, caption)
%DRAW_POWER_EMBEDDING  Colored oscillatory groups on one shared plane.
%
% shown selects group ids. 1 is oscillatory high power, 2 is the other
% high-power neurons, and 3 is oscillatory low power. Axis limits always
% use every row of xy, including groups that are not drawn.

open_full_figure(figure_name)
hold on
draw_order = [2 1 3];
draw_order = draw_order(ismember(draw_order, shown));
for g = draw_order
    draw_group_density(xy, xy(group_id == g, :), colors(g, :))
end
handles = gobjects(1, 0);
legend_labels = {};
for g = draw_order
    pts = xy(group_id == g, :);
    pts = pts(all(isfinite(pts), 2), :);
    if isempty(pts)
        continue
    end
    if g == 3
        marker = '^';
        marker_size = 28;
    else
        marker = 'o';
        marker_size = 18;
    end
    handles(end + 1) = scatter(pts(:, 1), pts(:, 2), marker_size, colors(g, :), ...
        'filled', 'Marker', marker, 'MarkerFaceAlpha', 0.85, ...
        'MarkerEdgeColor', 'none', 'DisplayName', labels{g}); %#ok<AGROW>
    legend_labels{end + 1} = labels{g}; %#ok<AGROW>
end
limits = embedding_axis_limits(xy);
if ~isempty(limits)
    xlim(limits(1:2))
    ylim(limits(3:4))
end
xlabel(x_label, 'Interpreter', 'none')
ylabel(y_label, 'Interpreter', 'none')
if ~isempty(handles)
    legend(handles, legend_labels, 'Location', 'best', 'Box', 'off', 'Interpreter', 'none')
end
grid on
box off
sgtitle(caption, 'Interpreter', 'none')
save_script_figure(folder)
end


function draw_group_density(all_xy, group_xy, face_color)
%DRAW_GROUP_DENSITY  Filled density bands for one group on the current axes.

all_xy = all_xy(all(isfinite(all_xy), 2), :);
group_xy = group_xy(all(isfinite(group_xy), 2), :);
if size(group_xy, 1) < 3 || isempty(all_xy)
    return
end
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
if isempty(density)
    return
end
density_grid = reshape(density, size(xg));
if ~any(density_grid(:) > 0)
    return
end
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
if numel(thresholds) < 2
    return
end
sorted_thresholds = sort(thresholds(:));
n_thresholds = numel(sorted_thresholds);
for band_n = (n_thresholds - 1):-1:1
    contourf(x1, x2, density_grid, ...
        sorted_thresholds([band_n band_n + 1]), ...
        'LineWidth', 0.5, 'FaceColor', face_color, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.28 + 0.42 * band_n / n_thresholds, ...
        'HandleVisibility', 'off');
end
end


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


function sample = two_way_sample(feature_table, neuron_area, neuron_type, feature_name)
%TWO_WAY_SAMPLE  Neurons kept for an area-by-cell-type ANOVA.

values = feature_column(feature_table, feature_name);
n = numel(values);
areas = cell(n, 1);
types = cell(n, 1);
keep = isfinite(values);
kept_types = {'Trough', 'Peak', 'Unlocked'};
for i = 1:n
    if i <= numel(neuron_area)
        areas{i} = area_display_name(neuron_area{i});
    else
        areas{i} = 'area unknown';
    end
    if i <= numel(neuron_type)
        types{i} = display_lock_label(neuron_type{i});
    else
        types{i} = 'Non-classified';
    end
    if ~keep(i)
        continue
    end
    if strcmpi(areas{i}, 'VLPAG') || strcmp(areas{i}, 'area unknown')
        keep(i) = false;
    elseif ~any(strcmp(types{i}, kept_types))
        keep(i) = false;
    end
end
sample.y = values(keep);
sample.area = areas(keep);
sample.type = types(keep);
sample.feature = feature_name;
preferred = {'SC', 'LPAG', 'DLPAG', 'DR'};
sample.area_names = {};
for k = 1:numel(preferred)
    if any(strcmp(sample.area, preferred{k}))
        sample.area_names{end + 1} = preferred{k}; %#ok<AGROW>
    end
end
for i = 1:numel(sample.area)
    if ~any(strcmp(sample.area_names, sample.area{i}))
        sample.area_names{end + 1} = sample.area{i}; %#ok<AGROW>
    end
end
sample.type_names = {};
for k = 1:numel(kept_types)
    if any(strcmp(sample.type, kept_types{k}))
        sample.type_names{end + 1} = kept_types{k}; %#ok<AGROW>
    end
end
end


function [area_idx, type_idx] = two_way_indices(sample)
area_idx = zeros(numel(sample.y), 1);
type_idx = zeros(numel(sample.y), 1);
for i = 1:numel(sample.y)
    area_idx(i) = find(strcmp(sample.area_names, sample.area{i}), 1);
    type_idx(i) = find(strcmp(sample.type_names, sample.type{i}), 1);
end
end


function model = two_way_anova(y, area_idx, type_idx, n_area, n_type, with_interaction)
%TWO_WAY_ANOVA  Area and cell type, with or without their interaction.
%
% Effect coding. With an interaction, each term is tested after the
% others. Without it, each main effect is tested after the other one.

model.ok = false;
model.note = '';
model.n = numel(y);
model.terms = empty_anova_term();
model.terms(1) = [];
if n_area < 2 || n_type < 2 || numel(y) < 4
    model.note = 'Fewer than two areas or two cell types had enough neurons.';
    return
end
y = y(:);
area_idx = area_idx(:);
type_idx = type_idx(:);
X = twoway_design(area_idx, type_idx, n_area, n_type, with_interaction);
n_a = n_area - 1;
n_t = n_type - 1;
area_cols = 2:(1 + n_a);
type_cols = (2 + n_a):(1 + n_a + n_t);
term_names = {'Area', 'Cell type'};
term_cols = {area_cols, type_cols};
if with_interaction
    term_names{end + 1} = 'Interaction'; %#ok<AGROW>
    term_cols{end + 1} = (2 + n_a + n_t):size(X, 2); %#ok<AGROW>
end
[rss_full, rank_full] = ols_rss(X, y);
df_error = numel(y) - rank_full;
if rank_full < size(X, 2)
    model.note = 'Some area and cell-type combinations are missing, so the model rank is reduced.';
end
if df_error < 1
    model.note = 'Too few neurons for the residual degrees of freedom.';
    return
end
model.ok = true;
model.terms = repmat(empty_anova_term(), numel(term_names), 1);
for k = 1:numel(term_names)
    keep = true(1, size(X, 2));
    keep(term_cols{k}) = false;
    [rss_reduced, rank_reduced] = ols_rss(X(:, keep), y);
    df1 = rank_full - rank_reduced;
    ss = rss_reduced - rss_full;
    scale = max([1, abs(rss_reduced), abs(rss_full)]);
    if ss < 0 && ss > -1e-8 * scale
        ss = 0;
    end
    model.terms(k).name = term_names{k};
    model.terms(k).df1 = df1;
    model.terms(k).df2 = df_error;
    if df1 < 1 || ss < 0
        continue
    end
    if rss_full <= 1e-10 * scale
        if ss > 1e-10 * scale
            model.terms(k).F = inf;
            model.terms(k).p = 0;
            model.terms(k).eta = 1;
        else
            model.terms(k).F = 0;
            model.terms(k).p = 1;
            model.terms(k).eta = 0;
        end
        continue
    end
    model.terms(k).F = (ss / df1) / (rss_full / df_error);
    model.terms(k).p = f_survival(model.terms(k).F, df1, df_error);
    model.terms(k).eta = ss / (ss + rss_full);
end
end


function s = empty_anova_term()
s = struct('name', '', 'F', nan, 'df1', nan, 'df2', nan, ...
    'p', nan, 'holm_p', nan, 'eta', nan);
end


function X = twoway_design(area_idx, type_idx, n_area, n_type, with_interaction)
n = numel(area_idx);
n_a = n_area - 1;
n_t = n_type - 1;
n_col = 1 + n_a + n_t;
if with_interaction
    n_col = n_col + n_a * n_t;
end
X = zeros(n, n_col);
X(:, 1) = 1;
for a = 1:n_a
    code = zeros(n, 1);
    code(area_idx == a) = 1;
    code(area_idx == n_area) = -1;
    X(:, 1 + a) = code;
end
for t = 1:n_t
    code = zeros(n, 1);
    code(type_idx == t) = 1;
    code(type_idx == n_type) = -1;
    X(:, 1 + n_a + t) = code;
end
if with_interaction
    col = 2 + n_a + n_t;
    for a = 1:n_a
        for t = 1:n_t
            X(:, col) = X(:, 1 + a) .* X(:, 1 + n_a + t);
            col = col + 1;
        end
    end
end
end


function [rss, xrank] = ols_rss(X, y)
y = y(:);
if isempty(X) || size(X, 2) < 1
    rss = sum((y - mean(y)) .^ 2);
    xrank = 0;
    return
end
beta = X \ y;
rss = sum((y - X * beta) .^ 2);
singular = svd(X);
if isempty(singular)
    xrank = 0;
    return
end
tol = max(size(X)) * eps(max(singular));
xrank = sum(singular > tol);
end


function models = assign_two_way_holm(models)
raw = [];
which_model = [];
which_term = [];
for j = 1:numel(models)
    for k = 1:numel(models{j}.terms)
        raw(end + 1) = models{j}.terms(k).p; %#ok<AGROW>
        which_model(end + 1) = j; %#ok<AGROW>
        which_term(end + 1) = k; %#ok<AGROW>
    end
end
adjusted = holm_sidak_adjust(raw);
for i = 1:numel(raw)
    models{which_model(i)}.terms(which_term(i)).holm_p = adjusted(i);
end
end


function draw_two_way_panels(samples, models, folder, figure_name, caption)
open_full_figure(figure_name)
type_gap = 0.24;
for j = 1:numel(samples)
    subplot(2, 3, j)
    hold on
    sample = samples{j};
    model = models{j};
    n_area = numel(sample.area_names);
    n_type = numel(sample.type_names);
    offsets = ((1:n_type) - (n_type + 1) / 2) * type_gap;
    for a = 1:n_area
        for t = 1:n_type
            mask = strcmp(sample.area, sample.area_names{a}) ...
                & strcmp(sample.type, sample.type_names{t});
            y = sample.y(mask);
            if isempty(y)
                continue
            end
            if numel(y) == 1
                x = a + offsets(t);
            else
                x = a + offsets(t) + linspace(-0.05, 0.05, numel(y));
            end
            scatter(x, y, 12, [0 0 0], 'filled', ...
                'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor', 'none')
            bar(a + offsets(t), mean(y), 0.2, 'FaceColor', lock_bar_color(sample.type_names{t}), ...
                'FaceAlpha', 0.45, 'EdgeColor', 'none')
            hold on
        end
    end
    if n_area > 0
        set(gca, 'XTick', 1:n_area, 'XTickLabel', sample.area_names)
        xlim([0.45, n_area + 0.55])
    end
    try
        set(gca, 'TickLabelInterpreter', 'none')
    catch
    end
    title(two_way_title(sample.feature, model), 'Interpreter', 'none')
    grid on
    box off
end
sgtitle(caption, 'Interpreter', 'none')
save_script_figure(folder)
end


function lines = two_way_title(feature_name, model)
if ~model.ok || isempty(model.terms)
    lines = {feature_name, 'Too few neurons'};
    return
end
lines = cell(1, numel(model.terms) + 1);
lines{1} = feature_name;
for k = 1:numel(model.terms)
    if ~isfinite(model.terms(k).holm_p)
        p_text = 'not tested';
    else
        p_text = sprintf('%.3g', model.terms(k).holm_p);
    end
    lines{k + 1} = sprintf('%s Holm-Sidak p = %s', model.terms(k).name, p_text);
end
end


function text = two_way_report(samples, models, intro)
lines = intro;
for j = 1:numel(samples)
    model = models{j};
    if ~model.ok
        lines{end + 1} = sprintf('%s: %s', samples{j}.feature, model.note); %#ok<AGROW>
        continue
    end
    significant = false(1, numel(model.terms));
    for k = 1:numel(model.terms)
        significant(k) = isfinite(model.terms(k).holm_p) && model.terms(k).holm_p < 0.05;
    end
    if ~any(significant)
        lines{end + 1} = sprintf('%s   n = %d   no Holm-Sidak significant term', ...
            samples{j}.feature, model.n); %#ok<AGROW>
        continue
    end
    lines{end + 1} = sprintf('%s   n = %d', samples{j}.feature, model.n); %#ok<AGROW>
    if ~isempty(model.note)
        lines{end + 1} = ['  ', model.note]; %#ok<AGROW>
    end
    for k = 1:numel(model.terms)
        if ~significant(k)
            continue
        end
        term = model.terms(k);
        lines{end + 1} = sprintf('  %-14s F(%s, %s) = %s   partial eta2 = %s   raw p = %s   Holm-Sidak p = %s', ...
            term.name, stat_num(term.df1, '%.0f'), stat_num(term.df2, '%.0f'), ...
            stat_num(term.F, '%.4g'), stat_num(term.eta, '%.4g'), ...
            stat_num(term.p, '%.4g'), stat_num(term.holm_p, '%.4g')); %#ok<AGROW>
    end
end
text = strjoin(lines, newline);
end


function [X, y] = lock_classifier_data(feature_table, neuron_type, feature_names, osc_rows)
%LOCK_CLASSIFIER_DATA  Selected features for trough, peak, or unlocked neurons.

n = size(feature_table, 1);
X = nan(n, numel(feature_names));
for j = 1:numel(feature_names)
    column = feature_column(feature_table, feature_names{j});
    n_copy = min(n, numel(column));
    X(1:n_copy, j) = column(1:n_copy);
end
class_names = {'Trough', 'Peak', 'Unlocked'};
y = zeros(n, 1);
for i = 1:min(n, numel(neuron_type))
    hit = find(strcmp(class_names, display_lock_label(neuron_type{i})), 1);
    if ~isempty(hit)
        y(i) = hit;
    end
end
osc_keep = false(n, 1);
osc_rows = osc_rows(:);
osc_rows = osc_rows(osc_rows >= 1 & osc_rows <= n);
osc_keep(round(osc_rows)) = true;
keep = osc_keep & y > 0 & all(isfinite(X), 2);
X = X(keep, :);
y = y(keep);
end


function methods = classifier_method_list(method)
allowed = {'lda', 'naivebayes', 'svm', 'tree', 'logit'};
if isstring(method)
    method = cellstr(method);
elseif ischar(method)
    method = {method};
end
method = method(:)';
if numel(method) == 1 && strcmpi(method{1}, 'all')
    methods = allowed;
    return
end
methods = cell(1, numel(method));
for i = 1:numel(method)
    hit = strcmpi(allowed, method{i});
    if ~any(hit)
        error('Plot_acg_population_logit:classifier', ...
            'classifier_method must be all, lda, naivebayes, svm, tree, logit, or a list of those names.');
    end
    methods{i} = allowed{hit};
end
end


function text = classifier_method_text(method)
if iscell(method)
    text = strjoin(method, ', ');
    return
end
text = char(method);
end


function result = classifier_crossval(X, y, class_names, method)
%CLASSIFIER_CROSSVAL  Stratified cross-validated accuracy for one method.

result.ok = false;
result.note = '';
result.method = method;
result.method_label = classifier_method_label(method);
result.accuracy = nan;
result.neuron_accuracy = nan;
result.chance = nan;
result.n_folds = 0;
result.n = numel(y);
result.class_names = class_names;
result.counts = zeros(numel(class_names));
result.recall = nan(1, numel(class_names));
result.class_n = zeros(1, numel(class_names));
if ~classifier_method_available(method)
    result.note = 'This option needs the Statistics and Machine Learning Toolbox.';
    return
end
if isempty(y) || size(X, 2) < 1
    result.note = 'No neuron had a lock type and finite values on every selected feature.';
    return
end
for c = 1:numel(class_names)
    result.class_n(c) = sum(y == c);
end
present = find(result.class_n > 0);
if numel(present) < 2
    result.note = 'Fewer than two lock types had neurons with every selected feature.';
    return
end
result.n_folds = min(10, min(result.class_n(present)));
if result.n_folds < 2
    result.note = 'A lock type has fewer than two neurons, so cross-validation was not run.';
    return
end
folds = stratified_folds(y, result.n_folds);
prediction = nan(size(y));
for fold = 1:result.n_folds
    test = folds == fold;
    try
        predicted = classifier_predict(method, X(~test, :), y(~test), X(test, :));
        prediction(test) = predicted(:);
    catch err
        result.note = err.message;
        return
    end
end
result.ok = true;
result.neuron_accuracy = mean(prediction == y);
for c = 1:numel(class_names)
    for p = 1:numel(class_names)
        result.counts(c, p) = sum(y == c & prediction == p);
    end
    if result.class_n(c) > 0
        result.recall(c) = result.counts(c, c) / result.class_n(c);
    end
end
present_recall = result.recall(result.class_n > 0);
result.accuracy = mean(present_recall);
result.chance = 1 / numel(present_recall);
end


function folds = stratified_folds(y, n_folds)
folds = zeros(size(y));
classes = unique(y(:))';
for c = classes
    idx = find(y == c);
    for i = 1:numel(idx)
        folds(idx(i)) = mod(i - 1, n_folds) + 1;
    end
end
end


function label = classifier_method_label(method)
switch lower(method)
    case 'lda'
        label = 'Linear discriminant';
    case 'naivebayes'
        label = 'Gaussian naive Bayes';
    case 'svm'
        label = 'Medium Gaussian SVM';
    case 'tree'
        label = 'Classification tree';
    case 'logit'
        label = 'Multinomial logit';
    otherwise
        label = method;
end
end


function available = classifier_method_available(method)
switch lower(method)
    case {'lda', 'naivebayes', 'logit'}
        available = true;
    case 'svm'
        available = exist('fitcecoc', 'file') == 2 && exist('templateSVM', 'file') == 2;
    case 'tree'
        available = exist('fitctree', 'file') == 2;
    otherwise
        available = false;
end
end


function prediction = classifier_predict(method, X_train, y_train, X_test, weighting)
%CLASSIFIER_PREDICT  One fold. Features are standardized on the training fold.

if nargin < 5 || isempty(weighting)
    weighting = 'weighted';
end
n_class = max([y_train(:); 1]);
[prediction, ~, ~] = fold_predictions(method, X_train, y_train, X_test, weighting, n_class);
end


function weights = inverse_frequency_weights(y)
%INVERSE_FREQUENCY_WEIGHTS  Equal total weight per lock type, mean weight 1.

y = y(:);
weights = ones(size(y));
classes = unique(y);
n_class = numel(classes);
n = numel(y);
if n_class < 1 || n < 1
    return
end
for c = 1:n_class
    members = y == classes(c);
    weights(members) = n / (n_class * sum(members));
end
end


function shuffled = shared_label_shuffles(y, n_perm)
%SHARED_LABEL_SHUFFLES  One label permutation list for every model.

y = y(:);
shuffled = zeros(n_perm, numel(y));
for shuffle = 1:n_perm
    shuffled(shuffle, :) = y(randperm(numel(y))).';
end
end


function [populations, loaded] = try_load_lock_model_cache(file, version, seed, n_perm, methods, class_names, feature_names, feature_table, neuron_type, populations)
%TRY_LOAD_LOCK_MODEL_CACHE  Use a saved fit only when the conditions match.

loaded = false;
if exist(file, 'file') ~= 2
    fprintf('No saved lock-type model at %s. The fit will run.\n', file)
    return
end
stored = load(file, 'lock_model_cache');
if ~isfield(stored, 'lock_model_cache')
    fprintf(['%s has no lock_model_cache variable. The fit will run.\n'], file)
    return
end
reason = lock_model_cache_mismatch(stored.lock_model_cache, version, seed, ...
    n_perm, methods, class_names, feature_names, feature_table, neuron_type, populations);
if ~isempty(reason)
    fprintf('Saved lock-type model was not used: %s The fit will run.\n', reason)
    return
end
populations = stored.lock_model_cache.populations;
loaded = true;
fprintf(['Loaded lock-type models from %s (code version %d, %d shuffles). ', ...
    'The fit was not repeated.\n'], file, version, n_perm)
end


function reason = lock_model_cache_mismatch(cache, version, seed, n_perm, methods, class_names, feature_names, feature_table, neuron_type, populations)
reason = '';
if ~isstruct(cache) || ~isfield(cache, 'version') || cache.version ~= version
    reason = 'the code version differs.';
    return
end
if ~isfield(cache, 'seed') || cache.seed ~= seed
    reason = 'the random seed differs.';
    return
end
if ~isfield(cache, 'n_permutations') || cache.n_permutations ~= n_perm
    reason = 'the number of shuffles differs.';
    return
end
if ~isfield(cache, 'methods') || ~isequal(cache.methods, methods)
    reason = 'the classifier list differs.';
    return
end
if ~isfield(cache, 'class_names') || ~isequal(cache.class_names, class_names)
    reason = 'the lock-type names differ.';
    return
end
if ~isfield(cache, 'feature_names') || ~isequal(cache.feature_names, feature_names)
    reason = 'the feature list differs.';
    return
end
if ~isfield(cache, 'populations') || numel(cache.populations) ~= numel(populations)
    reason = 'the neuron groups differ.';
    return
end
for p = 1:numel(populations)
    saved = cache.populations{p};
    if ~isstruct(saved) || ~isfield(saved, 'name') || ~strcmp(saved.name, populations{p}.name)
        reason = 'a neuron group name differs.';
        return
    end
    [X, y] = lock_classifier_data(feature_table, neuron_type, feature_names, populations{p}.rows);
    if ~isfield(saved, 'X') || ~isfield(saved, 'y') || ~isequal(saved.X, X) || ~isequal(saved.y, y)
        reason = 'the neurons or their feature values differ.';
        return
    end
    if ~isfield(saved, 'grid') || ~lock_model_grid_complete(saved.grid, n_perm)
        reason = 'the saved fit is incomplete.';
        return
    end
end
end


function complete = lock_model_grid_complete(grid, n_perm)
complete = iscell(grid) && isequal(size(grid), [3, 3]);
if ~complete
    return
end
for i = 1:numel(grid)
    item = grid{i};
    if ~isstruct(item) || ~isfield(item, 'ok')
        complete = false;
        return
    end
    if item.ok && (~isfield(item, 'n_perm') || item.n_perm ~= n_perm)
        complete = false;
        return
    end
end
end


function write_lock_model_cache(file, version, seed, n_perm, methods, class_names, feature_names, populations)
lock_model_cache.version = version;
lock_model_cache.seed = seed;
lock_model_cache.n_permutations = n_perm;
lock_model_cache.methods = methods;
lock_model_cache.class_names = class_names;
lock_model_cache.feature_names = feature_names;
lock_model_cache.populations = populations;
save(file, 'lock_model_cache')
fprintf('Saved lock-type models to %s (code version %d).\n', file, version)
end


function report_lock_model_population(population, class_names)
grid = population.grid;
for row = 1:size(grid, 1)
    for col = 1:size(grid, 2)
        result = grid{row, col};
        if result.ok
            fprintf(['%s, %s mean lock-type accuracy: %.1f%%. ', ...
                'Equal-type reference: %.1f%%. Neuron-weighted accuracy: %.1f%%. ', ...
                'Permutation p = %.3g.\n'], ...
                population.name, result.method_label, ...
                100 * result.accuracy, 100 * result.chance, ...
                100 * result.neuron_accuracy, result.perm_p)
            for c = 1:numel(class_names)
                fprintf('  %s accuracy: %.1f%%. Permutation p = %.3g.\n', ...
                    class_names{c}, 100 * result.recall(c), result.class_perm_p(c))
            end
        else
            fprintf('%s, %s: %s\n', population.name, result.method_label, result.note)
        end
    end
end
end


function results = classifier_weight_trio(X, y, class_names, method)
%CLASSIFIER_WEIGHT_TRIO  Weighted, unweighted, and joint decisions for one method.
%
% The joint decision mixes the two class probabilities. Each probability
% is scaled by that model's mean lock-type accuracy on the training fold.
% The test neuron's lock type is not an input.

row_names = {'Weighted', 'Unweighted', 'Joint decision'};
results = cell(3, 1);
for row = 1:3
    results{row} = empty_classifier_result(class_names, method, row_names{row});
end
if ~classifier_method_available(method)
    for row = 1:3
        results{row}.note = 'This option needs the Statistics and Machine Learning Toolbox.';
    end
    return
end
y = y(:);
n = numel(y);
n_class = numel(class_names);
if isempty(y) || size(X, 2) < 1
    for row = 1:3
        results{row}.note = 'No neuron had a lock type and finite values on every selected feature.';
    end
    return
end
for row = 1:3
    results{row}.n = n;
    for c = 1:n_class
        results{row}.class_n(c) = sum(y == c);
    end
end
present = find(results{1}.class_n > 0);
if numel(present) < 2
    for row = 1:3
        results{row}.note = 'Fewer than two lock types had neurons with every selected feature.';
    end
    return
end
n_folds = min(10, min(results{1}.class_n(present)));
if n_folds < 2
    for row = 1:3
        results{row}.note = 'A lock type has fewer than two neurons, so cross-validation was not run.';
    end
    return
end
folds = stratified_folds(y, n_folds);
prediction = nan(n, 3);
for fold = 1:n_folds
    test = folds == fold;
    train = ~test;
    try
        [pred_w, post_w, train_w] = fold_predictions(method, X(train, :), y(train), ...
            X(test, :), 'weighted', n_class);
        [pred_u, post_u, train_u] = fold_predictions(method, X(train, :), y(train), ...
            X(test, :), 'unweighted', n_class);
        acc_w = mean_lock_accuracy(y(train), train_w, n_class);
        acc_u = mean_lock_accuracy(y(train), train_u, n_class);
        prediction(test, 1) = pred_w;
        prediction(test, 2) = pred_u;
        prediction(test, 3) = mix_class_posteriors(post_w, post_u, acc_w, acc_u);
    catch err
        for row = 1:3
            results{row}.note = err.message;
        end
        return
    end
end
for row = 1:3
    results{row} = scored_prediction_result(results{row}, prediction(:, row), y, n_folds);
end
end


function prediction = mix_class_posteriors(post_w, post_u, acc_w, acc_u)
%MIX_CLASS_POSTERIORS  Scale each model's probabilities by its training accuracy.

if ~isfinite(acc_w) || acc_w < 0
    acc_w = 0;
end
if ~isfinite(acc_u) || acc_u < 0
    acc_u = 0;
end
if acc_w + acc_u <= 0
    acc_w = 1;
    acc_u = 1;
end
mixed = acc_w * post_w + acc_u * post_u;
[~, prediction] = max(mixed, [], 2);
prediction = prediction(:);
end


function accuracy = mean_lock_accuracy(y, prediction, n_class)
recalls = zeros(1, 0);
for c = 1:n_class
    members = y(:) == c;
    if any(members)
        recalls(end + 1) = mean(prediction(members) == c); %#ok<AGROW>
    end
end
if isempty(recalls)
    accuracy = nan;
    return
end
accuracy = mean(recalls);
end


function result = empty_classifier_result(class_names, method, row_name)
result.ok = false;
result.note = '';
result.method = method;
result.method_label = [row_name, ' ', classifier_method_label(method)];
result.accuracy = nan;
result.neuron_accuracy = nan;
result.chance = nan;
result.n_folds = 0;
result.n = 0;
result.class_names = class_names;
result.counts = zeros(numel(class_names));
result.recall = nan(1, numel(class_names));
result.class_n = zeros(1, numel(class_names));
result.n_perm = 0;
result.null_mean = nan;
result.null_p95 = nan;
result.perm_p = nan;
result.null_accuracy = zeros(0, 1);
result.null_recall = zeros(0, numel(class_names));
result.class_null_mean = nan(1, numel(class_names));
result.class_null_p95 = nan(1, numel(class_names));
result.class_perm_p = nan(1, numel(class_names));
result.perm_note = '';
end


function result = scored_prediction_result(result, prediction, y, n_folds)
prediction = prediction(:);
y = y(:);
if numel(prediction) ~= numel(y) || any(~isfinite(prediction))
    result.note = 'The classifier did not return a class for every neuron.';
    return
end
n_class = numel(result.class_names);
result.ok = true;
result.n_folds = n_folds;
result.neuron_accuracy = mean(prediction == y);
for c = 1:n_class
    for p = 1:n_class
        result.counts(c, p) = sum(y == c & prediction == p);
    end
    if result.class_n(c) > 0
        result.recall(c) = result.counts(c, c) / result.class_n(c);
    end
end
present_recall = result.recall(result.class_n > 0);
result.accuracy = mean(present_recall);
result.chance = 1 / numel(present_recall);
end


function result = apply_shared_null(result, null_accuracy, null_recall)
%APPLY_SHARED_NULL  Permutation p from one shared list of shuffled labels.

result.n_perm = 0;
result.null_mean = nan;
result.null_p95 = nan;
result.perm_p = nan;
result.null_accuracy = zeros(0, 1);
n_class = numel(result.class_names);
result.null_recall = zeros(0, n_class);
result.class_null_mean = nan(1, n_class);
result.class_null_p95 = nan(1, n_class);
result.class_perm_p = nan(1, n_class);
result.perm_note = '';
if ~result.ok
    result.perm_note = 'The permutation test was skipped because the classifier did not run.';
    return
end
kept = isfinite(null_accuracy);
null_accuracy = null_accuracy(kept);
null_recall = null_recall(kept, :);
result.n_perm = numel(null_accuracy);
result.null_accuracy = null_accuracy;
result.null_recall = null_recall;
if isempty(null_accuracy)
    result.perm_note = 'Every label shuffle failed, so the permutation p was not computed.';
    return
end
result.null_mean = mean(null_accuracy);
result.null_p95 = upper_percentile(null_accuracy, 0.95);
result.perm_p = permutation_p(null_accuracy, result.accuracy);
for c = 1:n_class
    samples = null_recall(isfinite(null_recall(:, c)), c);
    if isempty(samples) || ~isfinite(result.recall(c))
        continue
    end
    result.class_null_mean(c) = mean(samples);
    result.class_null_p95(c) = upper_percentile(samples, 0.95);
    result.class_perm_p(c) = permutation_p(samples, result.recall(c));
end
end


function [test_pred, test_post, train_pred] = fold_predictions(method, X_train, y_train, X_test, weighting, n_class)
%FOLD_PREDICTIONS  Fit one model and score the test rows and the training rows.

[X_train, X_test] = standardize_fold(X_train, X_test);
y_train = y_train(:);
if strcmpi(weighting, 'unweighted')
    weights = ones(numel(y_train), 1);
else
    weights = inverse_frequency_weights(y_train);
end
switch lower(method)
    case 'lda'
        [test_pred, test_post, train_pred] = lda_two_sets( ...
            X_train, y_train, X_test, weights, n_class);
    case 'naivebayes'
        [test_pred, test_post, train_pred] = naive_bayes_two_sets( ...
            X_train, y_train, X_test, weights, n_class);
    case 'svm'
        [test_pred, test_post, train_pred] = svm_two_sets( ...
            X_train, y_train, X_test, weights, n_class);
    otherwise
        error('Plot_acg_population_logit:classifier', ...
            'The joint comparison runs svm, naivebayes, and lda.');
end
end


function [X_train, X_test] = standardize_fold(X_train, X_test)
center = mean(X_train, 1);
scale = std(X_train, 0, 1);
scale(~isfinite(scale) | scale < eps(max(abs(center) + 1))) = 1;
X_train = (X_train - center) ./ scale;
X_test = (X_test - center) ./ scale;
end


function [test_pred, test_post, train_pred] = lda_two_sets(X_train, y_train, X_test, weights, n_class)
%LDA_TWO_SETS  One linear discriminant, scored on the test and training rows.

if exist('fitcdiscr', 'file') == 2
    weighted = max(weights) - min(weights) > 1e-10;
    try
        if weighted
            model = fitcdiscr(X_train, y_train, 'DiscrimType', 'linear', 'Weights', weights);
        else
            model = fitcdiscr(X_train, y_train, 'DiscrimType', 'linear');
        end
    catch
        if weighted
            model = fitcdiscr(X_train, y_train, 'DiscrimType', 'pseudoLinear', 'Prior', 'uniform');
        else
            model = fitcdiscr(X_train, y_train, 'DiscrimType', 'pseudoLinear');
        end
    end
    [test_pred, test_post] = posterior_from_model(model, X_test, n_class);
    [train_pred, ~] = posterior_from_model(model, X_train, n_class);
    return
end
[test_pred, test_post, train_pred] = lda_hand_rolled(X_train, y_train, X_test, weights, n_class);
end


function [test_pred, test_post, train_pred] = lda_hand_rolled(X_train, y_train, X_test, weights, n_class)
classes = unique(y_train(:))';
n_present = numel(classes);
n_feature = size(X_train, 2);
means = zeros(n_present, n_feature);
priors = zeros(n_present, 1);
covariance = zeros(n_feature, n_feature);
for c = 1:n_present
    rows = y_train == classes(c);
    priors(c) = sum(weights(rows)) / sum(weights);
    means(c, :) = mean(X_train(rows, :), 1);
    centered = X_train(rows, :) - means(c, :);
    class_df = max(sum(rows) - 1, 1);
    covariance = covariance + (centered' * centered) / class_df;
end
covariance = covariance / n_present;
ridge = 1e-6 * trace(covariance) / n_feature;
if ~isfinite(ridge) || ridge <= 0
    ridge = 1e-6;
end
covariance = covariance + ridge * eye(n_feature);
weighted_means = covariance \ means';
quadratic = 0.5 * sum(means' .* weighted_means, 1);
[test_pred, test_post] = prediction_from_scores( ...
    X_test * weighted_means - quadratic + log(max(priors, realmin)).', classes, n_class);
[train_pred, ~] = prediction_from_scores( ...
    X_train * weighted_means - quadratic + log(max(priors, realmin)).', classes, n_class);
end


function [test_pred, test_post, train_pred] = naive_bayes_two_sets(X_train, y_train, X_test, weights, n_class)
%NAIVE_BAYES_TWO_SETS  One Gaussian naive Bayes, scored on both row sets.

if exist('fitcnb', 'file') == 2
    weighted = max(weights) - min(weights) > 1e-10;
    try
        if weighted
            model = fitcnb(X_train, y_train, 'DistributionNames', 'normal', 'Weights', weights);
        else
            model = fitcnb(X_train, y_train, 'DistributionNames', 'normal');
        end
    catch
        if weighted
            model = fitcnb(X_train, y_train, 'DistributionNames', 'normal', 'Prior', 'uniform');
        else
            model = fitcnb(X_train, y_train, 'DistributionNames', 'normal');
        end
    end
    [test_pred, test_post] = posterior_from_model(model, X_test, n_class);
    [train_pred, ~] = posterior_from_model(model, X_train, n_class);
    return
end
classes = unique(y_train(:))';
n_present = numel(classes);
log_test = zeros(size(X_test, 1), n_present);
log_train = zeros(size(X_train, 1), n_present);
for c = 1:n_present
    rows = y_train == classes(c);
    mu = mean(X_train(rows, :), 1);
    sigma = var(X_train(rows, :), 0, 1);
    sigma(~isfinite(sigma) | sigma < 1e-6) = 1e-6;
    prior = sum(weights(rows)) / sum(weights);
    log_test(:, c) = log(max(prior, realmin)) ...
        - 0.5 * sum(((X_test - mu) ./ sqrt(sigma)) .^ 2, 2) ...
        - 0.5 * sum(log(2 * pi * sigma));
    log_train(:, c) = log(max(prior, realmin)) ...
        - 0.5 * sum(((X_train - mu) ./ sqrt(sigma)) .^ 2, 2) ...
        - 0.5 * sum(log(2 * pi * sigma));
end
[test_pred, test_post] = prediction_from_scores(log_test, classes, n_class);
[train_pred, ~] = prediction_from_scores(log_train, classes, n_class);
end


function [prediction, posterior] = naive_bayes_predict(X_train, y_train, X_test, weights, n_class)
%NAIVE_BAYES_PREDICT  Gaussian naive Bayes for one score matrix.

if nargin < 5 || isempty(n_class)
    n_class = max(y_train(:));
end
[prediction, posterior, ~] = naive_bayes_two_sets(X_train, y_train, X_test, weights, n_class);
end


function [test_pred, test_post, train_pred] = svm_two_sets(X_train, y_train, X_test, weights, n_class)
%SVM_TWO_SETS  One medium Gaussian SVM, scored on the test and training rows.

kernel_scale = sqrt(size(X_train, 2));
template = templateSVM('KernelFunction', 'gaussian', ...
    'KernelScale', kernel_scale, 'BoxConstraint', 1, 'Standardize', false);
model = fitcecoc(X_train, y_train, 'Learners', template, 'Weights', weights);
[test_pred, test_post] = svm_posterior(model, X_test, n_class);
[train_pred, ~] = svm_posterior(model, X_train, n_class);
end


function [prediction, posterior] = svm_posterior(model, X, n_class)
[prediction, neg_loss] = predict(model, X);
prediction = prediction(:);
posterior = align_posterior(softmax_rows(neg_loss), model.ClassNames(:), n_class);
end


function [prediction, posterior] = posterior_from_model(model, X, n_class)
[prediction, posterior] = predict(model, X);
prediction = prediction(:);
posterior = align_posterior(posterior, model.ClassNames(:), n_class);
end


function [prediction, posterior] = prediction_from_scores(scores, classes, n_class)
posterior = align_posterior(softmax_rows(scores), classes(:), n_class);
[~, prediction] = max(posterior, [], 2);
prediction = prediction(:);
end


function probability = softmax_rows(scores)
scores(~isfinite(scores)) = -inf;
scores = scores - max(scores, [], 2);
probability = exp(scores);
row_sum = sum(probability, 2);
row_sum(row_sum <= 0) = 1;
probability = probability ./ row_sum;
end


function aligned = align_posterior(posterior, classes, n_class)
aligned = zeros(size(posterior, 1), n_class);
classes = classes(:);
for c = 1:numel(classes)
    class_id = classes(c);
    if class_id >= 1 && class_id <= n_class
        aligned(:, class_id) = posterior(:, c);
    end
end
row_sum = sum(aligned, 2);
keep = row_sum > 0;
aligned(keep, :) = aligned(keep, :) ./ row_sum(keep);
end


function prediction = logit_predict(X_train, y_train, X_test, weights)
%LOGIT_PREDICT  Multinomial logistic regression.

if exist('mnrfit', 'file') == 2 && exist('mnrval', 'file') == 2 ...
        && max(weights) - min(weights) < 1e-10
    coefficients = mnrfit(X_train, y_train);
    probability = mnrval(coefficients, X_test);
    [~, pick] = max(probability, [], 2);
    classes = unique(y_train(:));
    prediction = classes(pick);
    prediction = prediction(:);
    return
end
prediction = logit_irls_predict(X_train, y_train, X_test, weights);
end


function prediction = logit_irls_predict(X_train, y_train, X_test, weights)
classes = unique(y_train(:));
n_class = numel(classes);
y_index = zeros(size(y_train));
for c = 1:n_class
    y_index(y_train == classes(c)) = c;
end
observation_weight = weights(:);
Z = [ones(size(X_train, 1), 1), X_train];
q = size(Z, 2);
n_coef = n_class - 1;
beta = zeros(q, n_coef);
penalty = 1e-4;
for iter = 1:40 %#ok<NASGU>
    eta = [Z * beta, zeros(size(X_train, 1), 1)];
    eta = eta - max(eta, [], 2);
    probability = exp(eta);
    probability = probability ./ sum(probability, 2);
    gradient = zeros(q * n_coef, 1);
    information = zeros(q * n_coef);
    for k = 1:n_coef
        gradient((k - 1) * q + (1:q)) = Z' * (observation_weight .* ...
            (double(y_index == k) - probability(:, k))) - penalty * beta(:, k);
        for m = 1:n_coef
            if k == m
                weight = probability(:, k) .* (1 - probability(:, k));
            else
                weight = -probability(:, k) .* probability(:, m);
            end
            information((k - 1) * q + (1:q), (m - 1) * q + (1:q)) = ...
                (Z .* (observation_weight .* weight))' * Z;
        end
    end
    information = information + penalty * eye(q * n_coef);
    step = information \ gradient;
    if ~all(isfinite(step))
        break
    end
    if max(abs(step)) > 5
        step = 5 * step / max(abs(step));
    end
    beta = beta + reshape(step, q, n_coef);
    if max(abs(step)) < 1e-8
        break
    end
end
eta = [([ones(size(X_test, 1), 1), X_test]) * beta, zeros(size(X_test, 1), 1)];
eta = eta - max(eta, [], 2);
probability = exp(eta);
probability = probability ./ sum(probability, 2);
[~, pick] = max(probability, [], 2);
prediction = classes(pick);
prediction = prediction(:);
end


function draw_lock_classifier(grid, folder, method, figure_name, population_name)
open_full_figure(figure_name)
[n_rows, n_cols] = size(grid);
chance = nan;
for row = 1:n_rows
    for col = 1:n_cols
        subplot(n_rows, n_cols, (row - 1) * n_cols + col)
        result = grid{row, col};
        if ~result.ok
            axis off
            title(result.method_label, 'Interpreter', 'none')
            if ~isempty(result.note)
                text(0.5, 0.5, result.note, 'HorizontalAlignment', 'center', ...
                    'Interpreter', 'none')
            end
            continue
        end
        imagesc(result.counts)
        set(gca, 'XTick', 1:numel(result.class_names), 'XTickLabel', result.class_names, ...
            'YTick', 1:numel(result.class_names), 'YTickLabel', result.class_names)
        try
            set(gca, 'TickLabelInterpreter', 'none')
        catch
        end
        colormap(gca, [linspace(0.97, 0.20, 64)', linspace(0.97, 0.45, 64)', linspace(0.97, 0.75, 64)'])
        for r = 1:size(result.counts, 1)
            for c = 1:size(result.counts, 2)
                text(c, r, sprintf('%d', result.counts(r, c)), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
            end
        end
        xlabel('Predicted lock type')
        ylabel('True lock type')
        title({result.method_label, sprintf('Mean type accuracy %.1f%%', 100 * result.accuracy)}, ...
            'Interpreter', 'none')
        box off
        if ~isfinite(chance)
            chance = result.chance;
        end
    end
end
if isfinite(chance)
    chance_line = sprintf('Equal-type reference %.1f%%. Rows are the true lock type.', ...
        100 * chance);
else
    chance_line = 'Rows are the true lock type.';
end
sgtitle({sprintf('%s. Columns: %s.', population_name, method), ...
    'Rows are weighted, unweighted, and the joint decision.', chance_line}, ...
    'Interpreter', 'none')
save_script_figure(folder)
end


function text = lock_classifier_report(feature_titles, feature_sources, method, n_perm, populations)
feature_lines = cell(1, numel(feature_titles));
for j = 1:numel(feature_titles)
    feature_lines{j} = sprintf('  %s, significant for %s.', feature_titles{j}, feature_sources{j});
end
lines = [{ ...
    'Classifier of lock type from the high-power features that differ by lock type.', ...
    'A feature is included when its one-way ANOVA Holm-Sidak p is below 0.05 for all classified neurons or for oscillatory neurons only.', ...
    'Holm-Sidak is computed separately inside each of those two comparisons.', ...
    'The selected features:'}, ...
    feature_lines, { ...
    'Naive Bayes, LDA, and the medium Gaussian SVM each have a weighted fit, an unweighted fit, and a joint decision.', ...
    'The weighted fit gives each neuron a weight equal to the inverse of its lock-type count, scaled so the mean weight is 1.', ...
    'The unweighted fit gives every neuron the same weight.', ...
    'The joint decision mixes the two class probabilities. Each probability is scaled by that model''s mean lock-type accuracy on the training fold.', ...
    'The mix uses the ACG features only. The lock type of the neuron being classified is not an input.', ...
    'Each classifier is run on every trough, peak, and unlocked neuron, and again on the oscillatory neurons only.', ...
    'Oscillatory neurons are cluster 1 from section 3.', ...
    'Non-classified cells are excluded.', ...
    'A neuron is used only when every selected feature is finite.', ...
    'Accuracy is stratified cross-validation. Each class is split across the folds in stored order.', ...
    'Within each training fold the features are standardized.', ...
    'The accuracy of the whole model is the unweighted mean of the lock-type accuracies.', ...
    'Neuron-weighted accuracy counts every neuron once and is listed beside that mean.', ...
    'The equal-type reference is one divided by the number of lock types. A constant guess at one type scores that value on the mean accuracy.', ...
    sprintf('Figure columns, left to right: %s.', method), ...
    'Figure rows, top to bottom: weighted, unweighted, joint decision.', ...
    'svm is a Gaussian-kernel SVM with medium scale sqrt(number of features) and box constraint 1.', ...
    'Three SVM classes are compared one against one. Its class probabilities are a softmax of the SVM scores.', ...
    'naivebayes is Gaussian naive Bayes. lda is a linear discriminant.', ...
    sprintf('The permutation test draws %d shuffled label lists once per population.', n_perm), ...
    'During the fit the command window reports each model and every 25th shuffle.', ...
    'A finished fit is saved as lock_type_model_cache.mat in the output folder.', ...
    'The same code version, features, neurons, labels, and permutation count load that file instead of fitting again.', ...
    'Every model in that population is tested on those same lists.', ...
    'The features stay fixed. The random stream is reset to 42 before the first population.', ...
    'The same shuffles give the p-value of each model and of each lock type.', ...
    'A lock-type accuracy is the fraction of neurons of that type that were classified as that type.', ...
    'p = (1 + number of shuffled accuracies at or above the real accuracy) / (1 + number of shuffles).', ...
    'A bar is filled when that permutation p is below 0.05. Otherwise only the border is drawn.'}];
for p = 1:numel(populations)
    lines{end + 1} = populations{p}.name; %#ok<AGROW>
    results = populations{p}.results;
    shared = [];
    for m = 1:numel(results)
        if results{m}.ok
            shared = results{m};
            break
        end
    end
    if ~isempty(shared)
        lines{end + 1} = sprintf('  Neurons: %d. Folds: %d.', shared.n, shared.n_folds); %#ok<AGROW>
        count_bits = cell(1, numel(shared.class_names));
        for c = 1:numel(shared.class_names)
            count_bits{c} = sprintf('%s %d', shared.class_names{c}, shared.class_n(c));
        end
        lines{end + 1} = ['  Class counts: ', strjoin(count_bits, ', ')]; %#ok<AGROW>
        lines{end + 1} = sprintf('  Equal-type reference: %.1f%%.', 100 * shared.chance); %#ok<AGROW>
    end
    for m = 1:numel(results)
        result = results{m};
        lines{end + 1} = ['  ', result.method_label]; %#ok<AGROW>
        if ~result.ok
            lines{end + 1} = ['    ', result.note]; %#ok<AGROW>
            continue
        end
        lines{end + 1} = sprintf('    Mean of the lock-type accuracies: %.1f%%.', 100 * result.accuracy); %#ok<AGROW>
        lines{end + 1} = sprintf('    Neuron-weighted accuracy: %.1f%%.', 100 * result.neuron_accuracy); %#ok<AGROW>
        if result.n_perm < 1 || ~isfinite(result.perm_p)
            lines{end + 1} = ['    ', result.perm_note]; %#ok<AGROW>
        else
            lines{end + 1} = sprintf('    Permutation shuffles: %d.', result.n_perm); %#ok<AGROW>
            lines{end + 1} = sprintf('    Whole model null mean: %.1f%%.', 100 * result.null_mean); %#ok<AGROW>
            lines{end + 1} = sprintf('    Whole model null 95th percentile: %.1f%%.', 100 * result.null_p95); %#ok<AGROW>
            lines{end + 1} = sprintf('    Whole model permutation p, on the mean of the lock-type accuracies: %.4g.', result.perm_p); %#ok<AGROW>
        end
        for c = 1:numel(result.class_names)
            lines{end + 1} = sprintf('    %-16s accuracy %s   null mean %s   null 95th %s   permutation p %s', ...
                result.class_names{c}, ...
                stat_num(100 * result.recall(c), '%.1f%%'), ...
                stat_num(100 * result.class_null_mean(c), '%.1f%%'), ...
                stat_num(100 * result.class_null_p95(c), '%.1f%%'), ...
                stat_num(result.class_perm_p(c), '%.4g')); %#ok<AGROW>
        end
        lines{end + 1} = '    Confusion counts, rows true and columns predicted:'; %#ok<AGROW>
        lines{end + 1} = labeled_number_matrix('True', result.class_names, ...
            result.class_names, result.counts, '%.0f'); %#ok<AGROW>
    end
end
text = strjoin(lines, newline);
end


function result = attach_accuracy_permutation(result, X, y, class_names, n_perm)
%ATTACH_ACCURACY_PERMUTATION  Label-shuffle null for the model and each lock type.

n_class = numel(result.class_names);
result.n_perm = 0;
result.null_mean = nan;
result.null_p95 = nan;
result.perm_p = nan;
result.null_accuracy = zeros(0, 1);
result.null_recall = zeros(0, n_class);
result.class_null_mean = nan(1, n_class);
result.class_null_p95 = nan(1, n_class);
result.class_perm_p = nan(1, n_class);
result.perm_note = '';
if ~result.ok
    result.perm_note = 'The permutation test was skipped because the classifier did not run.';
    return
end
if n_perm < 1
    result.perm_note = 'n_accuracy_permutations is below 1, so no shuffle was run.';
    return
end
null_accuracy = nan(n_perm, 1);
null_recall = nan(n_perm, n_class);
for shuffle = 1:n_perm
    shuffled = classifier_crossval(X, y(randperm(numel(y))), class_names, result.method);
    if shuffled.ok
        null_accuracy(shuffle) = shuffled.accuracy;
        null_recall(shuffle, :) = shuffled.recall(:).';
    end
end
kept = isfinite(null_accuracy);
null_accuracy = null_accuracy(kept);
null_recall = null_recall(kept, :);
result.n_perm = numel(null_accuracy);
result.null_accuracy = null_accuracy;
result.null_recall = null_recall;
if isempty(null_accuracy)
    result.perm_note = 'Every label shuffle failed, so the permutation p was not computed.';
    return
end
result.null_mean = mean(null_accuracy);
result.null_p95 = upper_percentile(null_accuracy, 0.95);
result.perm_p = permutation_p(null_accuracy, result.accuracy);
for c = 1:n_class
    samples = null_recall(isfinite(null_recall(:, c)), c);
    if isempty(samples) || ~isfinite(result.recall(c))
        continue
    end
    result.class_null_mean(c) = mean(samples);
    result.class_null_p95(c) = upper_percentile(samples, 0.95);
    result.class_perm_p(c) = permutation_p(samples, result.recall(c));
end
end


function p = permutation_p(null_values, observed)
p = (1 + sum(null_values >= observed - 1e-12)) / (1 + numel(null_values));
end


function value = upper_percentile(x, fraction)
x = sort(x(:));
if isempty(x)
    value = nan;
    return
end
index = max(1, ceil(fraction * numel(x)));
value = x(min(index, numel(x)));
end


function draw_accuracy_null(populations, folder, method)
for p = 1:numel(populations)
    open_full_figure(['Lock-type accuracy permutation, ', populations{p}.name])
    grid = populations{p}.grid;
    [n_rows, n_cols] = size(grid);
    for row = 1:n_rows
        for col = 1:n_cols
            subplot(n_rows, n_cols, (row - 1) * n_cols + col)
            result = grid{row, col};
            samples = result.null_accuracy;
            observed = result.accuracy;
            if ~result.ok || isempty(samples) || ~isfinite(observed)
                axis off
                title(result.method_label, 'Interpreter', 'none')
                if ~isempty(result.perm_note)
                    text(0.5, 0.5, result.perm_note, 'HorizontalAlignment', 'center', ...
                        'Interpreter', 'none')
                end
                continue
            end
            histogram(100 * samples, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
            hold on
            y_limit = ylim;
            plot(100 * [observed observed], y_limit, 'Color', [0.8 0.15 0.15], 'LineWidth', 1.6)
            xlabel('Mean of lock-type accuracies (%)')
            ylabel('Shuffles')
            title({sprintf('p = %.3g', result.perm_p), result.method_label}, ...
                'Interpreter', 'none')
            box off
        end
    end
    sgtitle({populations{p}.name, ...
        'Columns are SVM, naive Bayes, and LDA. Rows are weighted, unweighted, and the joint decision.', ...
        'All nine models use the same shuffled labels.', ...
        sprintf('Gray bars are label shuffles. The red line is the real mean accuracy. %s.', method)}, ...
        'Interpreter', 'none')
    save_script_figure(folder)
end
end


function draw_class_accuracy_bars(populations, folder)
for p = 1:numel(populations)
    open_full_figure(['Lock-type accuracy by cell type, ', populations{p}.name])
    grid = populations{p}.grid;
    [n_rows, n_cols] = size(grid);
    for row = 1:n_rows
        for col = 1:n_cols
            subplot(n_rows, n_cols, (row - 1) * n_cols + col)
            result = grid{row, col};
            if ~result.ok
                axis off
                title(result.method_label, 'Interpreter', 'none')
                if ~isempty(result.note)
                    text(0.5, 0.5, result.note, 'HorizontalAlignment', 'center', ...
                        'Interpreter', 'none')
                end
                continue
            end
            hold on
            y_top = 0;
            for c = 1:numel(result.class_names)
                if ~isfinite(result.recall(c))
                    continue
                end
                height = 100 * result.recall(c);
                color = lock_bar_color(result.class_names{c});
                significant = isfinite(result.class_perm_p(c)) && result.class_perm_p(c) < 0.05;
                if significant
                    bar(c, height, 0.72, 'FaceColor', color, 'EdgeColor', color, ...
                        'FaceAlpha', 0.9, 'LineWidth', 1.2)
                else
                    bar(c, height, 0.72, 'FaceColor', 'none', 'EdgeColor', color, ...
                        'LineWidth', 1.6)
                end
                if isfinite(result.class_perm_p(c))
                    text(c, height, sprintf('p = %.3g', result.class_perm_p(c)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                        'Interpreter', 'none')
                end
                y_top = max(y_top, height);
            end
            set(gca, 'XTick', 1:numel(result.class_names), 'XTickLabel', result.class_names)
            try
                set(gca, 'TickLabelInterpreter', 'none')
            catch
            end
            ylim([0, max(110, y_top + 15)])
            xlim([0.4, numel(result.class_names) + 0.6])
            ylabel('Accuracy (%)')
            title({result.method_label, sprintf('Mean p = %.3g', result.perm_p)}, ...
                'Interpreter', 'none')
            box off
        end
    end
    sgtitle({populations{p}.name, ...
        'Columns are SVM, naive Bayes, and LDA. Rows are weighted, unweighted, and the joint decision.', ...
        'Filled bars have permutation p < 0.05. Open bars do not.'}, ...
        'Interpreter', 'none')
    save_script_figure(folder)
end
end


function labels = lock_labels_for_rows(neuron_type, row_index)
labels = cell(numel(row_index), 1);
for i = 1:numel(row_index)
    row = row_index(i);
    if row >= 1 && row <= numel(neuron_type)
        labels{i} = display_lock_label(neuron_type{row});
    else
        labels{i} = 'Non-classified';
    end
end
end


function draw_lock_type_embedding(xy, lock_labels, folder, figure_name, x_label, y_label, caption)
names = {'Trough', 'Peak', 'Unlocked'};
open_full_figure(figure_name)
embedding_axes = gobjects(1, numel(names));
for g = 1:numel(names)
    embedding_axes(g) = subplot(1, numel(names), g);
    group_xy = xy(strcmp(lock_labels, names{g}), :);
    draw_iso_probability(xy, group_xy, lock_bar_color(names{g}), zeros(0, 2), [0 0 0])
    title(names{g}, 'Interpreter', 'none')
    xlabel(x_label, 'Interpreter', 'none')
    ylabel(y_label, 'Interpreter', 'none')
    grid on
    box off
end
limits = embedding_axis_limits(xy);
if ~isempty(limits)
    for g = 1:numel(names)
        set(embedding_axes(g), 'XLim', limits(1:2), 'YLim', limits(3:4))
    end
    linkaxes(embedding_axes, 'xy')
end
sgtitle(caption, 'Interpreter', 'none')
save_script_figure(folder)
end


function [groups, names] = values_by_lock(feature_table, neuron_type, feature_name, min_n, row_keep)
%VALUES_BY_LOCK  Feature values for trough, peak, and unlocked neurons.

values = feature_column(feature_table, feature_name);
if nargin < 5
    row_keep = [];
end
if nargin >= 5
    row_keep = row_keep(:);
    n = min(numel(values), numel(neuron_type));
    if isempty(row_keep)
        values = values([]);
        neuron_type = neuron_type([]);
    elseif islogical(row_keep)
        mask = false(n, 1);
        use = row_keep(1:min(numel(row_keep), n));
        mask(1:numel(use)) = use(:);
        values = values(mask);
        neuron_type = neuron_type(mask);
    else
        row_keep = unique(round(row_keep(row_keep >= 1 & row_keep <= n)), 'stable');
        values = values(row_keep);
        neuron_type = neuron_type(row_keep);
    end
end
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


function anova = lock_type_feature_anova(feature_table, neuron_type, feature_names, feature_titles, row_keep)
%LOCK_TYPE_FEATURE_ANOVA  One-way lock-type ANOVA for every high-power feature.

use_rows = nargin >= 5;
n_features = numel(feature_names);
anova = repmat(empty_area_anova(), n_features, 1);
for j = 1:n_features
    if use_rows
        [groups, names] = values_by_lock(feature_table, neuron_type, ...
            feature_names{j}, 2, row_keep);
    else
        [groups, names] = values_by_lock(feature_table, neuron_type, ...
            feature_names{j}, 2);
    end
    anova(j).feature = feature_titles{j};
    anova(j).areas = names;
    anova(j).groups = groups;
    stats = one_way_anova(groups);
    anova(j).n = stats.n;
    anova(j).mean = stats.mean;
    anova(j).F = stats.F;
    anova(j).df1 = stats.df1;
    anova(j).df2 = stats.df2;
    anova(j).raw_p = stats.p;
    anova(j).eta = stats.eta;
end
holm = holm_sidak_adjust([anova.raw_p]);
for j = 1:n_features
    anova(j).holm_p = holm(j);
    names = anova(j).areas;
    groups = anova(j).groups;
    if ~isfinite(anova(j).holm_p) || anova(j).holm_p >= 0.05 || numel(names) < 2
        continue
    end
    pair_list = empty_feature_stats();
    pair_list(1) = [];
    pair_index = zeros(0, 2);
    for a = 1:numel(names)
        for b = a + 1:numel(names)
            pair = compare_two_samples(groups{a}, groups{b}, 'unpaired');
            pair.feature = sprintf('%s vs %s', names{a}, names{b});
            pair.role = anova(j).feature;
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
    anova(j).pairs = pair_list;
    anova(j).pair_index = pair_index;
end
end


function draw_lock_feature_figure(anova, folder, n_cols, figure_name, title_lines)
open_full_figure(figure_name)
n_features = numel(anova);
tiledlayout(ceil(n_features / n_cols), n_cols, 'TileSpacing', 'compact')
for j = 1:n_features
    nexttile
    hold on
    names = anova(j).areas;
    groups = anova(j).groups;
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
        if ~isfinite(anova(j).mean(a))
            continue
        end
        bar(a, anova(j).mean(a), 0.55, 'FaceColor', lock_bar_color(names{a}), ...
            'FaceAlpha', 0.45, 'EdgeColor', 'none')
        hold on
    end
    y_top = nan;
    if ~isempty(anova(j).pairs)
        show = find([anova(j).pairs.holm_p] < 0.05);
        if ~isempty(show)
            pooled = [];
            for a = 1:numel(groups)
                pooled = [pooled; groups{a}(:)]; %#ok<AGROW>
            end
            y_top = draw_sig_brackets(anova(j).pair_index(show, :), ...
                [anova(j).pairs(show).holm_p], pooled);
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
    if ~isfinite(anova(j).holm_p)
        p_text = 'Too few neurons';
    else
        p_text = sprintf('Holm-Sidak p = %.3g', anova(j).holm_p);
    end
    title({anova(j).feature, p_text}, 'Interpreter', 'none')
    grid on
    box off
end
sgtitle(title_lines)
save_script_figure(folder)
end


function [names, titles, sources] = significant_lock_features(anova_all, anova_osc, feature_names, feature_titles)
%SIGNIFICANT_LOCK_FEATURES  Features with a lock-type difference in either set.

names = {};
titles = {};
sources = {};
for j = 1:numel(feature_names)
    in_all = isfinite(anova_all(j).holm_p) && anova_all(j).holm_p < 0.05;
    in_osc = isfinite(anova_osc(j).holm_p) && anova_osc(j).holm_p < 0.05;
    if ~(in_all || in_osc)
        continue
    end
    names{end + 1} = feature_names{j}; %#ok<AGROW>
    titles{end + 1} = feature_titles{j}; %#ok<AGROW>
    if in_all && in_osc
        sources{end + 1} = 'all classified neurons and oscillatory neurons'; %#ok<AGROW>
    elseif in_all
        sources{end + 1} = 'all classified neurons'; %#ok<AGROW>
    else
        sources{end + 1} = 'oscillatory neurons'; %#ok<AGROW>
    end
end
end


function color = lock_bar_color(name)
switch name
    case 'Trough'
        color = [0.20 0.45 0.80];
    case 'Peak'
        color = [0.85 0.25 0.20];
    case 'Unlocked'
        color = [0.18 0.62 0.18];
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


function path = save_oscillatory_label_table(folder, session, cluster_id, osc_member, complete_rows)
%SAVE_OSCILLATORY_LABEL_TABLE  Session and neuron id with the oscillatory label.

n = numel(osc_member);
if numel(session) ~= n || numel(cluster_id) ~= n
    error('Plot_acg_population_logit:cluster', ...
        'The session and neuron id must have one entry per population neuron.');
end
if exist(folder, 'dir') ~= 7
    mkdir(folder);
end
session_text = cell(n, 1);
for i = 1:n
    session_text{i} = value_text(session{i});
end
classification = repmat({'not classified'}, n, 1);
classified = false(n, 1);
classified(complete_rows) = true;
classification(classified & ~osc_member(:)) = {'non-oscillatory'};
classification(osc_member(:)) = {'oscillatory'};
labels = table((1:n)', session_text, double(cluster_id(:)), osc_member(:), classification, ...
    'VariableNames', {'population_row', 'session', 'neuron_id', 'oscillatory', 'classification'});
path = fullfile(folder, 'oscillatory_neuron_labels.csv');
writetable(labels, path);
fprintf('Saved %s\n', path)
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


function [power_high, power_low, power_all] = stratum_mean_delta(records, record_index)
%STRATUM_MEAN_DELTA  Mean channel delta power of the high, low, and all-run sets.
%
% score_high is the mean of the squared median delta-band envelopes of
% the runs at or above that neuron's 75th percentile. score_low is the
% same mean for the runs at or below the 25th percentile. score_all is
% that mean over every run with a finite score. Older checkpoints have
% no score_all, and those entries stay NaN.

n = numel(record_index);
power_high = nan(n, 1);
power_low = nan(n, 1);
power_all = nan(n, 1);
for i = 1:n
    record = records{record_index(i)};
    if ~isstruct(record) || ~isfield(record, 'power') || ~isstruct(record.power)
        continue
    end
    power_high(i) = scalar_power_field(record.power, 'score_high');
    power_low(i) = scalar_power_field(record.power, 'score_low');
    power_all(i) = scalar_power_field(record.power, 'score_all');
end
end


function value = scalar_power_field(power, field_name)
value = nan;
if ~isfield(power, field_name)
    return
end
raw = power.(field_name);
if isnumeric(raw) && isscalar(raw)
    value = double(raw);
end
end


function m = finite_mean(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
    return
end
m = mean(x);
end


function scale = choose_gap_scale(raw_gap, log_gap)
%CHOOSE_GAP_SCALE  Use the high-minus-low scale with smaller absolute skewness.

scale.skew_raw = moment_skewness(raw_gap);
scale.skew_log = moment_skewness(log_gap);
if isfinite(scale.skew_log) && abs(scale.skew_log) < abs(scale.skew_raw)
    scale.name = 'log';
    scale.reason = sprintf(['Absolute skewness is %.3f for the raw gap and %.3f for the log gap. ', ...
        'The log gap is less skewed, so the model uses log(high) minus log(low).'], ...
        abs(scale.skew_raw), abs(scale.skew_log));
    return
end
scale.name = 'raw';
scale.reason = sprintf(['Absolute skewness is %.3f for the raw gap and %.3f for the log gap. ', ...
    'The raw gap is less skewed, so the model uses high minus low.'], ...
    abs(scale.skew_raw), abs(scale.skew_log));
end


function g = moment_skewness(x)
x = x(:);
x = x(isfinite(x));
if numel(x) < 3
    g = nan;
    return
end
sigma = std(x, 0);
if ~isfinite(sigma) || sigma == 0
    g = 0;
    return
end
g = mean(((x - mean(x)) / sigma) .^ 3);
end


function fit = fit_power_feature_model(y, X, names, cluster_id, weights)
%FIT_POWER_FEATURE_MODEL  Weighted least squares with an intercept column.
%
% cluster_id groups rows that belong to one neuron. When it is present,
% the standard errors are cluster-robust, using the same G/(G-1) factor
% and G - p degrees of freedom as the locking-period oscillation fits.
% weights give each neuron total weight 1 when a neuron has two strata.

y = y(:);
n = size(X, 1);
if nargin < 5 || isempty(weights)
    weights = ones(n, 1);
end
weights = weights(:);
if numel(weights) ~= n || any(~isfinite(weights)) || any(weights < 0) || sum(weights) <= 0
    error('Plot_acg_population_logit:powerModel', ...
        'The weight vector must have one non-negative entry per row.');
end
if nargin < 4
    cluster_id = [];
end
A = [ones(n, 1), X];
fit.names = [{'Intercept'}, names];
fit.n = n;
fit.p = size(A, 2);
fit.df_resid = n - fit.p;
fit.df1 = fit.p - 1;
fit.clustered = ~isempty(cluster_id);
fit.n_clusters = nan;
fit.df = fit.df_resid;
sw = sqrt(weights);
Aw = A .* sw;
yw = y .* sw;
fit.beta = Aw \ yw;
fit.beta = fit.beta(:);
resid = y - A * fit.beta;
ss_res = sum(weights .* resid .^ 2);
y_bar = sum(weights .* y) / sum(weights);
ss_tot = sum(weights .* (y - y_bar) .^ 2);
if ss_tot > 0
    fit.rsq = 1 - ss_res / ss_tot;
else
    fit.rsq = nan;
end
if fit.df_resid > 0 && n > 1 && ss_tot > 0
    fit.adj_rsq = 1 - (ss_res / fit.df_resid) / (ss_tot / (n - 1));
    fit.sigma = sqrt(ss_res / fit.df_resid);
else
    fit.adj_rsq = nan;
    fit.sigma = nan;
end
fit.se = nan(fit.p, 1);
fit.t = nan(fit.p, 1);
fit.p_value = nan(fit.p, 1);
fit.F = nan;
fit.F_p = nan;
xtx = Aw' * Aw;
fit.collinear = ~isfinite(rcond(xtx)) || rcond(xtx) < 1e-8;
fit.ss_model = ss_tot - ss_res;
fit.resid_skew = moment_skewness(resid);
if fit.collinear
    return
end
if ~fit.clustered
    if fit.df > 0 && fit.sigma > 0
        covariance = (fit.sigma ^ 2) * inv(xtx);
        fit.se = sqrt(max(diag(covariance), 0));
        fit.t = fit.beta ./ fit.se;
        for i = 1:fit.p
            fit.p_value(i) = two_sided_t_p(fit.t(i), fit.df);
        end
    end
    if fit.df > 0 && fit.df1 > 0 && ss_res > 0
        fit.F = (fit.ss_model / fit.df1) / (ss_res / fit.df);
        fit.F_p = f_survival(fit.F, fit.df1, fit.df);
    end
    return
end
cluster_id = cluster_id(:);
if numel(cluster_id) ~= n
    error('Plot_acg_population_logit:powerModel', ...
        'The cluster id must have one entry per row.');
end
ids = unique(cluster_id);
G = numel(ids);
fit.n_clusters = G;
fit.df = G - fit.p;
bread = inv(xtx);
meat = zeros(fit.p, fit.p);
for i = 1:G
    in = cluster_id == ids(i);
    score = A(in, :)' * (weights(in) .* resid(in));
    meat = meat + score * score';
end
if G > 1
    meat = meat * (G / (G - 1));
end
V = bread * meat * bread;
V = (V + V') / 2;
fit.se = sqrt(max(real(diag(V)), 0));
fit.t = fit.beta ./ fit.se;
fit.t(~isfinite(fit.se) | fit.se <= 0) = nan;
if fit.df >= 1
    for i = 1:fit.p
        fit.p_value(i) = two_sided_t_p(fit.t(i), fit.df);
    end
end
if fit.df >= 1 && fit.df1 > 0
    contrast_rows = fit.df1;
    R = [zeros(contrast_rows, 1), eye(contrast_rows)];
    contrast = R * fit.beta;
    middle = R * V * R';
    if all(isfinite(contrast))
        wald = real(contrast' * pinv(middle) * contrast);
        fit.F = wald / fit.df1;
        fit.F_p = f_survival(fit.F, fit.df1, fit.df);
    end
end
end


function result = fit_within_neuron_gap(gap_raw, gap_log, gap_x, names, sample_name)
%FIT_WITHIN_NEURON_GAP  Power gap on predictor gaps, plus their correlation.

gap_raw = gap_raw(:);
gap_log = gap_log(:);
names = names(:)';
if size(gap_x, 1) ~= numel(gap_raw) || numel(gap_log) ~= numel(gap_raw)
    error('Plot_acg_population_logit:powerModel', ...
        'The power-gap rows and the predictor rows do not match.');
end
if size(gap_x, 1) < 3
    error('Plot_acg_population_logit:powerModel', ...
        'Too few neurons for the gap model (%s).', sample_name);
end
span = max(abs(gap_x(:)));
if isempty(span) || ~isfinite(span)
    span = 1;
end
gap_sd = std(gap_x, 0, 1);
keep = isfinite(gap_sd) & gap_sd > eps(span);
if ~any(keep)
    error('Plot_acg_population_logit:powerModel', ...
        'No predictor changes between the tails for %s.', sample_name);
end
dropped = names(~keep);
names = names(keep);
gap_x = gap_x(:, keep);
if isempty(dropped)
    result.dropped_line = 'Every listed predictor has some high-minus-low change.';
else
    result.dropped_line = ['Left out because the high-minus-low change has no variance: ', ...
        strjoin(dropped, ', '), '.'];
end
if size(gap_x, 1) < numel(names) + 2
    error('Plot_acg_population_logit:powerModel', ...
        'Too few neurons for the predictors that remain (%s).', sample_name);
end
result.scale = choose_gap_scale(gap_raw, gap_log);
if strcmp(result.scale.name, 'log')
    result.outcome = gap_log;
    result.outcome_label = 'log high minus log low delta power';
else
    result.outcome = gap_raw;
    result.outcome_label = 'high minus low delta power';
end
result.fit = fit_power_feature_model(result.outcome, gap_x, names);
result.predictors = gap_x;
result.names = names;
result.n = size(gap_x, 1);
result.correlation = corr(gap_x);
result.correlation_lines = correlation_table_lines(result.correlation, names);
result.sample_name = sample_name;
end


function print_gap_model(result)
fprintf('%s. n = %d.\n', result.sample_name, result.n)
fprintf('%s\n', result.dropped_line)
fprintf('Skewness of the raw gap: %.3f. Skewness of the log gap: %.3f.\n', ...
    result.scale.skew_raw, result.scale.skew_log)
fprintf('%s\n', result.scale.reason)
fprintf('Model outcome: %s. R^2 = %.3f. Adjusted R^2 = %.3f. F(%g,%g) = %.3f, p = %.3g.\n', ...
    result.outcome_label, result.fit.rsq, result.fit.adj_rsq, ...
    result.fit.df1, result.fit.df, result.fit.F, result.fit.F_p)
for i = 1:numel(result.fit.names)
    fprintf('  %-18s  coefficient %.4g  SE %.4g  t %.3f  p %.3g\n', ...
        result.fit.names{i}, result.fit.beta(i), result.fit.se(i), ...
        result.fit.t(i), result.fit.p_value(i))
end
fprintf('Predictor correlations:\n')
for i = 1:numel(result.correlation_lines)
    fprintf('%s\n', result.correlation_lines{i})
end
end


function lines = gap_fit_lines(result)
coeff = cell(numel(result.fit.names), 1);
for i = 1:numel(result.fit.names)
    coeff{i} = sprintf('%-18s %12.4g %12.4g %10.3f %10.3g', ...
        result.fit.names{i}, result.fit.beta(i), result.fit.se(i), ...
        result.fit.t(i), result.fit.p_value(i));
end
lines = [{ ...
    result.dropped_line, ...
    sprintf('Skewness of the raw gap: %.3f.', result.scale.skew_raw), ...
    sprintf('Skewness of the log gap: %.3f.', result.scale.skew_log), ...
    result.scale.reason, ...
    sprintf('Fitted outcome: %s.', result.outcome_label), ...
    sprintf('R^2 = %.3f. Adjusted R^2 = %.3f. F(%g, %g) = %.3f. p = %.3g.', ...
        result.fit.rsq, result.fit.adj_rsq, result.fit.df1, result.fit.df, ...
        result.fit.F, result.fit.F_p), ...
    sprintf('%-18s %12s %12s %10s %10s', 'Term', 'Coefficient', 'SE', 't', 'p')}, ...
    coeff(:)', { ...
    'Pearson correlation of the predictors:', ...
    }, result.correlation_lines(:)'];
if result.fit.collinear
    lines{end + 1} = 'The predictor columns are collinear, so the slope standard errors are not reported.';
end
end


function lines = correlation_table_lines(correlation, names)
names = names(:)';
n = numel(names);
width = max(cellfun(@numel, names));
width = max(width, 6);
header = blanks(width);
for j = 1:n
    header = [header, ' ', sprintf('%*s', width, names{j})]; %#ok<AGROW>
end
lines = cell(n + 2, 1);
lines{1} = header;
for i = 1:n
    row = sprintf('%*s', width, names{i});
    for j = 1:n
        row = [row, ' ', sprintf('%*.3f', width, correlation(i, j))]; %#ok<AGROW>
    end
    lines{i + 1} = row;
end
off = abs(correlation);
off(1:n + 1:end) = -inf;
[peak, idx] = max(off(:));
if isempty(idx) || ~isfinite(peak)
    lines{n + 2} = 'No off-diagonal correlation could be computed.';
    return
end
[row_i, col_i] = ind2sub([n, n], idx);
lines{n + 2} = sprintf( ...
    'Largest absolute correlation off the diagonal: %s and %s, r = %.3f.', ...
    names{row_i}, names{col_i}, correlation(row_i, col_i));
end


function draw_power_gap_scale(raw_gap, log_gap, scale, folder, figure_name)
if nargin < 5 || isempty(figure_name)
    figure_name = 'Delta power gap scale';
end
open_full_figure(figure_name)
subplot(1, 2, 1)
histogram(raw_gap, 'FaceColor', [0.25 0.35 0.55], 'EdgeColor', 'none')
xlabel('High minus low delta power')
ylabel('Neurons')
title(sprintf('Raw gap, skewness %.3f', scale.skew_raw), 'Interpreter', 'none')
box off
subplot(1, 2, 2)
histogram(log_gap, 'FaceColor', [0.75 0.45 0.25], 'EdgeColor', 'none')
xlabel('log high minus log low delta power')
ylabel('Neurons')
title(sprintf('Log gap, skewness %.3f', scale.skew_log), 'Interpreter', 'none')
box off
sgtitle({scale.reason, sprintf('The model uses the %s gap.', scale.name)}, ...
    'Interpreter', 'none')
save_script_figure(folder)
end


function draw_regressor_correlation(correlation, names, folder, figure_name)
open_full_figure(figure_name)
heatmap(names, names, correlation, 'ColorLimits', [-1 1])
title('Pearson correlation of the high-minus-low predictors')
save_script_figure(folder)
end


function draw_power_model_coefficients(fit, outcome_label, folder, figure_name)
if nargin < 4 || isempty(figure_name)
    figure_name = 'Within-neuron delta power gap';
end
open_full_figure(figure_name)
shown = 2:numel(fit.names);
hold on
for i = 1:numel(shown)
    k = shown(i);
    if strcmp(fit.names{k}, 'FiringRate')
        color = [0.75 0.42 0.18];
    else
        color = [0.25 0.35 0.55];
    end
    plot([i i], fit.beta(k) + fit.se(k) * [-1.96 1.96], 'Color', color, 'LineWidth', 1.4)
    plot(i, fit.beta(k), 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color)
end
plot([0.5, numel(shown) + 0.5], [0 0], 'Color', [0.5 0.5 0.5])
set(gca, 'XTick', 1:numel(shown), 'XTickLabel', fit.names(shown))
try
    set(gca, 'TickLabelInterpreter', 'none')
    xtickangle(30)
catch
end
ylabel(['Coefficient on ', outcome_label])
box off
sgtitle({sprintf('%s. R^2 = %.3f. F p = %.3g.', outcome_label, fit.rsq, fit.F_p), ...
    'Navy is a high-minus-low autocorrelogram feature.', ...
    'Orange is FiringRate, the high-minus-low firing rate in spikes/s.', ...
    'The intercept is the power gap when every predictor change is zero. It is not drawn.', ...
    'Whiskers are the coefficient plus or minus 1.96 standard errors.'}, ...
    'Interpreter', 'none')
save_script_figure(folder)
end


function path = stratum_summary_path(checkpoint, method, count)
if ~isempty(checkpoint)
    path = fullfile(fileparts(checkpoint), ...
        locking_summary_filename(false, method, count));
    return
end
path = fullfile(default_parsing_folder(), ...
    locking_summary_filename(false, method, count));
end


function [low, high] = population_stratum_rates(summary_path, session, cluster_id)
%POPULATION_STRATUM_RATES  Tail rates aligned to the population rows.
%
% Unmatched neurons, and neurons without a positive rate in both tails,
% are NaN.

n_pop = numel(cluster_id);
low = nan(n_pop, 1);
high = nan(n_pop, 1);
if isempty(summary_path) || exist(summary_path, 'file') ~= 2
    error('Plot_acg_population_logit:powerModel', ...
        'The locking summary was not found at %s.', summary_path);
end
try
    stored = load(summary_path, 'power_acg');
catch
    error('Plot_acg_population_logit:powerModel', ...
        'The locking summary could not be read.');
end
if ~isfield(stored, 'power_acg') || ~isfield(stored.power_acg, 'all_neurons')
    error('Plot_acg_population_logit:powerModel', ...
        'The summary has no power_acg.all_neurons.');
end
block = stored.power_acg.all_neurons;
if ~isfield(block, 'rate_low') || ~isfield(block, 'rate_high') ...
        || ~isfield(block, 'session') || ~isfield(block, 'cluster_id')
    error('Plot_acg_population_logit:powerModel', ...
        'The power summary has no session, cluster id, rate_low, and rate_high.');
end
rate_low = double(block.rate_low(:));
rate_high = double(block.rate_high(:));
n = min(numel(rate_low), numel(rate_high));
summary_keys = neuron_keys(block.session, block.cluster_id);
n = min(n, numel(summary_keys));
pop_keys = neuron_keys(session, cluster_id);
[found, loc] = ismember(pop_keys, summary_keys(1:n));
use = found(:) & loc(:) >= 1 & loc(:) <= n;
low(use) = rate_low(loc(use));
high(use) = rate_high(loc(use));
bad = ~(isfinite(low) & isfinite(high) & low > 0 & high > 0);
low(bad) = nan;
high(bad) = nan;
end


function report = feature_change_vs_rate(gap_x, names)
%FEATURE_CHANGE_VS_RATE  Relate each feature gap to the firing-rate gap.

names = names(:)';
rate_col = find(strcmp(names, 'FiringRate'), 1);
report.n = 0;
report.lines = {'FiringRate is not among the predictors, so the feature comparison was not run.'};
if isempty(rate_col)
    return
end
rate = gap_x(:, rate_col);
feature_names = names;
feature_names(rate_col) = [];
features = gap_x;
features(:, rate_col) = [];
n_features = numel(feature_names);
pearson = nan(n_features, 1);
pearson_p = nan(n_features, 1);
spearman = nan(n_features, 1);
spearman_p = nan(n_features, 1);
slope = nan(n_features, 1);
slope_se = nan(n_features, 1);
slope_p = nan(n_features, 1);
rsq = nan(n_features, 1);
for j = 1:n_features
    x = features(:, j);
    ok = isfinite(x) & isfinite(rate);
    n = sum(ok);
    pearson(j) = column_pearson(x(ok), rate(ok));
    pearson_p(j) = correlation_p(pearson(j), n);
    spearman(j) = column_pearson(average_rank(x(ok)), average_rank(rate(ok)));
    spearman_p(j) = correlation_p(spearman(j), n);
    if n >= 3
        fit = fit_power_feature_model(x(ok), rate(ok), {'FiringRate'});
        slope(j) = fit.beta(2);
        slope_se(j) = fit.se(2);
        slope_p(j) = fit.p_value(2);
        rsq(j) = fit.rsq;
    end
end
holm = holm_sidak_adjust(pearson_p);
report.n = sum(isfinite(rate));
report.names = feature_names;
report.rate = rate;
report.features = features;
report.pearson = pearson;
report.pearson_p = pearson_p;
report.holm = holm;
report.slope = slope;
report.slope_se = slope_se;
report.lines = cell(n_features + 1, 1);
report.lines{1} = sprintf('%-18s %8s %10s %10s %8s %10s %10s %10s %8s', ...
    'Feature', 'Pearson', 'p', 'Holm', 'Spearman', 'p', 'Slope', 'SE', 'R^2');
for j = 1:n_features
    report.lines{j + 1} = sprintf('%-18s %8.3f %10.3g %10.3g %8.3f %10.3g %10.3g %10.3g %8.3f', ...
        feature_names{j}, pearson(j), pearson_p(j), holm(j), ...
        spearman(j), spearman_p(j), slope(j), slope_se(j), rsq(j));
end
end


function r = column_pearson(a, b)
a = a(:) - mean(a);
b = b(:) - mean(b);
denom = sqrt(sum(a .^ 2) * sum(b .^ 2));
if ~(denom > 0)
    r = nan;
    return
end
r = sum(a .* b) / denom;
end


function p = correlation_p(r, n)
if n < 3 || ~isfinite(r)
    p = nan;
    return
end
if abs(r) >= 1
    p = 0;
    return
end
tstat = r * sqrt((n - 2) / (1 - r ^ 2));
p = two_sided_t_p(tstat, n - 2);
end


function ranks = average_rank(x)
x = x(:);
n = numel(x);
ranks = nan(n, 1);
if n == 0
    return
end
[sorted, order] = sort(x);
ranks_in_order = zeros(n, 1);
k = 1;
while k <= n
    j = k;
    while j < n && sorted(j + 1) == sorted(k)
        j = j + 1;
    end
    ranks_in_order(k:j) = (k + j) / 2;
    k = j + 1;
end
ranks(order) = ranks_in_order;
end


function draw_feature_rate_link(report, folder)
if ~isfield(report, 'pearson') || isempty(report.pearson)
    return
end
open_full_figure('Feature change versus firing-rate change')
hold on
for i = 1:numel(report.pearson)
    significant = isfinite(report.pearson_p(i)) && report.pearson_p(i) < 0.05;
    if significant
        color = [0.85 0.15 0.15];
    else
        color = [0 0 0];
    end
    plot(i, report.pearson(i), 'o', 'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, 'MarkerSize', 8)
    if isfinite(report.pearson_p(i))
        text(i, report.pearson(i), sprintf('p = %.3g', report.pearson_p(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none', 'Color', color)
    end
end
plot([0.5, numel(report.pearson) + 0.5], [0 0], 'Color', [0.5 0.5 0.5])
set(gca, 'XTick', 1:numel(report.names), 'XTickLabel', report.names)
try
    set(gca, 'TickLabelInterpreter', 'none')
    xtickangle(30)
catch
end
ylabel('Pearson r, feature change vs firing-rate change')
ylim([-1.25 1.35])
box off
sgtitle({'Each correlation uses one neuron with both power tails.', ...
    'The feature change and the firing-rate change are high minus low.', ...
    'Red marks a Pearson p below 0.05. The number above each dot is that p.'}, ...
    'Interpreter', 'none')
save_script_figure(folder)
end


function note = draw_ceteris_paribus(result, folder, figure_name)
%DRAW_CETERIS_PARIBUS  One panel per significant slope, other predictors at their means.

note = 'No slope has p below 0.05, so the ceteris paribus figure was not drawn.';
if ~isfield(result, 'fit') || ~isfield(result, 'predictors')
    return
end
fit = result.fit;
X = result.predictors;
y = result.outcome(:);
if size(X, 2) ~= numel(fit.names) - 1 || size(X, 1) ~= numel(y)
    return
end
significant = find(fit.p_value(2:end) < 0.05);
if isempty(significant)
    return
end
n_panels = numel(significant);
n_cols = min(3, n_panels);
n_rows = ceil(n_panels / n_cols);
open_full_figure(figure_name)
tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact')
predictor_means = mean(X, 1);
for s = 1:n_panels
    nexttile
    column = significant(s);
    [x_line, y_line, x_pts, y_pts] = ceteris_paribus_curve(y, X, fit.beta, column, predictor_means);
    hold on
    plot(x_pts, y_pts, 'o', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 5)
    plot(x_line, y_line, 'Color', [0.75 0.15 0.15], 'LineWidth', 1.6)
    xlabel(fit.names{column + 1}, 'Interpreter', 'none')
    ylabel(result.outcome_label, 'Interpreter', 'none')
    title(sprintf('p = %.3g', fit.p_value(column + 1)), 'Interpreter', 'none')
    box off
end
shown_names = fit.names(significant + 1);
sgtitle({figure_name, ...
    'Each panel varies one predictor. The others stay at their mean.', ...
    'Points are the power gap after that adjustment. The line is the model.'}, ...
    'Interpreter', 'none')
save_script_figure(folder)
note = sprintf(['Ceteris paribus panels, other predictors fixed at their mean: %s. ', ...
    'A panel is drawn when that slope has p below 0.05.'], strjoin(shown_names, ', '));
end


function [x_line, y_line, x_pts, y_pts] = ceteris_paribus_curve(y, X, beta, column, predictor_means)
%CETERIS_PARIBUS_CURVE  Outcome against one predictor, others fixed at their means.

y = y(:);
beta = beta(:);
n_predictors = size(X, 2);
others = true(1, n_predictors);
others(column) = false;
x_pts = X(:, column);
adjustment = zeros(numel(y), 1);
if any(others)
    other_beta = beta(find(others) + 1);
    adjustment = (X(:, others) - predictor_means(others)) * other_beta;
end
y_pts = y - adjustment;
x_line = linspace(min(x_pts), max(x_pts), 40).';
held = 0;
if any(others)
    held = predictor_means(others) * beta(find(others) + 1);
end
y_line = beta(1) + beta(column + 1) * x_line + held;
end


function note = draw_feature_rate_scatter(report, folder)
%DRAW_FEATURE_RATE_SCATTER  Significant feature changes against the firing-rate change.

note = 'No feature has a Pearson p below 0.05, so the scatter was not drawn.';
if ~isfield(report, 'pearson_p') || ~isfield(report, 'features') || isempty(report.features)
    return
end
significant = find(report.pearson_p(:).' < 0.05);
if isempty(significant)
    return
end
n_panels = numel(significant);
n_cols = min(3, n_panels);
n_rows = ceil(n_panels / n_cols);
open_full_figure('Significant feature changes versus firing-rate change')
tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact')
for s = 1:n_panels
    nexttile
    j = significant(s);
    x = report.rate;
    y = report.features(:, j);
    ok = isfinite(x) & isfinite(y);
    hold on
    plot(x(ok), y(ok), 'o', 'MarkerFaceColor', [0.25 0.25 0.25], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 5)
    if isfinite(report.slope(j))
        x_line = linspace(min(x(ok)), max(x(ok)), 40);
        intercept = mean(y(ok)) - report.slope(j) * mean(x(ok));
        plot(x_line, intercept + report.slope(j) * x_line, ...
            'Color', [0.85 0.15 0.15], 'LineWidth', 1.6)
    end
    xlabel('Firing-rate change (spikes/s)')
    ylabel(report.names{j}, 'Interpreter', 'none')
    title(sprintf('Pearson p = %.3g', report.pearson_p(j)), 'Interpreter', 'none')
    box off
end
sgtitle({'Significant feature changes against the firing-rate change.', ...
    'Each point is one neuron. The line is that feature change regressed on the firing-rate change.', ...
    'A panel is included when the Pearson p is below 0.05.'}, ...
    'Interpreter', 'none')
save_script_figure(folder)
note = sprintf('Scatter panels for Pearson p below 0.05: %s.', ...
    strjoin(report.names(significant), ', '));
end


function values = population_summary_field(summary_path, session, cluster_id, field_name)
%POPULATION_SUMMARY_FIELD  One saved power-summary column, aligned to the population.

n = min(numel(session), numel(cluster_id));
values = nan(n, 1);
if isempty(summary_path) || exist(summary_path, 'file') ~= 2
    error('Plot_acg_population_logit:areaNest', ...
        'The locking summary was not found at %s.', summary_path);
end
try
    stored = load(summary_path, 'power_acg');
catch
    error('Plot_acg_population_logit:areaNest', ...
        'The locking summary could not be read.');
end
if ~isfield(stored, 'power_acg') || ~isfield(stored.power_acg, 'all_neurons')
    error('Plot_acg_population_logit:areaNest', ...
        'The summary has no power_acg.all_neurons.');
end
block = stored.power_acg.all_neurons;
if ~isfield(block, field_name) || ~isfield(block, 'session') || ~isfield(block, 'cluster_id')
    error('Plot_acg_population_logit:areaNest', ...
        'The power summary has no %s aligned by session and cluster id.', field_name);
end
score = double(block.(field_name)(:));
summary_keys = neuron_keys(block.session, block.cluster_id);
n_score = min(numel(score), numel(summary_keys));
pop_keys = neuron_keys(session(1:n), cluster_id(1:n));
[found, loc] = ismember(pop_keys, summary_keys(1:n_score));
use = found(:) & loc(:) >= 1 & loc(:) <= n_score;
values(use) = score(loc(use));
if ~any(isfinite(values))
    error('Plot_acg_population_logit:areaNest', ...
        'No population neuron matched a finite %s in the power summary.', field_name);
end
end


function root = experiment_data_root()
root = ['\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\' ...
    'NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data'];
end


function folder = logit_result_folder()
%LOGIT_RESULT_FOLDER  Neuron table, population matrices, and the lock-type cache.
folder = fullfile(experiment_data_root(), 'Analysis results', 'acg population logit');
end


function [figure_folder, estimate_folder] = logit_folders(figure_folder)
%LOGIT_FOLDERS  Use the output_folder given by the caller for figures.
estimate_folder = logit_result_folder();
if isempty(figure_folder)
    error('Plot_acg_population_logit:outputFolder', ...
        ['Set output_folder to the folder where figures should be saved, ' ...
        'then run this section.']);
end
figure_folder = char(figure_folder);
if exist(figure_folder, 'dir') ~= 7
    mkdir(figure_folder);
end
if exist(estimate_folder, 'dir') ~= 7
    mkdir(estimate_folder);
end
end


function path = logit_prep_file(estimate_folder, figure_folder)
%LOGIT_PREP_FILE  Population matrices, including one saved with an older run.
path = fullfile(estimate_folder, 'acg_population_logit_prep.mat');
if exist(path, 'file') == 2
    return
end
candidates = { ...
    fullfile(default_parsing_folder(), 'acg_population_logit_prep.mat'), ...
    fullfile(figure_folder, 'acg_population_logit_prep.mat')};
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file') == 2
        path = candidates{i};
        return
    end
end
end


function path = lock_cache_file(estimate_folder, figure_folder)
%LOCK_CACHE_FILE  Lock-type fit, copied into the estimate folder when found elsewhere.
path = fullfile(estimate_folder, 'lock_type_model_cache.mat');
if exist(path, 'file') == 2
    return
end
candidates = { ...
    fullfile(figure_folder, 'lock_type_model_cache.mat'), ...
    fullfile(fileparts(mfilename('fullpath')), ...
        'acg_population_logit_figures', 'lock_type_model_cache.mat')};
for i = 1:numel(candidates)
    legacy = candidates{i};
    if exist(legacy, 'file') ~= 2 || strcmpi(legacy, path)
        continue
    end
    if exist(estimate_folder, 'dir') ~= 7
        mkdir(estimate_folder);
    end
    try
        copyfile(legacy, path);
        fprintf('Copied the saved lock-type models to %s\n', path)
        return
    catch
        fprintf('Using the saved lock-type models at %s\n', legacy)
        path = legacy;
        return
    end
end
end


function folder = default_parsing_folder()
folder = fullfile(experiment_data_root(), 'Analysis results', 'spike train parsing');
end


function channel = channels_for_neurons(parsing_folder, session, cluster_id)
%CHANNELS_FOR_NEURONS  Channel number from the parsing cluster tables.
%
% The stored channel number is used as-is. Neurons with the same session
% and the same number are treated as one channel.

n = min(numel(session), numel(cluster_id));
channel = nan(n, 1);
if isempty(parsing_folder) || exist(parsing_folder, 'dir') ~= 7
    error('Plot_acg_population_logit:areaNest', ...
        'The parsing folder was not found at %s.', parsing_folder);
end
files = parsing_mat_files(parsing_folder);
if isempty(files)
    error('Plot_acg_population_logit:areaNest', ...
        'No spike-train parsing file was found in %s.', parsing_folder);
end
map_keys = cell(0, 1);
map_channel = zeros(0, 1);
seen_session = {};
for f = 1:numel(files)
    try
        stored = load(files(f).path, 'parsing_struct');
    catch
        continue
    end
    if ~isfield(stored, 'parsing_struct')
        continue
    end
    blocks = stored.parsing_struct;
    for b = 1:numel(blocks)
        block = blocks(b);
        if ~isfield(block, 'animal_code') || ~isfield(block, 'cluster_info')
            continue
        end
        animal = strtrim(value_text(block.animal_code));
        if isempty(animal) || any(strcmp(seen_session, animal))
            continue
        end
        [id_col, ch_col] = cluster_channel_columns(block.cluster_info);
        if isempty(id_col)
            continue
        end
        seen_session{end + 1} = animal; %#ok<AGROW>
        n_row = min(numel(id_col), numel(ch_col));
        new_keys = neuron_keys(repmat({animal}, n_row, 1), id_col(1:n_row));
        map_keys = [map_keys; new_keys]; %#ok<AGROW>
        map_channel = [map_channel; double(ch_col(1:n_row))]; %#ok<AGROW>
    end
    clear stored
end
if isempty(map_keys)
    error('Plot_acg_population_logit:areaNest', ...
        'The parsing files have no cluster channel column.');
end
pop_keys = neuron_keys(session(1:n), cluster_id(1:n));
[found, loc] = ismember(pop_keys, map_keys);
use = found(:) & loc(:) >= 1;
channel(use) = map_channel(loc(use));
end


function files = parsing_mat_files(parsing_folder)
files = struct('path', {}, 'name', {});
combined = fullfile(parsing_folder, 'spike_train_parsing_structure.mat');
if exist(combined, 'file') == 2
    files(end + 1) = struct('path', combined, ...
        'name', 'spike_train_parsing_structure.mat'); %#ok<AGROW>
end
flist = dir(fullfile(parsing_folder, '*spike_train_parsing*.mat'));
skip = {'spike_train_parsing_structure.mat', 'spike_train_parsing_animal_names.mat', ...
    'spike_train_parsing_animal_names_v1.mat'};
for f = 1:numel(flist)
    if ismember(flist(f).name, skip)
        continue
    end
    files(end + 1) = struct( ...
        'path', fullfile(flist(f).folder, flist(f).name), ...
        'name', flist(f).name); %#ok<AGROW>
end
end


function [id_col, ch_col] = cluster_channel_columns(info)
id_col = [];
ch_col = [];
names = cluster_table_names(info);
if isempty(names)
    return
end
id_name = first_present(names, {'cluster_id', 'id'});
ch_name = first_present(names, {'ch', 'channel'});
if isempty(id_name) || isempty(ch_name)
    return
end
id_col = double(cluster_table_column(info, id_name));
ch_col = double(cluster_table_column(info, ch_name));
end


function names = cluster_table_names(info)
names = {};
if isstruct(info)
    names = fieldnames(info);
    return
end
if exist('istable', 'file') == 2
    try
        if istable(info)
            names = info.Properties.VariableNames;
        end
    catch
    end
end
end


function column = cluster_table_column(info, name)
if isstruct(info)
    column = info.(name);
else
    column = info.(name);
end
column = column(:);
end


function name = first_present(names, candidates)
name = '';
for i = 1:numel(candidates)
    if any(strcmp(names, candidates{i}))
        name = candidates{i};
        return
    end
end
end


function result = nested_area_power(power, area, session, channel, condition)
%NESTED_AREA_POWER  Area test after collapsing the nested repeated measures.
%
% Neurons on one channel are averaged. Those channel means, within one
% area in one session, are then averaged. Area is tested on that
% session-area table after a session term.

result.ok = false;
result.condition = condition;
if strcmp(condition, 'all')
    condition_label = 'all runs';
else
    condition_label = [condition, '-power tail'];
end
result.condition_label = condition_label;
result.lines = {['Condition: ', condition_label, '.']};
if strcmp(condition, 'all')
    result.lines{end + 1} = ['All-runs power is the mean delta power of every run with a finite score. ', ...
        'Those runs are not restricted to the outer power percentiles.'];
end
n = min([numel(power), numel(area), numel(session), numel(channel)]);
power = power(1:n);
area = area(1:n);
session = session(1:n);
channel = channel(1:n);
area = area_labels(area);
session = session_labels(session);
positive = isfinite(power) & power > 0;
scale = choose_level_scale(power(positive));
if strcmp(scale.name, 'log')
    analyzed = nan(size(power));
    analyzed(positive) = log(power(positive));
    scale_label = ['log ', condition_label, ' delta power'];
else
    analyzed = power;
    scale_label = [condition_label, ' delta power'];
end
result.scale_label = scale_label;
packed = collapse_to_session_area(analyzed, area, session, channel);
result.packed = packed;
result.lines = [result.lines, { ...
    'Individual neuron values are not the rows of the test. Neurons that share a session and a channel are averaged first.', ...
    'Those channel means that share a session and an area are then averaged. Each remaining row is one area in one session.', ...
    'The model has a separate intercept for each session, because sessions differ in overall power.', ...
    'Area is tested after those session intercepts. The error is the remaining variation of the area differences across sessions.', ...
    'The figure subtracts the fitted session intercept from each point.', ...
    'When the area test is below 0.05, every area pair is compared in the same model.', ...
    'Holm-Sidak adjusts those pairs. A line and stars are drawn for each pair that stays below 0.05.', ...
    sprintf('Neurons with a finite %s value and a known channel: %d.', condition_label, packed.n_neurons), ...
    sprintf('Channels: %d. Session-area rows: %d.', packed.n_channels, numel(packed.y)), ...
    sprintf('Neurons with no channel match: %d.', packed.n_unmatched), ...
    sprintf('Channels whose neurons were not all labeled with one area: %d. Those channels use the majority area.', ...
        packed.n_area_conflict), ...
    scale.reason, ...
    sprintf('Fitted values: %s.', scale_label)}];
if numel(packed.y) < 4 || numel(unique(packed.area)) < 2 || numel(unique(packed.session)) < 2
    result.lines{end + 1} = 'Too few session-area rows, areas, or sessions for the nested area test.';
    return
end
tested = test_session_area(packed.y, packed.area, packed.session, ...
    packed.row_channels, packed.row_neurons);
result.ok = tested.ok;
result.tested = tested;
result.lines = [result.lines, tested.lines];
end


function labels = area_labels(raw)
labels = cell(numel(raw), 1);
for i = 1:numel(raw)
    labels{i} = area_display_name(raw{i});
end
end


function labels = session_labels(raw)
labels = cell(numel(raw), 1);
for i = 1:numel(raw)
    labels{i} = strtrim(value_text(raw{i}));
end
end


function scale = choose_level_scale(y)
y = y(isfinite(y));
scale.skew_raw = moment_skewness(y);
scale.skew_log = nan;
scale.name = 'raw';
if numel(y) < 3 || any(y <= 0) || ~isfinite(scale.skew_raw)
    scale.reason = 'The power is not positive for every used neuron, so the nested test uses the power itself.';
    return
end
scale.skew_log = moment_skewness(log(y));
if isfinite(scale.skew_log) && abs(scale.skew_log) < abs(scale.skew_raw)
    scale.name = 'log';
    scale.reason = sprintf(['Absolute skewness is %.3f for the power and %.3f for its logarithm. ', ...
        'The logarithm is less skewed, so the nested test uses log power. ', ...
        'The logarithm is taken before the channel and session-area averages.'], ...
        abs(scale.skew_raw), abs(scale.skew_log));
    return
end
scale.reason = sprintf(['Absolute skewness is %.3f for the power and %.3f for its logarithm. ', ...
    'The power is less skewed, so the nested test uses the power itself.'], ...
    abs(scale.skew_raw), abs(scale.skew_log));
end


function packed = collapse_to_session_area(y, area, session, channel)
y = y(:);
channel = channel(:);
area = area(:);
session = session(:);
n = numel(y);
known = isfinite(y) & isfinite(channel) & ~strcmp(area, 'area unknown');
for i = 1:n
    if isempty(session{i})
        known(i) = false;
    end
end
packed.n_neurons = sum(known);
packed.n_unmatched = sum(isfinite(y) & ~isfinite(channel));
channel_key = cell(n, 1);
for i = 1:n
    if ~known(i)
        channel_key{i} = '';
        continue
    end
    channel_key{i} = [session{i}, '|', rate_id_text(channel(i))];
end
keys = unique(channel_key(known), 'stable');
n_channel = numel(keys);
channel_y = nan(n_channel, 1);
channel_area = cell(n_channel, 1);
channel_session = cell(n_channel, 1);
channel_n = zeros(n_channel, 1);
n_conflict = 0;
for c = 1:n_channel
    in = known & strcmp(channel_key, keys{c});
    channel_y(c) = mean(y(in));
    channel_n(c) = sum(in);
    channel_area{c} = majority_label(area(in));
    channel_session{c} = session{find(in, 1)};
    if numel(unique(area(in))) > 1
        n_conflict = n_conflict + 1;
    end
end
packed.n_channels = n_channel;
packed.n_area_conflict = n_conflict;
pair_key = cell(n_channel, 1);
for c = 1:n_channel
    pair_key{c} = [channel_session{c}, '|', channel_area{c}];
end
if n_channel == 0
    pairs = cell(0, 1);
else
    pairs = unique(pair_key, 'stable');
end
n_pair = numel(pairs);
packed.y = nan(n_pair, 1);
packed.area = cell(n_pair, 1);
packed.session = cell(n_pair, 1);
packed.row_channels = zeros(n_pair, 1);
packed.row_neurons = zeros(n_pair, 1);
for p = 1:n_pair
    in = strcmp(pair_key, pairs{p});
    packed.y(p) = mean(channel_y(in));
    packed.row_channels(p) = sum(in);
    packed.row_neurons(p) = sum(channel_n(in));
    hit = find(in, 1);
    packed.area{p} = channel_area{hit};
    packed.session{p} = channel_session{hit};
end
end


function label = majority_label(labels)
labels = labels(:);
levels = unique(labels);
counts = zeros(numel(levels), 1);
for i = 1:numel(levels)
    counts(i) = sum(strcmp(labels, levels{i}));
end
[~, pick] = max(counts);
label = levels{pick};
end


function tested = test_session_area(y, area, session, n_channels, n_neurons)
y = y(:);
area_levels = sort(unique(area));
session_levels = sort(unique(session));
n_area = numel(area_levels);
n_session = numel(session_levels);
X = session_area_design(area, area_levels, session, session_levels);
[rss_full, rank_full] = ols_rss(X, y);
n_a = n_area - 1;
keep = true(1, size(X, 2));
keep(2:(1 + n_a)) = false;
[rss_reduced, rank_reduced] = ols_rss(X(:, keep), y);
df1 = rank_full - rank_reduced;
df2 = numel(y) - rank_full;
tested.ok = df1 >= 1 && df2 >= 1;
tested.areas = area_levels;
tested.y = y;
tested.area = area;
tested.session = session;
tested.lines = { ...
    sprintf('Areas: %d. Sessions: %d. Session-area rows: %d.', n_area, n_session, numel(y)), ...
    sprintf('%-12s %8s %8s %8s %12s', 'Area', 'Sessions', 'Channels', 'Neurons', 'Raw mean')};
tested.lines{end + 1} = 'The raw mean still includes each session''s power level.';
for a = 1:n_area
    in = strcmp(area, area_levels{a});
    tested.lines{end + 1} = sprintf('%-12s %8d %8d %8d %12.4g', ...
        area_levels{a}, numel(unique(session(in))), sum(n_channels(in)), ...
        sum(n_neurons(in)), mean(y(in))); %#ok<AGROW>
end
if ~tested.ok
    tested.lines{end + 1} = 'The session and area terms do not leave a residual, so the area test was not computed.';
    return
end
ss = rss_reduced - rss_full;
scale = max([1, abs(rss_reduced), abs(rss_full)]);
if ss < 0 && ss > -1e-8 * scale
    ss = 0;
end
if rss_full <= 1e-10 * scale && ss > 1e-10 * scale
    F = inf;
    p = 0;
else
    F = (ss / df1) / (rss_full / df2);
    p = f_survival(F, df1, df2);
end
tested.F = F;
tested.df1 = df1;
tested.df2 = df2;
tested.p = p;
tested.lines{end + 1} = sprintf(['Area after the session intercepts: F(%g, %g) = %.3f, p = %.3g. ', ...
    'This is the nested test of area differences.'], df1, df2, F, p);
beta = X \ y;
n_s = n_session - 1;
session_cols = (2 + n_a):(1 + n_a + n_s);
tested.adjusted = y - X(:, session_cols) * beta(session_cols);
sigma2 = rss_full / df2;
xtx = X' * X;
if ~isfinite(rcond(xtx)) || rcond(xtx) < 1e-8
    C = pinv(xtx);
else
    C = inv(xtx);
end
C = sigma2 * C;
area_cov = C(2:(1 + n_a), 2:(1 + n_a));
area_beta = beta(2:(1 + n_a));
effect = [area_beta; -sum(area_beta)];
se = nan(n_area, 1);
for a = 1:n_a
    se(a) = sqrt(max(real(area_cov(a, a)), 0));
end
contrast = -ones(n_a, 1);
se(n_area) = sqrt(max(real(contrast' * area_cov * contrast), 0));
tstat = effect ./ se;
tstat(~isfinite(se) | se <= 0) = nan;
raw_p = nan(n_area, 1);
for a = 1:n_area
    raw_p(a) = two_sided_t_p(tstat(a), df2);
end
holm = holm_sidak_adjust(raw_p);
tested.effect = effect;
tested.area_level = beta(1) + effect;
tested.lines{end + 1} = ['Each area level is the intercept plus that area''s deviation. ', ...
    'The deviation is the area minus the unweighted mean of the areas, after the session intercepts.'];
tested.lines{end + 1} = sprintf('%-12s %12s %12s %12s %10s %10s %10s', ...
    'Area', 'Level', 'Deviation', 'SE', 't', 'p', 'Holm');
for a = 1:n_area
    tested.lines{end + 1} = sprintf('%-12s %12.4g %12.4g %12.4g %10.3f %10.3g %10.3g', ...
        area_levels{a}, tested.area_level(a), effect(a), se(a), tstat(a), raw_p(a), holm(a)); %#ok<AGROW>
end
n_pairs = n_area * (n_area - 1) / 2;
pair_index = zeros(n_pairs, 2);
pair_diff = nan(n_pairs, 1);
pair_se = nan(n_pairs, 1);
pair_raw = nan(n_pairs, 1);
k = 0;
for a = 1:n_area
    for b = a + 1:n_area
        k = k + 1;
        pair_index(k, :) = [a, b];
        contrast = area_pair_contrast(a, b, n_a);
        difference = contrast' * area_beta;
        pair_se_k = sqrt(max(real(contrast' * area_cov * contrast), 0));
        pair_diff(k) = difference;
        pair_se(k) = pair_se_k;
        if pair_se_k > 0
            pair_raw(k) = two_sided_t_p(difference / pair_se_k, df2);
        elseif difference ~= 0
            pair_raw(k) = 0;
        else
            pair_raw(k) = 1;
        end
    end
end
run_pairs = isfinite(p) && p < 0.05;
tested.pair_index = zeros(0, 2);
tested.pair_holm = zeros(0, 1);
if ~run_pairs
    tested.lines{end + 1} = 'Pairwise area comparisons were not run, because the area test was not below 0.05.';
    return
end
pair_holm = holm_sidak_adjust(pair_raw);
tested.pair_index = pair_index;
tested.pair_holm = pair_holm(:);
tested.lines{end + 1} = ['Pairwise differences are area levels from the same model. ', ...
    'Holm-Sidak adjusts those pairs. The figure draws a line and stars when the adjusted p is below 0.05.'];
tested.lines{end + 1} = sprintf('%-28s %12s %12s %10s %10s %10s', ...
    'Pair', 'Difference', 'SE', 't', 'p', 'Holm');
for k = 1:n_pairs
    if pair_se(k) > 0
        t_pair = pair_diff(k) / pair_se(k);
    elseif pair_diff(k) ~= 0
        t_pair = inf;
    else
        t_pair = 0;
    end
    tested.lines{end + 1} = sprintf('%-28s %12.4g %12.4g %10.3f %10.3g %10.3g', ...
        sprintf('%s vs %s', area_levels{pair_index(k, 1)}, area_levels{pair_index(k, 2)}), ...
        pair_diff(k), pair_se(k), t_pair, pair_raw(k), pair_holm(k)); %#ok<AGROW>
end
end


function contrast = area_pair_contrast(first_area, second_area, n_coded)
%AREA_PAIR_CONTRAST  Coefficients for the first area level minus the second.
%
% The last area is the omitted effect-coded level, equal to minus the sum
% of the coded area coefficients.

contrast = zeros(n_coded, 1);
if first_area <= n_coded
    contrast(first_area) = contrast(first_area) + 1;
else
    contrast = contrast - 1;
end
if second_area <= n_coded
    contrast(second_area) = contrast(second_area) - 1;
else
    contrast = contrast + 1;
end
end


function X = session_area_design(area, area_levels, session, session_levels)
n = numel(area);
n_a = numel(area_levels) - 1;
n_s = numel(session_levels) - 1;
X = zeros(n, 1 + n_a + n_s);
X(:, 1) = 1;
for a = 1:n_a
    code = zeros(n, 1);
    code(strcmp(area, area_levels{a})) = 1;
    code(strcmp(area, area_levels{end})) = -1;
    X(:, 1 + a) = code;
end
for s = 1:n_s
    code = zeros(n, 1);
    code(strcmp(session, session_levels{s})) = 1;
    code(strcmp(session, session_levels{end})) = -1;
    X(:, 1 + n_a + s) = code;
end
end


function draw_nested_area_power(result, folder)
open_full_figure(['Delta power by area, ', result.condition])
areas = result.tested.areas;
hold on
for a = 1:numel(areas)
    in = strcmp(result.tested.area, areas{a});
    values = result.tested.adjusted(in);
    x = a + violin_point_offset(numel(values));
    color = area_bar_color(areas{a});
    plot(x, values, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none')
    plot([a - 0.25, a + 0.25], result.tested.area_level(a) * [1 1], ...
        'Color', [0.1 0.1 0.1], 'LineWidth', 1.6)
end
y_top = nan;
if isfield(result.tested, 'pair_holm') && any(result.tested.pair_holm < 0.05)
    show = result.tested.pair_holm < 0.05;
    y_top = draw_sig_brackets(result.tested.pair_index(show, :), ...
        result.tested.pair_holm(show), result.tested.adjusted);
end
xlim([0.4, numel(areas) + 0.6])
set(gca, 'XTick', 1:numel(areas), 'XTickLabel', areas)
try
    set(gca, 'TickLabelInterpreter', 'none')
catch
end
ylabel(['Session-adjusted ', result.scale_label])
box off
if isfinite(y_top)
    yl = ylim;
    ylim([yl(1), max(yl(2), y_top)])
end
sgtitle({sprintf('%s. Area after session intercepts: F(%g, %g) = %.3f, p = %.3g.', ...
    result.condition_label, result.tested.df1, result.tested.df2, ...
    result.tested.F, result.tested.p), ...
    'Each point is one session, after that session intercept is subtracted.', ...
    'It starts as the mean of the channel means in that area.', ...
    'The line is the area level from the model.', ...
    'A connecting line and stars mark a pair whose Holm-Sidak p is below 0.05.'}, ...
    'Interpreter', 'none')
save_script_figure(folder)
end


function report = load_stratum_firing_rates(summary_path, session, cluster_id)
%LOAD_STRATUM_FIRING_RATES  Paired rates from the saved power summary.
%
% rate_low and rate_high are the duration-weighted firing rates of the
% runs in each power tail. model_fr is that rate for the runs that enter
% the all-runs autocorrelogram.

report.low = zeros(0, 1);
report.high = zeros(0, 1);
report.all_rate = zeros(0, 1);
report.lines = {['Summary: ', summary_path]};
if isempty(summary_path) || exist(summary_path, 'file') ~= 2
    report.lines{end + 1} = 'The locking summary was not found, so the stratum firing rates were not read.';
    return
end
try
    stored = load(summary_path, 'power_acg');
catch
    report.lines{end + 1} = 'The locking summary could not be read.';
    return
end
if ~isfield(stored, 'power_acg') || ~isfield(stored.power_acg, 'all_neurons')
    report.lines{end + 1} = 'The summary has no power_acg.all_neurons.';
    return
end
block = stored.power_acg.all_neurons;
if ~isfield(block, 'rate_low') || ~isfield(block, 'rate_high')
    report.lines{end + 1} = 'The power summary has no rate_low and rate_high.';
    return
end
low = double(block.rate_low(:));
high = double(block.rate_high(:));
if isfield(block, 'model_fr') && numel(block.model_fr) == numel(low)
    all_rate = double(block.model_fr(:));
    all_name = 'all-runs autocorrelogram';
    all_detail = ['model_fr is the spike count of the runs that enter the all-runs autocorrelogram, ', ...
        'divided by the summed duration of those runs.'];
elseif isfield(block, 'rate_all')
    all_rate = double(block.rate_all(:));
    all_name = 'all firing runs';
    all_detail = ['rate_all is the spike count of every firing run, divided by the summed duration. ', ...
        'The summary has no model_fr for the all-runs autocorrelogram alone.'];
else
    report.lines{end + 1} = 'The power summary has no model_fr or rate_all.';
    return
end
n = min([numel(low), numel(high), numel(all_rate)]);
low = low(1:n);
high = high(1:n);
all_rate = all_rate(1:n);
if isempty(session)
    keep = true(n, 1);
    match_line = sprintf(['Neurons in the power summary: %d. Population identifiers were not ', ...
        'in the workspace, so every summary neuron is included.'], n);
elseif ~isfield(block, 'session') || ~isfield(block, 'cluster_id')
    keep = true(n, 1);
    match_line = sprintf(['Neurons in the power summary: %d. The summary has no session and ', ...
        'cluster id, so every summary neuron is included.'], n);
else
    summary_keys = neuron_keys(block.session, block.cluster_id);
    if numel(summary_keys) < n
        keep = false(n, 1);
    else
        keep = ismember(summary_keys(1:n), neuron_keys(session, cluster_id));
    end
    match_line = sprintf('Neurons in the power summary: %d. Matched to this population: %d.', ...
        n, sum(keep));
end
ok = keep & isfinite(low) & isfinite(high) & isfinite(all_rate) ...
    & low > 0 & high > 0 & all_rate > 0;
report.low = low(ok);
report.high = high(ok);
report.all_rate = all_rate(ok);
report.lines = [report.lines, { ...
    'These rates are already stored on the power-split summary. They are not estimated again.', ...
    'Low-power rate: spikes in runs at or below the neuron''s 25th percentile of delta power, divided by the summed duration of those runs.', ...
    'High-power rate: the same quantity for runs at or above the 75th percentile.', ...
    all_detail, ...
    match_line, ...
    sprintf('Neurons with a positive rate in the low-power tail, the high-power tail, and the %s: %d.', ...
        all_name, numel(report.low))}];
if numel(report.low) < 2
    report.lines{end + 1} = 'Fewer than two neurons have all three rates, so the comparison stops here.';
    return
end
report.lines = [report.lines, { ...
    sprintf('Median low-power rate: %.4g sp/s. Median high-power rate: %.4g sp/s. Median %s rate: %.4g sp/s.', ...
        median(report.low), median(report.high), all_name, median(report.all_rate)), ...
    sprintf('Mean low-power rate: %.4g sp/s. Mean high-power rate: %.4g sp/s. Mean %s rate: %.4g sp/s.', ...
        mean(report.low), mean(report.high), all_name, mean(report.all_rate)), ...
    'Four paired tests are reported for each comparison. All four use the same neurons.', ...
    'The paired t-test asks whether the mean of the log ratio is zero. That is equal geometric means.', ...
    'The Wilcoxon signed-rank test asks whether those paired log differences are symmetric around zero.', ...
    'The sign test asks whether each rate is the higher one equally often. Ties are left out. Its p-value is the exact two-sided binomial probability.', ...
    'The Kolmogorov-Smirnov test asks whether the two firing-rate distributions match. Its p-value comes from swapping the two rates inside each neuron, 4999 times, with random seed 1. The reported p is one plus the number of swaps at least as large, divided by one plus the number of swaps.'}, ...
    ratio_lines('High power', 'low power', report.high, report.low), ...
    ratio_lines(all_name, 'low power', report.all_rate, report.low), ...
    ratio_lines(all_name, 'high power', report.all_rate, report.high)];
end


function lines = ratio_lines(name_a, name_b, a, b)
ratio = a(:) ./ b(:);
tests = paired_rate_tests(a, b);
within_two = mean(ratio >= 0.5 & ratio <= 2);
within_one_half = mean(ratio >= (2 / 3) & ratio <= 1.5);
lines = { ...
    sprintf(['%s / %s: median %.3g. Geometric mean %.3g. ', ...
        'Within a factor of 2: %.1f%%. Within a factor of 1.5: %.1f%%.'], ...
        name_a, name_b, median(ratio), exp(mean(log(ratio))), ...
        100 * within_two, 100 * within_one_half), ...
    sprintf('Paired t, log(%s) minus log(%s): mean %.4g, t %.3f, df %.0f, p %.3g.', ...
        name_a, name_b, tests.t.difference, tests.t.t, tests.t.df, tests.t.raw_p), ...
    sprintf('Wilcoxon signed-rank on that log difference: p %.3g.', tests.signrank_p), ...
    sprintf('Sign test: %s higher in %d neurons, %s higher in %d, ties %d, p %.3g.', ...
        name_a, tests.n_pos, name_b, tests.n_neg, tests.n_tie, tests.sign_p), ...
    sprintf('Paired permutation Kolmogorov-Smirnov: D %.3g, p %.3g.', ...
        tests.ks_d, tests.ks_p)};
end


function keys = neuron_keys(session, cluster_id)
if ischar(session)
    session = cellstr(session);
elseif isstring(session)
    session = cellstr(session);
end
session = session(:);
cluster_id = double(cluster_id(:));
n = min(numel(session), numel(cluster_id));
keys = cell(n, 1);
for i = 1:n
    label = session{i};
    if isstring(label)
        label = char(label);
    end
    if ~ischar(label)
        label = '';
    end
    keys{i} = [strtrim(label), '|', rate_id_text(cluster_id(i))];
end
end


function text = rate_id_text(value)
if isfinite(value) && value == round(value)
    text = sprintf('%.0f', value);
else
    text = sprintf('%.12g', value);
end
end


function tests = paired_rate_tests(a, b)
%PAIRED_RATE_TESTS  Location and distribution tests for one paired rate comparison.
%
% The t-test and the signed-rank test use log(a) minus log(b). The sign
% test uses the sign of a minus b. The Kolmogorov-Smirnov p-value swaps
% the two rates inside each neuron.

a = a(:);
b = b(:);
tests.t = compare_two_samples(log(a), log(b), 'paired');
tests.signrank_p = paired_log_signrank_p(a, b);
[tests.sign_p, tests.n_pos, tests.n_neg, tests.n_tie] = paired_sign_p(a, b);
[tests.ks_d, tests.ks_p] = paired_exchange_ks(a, b);
end


function p = paired_log_signrank_p(a, b)
log_a = log(a(:));
log_b = log(b(:));
if exist('paired_signrank_p', 'file') == 2
    try
        p = paired_signrank_p(log_a, log_b, 'both');
        if isfinite(p)
            return
        end
    catch
    end
end
p = signed_rank_normal_p(log_a - log_b);
end


function p = signed_rank_normal_p(d)
d = d(isfinite(d) & d ~= 0);
n = numel(d);
if n < 2
    p = nan;
    return
end
magnitude = abs(d);
[~, order] = sort(magnitude);
ranks = zeros(n, 1);
tie_cubic = 0;
k = 1;
while k <= n
    j = k;
    while j < n && magnitude(order(j + 1)) == magnitude(order(k))
        j = j + 1;
    end
    ranks(order(k:j)) = (k + j) / 2;
    t = j - k + 1;
    if t > 1
        tie_cubic = tie_cubic + t ^ 3 - t;
    end
    k = j + 1;
end
w_pos = sum(ranks(d > 0));
mu = n * (n + 1) / 4;
sigma2 = n * (n + 1) * (2 * n + 1) / 24 - tie_cubic / 48;
sigma = sqrt(max(sigma2, 0));
if sigma <= 0
    p = nan;
    return
end
delta = w_pos - mu;
z = (delta - 0.5 * sign(delta)) / sigma;
if delta == 0
    z = 0;
end
p = min(1, erfc(abs(z) / sqrt(2)));
end


function [p, n_pos, n_neg, n_tie] = paired_sign_p(a, b)
d = a(:) - b(:);
d = d(isfinite(d));
n_pos = sum(d > 0);
n_neg = sum(d < 0);
n_tie = sum(d == 0);
n = n_pos + n_neg;
if n < 1
    p = nan;
    return
end
p = min(1, 2 * binomial_cdf(min(n_pos, n_neg), n, 0.5));
end


function cdf = binomial_cdf(k, n, probability)
if k < 0
    cdf = 0;
    return
end
if k >= n
    cdf = 1;
    return
end
cdf = 1 - betainc(probability, k + 1, n - k);
end


function [D, p] = paired_exchange_ks(a, b)
a = a(:);
b = b(:);
ok = isfinite(a) & isfinite(b);
a = a(ok);
b = b(ok);
n = numel(a);
D = two_sample_ks(a, b);
n_perm = 4999;
if n < 2 || D == 0
    p = 1;
    return
end
previous = rng;
rng(1, 'twister');
extreme = 0;
for i = 1:n_perm
    swap = rand(n, 1) < 0.5;
    perm_a = a;
    perm_b = b;
    perm_a(swap) = b(swap);
    perm_b(swap) = a(swap);
    if two_sample_ks(perm_a, perm_b) + eps(D) >= D
        extreme = extreme + 1;
    end
end
rng(previous);
p = (extreme + 1) / (n_perm + 1);
end


function D = two_sample_ks(a, b)
sa = sort(a(:));
sb = sort(b(:));
n = numel(sa);
m = numel(sb);
if n == 0 || m == 0
    D = nan;
    return
end
i = 1;
j = 1;
D = 0;
while i <= n || j <= m
    if j > m || (i <= n && sa(i) <= sb(j))
        x = sa(i);
    else
        x = sb(j);
    end
    while i <= n && sa(i) == x
        i = i + 1;
    end
    while j <= m && sb(j) == x
        j = j + 1;
    end
    D = max(D, abs((i - 1) / n - (j - 1) / m));
end
end


function draw_stratum_rate_figure(low, high, folder)
open_full_figure('Power stratum firing rates')
low = low(:);
high = high(:);
ratio = high ./ low;
tests = paired_rate_tests(high, low);
subplot(1, 2, 1)
hold on
span = [min([low; high]), max([low; high])];
if ~(span(2) > span(1))
    span = span(1) * [0.5 2];
end
plot(span, span, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1)
plot(low, high, 'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.25 0.35 0.55], ...
    'MarkerEdgeColor', 'none')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlim(span)
ylim(span)
xlabel('Low-power firing rate (spikes/s)')
ylabel('High-power firing rate (spikes/s)')
title(sprintf('Median high/low = %.3g', median(ratio)), 'Interpreter', 'none')
box off
subplot(1, 2, 2)
histogram(log2(ratio), 'FaceColor', [0.25 0.35 0.55], 'EdgeColor', 'none')
hold on
yl = ylim;
plot([0 0], yl, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1)
xlabel('log2(high-power rate / low-power rate)')
ylabel('Neurons')
title({sprintf('Paired t p %.3g. Signed-rank p %.3g.', tests.t.raw_p, tests.signrank_p), ...
    sprintf('Sign test p %.3g. Permutation KS p %.3g.', tests.sign_p, tests.ks_p)}, ...
    'Interpreter', 'none')
box off
sgtitle({'Each point is one neuron. The dashed line is equal firing rate.', ...
    'The four p-values test whether the low-power and high-power rates match.'}, ...
    'Interpreter', 'none')
save_script_figure(folder)
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
%LOAD_EXAMPLE_ACG_BANK  Mirrored high-power and low-power autocorrelograms.

bank = empty_acg_bank('The locking summary file was not found, so the original autocorrelograms are not drawn.');
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
    low_name = '';
    if has_neuron_curve(stored.power_acg, 'low')
        low_name = 'low';
    end
    note = ['High-power and low-power sum-to-1 autocorrelograms, positive lags, both sides averaged.'];
    if isempty(low_name)
        note = [note, ' The summary has no low-power curve.'];
    end
    bank = bank_from_curves(stored.power_acg.all_neurons, stored.power_acg.centers, ...
        'high', low_name, note);
    return
end
try
    stored = load(summary_path, 'locking_acg');
catch
    bank.note = 'The locking summary file could not be read.';
    return
end
if isfield(stored, 'locking_acg') && has_neuron_curve(stored.locking_acg, 'all')
    bank = bank_from_curves(stored.locking_acg.all_neurons, stored.locking_acg.centers, ...
        'all', '', ...
        ['The power-stratum curves were not in the summary. ', ...
        'The high-power row uses the all-runs autocorrelogram, positive lags, both sides averaged.']);
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
bank.high = zeros(0, 0);
bank.low = zeros(0, 0);
bank.session = {};
bank.cluster_id = zeros(0, 1);
bank.note = note;
end


function bank = bank_from_curves(block, centers, high_name, low_name, note)
bank = empty_acg_bank(note);
centers = double(centers(:).');
high_curves = double(block.(high_name));
if size(high_curves, 2) ~= numel(centers)
    bank.note = 'The autocorrelogram width does not match its lag axis.';
    return
end
[bank.high, bank.lag] = mirror_acg_matrix(high_curves, centers);
bank.low = nan(size(bank.high));
if ~isempty(low_name) && isfield(block, low_name)
    low_curves = double(block.(low_name));
    if size(low_curves, 2) == numel(centers) && size(low_curves, 1) == size(high_curves, 1)
        [bank.low, ~] = mirror_acg_matrix(low_curves, centers);
    end
end
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


function [acg_lag, acg_high, acg_low, found_high, found_low] = example_stratum_curves(bank, session, cluster_id)
acg_lag = bank.lag;
acg_high = nan(1, numel(bank.lag));
acg_low = acg_high;
found_high = false;
found_low = false;
row = find_example_neuron(bank, session, cluster_id);
if isempty(row) || row > size(bank.high, 1)
    return
end
acg_high = bank.high(row, :);
found_high = any(isfinite(acg_high));
if row <= size(bank.low, 1)
    acg_low = bank.low(row, :);
    found_low = any(isfinite(acg_low));
end
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


function y = smooth_example_trace(x, y, width_ms)
%SMOOTH_EXAMPLE_TRACE  Gaussian smooth. width_ms is the standard deviation.

y = double(y(:)).';
x = double(x(:));
n = min(numel(x), numel(y));
x = x(1:n);
y = y(1:n);
if n < 3 || ~isfinite(width_ms) || width_ms <= 0
    return
end
bin_s = median(diff(x));
if ~isfinite(bin_s) || bin_s <= 0
    return
end
sigma_bins = (width_ms / 1000) / bin_s;
if sigma_bins < 0.05
    return
end
half_width = max(1, ceil(4 * sigma_bins));
t = -half_width:half_width;
kernel = exp(-0.5 * (t / sigma_bins) .^ 2);
kernel = kernel / sum(kernel);
finite = isfinite(y);
filled = y;
filled(~finite) = 0;
numer = conv(filled, kernel, 'same');
denom = conv(double(finite), kernel, 'same');
out = y;
use = denom > 0;
out(use) = numer(use) ./ denom(use);
out(~finite) = nan;
y = out;
end


function y_lim = example_acg_ylim(acg_lag, acg_y, lag_limits)
%EXAMPLE_ACG_YLIM  Scale one or more autocorrelograms without the near-zero peak.
%
% acg_y has one curve per column. Every column uses acg_lag.

acg_lag = double(acg_lag(:));
acg_y = example_acg_columns(acg_lag, acg_y);
in_window = acg_lag >= lag_limits(1) & acg_lag <= lag_limits(2);
tail = in_window & acg_lag >= 0.05;
values = acg_y(tail, :);
values = values(isfinite(values(:)));
if numel(values) < 5
    shown = acg_y(in_window, :);
    values = shown(isfinite(shown(:)));
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


function curves = example_acg_columns(acg_lag, acg_y)
%EXAMPLE_ACG_COLUMNS  Curves as columns aligned to acg_lag.

n = numel(acg_lag);
acg_y = double(acg_y);
if isvector(acg_y)
    flat = acg_y(:);
    if numel(flat) == n
        curves = flat;
        return
    end
    if n > 0 && mod(numel(flat), n) == 0
        curves = reshape(flat, n, []);
        return
    end
    curves = nan(n, 1);
    n_copy = min(n, numel(flat));
    curves(1:n_copy) = flat(1:n_copy);
    return
end
if size(acg_y, 1) == n
    curves = acg_y;
    return
end
if size(acg_y, 2) == n
    curves = acg_y.';
    return
end
curves = nan(n, size(acg_y, 2));
n_copy = min(n, size(acg_y, 1));
curves(1:n_copy, :) = acg_y(1:n_copy, :);
end


function text = example_neuron_summary(index, session, cluster_id, area_label, ...
    lock_label, modulation, found_high, found_low, feature_table, feature_table_low, row, names)
header = sprintf('Neuron %d   session %s   id %g   area %s   lock %s   modulation %.4g', ...
    index, session, cluster_id, area_label, lock_label, modulation);
lines = {header};
if ~found_high
    lines{end + 1} = 'High-power autocorrelogram was not matched in the summary file.'; %#ok<AGROW>
end
if ~found_low
    lines{end + 1} = 'Low-power autocorrelogram was not matched in the summary file.'; %#ok<AGROW>
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


function [cycle_n, cycle_r, cycle_p] = draw_power_logit_pair(logit, rows, columns, titles, ...
    group_name, n_group, n_population, lag, color_limit, color_map, ...
    first_peak_lag, first_trough_lag, second_peak_lag, second_trough_lag, ...
    scale_name, order_note, folder)
%DRAW_POWER_LOGIT_PAIR  High-power and low-power logit image above the mean.

n_cols = numel(columns);
n_rows = 2;
cycle_n = nan(1, n_cols);
cycle_r = nan(1, n_cols);
cycle_p = nan(1, n_cols);
y_limit = color_limit;
for c = 1:n_cols
    column_mean = mean(logit.sum.(columns{c})(rows, :), 1, 'omitnan');
    column_peak = max(abs(column_mean(isfinite(column_mean))));
    if isfinite(column_peak)
        y_limit = max(y_limit, column_peak);
    end
end
y_limit = y_limit * 1.08;
open_full_figure(['ACG population logit, ', group_name, ', ', scale_name])
for c = 1:n_cols
    this_logit = logit.sum.(columns{c})(rows, :);
    n = size(this_logit, 1);
    subplot(n_rows, n_cols, c)
    if n < 1
        title({titles{c}, sprintf('%s/n = 0/%d', group_name, n_population)})
        axis off
    elseif n == 1
        imagesc(lag, [0.5 1.5], this_logit)
        ylim([0.5 1.5])
        axis tight
        clim([-color_limit color_limit])
        colormap(gca, color_map)
        ylabel('Neuron')
        title({titles{c}, sprintf('%s/n = %d/%d', group_name, n_group, n_population)})
        box off
    else
        imagesc(lag, 1:n, this_logit)
        set(gca, 'YTick', [1 n])
        axis tight
        clim([-color_limit color_limit])
        colormap(gca, color_map)
        ylabel('Neuron')
        title({titles{c}, sprintf('%s/n = %d/%d', group_name, n_group, n_population)})
        box off
    end
    if c == n_cols && n >= 1
        colorbar
    end
    if n >= 1
        first_diff = mean(this_logit(:, first_peak_lag), 2, 'omitnan') ...
            - mean(this_logit(:, first_trough_lag), 2, 'omitnan');
        second_diff = mean(this_logit(:, second_peak_lag), 2, 'omitnan') ...
            - mean(this_logit(:, second_trough_lag), 2, 'omitnan');
        finite_pair = isfinite(first_diff) & isfinite(second_diff);
        cycle_n(c) = sum(finite_pair);
        if sum(finite_pair) >= 3
            [r_mat, p_mat] = corrcoef(first_diff(finite_pair), second_diff(finite_pair));
            cycle_r(c) = r_mat(1, 2);
            cycle_p(c) = p_mat(1, 2);
        end
    end
    subplot(n_rows, n_cols, n_cols + c)
    hold on
    plot(lag([1 end]), [0 0], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
    if n >= 1
        plot(lag, mean(this_logit, 1, 'omitnan'), 'k', 'LineWidth', 2.4)
    end
    xlabel('Lag (s)')
    ylabel('Mean logit')
    ylim([-y_limit y_limit])
    box off
end
figure_group = group_name;
figure_group(1) = upper(figure_group(1));
sgtitle({[figure_group, ' neurons, ', scale_name], order_note}, 'Interpreter', 'none')
save_script_figure(folder)
end


function lines = cycle_summary_lines(titles, cycle_n, cycle_r, cycle_p)
lines = cell(1, numel(titles));
for c = 1:numel(titles)
    if ~isfinite(cycle_r(c))
        lines{c} = sprintf('%-32s n = %.0f    too few neurons for r', ...
            titles{c}, cycle_n(c));
    else
        lines{c} = sprintf('%-32s n = %.0f    r = %.4g    uncorrected p = %.4g', ...
            titles{c}, cycle_n(c), cycle_r(c), cycle_p(c));
    end
end
end


function write_section_4_summary(folder, column_title, cycle_n, cycle_r, cycle_p, ...
    lag, high_vs_low_h, high_vs_low_p, n_osc, ...
    n_osc_high, n_other_high, n_osc_low_pca, n_osc_low_tsne, n_omitted, ...
    embed_perplexity, paired_stats, lock_type_label, type_observed, type_total, type_percent)

cycle_lines = cycle_summary_lines(column_title, cycle_n, cycle_r, cycle_p);
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
    'The logit figure is two columns, high power and low power, and two rows, the neuron image and the mean logit.', ...
    'Pearson correlation of the first cycle (peak minus trough) with the second cycle is reported here and is not drawn.', ...
    'The p-values are uncorrected and are not in the Holm-Sidak family.', ...
    cycle_lines))
write_summary_section(folder, '4. High power minus low power across lag', join_summary_lines( ...
    'Each lag is a paired t-test of the same oscillatory neurons, high power minus low power.', ...
    'Red points on the figure mark uncorrected p < 0.05.', ...
    'These lag tests are not part of the Holm-Sidak correction used for the feature panels.', ...
    lag_line))
write_summary_section(folder, '4. Low-power points in the PCA and t-SNE', join_summary_lines( ...
    'One PCA and one t-SNE are fit on the high-power neurons plus the oscillatory neurons at low power.', ...
    'The oscillatory and non-oscillatory figures, and the lock-type figures, use that shared plane.', ...
    'One PCA figure and one t-SNE figure draw only the high-power neurons. Another pair draws those same points and adds the low-power oscillatory neurons as triangles.', ...
    'Orange is oscillatory at high power. Purple is non-oscillatory at high power. Blue triangles are oscillatory at low power.', ...
    'Green is reserved for unlocked neurons on the lock-type figures.', ...
    sprintf('Oscillatory high power: %d. Other high-power neurons: %d. Oscillatory low power included: %d. Omitted for missing low-power features: %d.', ...
        n_osc_high, n_other_high, n_osc_low_pca, n_omitted), ...
    sprintf('Shared t-SNE perplexity: %g. Low-power rows in that t-SNE: %d.', ...
        embed_perplexity, n_osc_low_tsne)))
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
