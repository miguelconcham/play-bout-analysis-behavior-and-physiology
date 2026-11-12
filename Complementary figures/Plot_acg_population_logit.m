% Plot_acg_population_logit
% All neurons together, on a logit scale, from a saved ACG shuffle.
%
% Overlap is the standard method. Section 1 loads
% locking_period_acg_shuffle_n500.mat and builds the matrices and
% feature_table. Section 2 runs PCA, t-SNE, and the two-group split.
% Section 3 keeps k-means cluster 1, the higher-modulation group.
% Section 4 only draws. After section 3 has run, edit section 4 and run
% that section alone.
%
% Every matrix shares one row order. Row i of imagesc_logit.sum.all is
% the same neuron as row i of imagesc_logit.sum.power_high, power_low,
% and the rate-ratio matrices. neuron_type{i} is that neuron's lock group:
% 'trough', 'peak', 'nonentrained' (unlocked), or '' when it was not
% classified. sort_score.* is the mean percentile from 200 to 400 ms.
% Row 1 is the highest all-runs sum-to-1 score.
%
% Section 2 builds ACG features from the high-power sum-to-1 logit and
% splits the neurons with k-means. Cluster 1 is the group with the higher
% modulation. shown_logit is that group, strongest modulation first.
% imagesc_osc_sum_all is the same sum-to-1 image under a plain name.
% osc_neuron_type is its lock group. feature_table has one row per neuron.

%% 1. Load the checkpoint and build the imagesc matrices
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
    'neuron_type', 'neuron_record_index', 'feature_table', ...
    'column_id', 'column_title', 'scale_id', 'scale_title', ...
    'n_shuffles', 'sort_window_s', 'checkpoint', 'prep_file');
disp(['Section 1 saved the matrices to ', prep_file])
disp('imagesc_sum_all is neurons by lag. neuron_type is the lock group of each row.')
disp('feature_table is one row of ACG features per neuron. Run section 2 to classify.')

%% 2. ACG features, PCA, t-SNE, and the two-group split
% Re-run this section after changing the feature list. Section 1 does not
% need to run again while feature_table is still in the workspace.
% Cluster 1 is the k-means group with the higher modulation. Section 3
% uses that group. The PCA and t-SNE figures are there to see whether the
% same split is visible. They do not choose the neurons. The boxcharts
% after the split also show frequency, amplitude, and the peak-to-peak
% drop. Those extra columns do not enter the split.
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
figure('Name', 'ACG feature correlations', 'Color', 'w')
heatmap(feature_names, feature_names, feature_correlation, 'ColorLimits', [-1 1])
title('Correlation between ACG features')

feature_mean = mean(X, 1);
feature_sigma = std(X, 0, 1);
Xz = (X - feature_mean) ./ feature_sigma;

% pca_mu is the centering of this PCA. Section 4 uses it to place the
% low-power features in the same plane. t-SNE has no equivalent transform.
[pca_coeff, pca_score, pca_latent, ~, pca_explained, pca_mu] = pca(Xz);
figure('Name', 'PCA variance explained', 'Color', 'w')
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

figure('Name', 'PCA population structure', 'Color', 'w')
scatter(pca_score(:, 1), pca_score(:, 2), 35, 'filled', 'MarkerFaceAlpha', 0.7)
xlabel(sprintf('PC1 (%.1f%%)', pca_explained(1)))
ylabel(sprintf('PC2 (%.1f%%)', pca_explained(2)))
title('ACG features: PCA')
grid on
axis square

figure('Name', 'PCA feature loadings', 'Color', 'w')
biplot(pca_coeff(:, 1:2), 'Scores', pca_score(:, 1:2), 'VarLabels', feature_names)
title('PCA loadings and neuron scores')

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

figure('Name', 'ACG t-SNE', 'Color', 'w')
for k = 1:numel(perplexity_values)
    subplot(1, numel(perplexity_values), k)
    scatter(Y_tsne{k}(:, 1), Y_tsne{k}(:, 2), 30, 'filled', 'MarkerFaceAlpha', 0.7)
    xlabel('t-SNE 1')
    ylabel('t-SNE 2')
    title(sprintf('Perplexity = %d', perplexity_values(k)))
    axis square
    grid on
end

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

figure('Name', 'Comparison of exploratory partitions', 'Color', 'w')
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

figure('Name', 'Silhouette comparison', 'Color', 'w')
subplot(1, 2, 1)
silhouette(Xz, idx_kmeans, 'euclidean')
title('K-means')
subplot(1, 2, 2)
silhouette(Xz, idx_ward, 'euclidean')
title('Ward')

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
    figure('Name', 'Feature distributions by k-means group', 'Color', 'w')
    n_descriptor_cols = 5;
    tiledlayout(ceil(numel(shown_titles) / n_descriptor_cols), n_descriptor_cols, ...
        'TileSpacing', 'compact')
    for j = 1:numel(shown_titles)
        nexttile
        boxchart(categorical(idx_kmeans), descriptor_values(:, j))
        title({shown_titles{j}, shown_roles{j}}, 'Interpreter', 'none')
        xlabel('K-means group')
        grid on
    end
    sgtitle({'Cluster 1 is the higher modulation.', ...
        'Peak1 - Peak2 is the drop from the first peak to the second.'})
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
    'kmeans_centers', kmeans_centers, ...
    'modulation_score', modulation_score, ...
    'valid_rows', valid_rows);
disp('Run section 3 to keep k-means cluster 1.')

%% 3. Choose the oscillatory neurons
% Cluster 1 is the higher-modulation k-means group. Its rows are the
% complete feature rows, not yet the image rows. The strongest modulation
% is placed first.
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
if ~exist('idx_kmeans', 'var') || ~exist('complete_rows', 'var') ...
        || ~exist('modulation_score', 'var') || ~exist('neuron_type', 'var')
    error('Plot_acg_population_logit:cluster', ...
        'The k-means split is not in the workspace. Run section 2 first.');
end

osc_score = nan(numel(neuron_type), 1);
osc_score(complete_rows) = modulation_score;
osc_member = false(numel(neuron_type), 1);
osc_member(complete_rows(idx_kmeans == 1)) = true;
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

%% 4. Figures
% Re-run this section after editing the figure. Section 3 does not need
% to run again while shown_logit is still in the workspace. The low-power
% PCA figure and the paired feature figure also need the section 2 PCA.
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

for s = 1
    figure('Name', ['ACG population logit, ', scale_title{s}], ...
        'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.06 0.08 0.88 0.84])
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
        if sum(finite_pair) < 3
            plot(first_diff, second_diff, '.k')
            title('Too few neurons')
        else
            cycle_fit = fitlm(first_diff(finite_pair), second_diff(finite_pair));
            plot(cycle_fit)
            legend off
            [r_mat, p_mat] = corrcoef(first_diff(finite_pair), second_diff(finite_pair));
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
        'K-means group with the higher modulation, strongest first.'})
end

figure('Name', 'ACG population logit, high minus low power', ...
    'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.18 0.16 0.52 0.68])
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

% t-SNE cannot place a new point in an existing map. The low-power
% features are projected into the PCA built from the high-power features.
if ~exist('pca_coeff', 'var') || ~exist('pca_mu', 'var') ...
        || ~exist('pca_score', 'var') || ~exist('pca_explained', 'var') ...
        || ~exist('feature_mean', 'var') || ~exist('feature_sigma', 'var') ...
        || ~exist('feature_names', 'var') || ~exist('complete_rows', 'var') ...
        || ~exist('feature_table', 'var') || ~exist('imagesc_sum_power_low', 'var') ...
        || ~exist('osc_rows', 'var')
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

figure('Name', 'Low-power ACGs in the high-power PCA', ...
    'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.16 0.16 0.46 0.62])
hold on
high_cloud = scatter(pca_score(:, 1), pca_score(:, 2), 18, [0.75 0.75 0.75], ...
    'filled', 'MarkerFaceAlpha', 0.35);
paired_line = gobjects(0);
if any(finite_projection)
    paired_lines = plot( ...
        [high_xy(finite_projection, 1).'; low_xy(finite_projection, 1).'], ...
        [high_xy(finite_projection, 2).'; low_xy(finite_projection, 2).'], ...
        '-', 'Color', [0.55 0.55 0.55]);
    paired_line = paired_lines(1);
    if numel(paired_lines) > 1
        set(paired_lines(2:end), 'HandleVisibility', 'off')
    end
end
high_points = scatter(high_xy(finite_projection, 1), high_xy(finite_projection, 2), ...
    36, [0.85 0.25 0.20], 'filled');
low_points = scatter(low_xy(finite_projection, 1), low_xy(finite_projection, 2), ...
    36, [0.20 0.45 0.80], 'filled');
xlabel(sprintf('PC1 (%.1f%%)', pca_explained(1)))
ylabel(sprintf('PC2 (%.1f%%)', pca_explained(2)))
title({'Oscillatory neurons in the high-power PCA', ...
    'Red is high power. Blue is the same neuron at low power.'})
if isempty(paired_line)
    legend([high_cloud, high_points, low_points], ...
        {'All neurons, high power', 'Oscillatory, high power', ...
        'Oscillatory, low power'}, 'Location', 'best', 'Box', 'off')
else
    legend([high_cloud, paired_line, high_points, low_points], ...
        {'All neurons, high power', 'Same neuron', 'Oscillatory, high power', ...
        'Oscillatory, low power'}, 'Location', 'best', 'Box', 'off')
end
grid on
axis square
box off

% Same columns as the cluster boxcharts. High values come from the
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
    figure('Name', 'High versus low power features, oscillatory neurons', ...
        'Color', 'w')
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
            boxchart(categorical( ...
                [ones(sum(pair), 1); 2 * ones(sum(pair), 1)], ...
                [1 2], {'High power', 'Low power'}), ...
                [high_vals(pair); low_vals(pair)])
            plot([1 2], [high_vals(pair), low_vals(pair)].', ...
                '-', 'Color', [0.65 0.65 0.65])
        end
        if sum(pair) < 2
            title({compare_titles{j}, 'Too few neurons'}, 'Interpreter', 'none')
        else
            [~, feature_p] = ttest(high_vals(pair), low_vals(pair));
            title({compare_titles{j}, sprintf('uncorrected p = %.3g', feature_p)}, ...
                'Interpreter', 'none')
        end
        xlim([0.5 2.5])
        grid on
    end
    sgtitle({'Same oscillatory neurons at high and low power.', ...
        'Each p is a paired t-test and is not corrected across features.'})
end

type_color = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.45 0.45 0.45; 0.75 0.75 0.75];
figure('Name', 'Oscillatory cells within each lock type', ...
    'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.22 0.22 0.48 0.46])
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
