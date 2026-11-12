% Plot_acg_population_logit
% All neurons together, on a logit scale, from a saved ACG shuffle.
%
% Overlap is the standard method. Section 1 loads
% locking_period_acg_shuffle_n500.mat and builds the matrices. Section 2
% chooses the oscillatory neurons. Section 3 only draws. After section 2
% has run, edit section 3 and run that section alone.
%
% Every matrix shares one row order. Row i of imagesc_logit.sum.all is
% the same neuron as row i of imagesc_logit.sum.power_high, power_low,
% and the rate-ratio matrices. neuron_type{i} is that neuron's lock group:
% 'trough', 'peak', 'nonentrained' (unlocked), or '' when it was not
% classified. sort_score.* is the mean percentile from 200 to 400 ms.
% Row 1 is the highest all-runs sum-to-1 score.
%
% Section 2 keeps a neuron when, on the high 25% power sum-to-1 logit,
% the mean in the first peak window exceeds the mean in the first trough
% window by more than osc_difference_threshold. shown_logit is that
% subset, strongest difference first. imagesc_osc_sum_all is the same
% sum-to-1 image under a plain name. osc_neuron_type is its lock group.

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
save(prep_file, 'lag', 'imagesc_logit', 'percentile', 'sort_score', ...
    'imagesc_sum_all', 'imagesc_sum_power_high', 'imagesc_sum_power_low', ...
    'neuron_type', 'neuron_record_index', ...
    'column_id', 'column_title', 'scale_id', 'scale_title', ...
    'n_shuffles', 'sort_window_s', 'checkpoint', 'prep_file');
disp(['Section 1 saved the matrices to ', prep_file])
disp('imagesc_sum_all is neurons by lag. neuron_type is the lock group of each row.')
disp('Run section 2 to choose the oscillatory neurons.')

%% 2. Choose the oscillatory neurons
% Re-run this section after changing the windows or the threshold.
% Section 1 does not need to run again while its variables are loaded.
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
if ~exist('imagesc_logit', 'var') || ~exist('neuron_type', 'var')
    data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\' ...
        'NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data'];
    parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
    prep_file = fullfile(parsing_folder, 'acg_population_logit_prep.mat');
    load(prep_file)
end
if ~exist('neuron_type', 'var')
    error('Plot_acg_population_logit:neuronType', ...
        'The prep file has no neuron type. Run section 1 again.');
end

% Open intervals, in seconds. Row 1 is the first cycle. Row 2 is the second.
trough_windows_s = [0.05 0.15; 0.35 0.45];
peak_windows_s = [0.25 0.35; 0.45 0.55];
% High-power sum-to-1 logit: mean(first peak) - mean(first trough) must
% exceed this value. The same neurons are then used for every column.
osc_difference_threshold = 0;

[osc_rows, osc_score] = oscillatory_row_index( ...
    imagesc_logit.sum.power_high, lag, peak_windows_s(1, :), trough_windows_s(1, :), ...
    osc_difference_threshold);
osc_neuron_type = neuron_type(osc_rows);
[lock_type_count, lock_type_proportion, lock_type_label] = lock_type_proportion(osc_neuron_type);
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
% High-power minus low-power logit for the oscillatory neurons.
high_minus_low = shown_logit.sum.power_high - shown_logit.sum.power_low;
% Plain names for the oscillatory sum-to-1 images. Row 1 is the strongest
% first-cycle peak minus trough. The same rows are neuron_type(osc_rows).
imagesc_osc_sum_all = shown_logit.sum.all;
imagesc_osc_sum_power_high = shown_logit.sum.power_high;
imagesc_osc_sum_power_low = shown_logit.sum.power_low;

disp(sprintf('Oscillatory neurons: %d of %d. shown_logit uses that row order.', ...
    n_osc, n_before_osc))
disp('Run section 3 to draw.')

%% 3. Figures
% Re-run this section after editing the figure. Section 2 does not need
% to run again while shown_logit is still in the workspace.
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
if ~exist('shown_logit', 'var')
    error('Plot_acg_population_logit:shown', ...
        'shown_logit is not in the workspace. Run section 2 first.');
end

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

for s = 1:numel(scale_id)
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
        plot(first_diff, second_diff, '.k')
        xlabel('First cycle, peak - trough')
        ylabel('Second cycle, peak - trough')
        finite_pair = isfinite(first_diff) & isfinite(second_diff);
        if sum(finite_pair) < 3
            title('Too few neurons')
        else
            [r_mat, p_mat] = corrcoef(first_diff(finite_pair), second_diff(finite_pair));
            title(sprintf('r = %.2f, p = %.3g', r_mat(1, 2), p_mat(1, 2)))
        end
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
        'Sorted by high-power first peak minus first trough, shown as logit.'})
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

type_color = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.45 0.45 0.45; 0.75 0.75 0.75];
figure('Name', 'Oscillatory neuron lock types', ...
    'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.28 0.28 0.38 0.42])
type_percent = 100 * lock_type_proportion;
if n_osc < 1
    type_percent = zeros(1, 4);
end
type_bar = bar(type_percent);
type_bar.FaceColor = 'flat';
type_bar.CData = type_color;
set(gca, 'XTick', 1:4, 'XTickLabel', lock_type_label)
ylabel('Oscillatory neurons (%)')
ylim([0 100])
if n_osc < 1
    title('No oscillatory neurons')
else
    title(sprintf('Lock type among %d oscillatory neurons', n_osc))
    for k = 1:4
        text(k, min(type_percent(k) + 3, 96), sprintf('n = %d', lock_type_count(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    end
end
box off
