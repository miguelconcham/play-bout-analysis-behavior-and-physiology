function Plot_locking_period_acg_shuffle(checkpoint)
%PLOT_LOCKING_PERIOD_ACG_SHUFFLE  Figures for a saved ACG shuffle checkpoint.
%
% This is the plotting half of Analyze_locking_period_acg_shuffle, taken
% from the eligible-spike version of that script. Overlap stays the
% standard method. With no argument this opens
% locking_period_acg_shuffle_n500.mat. Change acg_method below to
% 'eligible' to open the eligible-spike checkpoint instead, or pass a
% checkpoint path.
%
% The neuron percentile figure has four rows in each cell-type column.
% The top row is an imagesc with one autocorrelogram percentile per row.
% Under it are the same curves, with their mean drawn bold. The third row
% is the percent of neurons at or above the 95th percentile, drawn on top
% of the percent at or below the 5th. The bottom row is the mean across
% neurons of log10(p/(1-p)), after percentiles of 0 and 100 are pulled
% inside (0, 1).

if nargin < 1
    checkpoint = '';
end
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
addpath(this_dir);
addpath(repo_root);
if exist('acg_delta_peaks', 'file') ~= 2
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

acg_method = 'overlap';
n_shuffles = 500;
if isempty(checkpoint)
    data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\' ...
        'NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data'];
    parsing_folder = fullfile(data_root, 'Analysis results', 'spike train parsing');
    checkpoint = fullfile(parsing_folder, shuffle_checkpoint_filename(n_shuffles, acg_method));
end
if exist(checkpoint, 'file') ~= 2
    error('Plot_locking_period_acg_shuffle:checkpoint', ...
        'No shuffle checkpoint at %s.', checkpoint);
end

stored = load(checkpoint, 'shuffle_records', 'parameters', 'grand_average');
if ~isfield(stored, 'shuffle_records') || ~isfield(stored, 'parameters') ...
        || ~isfield(stored, 'grand_average')
    error('Plot_locking_period_acg_shuffle:checkpoint', ...
        '%s does not contain shuffle_records, parameters, and grand_average.', checkpoint);
end
records = stored.shuffle_records;
if ~iscell(records)
    records = num2cell(records(:));
end
records = records(:);
parameters = stored.parameters;
required = {'n_shuffles', 'smooth_n', 'acg_lag', 'bin_w', 'n_bins'};
for i = 1:numel(required)
    if ~isfield(parameters, required{i}) || isempty(parameters.(required{i}))
        error('Plot_locking_period_acg_shuffle:checkpoint', ...
            '%s has no parameters.%s.', checkpoint, required{i});
    end
end
grand = stored.grand_average;
method = checkpoint_acg_method(parameters);
centers = centers_from_parameters(parameters);
disp(['Plotting ', method, ' shuffle checkpoint ', checkpoint])

group_order = {'trough', 'peak', 'nonentrained'};
group_title = cellfun(@lock_group_title, group_order, 'UniformOutput', false);
if isfield(grand, 'group_order') && ~isempty(grand.group_order)
    grand_order = grand.group_order;
    if isstring(grand_order)
        grand_order = cellstr(grand_order);
    end
else
    grand_order = {'trough', 'peak', 'delta', 'nonentrained'};
end
grand_title = cellfun(@lock_group_title, grand_order, 'UniformOutput', false);

plot_shuffle_exception(records, group_order, group_title, 'locking', ...
    'Mean PPC of the runs');
plot_shuffle_exception(records, group_order, group_title, 'power', ...
    'Mean delta power of the runs');
plot_all_run_exception(records, group_order, group_title);
report_grand_shuffle(grand, centers, parameters.smooth_n, grand_order, grand_title, method);
plot_percentile_curves(records, grand, centers, parameters.smooth_n, ...
    grand_order, grand_title, parameters.n_shuffles);
end


function method = checkpoint_acg_method(parameters)
% Checkpoints saved before the method selector are overlap.
method = 'overlap';
if isstruct(parameters) && isfield(parameters, 'acg_method') && ~isempty(parameters.acg_method)
    method = lower(strtrim(char(parameters.acg_method)));
end
if ~any(strcmp(method, {'overlap', 'eligible'}))
    error('Plot_locking_period_acg_shuffle:method', ...
        'acg_method must be ''overlap'' or ''eligible''.');
end
end


function centers = centers_from_parameters(parameters)
edges = (-parameters.acg_lag):parameters.bin_w:parameters.acg_lag;
centers = edges(1:end - 1) + parameters.bin_w / 2;
if numel(centers) ~= parameters.n_bins
    error('Plot_locking_period_acg_shuffle:axis', ...
        'Checkpoint bin count %d does not match the lag axis (%d).', ...
        parameters.n_bins, numel(centers));
end
end


function title_text = lock_group_title(name)
switch char(name)
    case 'trough'
        title_text = 'Trough-locked';
    case 'peak'
        title_text = 'Peak-locked';
    case 'delta'
        title_text = 'Delta-locked';
    case 'nonentrained'
        title_text = 'Non-entrained';
    otherwise
        title_text = char(name);
end
end


function plot_shuffle_exception(records, group_order, group_title, split_name, score_label)
if isempty(records)
    return
end
if strcmp(split_name, 'locking')
    split_label = 'locking';
else
    split_label = 'delta power';
end
measures = {'first_p', 'second_p'};
measure_labels = {'First-peak exceptionality', 'Second-cycle exceptionality'};
colors = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.45 0.45 0.45];
figure('Name', ['ACG peak shuffle, ', split_label], 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.08 0.08 0.84 0.8])
for m = 1:numel(measures)
    [p_low, p_high, ok, x_low, x_high, groups] = paired_shuffle_values( ...
        records, split_name, measures{m});
    y_low = -log10(p_low);
    y_high = -log10(p_high);
    subplot(2, 2, 2 * m - 1)
    plot_paired_exception(groups(ok), y_low(ok), y_high(ok), ...
        group_order, colors, measure_labels{m})
    subplot(2, 2, 2 * m)
    plot_exception_against_score(groups(ok), x_low(ok), x_high(ok), ...
        y_low(ok), y_high(ok), group_order, colors, score_label, measure_labels{m})
    report_exception_fraction(groups, p_low, p_high, ok, group_order, group_title, ...
        split_label, measure_labels{m});
end
hold on
handles = gobjects(1, numel(group_order));
for g = 1:numel(group_order)
    handles(g) = plot(nan, nan, '-', 'Color', colors(g, :), 'LineWidth', 2.5);
end
legend(handles, group_title, 'Location', 'eastoutside')
end


function plot_all_run_exception(records, group_order, group_title)
if isempty(records)
    return
end
measures = {'first_p', 'second_p'};
measure_labels = {'First-peak exceptionality', 'Second-cycle exceptionality'};
colors = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.45 0.45 0.45];
n = numel(records);
groups = cell(n, 1);
for i = 1:n
    groups{i} = records{i}.lock_group;
end
figure('Name', 'ACG peak shuffle, all firing runs', 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.12 0.15 0.7 0.55])
for m = 1:numel(measures)
    p = nan(n, 1);
    for i = 1:n
        if ~isfield(records{i}, 'all')
            continue
        end
        p(i) = records{i}.all.(measures{m});
        if strcmp(measures{m}, 'second_p') && ~records{i}.all.second_measurable
            p(i) = NaN;
        end
    end
    ok = isfinite(p) & p > 0;
    subplot(1, 2, m)
    hold on
    plot([0.5, numel(group_order) + 0.5], [-log10(0.05), -log10(0.05)], '--', ...
        'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
    for g = 1:numel(group_order)
        idx = ok & strcmp(groups, group_order{g});
        if ~any(idx)
            continue
        end
        plot(g * ones(sum(idx), 1), -log10(p(idx)), 'o', 'Color', colors(g, :), ...
            'MarkerFaceColor', colors(g, :), 'MarkerSize', 5, 'HandleVisibility', 'off')
        disp(sprintf('  %s, all firing runs, %s: p<0.05 = %d/%d', ...
            group_title{g}, measure_labels{m}, sum(p(idx) < 0.05), sum(idx)))
    end
    xlim([0.5, numel(group_order) + 0.5])
    set(gca, 'XTick', 1:numel(group_order), 'XTickLabel', group_title)
    ylabel([measure_labels{m}, ', -log10 p'])
    box off
end
sgtitle('All firing runs against shuffles of every firing run')
end


function [p_low, p_high, ok, x_low, x_high, groups] = paired_shuffle_values(records, split_name, measure_name)
n = numel(records);
p_low = nan(n, 1);
p_high = nan(n, 1);
x_low = nan(n, 1);
x_high = nan(n, 1);
groups = cell(n, 1);
for i = 1:n
    groups{i} = records{i}.lock_group;
    block = records{i}.(split_name);
    x_low(i) = block.score_low;
    x_high(i) = block.score_high;
    p_low(i) = block.low.(measure_name);
    p_high(i) = block.high.(measure_name);
    if strcmp(measure_name, 'second_p')
        if ~block.low.second_measurable || ~block.high.second_measurable
            p_low(i) = NaN;
            p_high(i) = NaN;
        end
    end
end
ok = isfinite(p_low) & isfinite(p_high) & p_low > 0 & p_high > 0;
end


function plot_paired_exception(lock_group, y_low, y_high, group_order, colors, y_label)
hold on
threshold = -log10(0.05);
plot([0.6 2.4], [threshold threshold], '--', 'Color', [0.45 0.45 0.45], ...
    'HandleVisibility', 'off')
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    plot([ones(sum(idx), 1), 2 * ones(sum(idx), 1)]', [y_low(idx), y_high(idx)]', '-', ...
        'Color', [colors(g, :), 0.25], 'HandleVisibility', 'off')
    plot([1 2], [mean(y_low(idx), 'omitnan'), mean(y_high(idx), 'omitnan')], '-', ...
        'Color', colors(g, :), 'LineWidth', 2.5, 'HandleVisibility', 'off')
end
xlim([0.6 2.4])
set(gca, 'XTick', [1 2], 'XTickLabel', {'Low 25%', 'High 25%'})
ylabel([y_label, ', -log10 p'])
title(sprintf('Paired, n = %d. Dashed line is p = 0.05', numel(y_low)))
box off
end


function plot_exception_against_score(lock_group, x_low, x_high, y_low, y_high, group_order, colors, x_label, y_label)
hold on
span = [x_low(:); x_high(:)];
span = span(isfinite(span));
if numel(span) >= 2
    plot([min(span) max(span)], [-log10(0.05) -log10(0.05)], '--', ...
        'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
end
for g = 1:numel(group_order)
    idx = strcmp(lock_group, group_order{g});
    if ~any(idx)
        continue
    end
    plot([x_low(idx), x_high(idx)]', [y_low(idx), y_high(idx)]', '-', ...
        'Color', [colors(g, :), 0.3], 'HandleVisibility', 'off')
end
xlabel(x_label)
ylabel([y_label, ', -log10 p'])
title('Each line is one neuron')
box off
end


function report_exception_fraction(groups, p_low, p_high, ok, group_order, group_title, split_label, measure_label)
for g = 1:numel(group_order)
    idx = ok & strcmp(groups, group_order{g});
    n = sum(idx);
    if n < 1
        continue
    end
    disp(sprintf('  %s, %s, %s: low p<0.05 = %d/%d, high p<0.05 = %d/%d', ...
        group_title{g}, split_label, measure_label, ...
        sum(p_low(idx) < 0.05), n, sum(p_high(idx) < 0.05), n))
end
end


function report_grand_shuffle(grand, centers, smooth_n, group_order, group_title, method)
scales = {'sum_to_one', 'rate_ratio'};
scale_labels = {'sum-to-1', 'rate ratio (Poisson = 1)'};
y_labels = {'Smoothed mean ACG', 'Ratio to Poisson rate'};
splits = {'locking', 'power'};
split_labels = {'locking', 'delta power'};
tails = {'low', 'high'};
tail_labels = {'Low 25%', 'High 25%'};
disp(' ')
disp('Grand-average shuffle. Peaks are scored on the Gaussian-smoothed mean autocorrelogram.')
disp('Sum-to-1 gives every neuron the same weight.')
if strcmp(method, 'eligible')
    disp('Rate ratio divides each eligible-spike curve by its own firing rate.')
else
    disp('Rate ratio divides each curve by its own firing rate squared.')
end
disp('The all-runs curve pools every firing run. Its null shuffles those runs, including the ones between the tails.')
for k = 1:numel(scales)
    scale = grand.(scales{k});
    for s = 1:numel(splits)
        figure('Name', ['Grand average shuffle, ', split_labels{s}, ', ', scale_labels{k}], ...
            'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.04 0.08 0.92 0.8])
        for t = 1:numel(tails)
            for g = 1:numel(group_order)
                block = scale.(splits{s}).(tails{t}).(group_order{g});
                subplot(numel(tails), numel(group_order), (t - 1) * numel(group_order) + g)
                result = grand_peak_result(block, centers, smooth_n);
                plot_grand_panel(result, [group_title{g}, ', ', tail_labels{t}], y_labels{k})
                disp(sprintf(['  %s, %s, %s, %s: n = %d, first prom = %.3g, shuffle p = %.3g, ' ...
                    'second prom = %.3g, shuffle p = %.3g'], ...
                    scale_labels{k}, group_title{g}, split_labels{s}, tail_labels{t}, result.n, ...
                    result.first_prominence, result.first_p, ...
                    result.second_prominence, result.second_p))
            end
        end
        sgtitle({['Grand mean versus shuffle, ', split_labels{s}, ', ', scale_labels{k}], ...
            'Black is the observed mean. Grey is the central 95% of shuffled means.'})
    end
    if ~isfield(scale, 'all')
        continue
    end
    figure('Name', ['Grand average shuffle, all firing runs, ', scale_labels{k}], ...
        'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.04 0.2 0.92 0.55])
    for g = 1:numel(group_order)
        block = scale.all.(group_order{g});
        subplot(1, numel(group_order), g)
        result = grand_peak_result(block, centers, smooth_n);
        plot_grand_panel(result, [group_title{g}, ', all firing runs'], y_labels{k})
        disp(sprintf(['  %s, %s, all firing runs: n = %d, first prom = %.3g, shuffle p = %.3g, ' ...
            'second prom = %.3g, shuffle p = %.3g'], ...
            scale_labels{k}, group_title{g}, result.n, ...
            result.first_prominence, result.first_p, ...
            result.second_prominence, result.second_p))
    end
    sgtitle({['Grand mean versus shuffle, all firing runs, ', scale_labels{k}], ...
        'Black is the observed mean. Grey is the central 95% of shuffled means.'})
end
end


function result = grand_peak_result(block, centers, smooth_n)
result.n = block.n_neurons;
result.first_prominence = NaN;
result.first_p = NaN;
result.second_prominence = NaN;
result.second_p = NaN;
result.first_lag = NaN;
result.second_lag = NaN;
result.fold_t = centers(centers > 0);
result.observed = nan(size(result.fold_t));
result.null_lo = result.observed;
result.null_hi = result.observed;
result.null_median = result.observed;
result.percentile = result.observed;
if block.n_neurons < 1
    return
end
mu = block.sum_obs ./ block.n_obs;
mu(block.n_obs == 0) = NaN;
obs = acg_delta_peaks(mu, centers, smooth_n);
result.fold_t = obs.fold_t;
result.observed = obs.fold_y;
result.first_lag = obs.first_lag;
result.second_lag = obs.second_lag;
result.first_prominence = obs.first_prominence;
result.second_prominence = obs.second_prominence;
n_s = size(block.sum_null, 1);
null_fold = nan(n_s, numel(obs.fold_t));
null_first = nan(n_s, 1);
null_second = nan(n_s, 1);
for s = 1:n_s
    sample = block.sum_null(s, :) ./ block.n_null(s, :);
    sample(block.n_null(s, :) == 0) = NaN;
    pk = acg_delta_peaks(sample, centers, smooth_n);
    null_first(s) = null_first_value(pk);
    null_second(s) = null_second_value(pk);
    if numel(pk.fold_y) == numel(obs.fold_t)
        null_fold(s, :) = pk.fold_y;
    end
end
result.first_p = exception_p(obs.first_prominence, null_first);
if obs.second_measurable
    result.second_p = exception_p(obs.second_prominence, null_second);
end
result.null_lo = prctile(null_fold, 2.5, 1);
result.null_hi = prctile(null_fold, 97.5, 1);
result.null_median = median(null_fold, 1, 'omitnan');
result.percentile = amplitude_percentile(result.observed, null_fold);
end


function v = null_first_value(peaks)
if ~peaks.first_measurable
    v = NaN;
else
    v = peaks.first_prominence;
end
end


function v = null_second_value(peaks)
if ~peaks.first_measurable
    v = NaN;
elseif peaks.second_measurable
    v = peaks.second_prominence;
else
    v = 0;
end
end


function p = exception_p(observed, nulls)
nulls = nulls(isfinite(nulls));
if ~isfinite(observed) || isempty(nulls)
    p = NaN;
    return
end
p = (1 + sum(nulls >= observed)) / (numel(nulls) + 1);
end


function pctl = amplitude_percentile(observed, null_rows)
% 100 means the observed value is above every shuffle. 0 means it is below
% every shuffle. A tie counts as half, so a high autocorrelogram stays high.
observed = double(observed(:).');
pctl = nan(size(observed));
if isempty(null_rows)
    return
end
finite = isfinite(null_rows);
n = sum(finite, 1);
below = sum(finite & null_rows < observed, 1);
equal = sum(finite & null_rows == observed, 1);
pctl = 100 * (below + 0.5 * equal) ./ n;
pctl(n == 0 | ~isfinite(observed)) = NaN;
end


function plot_grand_panel(result, panel_title, y_label)
if result.n < 1
    title({panel_title, 'n = 0'})
    axis off
    return
end
t = result.fold_t;
hold on
fill([t, fliplr(t)], [result.null_lo, fliplr(result.null_hi)], ...
    [0.75 0.75 0.75], 'EdgeColor', 'none', 'HandleVisibility', 'off')
plot(t, result.null_median, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1)
plot(t, result.observed, 'k', 'LineWidth', 1.6)
mark_grand_peak(result.first_lag, result.observed, result.fold_t)
mark_grand_peak(result.second_lag, result.observed, result.fold_t)
xlabel('Lag (s)')
ylabel(y_label)
title({panel_title, sprintf('n = %d\nfirst p = %.3g, second p = %.3g', ...
    result.n, result.first_p, result.second_p)})
box off
end


function plot_percentile_curves(records, grand, centers, smooth_n, group_order, group_title, n_shuffles)
scales = {'sum_pctl', 'ratio_pctl'};
scale_names = {'sum_to_one', 'rate_ratio'};
scale_labels = {'sum-to-1', 'rate ratio (Poisson = 1)'};
splits = {'locking', 'power'};
split_labels = {'locking', 'delta power'};
tails = {'low', 'high'};
tail_labels = {'Low 25%', 'High 25%'};
row_note = ['Top to bottom: each neuron, individual curves with the mean in bold, ' ...
    'percent of neurons at or above 95 (red) and at or below 5 (blue), mean logit.'];
for k = 1:numel(scales)
    for s = 1:numel(splits)
        figure('Name', ['Neuron shuffle percentile, ', split_labels{s}, ', ', scale_labels{k}], ...
            'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.02 0.02 0.96 0.96])
        for t = 1:numel(tails)
            for g = 1:numel(group_order)
                chosen = records_in_group(records, group_order{g});
                plot_percentile_column(chosen, splits{s}, tails{t}, scales{k}, ...
                    [group_title{g}, ', ', tail_labels{t}], n_shuffles, ...
                    4 * numel(tails), numel(group_order), 4 * (t - 1), g)
            end
        end
        sgtitle({['Each neuron against its own shuffles, ', split_labels{s}, ', ', scale_labels{k}], ...
            row_note})
    end
    for s = 1:numel(splits)
        figure('Name', ['Grand-average shuffle percentile, ', split_labels{s}, ', ', scale_labels{k}], ...
            'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.04 0.06 0.92 0.84])
        for t = 1:numel(tails)
            for g = 1:numel(group_order)
                block = grand.(scale_names{k}).(splits{s}).(tails{t}).(group_order{g});
                subplot(numel(tails), numel(group_order), (t - 1) * numel(group_order) + g)
                result = grand_peak_result(block, centers, smooth_n);
                plot_grand_percentile(result, [group_title{g}, ', ', tail_labels{t}])
            end
        end
        sgtitle({['Grand mean against shuffled means, ', split_labels{s}, ', ', scale_labels{k}], ...
            '100 means the mean autocorrelogram is above every shuffled mean.'})
    end
end
plot_all_run_percentiles(records, grand, centers, smooth_n, group_order, group_title, n_shuffles);
end


function plot_all_run_percentiles(records, grand, centers, smooth_n, group_order, group_title, n_shuffles)
scales = {'sum_pctl', 'ratio_pctl'};
scale_names = {'sum_to_one', 'rate_ratio'};
scale_labels = {'sum-to-1', 'rate ratio (Poisson = 1)'};
row_note = ['Top to bottom: each neuron, individual curves with the mean in bold, ' ...
    'percent of neurons at or above 95 (red) and at or below 5 (blue), mean logit.'];
for k = 1:numel(scales)
    figure('Name', ['Neuron shuffle percentile, all firing runs, ', scale_labels{k}], ...
        'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.02 0.04 0.96 0.92])
    for g = 1:numel(group_order)
        chosen = records_in_group(records, group_order{g});
        plot_percentile_column(chosen, 'all', '', scales{k}, ...
            [group_title{g}, ', all firing runs'], n_shuffles, ...
            4, numel(group_order), 0, g)
    end
    sgtitle({['Each neuron against its own shuffles, all firing runs, ', scale_labels{k}], ...
        row_note})
    if ~isfield(grand.(scale_names{k}), 'all')
        continue
    end
    figure('Name', ['Grand-average shuffle percentile, all firing runs, ', scale_labels{k}], ...
        'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.04 0.18 0.92 0.5])
    for g = 1:numel(group_order)
        block = grand.(scale_names{k}).all.(group_order{g});
        subplot(1, numel(group_order), g)
        result = grand_peak_result(block, centers, smooth_n);
        plot_grand_percentile(result, [group_title{g}, ', all firing runs'])
    end
    sgtitle({['Grand mean against shuffled means, all firing runs, ', scale_labels{k}], ...
        '100 means the mean autocorrelogram is above every shuffled mean.'})
end
end


function plot_percentile_column(records, split_name, tail_name, field_name, panel_title, n_shuffles, n_rows, n_cols, row_offset, col)
[lag, stack] = percentile_curve_stack(records, split_name, tail_name, field_name);
n = size(stack, 1);
if n < 1
    for row = 1:4
        subplot(n_rows, n_cols, (row_offset + row - 1) * n_cols + col)
        if row == 1
            title({panel_title, 'n = 0'})
        end
        axis off
    end
    return
end
draw_percentile_image(lag, stack, panel_title, n_rows, n_cols, row_offset + 1, col)
draw_percentile_lines(lag, stack, n_rows, n_cols, row_offset + 2, col)
draw_percentile_extreme_rates(lag, stack, n_rows, n_cols, row_offset + 3, col)
draw_percentile_logit(lag, stack, n_shuffles, n_rows, n_cols, row_offset + 4, col)
end


function draw_percentile_image(lag, stack, panel_title, n_rows, n_cols, row, col)
n = size(stack, 1);
subplot(n_rows, n_cols, (row - 1) * n_cols + col)
if n == 1
    imagesc(lag, [0.5 1.5], stack)
else
    imagesc(lag, 1:n, stack)
end
axis tight
clim([0 100])
if n == 1
    ylim([0.5 1.5])
end
colormap(gca, percentile_colormap())
if n > 1
    set(gca, 'YTick', [1 n])
end
ylabel('Neuron')
title({panel_title, sprintf('n = %d', n)})
box off
if col == n_cols
    colorbar
end
end


function draw_percentile_lines(lag, stack, n_rows, n_cols, row, col)
subplot(n_rows, n_cols, (row - 1) * n_cols + col)
hold on
for i = 1:size(stack, 1)
    plot(lag, stack(i, :), 'Color', [0.25 0.45 0.75 0.28], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off')
end
plot(lag, mean(stack, 1, 'omitnan'), 'k', 'LineWidth', 2.4)
plot(lag([1 end]), [50 50], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
plot(lag([1 end]), [5 5], '--', 'Color', [0.25 0.45 0.75], 'HandleVisibility', 'off')
plot(lag([1 end]), [95 95], '--', 'Color', [0.7 0.25 0.25], 'HandleVisibility', 'off')
ylim([0 100])
ylabel('Shuffle percentile')
box off
end


function draw_percentile_extreme_rates(lag, stack, n_rows, n_cols, row, col)
[high_percent, low_percent] = percentile_tail_fractions(stack, 95, 5);
subplot(n_rows, n_cols, (row - 1) * n_cols + col)
hold on
plot(lag, high_percent, 'Color', [0.7 0.25 0.25], 'LineWidth', 1.8)
plot(lag, low_percent, 'Color', [0.25 0.45 0.75], 'LineWidth', 1.8)
ylim([0 100])
ylabel('Neurons (%)')
if col == 1
    legend({'Percentile >= 95', 'Percentile <= 5'}, 'Location', 'northeast', 'Box', 'off')
end
box off
end


function draw_percentile_logit(lag, stack, n_shuffles, n_rows, n_cols, row, col)
logit = correct_percentile_logit(stack, n_shuffles);
subplot(n_rows, n_cols, (row - 1) * n_cols + col)
hold on
plot(lag([1 end]), [0 0], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
plot(lag, mean(logit, 1, 'omitnan'), 'k', 'LineWidth', 2.4)
xlabel('Lag (s)')
ylabel('Mean logit')
box off
end


function map = percentile_colormap()
n = 256;
half = n / 2;
low = [linspace(0.15, 0.96, half)', linspace(0.28, 0.96, half)', linspace(0.70, 0.96, half)'];
high = [linspace(0.96, 0.70, half)', linspace(0.96, 0.18, half)', linspace(0.96, 0.18, half)'];
map = [low; high];
end


function chosen = records_in_group(records, group_name)
chosen = {};
for i = 1:numel(records)
    group = records{i}.lock_group;
    if strcmp(group_name, 'delta')
        keep = any(strcmp(group, {'trough', 'peak'}));
    else
        keep = strcmp(group, group_name);
    end
    if keep
        chosen{end + 1, 1} = records{i}; %#ok<AGROW>
    end
end
end


function plot_grand_percentile(result, panel_title)
if result.n < 1 || isempty(result.fold_t)
    title({panel_title, 'n = 0'})
    axis off
    return
end
t = result.fold_t;
hold on
plot(t([1 end]), [50 50], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
plot(t([1 end]), [95 95], '--', 'Color', [0.7 0.25 0.25], 'HandleVisibility', 'off')
plot(t, result.percentile, 'k', 'LineWidth', 1.6)
mark_grand_peak(result.first_lag, result.percentile, t)
mark_grand_peak(result.second_lag, result.percentile, t)
ylim([0 100])
xlabel('Lag (s)')
ylabel('Shuffle percentile')
title({panel_title, sprintf('n = %d', result.n)})
box off
end


function mark_grand_peak(lag, y, t)
if ~isfinite(lag) || isempty(t)
    return
end
[~, idx] = min(abs(t - lag));
plot(t(idx), y(idx), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 5, 'HandleVisibility', 'off')
end
