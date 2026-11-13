%% Analyze_session_acg_shuffle
% Figures for the folder written by Estimate_session_acg_shuffle.
%
% Run this after the estimate. It reads grand_average.mat and the neuron
% files in the session ACG folder. Spike times are not read, and the
% shuffle is not run again.
%
% The grand-mean figure is the saved sum-to-1 autocorrelogram. The gray
% band is the 2.5% to 97.5% of the shuffle means, and the dashed line is
% their median. The percentile figure is the folded, smoothed grand mean
% as a percentile of the folded, smoothed shuffle means. 100 means that
% lag is above every shuffle mean.
%
% Neuron figures use the percentile curve already stored on each neuron.
% Row 1 of each image is the neuron with the highest mean percentile from
% 200 to 500 ms. A tail that was not estimated for that neuron is left out.
% Delta-locked is trough-locked and peak-locked together.
%
% Each show_* pair is [save, keep open]. [false false] skips that figure.
% figure_folder empty keeps the open figures on screen only. Set
% rebuild_grand true to rebuild grand_average.mat from the saved nulls
% before drawing. A missing grand_average.mat is built either way.

%% Paths and parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

save_root = fullfile(data_root, 'Analysis results', 'session acg shuffle');
figure_folder = '';
rebuild_grand = false;
sort_window_s = [0.2 0.5];
show_grand_mean = [false true];
show_grand_percentile = [false true];
show_power_overlay = [false true];
show_oscillation_counts = [false true];
show_neuron_percentiles = [false true];

show_grand_mean = figure_flag(show_grand_mean, 'show_grand_mean');
show_grand_percentile = figure_flag(show_grand_percentile, 'show_grand_percentile');
show_power_overlay = figure_flag(show_power_overlay, 'show_power_overlay');
show_oscillation_counts = figure_flag(show_oscillation_counts, 'show_oscillation_counts');
show_neuron_percentiles = figure_flag(show_neuron_percentiles, 'show_neuron_percentiles');

%% Grand average
grand_file = fullfile(save_root, 'grand_average.mat');
rebuilt = false;
if exist(grand_file, 'file') ~= 2
    session_acg_log('', 'grand_average.mat is missing. Building it from the saved nulls.');
    session_acg_grand_average(save_root);
    rebuilt = true;
elseif rebuild_grand
    session_acg_log('', 'Rebuilding the grand average from the saved nulls.');
    session_acg_grand_average(save_root);
    rebuilt = true;
end
if exist(grand_file, 'file') ~= 2
    error('session_acg:plot', ...
        'No grand average at %s. Run Estimate_session_acg_shuffle first.', grand_file);
end
stored = load(grand_file, 'grand');
if ~isfield(stored, 'grand')
    error('session_acg:plot', '%s has no grand average.', grand_file);
end
grand = stored.grand;
clear stored
assert_grand(grand);
if ~rebuilt && isfield(grand, 'report') && (ischar(grand.report) || isstring(grand.report)) ...
        && strlength(string(grand.report)) > 0
    fprintf('%s', grand.report);
    if ~endsWith(string(grand.report), newline)
        fprintf('\n');
    end
end

curve_names = session_acg_names('curves');
group_names = session_acg_names('groups');
curve_title = session_curve_titles();
curve_short = session_curve_short();
group_title = session_group_titles();
centers = grand_centers(grand);
n_shuffles = grand_shuffle_count(grand);
session_acg_log('', sprintf('Plotting %s. Method %s. Shuffles %g.', ...
    save_root, grand_method(grand), n_shuffles));

%% Grand-mean figures
if figure_requested(show_grand_mean, figure_folder)
    plot_grand_means(grand, centers, curve_names, group_names, curve_title, group_title, ...
        show_grand_mean, figure_folder);
end
if figure_requested(show_grand_percentile, figure_folder)
    plot_grand_percentiles(grand, curve_names, group_names, curve_title, group_title, ...
        show_grand_percentile, figure_folder);
end
if figure_requested(show_power_overlay, figure_folder)
    plot_power_overlay(grand, centers, group_names, group_title, ...
        show_power_overlay, figure_folder);
end
if figure_requested(show_oscillation_counts, figure_folder)
    plot_oscillation_counts(grand, curve_names, curve_short, group_names, group_title, ...
        show_oscillation_counts, figure_folder);
end

%% Neuron percentile figures
if figure_requested(show_neuron_percentiles, figure_folder)
    session_acg_log('', ['Reading neuron files for the percentile figures. ' ...
        'Each file is read in full; only the saved percentile curve is kept.']);
    records = collect_session_acg_percentiles(save_root);
    plot_neuron_percentiles(records, curve_names, group_names, curve_title, group_title, ...
        n_shuffles, sort_window_s, show_neuron_percentiles, figure_folder);
end
session_acg_log('', 'Session ACG plotting finished.');


function flag = figure_flag(flag, name)
if ~(isnumeric(flag) || islogical(flag)) || numel(flag) ~= 2
    error('session_acg:plot', '%s must be [save, keep open].', name);
end
flag = logical(flag(:))';
end


function requested = figure_requested(flag, folder)
requested = flag(2) || (flag(1) && ~isempty(folder));
end


function assert_grand(grand)
if ~isstruct(grand) || ~isfield(grand, 'groups') || ~isfield(grand, 'centers')
    error('session_acg:plot', ...
        'grand_average.mat must contain grand.groups and grand.centers.');
end
end


function centers = grand_centers(grand)
centers = double(grand.centers(:))';
if numel(centers) < 2
    error('session_acg:plot', 'The grand-average lag axis has fewer than 2 bins.');
end
end


function n_shuffles = grand_shuffle_count(grand)
n_shuffles = NaN;
if isfield(grand, 'n_shuffles') && isscalar(grand.n_shuffles)
    n_shuffles = double(grand.n_shuffles);
end
if ~(isfinite(n_shuffles) && n_shuffles >= 1) && isfield(grand, 'settings') ...
        && isfield(grand.settings, 'n_shuffles')
    n_shuffles = double(grand.settings.n_shuffles);
end
if ~(isfinite(n_shuffles) && n_shuffles >= 1)
    error('session_acg:plot', 'The grand average has no shuffle count.');
end
end


function method = grand_method(grand)
method = 'saved';
if isfield(grand, 'settings') && isfield(grand.settings, 'acg_method') && ~isempty(grand.settings.acg_method)
    method = char(string(grand.settings.acg_method));
end
end


function note = areas_note(grand)
note = '';
if ~isfield(grand, 'settings') || ~isfield(grand.settings, 'areas_to_keep') || isempty(grand.settings.areas_to_keep)
    return
end
note = sprintf('Areas: %s.', strjoin(cellstr(grand.settings.areas_to_keep(:)'), ', '));
end


function titles = session_curve_titles()
titles.session = 'Whole-session spikes';
titles.fooof_low = 'Low FOOOF-delta windows';
titles.fooof_high = 'High FOOOF-delta windows';
titles.hilbert_low = 'Low Hilbert-delta windows';
titles.hilbert_high = 'High Hilbert-delta windows';
end


function titles = session_curve_short()
titles.session = 'Session';
titles.fooof_low = 'FOOOF low';
titles.fooof_high = 'FOOOF high';
titles.hilbert_low = 'Hilbert low';
titles.hilbert_high = 'Hilbert high';
end


function lines = figure_title(lines, grand)
note = areas_note(grand);
if ~isempty(note)
    lines{end + 1} = note; %#ok<AGROW>
end
end


function titles = session_group_titles()
titles.trough = 'Trough-locked';
titles.peak = 'Peak-locked';
titles.delta = 'Delta-locked';
titles.nonentrained = 'Non-entrained';
end


function counts = grand_counts(grand, group_name, curve_name)
counts = struct('n', 0, 'n_oscillates', 0);
if ~isfield(grand, 'counts') || ~isfield(grand.counts, group_name) ...
        || ~isfield(grand.counts.(group_name), curve_name)
    return
end
stored = grand.counts.(group_name).(curve_name);
if isfield(stored, 'n') && isscalar(stored.n)
    counts.n = double(stored.n);
end
if isfield(stored, 'n_oscillates') && isscalar(stored.n_oscillates)
    counts.n_oscillates = double(stored.n_oscillates);
end
end


function [block, result] = grand_curve(grand, group_name, curve_name)
block = [];
result = [];
if ~isfield(grand.groups, group_name) || ~isfield(grand.groups.(group_name), curve_name)
    return
end
block = grand.groups.(group_name).(curve_name);
if isfield(block, 'result')
    result = block.result;
end
end


function heading = panel_heading(group_title, counts, result)
heading = {group_title, peak_detail(counts, result)};
end


function detail = peak_detail(counts, result)
if counts.n < 1
    detail = 'n = 0';
    return
end
if isempty(result) || ~isfield(result, 'first_p') || ~isscalar(result.first_p) || ~isfinite(result.first_p)
    detail = sprintf('n = %d, no delta-peak test', counts.n);
    return
end
if isfield(result, 'oscillates') && isscalar(result.oscillates) && result.oscillates
    call = 'oscillates';
else
    call = 'not above the null';
end
detail = sprintf('n = %d, p = %.3g, %s', counts.n, result.first_p, call);
end


function plot_grand_means(grand, centers, curve_names, group_names, curve_title, group_title, flag, folder)
fig = figure('Name', 'Session ACG grand mean', 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.02 0.04 0.96 0.92]);
n_rows = numel(curve_names);
n_cols = numel(group_names);
for c = 1:n_rows
    for g = 1:n_cols
        subplot(n_rows, n_cols, (c - 1) * n_cols + g)
        counts = grand_counts(grand, group_names{g}, curve_names{c});
        [block, result] = grand_curve(grand, group_names{g}, curve_names{c});
        draw_mean_panel(block, centers, counts, result, ...
            [{curve_title.(curve_names{c})}, panel_heading(group_title.(group_names{g}), counts, result)]);
        if g == 1 && strcmp(get(gca, 'Visible'), 'on')
            ylabel('Sum-to-1')
        end
    end
end
sgtitle(figure_title({'Grand mean autocorrelogram', ...
    'Gray band: 2.5% to 97.5% of the shuffle means. Dashed line: their median. Black line: mean autocorrelogram.', ...
    'Dotted lines mark 200 and 500 ms. A neuron oscillates when its delta-peak prominence beats its own null.'}, ...
    grand), 'Interpreter', 'none')
finish_figure(fig, flag, folder);
end


function draw_mean_panel(block, centers, counts, result, heading)
if counts.n < 1 || isempty(block) || ~isfield(block, 'mean') || numel(block.mean) ~= numel(centers)
    title(heading, 'Interpreter', 'none')
    axis off
    return
end
observed = double(block.mean(:))';
if ~any(isfinite(observed))
    title(heading, 'Interpreter', 'none')
    axis off
    return
end
hold on
if isfield(block, 'null_mean') && ~isempty(block.null_mean) && size(block.null_mean, 2) == numel(centers)
    band = session_acg_null_band(block.null_mean);
    fill_band(centers, band.low, band.high, [0.75 0.75 0.75])
    plot(centers, band.median, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1)
end
plot(centers, observed, 'k', 'LineWidth', 1.6)
mark_signed_peaks(result, observed, centers)
guide_lags(centers, [-0.5 -0.2 0.2 0.5])
xlim(centers([1 end]))
xlabel('Lag (s)')
title(heading, 'Interpreter', 'none')
set(gca, 'FontSize', 8)
box off
end


function plot_grand_percentiles(grand, curve_names, group_names, curve_title, group_title, flag, folder)
fig = figure('Name', 'Session ACG grand percentile', 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.02 0.04 0.96 0.92]);
n_rows = numel(curve_names);
n_cols = numel(group_names);
for c = 1:n_rows
    for g = 1:n_cols
        subplot(n_rows, n_cols, (c - 1) * n_cols + g)
        counts = grand_counts(grand, group_names{g}, curve_names{c});
        [~, result] = grand_curve(grand, group_names{g}, curve_names{c});
        draw_percentile_panel(result, counts, ...
            [{curve_title.(curve_names{c})}, panel_heading(group_title.(group_names{g}), counts, result)]);
        if g == 1 && strcmp(get(gca, 'Visible'), 'on')
            ylabel('Shuffle percentile')
        end
    end
end
sgtitle(figure_title({'Grand mean as a shuffle percentile', ...
    'The folded, smoothed mean against the folded, smoothed shuffle means. 100 is above every shuffle.', ...
    'Dotted lines mark 200 and 500 ms, and the 50th and 95th percentiles.'}, ...
    grand), 'Interpreter', 'none')
finish_figure(fig, flag, folder);
end


function draw_percentile_panel(result, counts, heading)
if counts.n < 1 || isempty(result) || ~isfield(result, 'pctl') || ~isfield(result, 'pctl_lag')
    title(heading, 'Interpreter', 'none')
    axis off
    return
end
lag = double(result.pctl_lag(:))';
y = double(result.pctl(:))';
if isempty(lag) || numel(lag) ~= numel(y) || ~any(isfinite(y))
    title(heading, 'Interpreter', 'none')
    axis off
    return
end
hold on
plot(lag([1 end]), [50 50], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
plot(lag([1 end]), [95 95], '--', 'Color', [0.7 0.25 0.25], 'HandleVisibility', 'off')
plot(lag, y, 'k', 'LineWidth', 1.6)
mark_folded_peak(result, 'first', y, lag)
mark_folded_peak(result, 'second', y, lag)
xlim(lag([1 end]))
ylim([0 100])
guide_lags(lag, [0.2 0.5])
ylim([0 100])
xlabel('Lag (s)')
title(heading, 'Interpreter', 'none')
set(gca, 'FontSize', 8)
box off
end


function plot_power_overlay(grand, centers, group_names, group_title, flag, folder)
pairs = { ...
    'fooof_low', 'fooof_high', 'FOOOF delta'; ...
    'hilbert_low', 'hilbert_high', 'Hilbert delta'};
low_color = [0.20 0.45 0.80];
high_color = [0.85 0.25 0.20];
fig = figure('Name', 'Session ACG low versus high delta', 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.06 0.12 0.88 0.72]);
n_cols = numel(group_names);
for row = 1:size(pairs, 1)
    for g = 1:n_cols
        subplot(size(pairs, 1), n_cols, (row - 1) * n_cols + g)
        draw_overlay_panel(grand, centers, group_names{g}, pairs{row, 1}, pairs{row, 2}, ...
            group_title.(group_names{g}), low_color, high_color, g == 1)
        if g == 1
            ylabel(pairs{row, 3})
        end
    end
end
sgtitle(figure_title({'Low 25% versus high 25% delta power', ...
    'Solid lines are the mean autocorrelogram. Dashed lines are the median shuffle. Blue is low and red is high.'}, ...
    grand), 'Interpreter', 'none')
finish_figure(fig, flag, folder);
end


function draw_overlay_panel(grand, centers, group_name, low_name, high_name, panel_title, low_color, high_color, show_legend)
low_counts = grand_counts(grand, group_name, low_name);
high_counts = grand_counts(grand, group_name, high_name);
[low_block, low_result] = grand_curve(grand, group_name, low_name);
[high_block, high_result] = grand_curve(grand, group_name, high_name);
low_row = curve_row(low_block, centers);
high_row = curve_row(high_block, centers);
if isempty(low_row) && isempty(high_row)
    title({panel_title, 'n = 0'}, 'Interpreter', 'none')
    axis off
    return
end
hold on
low_handle = plot_overlay_curve(low_block, low_row, centers, low_result, low_color);
high_handle = plot_overlay_curve(high_block, high_row, centers, high_result, high_color);
guide_lags(centers, [-0.5 -0.2 0.2 0.5])
xlim(centers([1 end]))
xlabel('Lag (s)')
title({panel_title, sprintf('low n = %d, high n = %d', low_counts.n, high_counts.n)}, ...
    'Interpreter', 'none')
box off
if show_legend && ~isempty(low_handle) && ~isempty(high_handle)
    legend([low_handle high_handle], {'Low 25%', 'High 25%'}, ...
        'Location', 'northeast', 'Box', 'off')
elseif show_legend && ~isempty(low_handle)
    legend(low_handle, {'Low 25%'}, 'Location', 'northeast', 'Box', 'off')
elseif show_legend && ~isempty(high_handle)
    legend(high_handle, {'High 25%'}, 'Location', 'northeast', 'Box', 'off')
end
set(gca, 'FontSize', 8)
end


function handle = plot_overlay_curve(block, observed, centers, result, color)
handle = gobjects(0);
if isempty(observed)
    return
end
if isfield(block, 'null_mean') && ~isempty(block.null_mean) && size(block.null_mean, 2) == numel(centers)
    band = session_acg_null_band(block.null_mean);
    plot(centers, band.median, '--', 'Color', color, 'LineWidth', 1, 'HandleVisibility', 'off')
end
handle = plot(centers, observed, 'Color', color, 'LineWidth', 1.6);
mark_signed_peaks(result, observed, centers)
end


function row = curve_row(block, centers)
row = [];
if isempty(block) || ~isfield(block, 'mean') || numel(block.mean) ~= numel(centers)
    return
end
row = double(block.mean(:))';
if ~any(isfinite(row))
    row = [];
end
end


function plot_oscillation_counts(grand, curve_names, curve_short, group_names, group_title, flag, folder)
percent = nan(numel(curve_names), numel(group_names));
for c = 1:numel(curve_names)
    for g = 1:numel(group_names)
        counts = grand_counts(grand, group_names{g}, curve_names{c});
        if counts.n > 0
            percent(c, g) = 100 * counts.n_oscillates / counts.n;
        end
    end
end
short_curve = cell(1, numel(curve_names));
for c = 1:numel(curve_names)
    short_curve{c} = curve_short.(curve_names{c});
end
legend_text = cell(1, numel(group_names));
for g = 1:numel(group_names)
    legend_text{g} = group_title.(group_names{g});
end
colors = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.55 0.35 0.70; 0.45 0.45 0.45];
alpha = NaN;
if isfield(grand, 'settings') && isfield(grand.settings, 'oscillation_alpha')
    alpha = double(grand.settings.oscillation_alpha);
end
fig = figure('Name', 'Session ACG delta oscillation', 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.12 0.16 0.7 0.62]);
bars = bar(percent);
for g = 1:min(numel(bars), size(colors, 1))
    bars(g).FaceColor = colors(g, :);
end
set(gca, 'XTick', 1:numel(curve_names), 'XTickLabel', short_curve)
ylabel(sprintf('Neurons with a delta oscillation (%%), alpha %g', alpha))
ylim([0 100])
legend(legend_text, 'Location', 'eastoutside', 'Box', 'off')
title('Share of neurons whose delta-peak prominence beats their own null', 'Interpreter', 'none')
box off
finish_figure(fig, flag, folder);
end


function plot_neuron_percentiles(records, curve_names, group_names, curve_title, group_title, n_shuffles, sort_window_s, flag, folder)
if isempty(records)
    session_acg_log('', 'No neuron files, so the neuron percentile figures were skipped.');
    return
end
row_note = ['Top to bottom: each neuron, individual curves with the mean in bold, ' ...
    'percent of neurons at or above 95 (red) and at or below 5 (blue), mean logit.'];
sort_note = sprintf('Row 1 has the highest mean percentile from %.0f to %.0f ms.', ...
    1000 * sort_window_s(1), 1000 * sort_window_s(2));
for c = 1:numel(curve_names)
    fig = figure('Name', ['Session ACG neurons, ' curve_title.(curve_names{c})], 'Color', 'w', ...
        'Units', 'normalized', 'OuterPosition', [0.02 0.04 0.96 0.92]);
    for g = 1:numel(group_names)
        packed = session_acg_percentile_matrix(records, group_names{g}, curve_names{c}, ...
            n_shuffles, sort_window_s);
        plot_percentile_column(packed, group_title.(group_names{g}), 4, numel(group_names), g)
    end
    sgtitle({['Each neuron against its own shuffles, ' curve_title.(curve_names{c})], ...
        row_note, sort_note}, 'Interpreter', 'none')
    finish_figure(fig, flag, folder);
end
fig = figure('Name', 'Session ACG neurons, all lock groups', 'Color', 'w', ...
    'Units', 'normalized', 'OuterPosition', [0.02 0.04 0.96 0.92]);
for c = 1:numel(curve_names)
    packed = session_acg_percentile_matrix(records, 'all', curve_names{c}, ...
        n_shuffles, sort_window_s);
    plot_percentile_column(packed, curve_title.(curve_names{c}), 4, numel(curve_names), c)
end
sgtitle({'All lock groups together', row_note, sort_note}, 'Interpreter', 'none')
finish_figure(fig, flag, folder);
end


function plot_percentile_column(packed, panel_title, n_rows, n_cols, col)
heading = percentile_heading(panel_title, packed);
if packed.n < 1
    for row = 1:4
        subplot(n_rows, n_cols, (row - 1) * n_cols + col)
        if row == 1
            title(heading, 'Interpreter', 'none')
        end
        axis off
    end
    return
end
draw_percentile_image(packed, heading, n_rows, n_cols, 1, col)
draw_percentile_lines(packed, n_rows, n_cols, 2, col)
draw_percentile_extreme_rates(packed, n_rows, n_cols, 3, col)
draw_percentile_logit(packed, n_rows, n_cols, 4, col)
end


function heading = percentile_heading(panel_title, packed)
if packed.n_dropped > 0
    detail = sprintf('n = %d, %d other axis left out', packed.n, packed.n_dropped);
else
    detail = sprintf('n = %d', packed.n);
end
heading = {panel_title, detail};
end


function draw_percentile_image(packed, heading, n_rows, n_cols, row, col)
lag = packed.lag;
stack = packed.percentile;
n = packed.n;
subplot(n_rows, n_cols, (row - 1) * n_cols + col)
if n == 1
    imagesc(lag, [0.5 1.5], stack)
    ylim([0.5 1.5])
else
    imagesc(lag, 1:n, stack)
end
axis tight
clim([0 100])
colormap(gca, percentile_colormap())
if n > 1
    set(gca, 'YTick', [1 n])
end
ylabel('Neuron')
title(heading, 'Interpreter', 'none')
box off
if col == n_cols
    colorbar
end
end


function draw_percentile_lines(packed, n_rows, n_cols, row, col)
lag = packed.lag;
stack = packed.percentile;
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


function draw_percentile_extreme_rates(packed, n_rows, n_cols, row, col)
[high_percent, low_percent] = percentile_tail_fractions(packed.percentile, 95, 5);
subplot(n_rows, n_cols, (row - 1) * n_cols + col)
hold on
plot(packed.lag, high_percent, 'Color', [0.7 0.25 0.25], 'LineWidth', 1.8)
plot(packed.lag, low_percent, 'Color', [0.25 0.45 0.75], 'LineWidth', 1.8)
ylim([0 100])
ylabel('Neurons (%)')
if col == 1
    legend({'Percentile >= 95', 'Percentile <= 5'}, 'Location', 'northeast', 'Box', 'off')
end
box off
end


function draw_percentile_logit(packed, n_rows, n_cols, row, col)
subplot(n_rows, n_cols, (row - 1) * n_cols + col)
hold on
plot(packed.lag([1 end]), [0 0], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
plot(packed.lag, mean(packed.logit, 1, 'omitnan'), 'k', 'LineWidth', 2.4)
xlabel('Lag (s)')
ylabel('Mean logit')
box off
end


function fill_band(lag, low, high, color)
lag = double(lag(:))';
low = double(low(:))';
high = double(high(:))';
ok = isfinite(lag) & isfinite(low) & isfinite(high);
breaks = diff([false, ok, false]);
starts = find(breaks == 1);
stops = find(breaks == -1) - 1;
for k = 1:numel(starts)
    idx = starts(k):stops(k);
    if numel(idx) < 2
        continue
    end
    fill([lag(idx), fliplr(lag(idx))], [low(idx), fliplr(high(idx))], color, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off')
end
end


function guide_lags(lag, marks)
if isempty(lag)
    return
end
yl = ylim;
for i = 1:numel(marks)
    if marks(i) < lag(1) || marks(i) > lag(end)
        continue
    end
    plot([marks(i) marks(i)], yl, ':', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off')
end
ylim(yl)
end


function mark_signed_peaks(result, observed, centers)
if isempty(result)
    return
end
mark_one_signed_peak(result, 'first', observed, centers)
mark_one_signed_peak(result, 'second', observed, centers)
end


function mark_one_signed_peak(result, which, observed, centers)
present_name = [which '_present'];
lag_name = [which '_lag'];
if ~isfield(result, present_name) || ~isfield(result, lag_name) || ~result.(present_name)
    return
end
lag = double(result.(lag_name));
if ~isscalar(lag) || ~isfinite(lag)
    return
end
for sign = [-1 1]
    [~, index] = min(abs(centers - sign * lag));
    if isfinite(observed(index))
        plot(centers(index), observed(index), 'o', 'Color', 'k', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 4, 'HandleVisibility', 'off')
    end
end
end


function mark_folded_peak(result, which, y, lag)
present_name = [which '_present'];
lag_name = [which '_lag'];
if isempty(result) || ~isfield(result, present_name) || ~isfield(result, lag_name) || ~result.(present_name)
    return
end
peak_lag = double(result.(lag_name));
if ~isscalar(peak_lag) || ~isfinite(peak_lag) || isempty(lag)
    return
end
[~, index] = min(abs(lag - peak_lag));
if ~isfinite(y(index))
    return
end
plot(lag(index), y(index), 'o', 'Color', 'k', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 4, 'HandleVisibility', 'off')
end


function map = percentile_colormap()
n = 256;
half = n / 2;
low = [linspace(0.15, 0.96, half)', linspace(0.28, 0.96, half)', linspace(0.70, 0.96, half)'];
high = [linspace(0.96, 0.70, half)', linspace(0.96, 0.18, half)', linspace(0.96, 0.18, half)'];
map = [low; high];
end


function finish_figure(fig, flag, folder)
do_save = flag(1) && ~isempty(folder);
keep_open = flag(2);
if do_save
    if exist(folder, 'dir') ~= 7
        mkdir(folder);
    end
    stamp = char(datetime('now', 'TimeZone', 'Europe/Berlin', 'Format', 'yyyy-MM-dd'));
    day_folder = fullfile(folder, ['SessionACG_' stamp]);
    if exist(day_folder, 'dir') ~= 7
        mkdir(day_folder);
    end
    name = regexprep(strtrim(fig.Name), '[^\w]+', '_');
    name = regexprep(name, '^_|_$', '');
    if isempty(name)
        name = sprintf('figure_%d', fig.Number);
    end
    set(fig, 'PaperPositionMode', 'auto')
    print(fig, fullfile(day_folder, [name, '.svg']), '-dsvg')
    print(fig, fullfile(day_folder, [name, '.jpg']), '-djpeg', '-r200')
    session_acg_log('', sprintf('Saved %s', fullfile(day_folder, [name, '.svg'])));
end
if ~keep_open
    close(fig)
end
drawnow
end
