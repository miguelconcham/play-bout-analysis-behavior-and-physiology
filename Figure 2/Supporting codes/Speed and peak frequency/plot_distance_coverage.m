function plot_distance_coverage(report)
%PLOT_DISTANCE_COVERAGE  Session counts and play-minus-non-play power.
%
% The count panels are the number of electrode-1 sessions at each
% inter-animal distance. Positive distance z is farther than that
% session's own mean.

if isempty(report.animals)
    return
end

options = report.options;
x_label = 'Inter-animal distance (session z)';
y_label = 'Play - non-play (log10 amplitude)';

figure('Name', 'Distance coverage')
subplot(2, 3, 1)
plot_diff_curve(report.curve_x, report.diff_all, x_label, y_label, ...
    'All electrode-1 sessions', options.high_z)
subplot(2, 3, 2)
plot_diff_curve(report.curve_x, report.diff_far_curve, x_label, y_label, ...
    sprintf('Enough overlap at z >= %.1f', options.high_z), options.high_z)
subplot(2, 3, 3)
plot_diff_curve(report.curve_x, report.diff_close_curve, x_label, y_label, ...
    sprintf('Enough overlap at z <= %.1f', -options.high_z), options.high_z)
subplot(2, 3, 4)
plot_counts(report.curve_x, report.n_nonplay, report.n_both, ...
    options.min_bin_play_s, options.min_bin_nonplay_s, ...
    'All electrode-1 sessions', options.high_z)
subplot(2, 3, 5)
plot_counts(report.curve_x, report.n_nonplay_far, report.n_both_far, ...
    options.min_bin_play_s, options.min_bin_nonplay_s, ...
    sprintf('Enough overlap at z >= %.1f', options.high_z), options.high_z)
subplot(2, 3, 6)
plot_counts(report.curve_x, report.n_nonplay_close, report.n_both_close, ...
    options.min_bin_play_s, options.min_bin_nonplay_s, ...
    sprintf('Enough overlap at z <= %.1f', -options.high_z), options.high_z)

figure('Name', 'Distance sessions')
subplot(2, 2, 1)
plot_duration_bars(report.nonplay_far_s, report.enough_far, ...
    report.animals, report.partners, options.min_nonplay_s, ...
    sprintf('Non-play frame span at distance z >= %.1f (s)', options.high_z), ...
    'Black: play and non-play both pass')
subplot(2, 2, 2)
plot_duration_bars(report.nonplay_close_s, report.enough_close, ...
    report.animals, report.partners, options.min_nonplay_s, ...
    sprintf('Non-play frame span at distance z <= %.1f (s)', -options.high_z), ...
    'Black: play and non-play both pass')
subplot(2, 2, 3)
plot_supported_bins(report.far_bin_lo, report.diff_far_bin, ...
    report.far_bin_lo >= options.high_z, options.high_z, ...
    'Start of farthest supported distance bin (z)', y_label, ...
    'Bin starts below the far cut', 'Bin starts at the far cut')
subplot(2, 2, 4)
plot_supported_bins(report.close_bin_lo, report.diff_close_bin, ...
    report.close_bin_hi <= -options.high_z, -options.high_z, ...
    'Start of closest supported distance bin (z)', y_label, ...
    'Bin extends past the close cut', 'Bin stays at the close cut')
end

function plot_diff_curve(x, diff_mat, x_label, y_label, plot_title, high_z)
hold on
mean_diff = mean(diff_mat, 1, 'omitnan');
if isempty(x) || ~any(isfinite(mean_diff))
    if ~isempty(x)
        xlim([x(1) x(end)])
    end
    ylim([-1 1])
    yline(0, 'k:')
    mark_distance_cuts(high_z)
    title(plot_title)
    xlabel(x_label)
    ylabel(y_label)
    text(mean(xlim), 0, 'No bin has 3 sessions', 'HorizontalAlignment', 'center')
    return
end
p_val = nan(1, numel(x));
ci = nan(2, numel(x));
for b = 1:numel(x)
    y = diff_mat(:, b);
    y = y(isfinite(y));
    if numel(y) < 3
        continue
    end
    [~, p_val(b), ci(:, b)] = ttest(y);
end
finite_ci = ~any(isnan(ci), 1);
if any(finite_ci)
    fill([x(finite_ci) fliplr(x(finite_ci))], ...
        [ci(1, finite_ci) fliplr(ci(2, finite_ci))], ...
        'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
end
plot(x, mean_diff, 'k', 'LineWidth', 1.5)
yline(0, 'k:')
axis tight
y_lim = ylim;
if all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
    fill_lo = y_lim(1) + 0.9 * (y_lim(2) - y_lim(1));
    fill_hi = y_lim(2);
    sig = p_val < 0.05;
    sig(isnan(p_val)) = false;
    start_end = [find(diff([0 sig 0]) == 1)' (find(diff([0 sig 0]) == -1) - 1)'];
    for j = 1:size(start_end, 1)
        fill(x(start_end(j, [1 2 2 1])), [fill_lo fill_lo fill_hi fill_hi], ...
            'k', 'EdgeColor', 'none', 'FaceAlpha', .35)
    end
end
mark_distance_cuts(high_z)
xlabel(x_label)
ylabel(y_label)
title(plot_title)
end

function plot_counts(x, n_nonplay, n_both, min_play_s, min_nonplay_s, plot_title, high_z)
hold on
plot(x, n_nonplay, 'Color', [0.65 0.65 0.65], 'LineWidth', 1.25)
plot(x, n_both, 'k', 'LineWidth', 1.5)
mark_distance_cuts(high_z)
xlabel('Inter-animal distance (session z)')
ylabel('Sessions')
title(plot_title)
legend({sprintf('Non-play >= %.1f s', min_nonplay_s), ...
    sprintf('Play >= %.1f s and non-play >= %.1f s', min_play_s, min_nonplay_s)}, ...
    'Location', 'best')
axis tight
end

function plot_duration_bars(duration_s, enough, animals, partners, threshold_s, x_label, plot_title)
[~, order] = sort(duration_s, 'descend');
n_sess = numel(duration_s);
labels = cell(n_sess, 1);
for s = 1:n_sess
    labels{s} = sprintf('%s P%d', animals{order(s)}, partners(order(s)));
end
bar_colors = repmat([0.75 0.75 0.75], n_sess, 1);
if any(enough)
    bar_colors(enough(order), :) = repmat([0 0 0], nnz(enough), 1);
end
duration_bars = barh(duration_s(order));
duration_bars.FaceColor = 'flat';
duration_bars.CData = bar_colors;
xline(threshold_s, 'k:')
set(gca, 'YTick', 1:n_sess, 'YTickLabel', labels, 'YDir', 'reverse')
xlabel(x_label)
title(plot_title)
end

function plot_supported_bins(bin_lo, diff_val, at_cut, cut_z, x_label, y_label, below_label, above_label)
hold on
has_bin = isfinite(bin_lo) & isfinite(diff_val);
at_cut = at_cut(:);
below = has_bin & ~at_cut;
above = has_bin & at_cut;
scatter(bin_lo(below), diff_val(below), 28, [0.6 0.6 0.6], 'filled')
scatter(bin_lo(above), diff_val(above), 36, 'k', 'filled')
yline(0, 'k:')
xline(cut_z, 'k:')
xlabel(x_label)
ylabel(y_label)
title('Each session at its extreme supported distance')
if any(below) && any(above)
    legend({below_label, above_label}, 'Location', 'best')
end
end

function mark_distance_cuts(high_z)
xline(-high_z, 'k:')
xline(high_z, 'k:')
end
