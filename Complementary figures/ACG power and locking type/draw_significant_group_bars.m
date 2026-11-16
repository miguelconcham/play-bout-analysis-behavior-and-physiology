function fig = draw_significant_group_bars(panels, method_label, group_names, figure_name, title_lines)
%DRAW_SIGNIFICANT_GROUP_BARS  Two-row summary of significant feature tests.
%
% The top axes are the oscillatory features. The bottom axes are the
% non-oscillatory features. Each feature has two adjacent mean bars.
% The neurons are drawn first, so they sit behind the bars. There is no
% violin. Stars are the Holm-Sidak labels stored on the panels.

if nargin < 2 || isempty(method_label)
    method_label = 'cluster';
end
if nargin < 3 || isempty(group_names)
    group_names = {'Oscillatory', 'Non-oscillatory'};
end
if nargin < 4 || isempty(figure_name)
    figure_name = 'Unpaired t summary oscillatory versus non-oscillatory';
end
if nargin < 5 || isempty(title_lines)
    title_lines = { ...
        ['Unpaired t-test. Oscillatory neurons are ', method_label, ' cluster 1.'], ...
        'Only features with Holm-Sidak p < 0.05 in that family. Bars are group means of the within-feature z-score.', ...
        '* p < 0.05, ** p < 0.01, *** p < 0.001.'};
end
fig = figure('Name', figure_name, 'Color', 'w');
set(fig, 'Units', 'pixels', 'Position', [80 80 1400 900])
% Keep each axes in its own variable. A numeric array cannot hold a
% MATLAB axes object, and that assignment stopped the figure before it
% was saved.
positions = [0.08 0.54 0.88 0.30; 0.08 0.12 0.88 0.30];
axes_list = cell(1, 2);
for g = 1:2
    ax = subplot(2, 1, g);
    draw_one_group_panel(panels(g), group_names)
    set(ax, 'Tag', 'wide-bars')
    axes_list{g} = ax;
end
for g = 1:2
    set(axes_list{g}, 'Position', positions(g, :))
end
annotation(fig, 'textbox', [0.06 0.88 0.88 0.11], 'String', title_lines, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Interpreter', 'none', ...
    'FontSize', 11, 'FitBoxToText', 'off')
end


function draw_one_group_panel(panel, group_names)
titles = panel.titles;
n = numel(titles);
if n == 0
    axis off
    text(0.5, 0.5, [panel.family, ' features: no Holm-Sidak p below 0.05.'], ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Interpreter', 'none')
    return
end
hold on
osc_color = [0.85 0.33 0.10];
non_color = [0.12 0.45 0.70];
feature_step = 1.20;
within = 0.38;
bar_width = 0.32;
x1 = zeros(1, n);
x2 = zeros(1, n);
center = zeros(1, n);
for k = 1:n
    x1(k) = (k - 1) * feature_step + 1;
    x2(k) = x1(k) + within;
    center(k) = (x1(k) + x2(k)) / 2;
end
xlim([x1(1) - 0.55, x2(end) + 0.55])
plot(xlim, [0 0], ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off')
for k = 1:n
    scatter_group(x1(k), panel.z1{k}, osc_color)
    scatter_group(x2(k), panel.z2{k}, non_color)
end
for k = 1:n
    if k == 1
        vis = 'on';
        osc_name = group_names{1};
        non_name = group_names{2};
    else
        vis = 'off';
        osc_name = '';
        non_name = '';
    end
    draw_mean_bar(x1(k), panel.mean1(k), bar_width, osc_color, vis, osc_name)
    draw_mean_bar(x2(k), panel.mean2(k), bar_width, non_color, vis, non_name)
end
tops = zeros(1, n);
bottoms = zeros(1, n);
for k = 1:n
    cloud = [panel.z1{k}(:); panel.z2{k}(:); panel.mean1(k); panel.mean2(k)];
    cloud = cloud(isfinite(cloud));
    if isempty(cloud)
        tops(k) = 0;
        bottoms(k) = 0;
    else
        tops(k) = max(cloud);
        bottoms(k) = min(cloud);
    end
end
span = max(tops) - min(bottoms);
if ~isfinite(span) || span <= 0
    span = 1;
end
gap = 0.08 * span;
for k = 1:n
    if isempty(panel.stars{k})
        continue
    end
    text(center(k), tops(k) + gap, panel.stars{k}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 13, 'Interpreter', 'none')
end
low = min([bottoms, 0]);
high = max([tops, 0]);
ylim([low - 0.10 * span, high + gap + 0.28 * span])
set(gca, 'XTick', center, 'XTickLabel', titles, ...
    'TickLabelInterpreter', 'none', 'FontSize', 9)
xtickangle(30)
xlim([x1(1) - 0.55, x2(end) + 0.55])
ylabel('Mean z-score within feature')
title([panel.family, ' features'], 'Interpreter', 'none')
box off
legend('Location', 'northeast', 'Interpreter', 'none')
end


function scatter_group(x_center, z, color)
z = z(:);
z = z(isfinite(z));
if isempty(z)
    return
end
offset = bar_point_offset(numel(z));
scatter(x_center + offset, z, 12, color, 'filled', ...
    'MarkerFaceAlpha', 0.45, 'MarkerEdgeColor', 'none', ...
    'HandleVisibility', 'off')
end


function draw_mean_bar(x_center, mean_value, width, face_color, visibility, display_name)
if ~isfinite(mean_value)
    return
end
half = width / 2;
x = [x_center - half, x_center + half, x_center + half, x_center - half];
y = [0, 0, mean_value, mean_value];
fill(x, y, face_color, 'FaceAlpha', 0.45, 'EdgeColor', 'k', ...
    'LineWidth', 0.9, 'HandleVisibility', visibility, ...
    'DisplayName', display_name)
end


function offset = bar_point_offset(n_points)
if n_points <= 1
    offset = zeros(n_points, 1);
    return
end
offset = linspace(-0.11, 0.11, n_points).';
end
