function draw_ccg_columns(figure_folder, file_name, reference_title, columns, centers, x_lim)
%DRAW_CCG_COLUMNS  One reference group: a column per other group, three rows.
%
% columns(c).title is the target group. columns(c).top and
% columns(c).middle are images, rows by lag. columns(c).mean, lo, and hi
% are the bottom trace and its 95% interval. The two images share one
% color scale. x_lim is the lag window in seconds, [-0.5 0.5] when omitted.
% The bottom y-axis is set from the mean curve and its interval inside
% that window. Zero on the bottom row is independence after normalization.

if nargin < 6 || isempty(x_lim)
    x_lim = [-0.5 0.5];
end
n_columns = numel(columns);
if n_columns < 1
    return
end
fig = figure('Name', reference_title, 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.02 0.06 0.96 0.86]);
limit = color_limit(columns, centers, x_lim);
for c = 1:n_columns
    draw_image(centers, columns(c).top, 3, n_columns, c, limit, x_lim);
    ylabel(columns(c).top_label);
    title(columns(c).title);
    draw_image(centers, columns(c).middle, 3, n_columns, n_columns + c, limit, x_lim);
    ylabel(columns(c).middle_label);
    subplot(3, n_columns, 2 * n_columns + c);
    hold on
    draw_interval(centers, columns(c).lo, columns(c).hi);
    if any(isfinite(columns(c).mean))
        plot(centers, columns(c).mean, 'k', 'LineWidth', 1.4);
    end
    plot(x_lim, [0 0], 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
    xlim(x_lim);
    ylim(ccg_mean_ylim(columns(c).mean, columns(c).lo, columns(c).hi, centers, x_lim));
    xlabel('Lag (s)');
    if c == 1
        ylabel('Normalized excess');
    end
    if isfield(columns(c), 'note') && ~isempty(columns(c).note)
        text(0.02, 0.92, columns(c).note, 'Units', 'normalized', 'FontSize', 8);
    end
end
if exist('sgtitle', 'file') == 2 || exist('sgtitle', 'builtin') == 5
    sgtitle(reference_title);
end
if isempty(figure_folder)
    return
end
if exist(figure_folder, 'dir') ~= 7
    mkdir(figure_folder);
end
print(fig, fullfile(figure_folder, [file_name, '.svg']), '-dsvg');
print(fig, fullfile(figure_folder, [file_name, '.jpg']), '-djpeg', '-r200');
end


function draw_image(centers, image, n_rows, n_columns, slot, limit, x_lim)
subplot(n_rows, n_columns, slot);
if isempty(image) || ~any(isfinite(image(:)))
    axis off
    return
end
imagesc(centers, 1:size(image, 1), image);
set(gca, 'Color', [0.93 0.93 0.93]);
caxis([-limit, limit]);
colormap(gca, blue_red_map());
xlim(x_lim);
end


function draw_interval(centers, lo, hi)
covered = isfinite(lo) & isfinite(hi);
if sum(covered) < 2
    return
end
fill([centers(covered), fliplr(centers(covered))], ...
    [lo(covered), fliplr(hi(covered))], [0.2 0.2 0.2], ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end


function limit = color_limit(columns, centers, x_lim)
shown = centers >= x_lim(1) & centers <= x_lim(2);
values = zeros(0, 1);
for c = 1:numel(columns)
    top = visible_values(columns(c).top, shown);
    middle = visible_values(columns(c).middle, shown);
    values = [values; top(:); middle(:)]; %#ok<AGROW>
end
if isempty(values)
    limit = 1;
    return
end
limit = prctile(abs(values), 98);
if ~isfinite(limit) || limit <= 0
    limit = max(abs(values));
end
if ~isfinite(limit) || limit <= 0
    limit = 1;
end
end


function values = visible_values(image, shown)
if isempty(image) || ~any(shown)
    values = [];
    return
end
if size(image, 2) ~= numel(shown)
    shown = true(1, size(image, 2));
end
block = image(:, shown);
values = block(isfinite(block));
end


function colors = blue_red_map()
n = 256;
half = n / 2;
blue = [linspace(0.15, 1, half)', linspace(0.25, 1, half)', linspace(0.70, 1, half)'];
red = [linspace(1, 0.70, half)', linspace(1, 0.15, half)', linspace(1, 0.15, half)'];
colors = [blue; red];
end
