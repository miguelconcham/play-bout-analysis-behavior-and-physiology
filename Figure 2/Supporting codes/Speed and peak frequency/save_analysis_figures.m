function saved = save_analysis_figures(figure_folder, figures_before)
%SAVE_ANALYSIS_FIGURES  Write figures opened during this script.
%
% figures_before is the list from findall before the script drew anything.
% Each new figure is saved as PNG, SVG, and a MATLAB .fig file. The name
% comes from the figure Name.

if exist(figure_folder, 'dir') ~= 7
    mkdir(figure_folder);
end
if nargin < 2
    figures_before = [];
end

open_figures = findall(0, 'Type', 'figure');
open_figures = open_figures(:);
keep = false(numel(open_figures), 1);
for k = 1:numel(open_figures)
    keep(k) = ~any(open_figures(k) == figures_before);
end
open_figures = open_figures(keep);
numbers = zeros(numel(open_figures), 1);
for k = 1:numel(open_figures)
    numbers(k) = get(open_figures(k), 'Number');
end
[~, order] = sort(numbers);
open_figures = open_figures(order);

saved = cell(numel(open_figures), 1);
used_names = cell(numel(open_figures), 1);
for k = 1:numel(open_figures)
    fig = open_figures(k);
    base = analysis_figure_basename(get(fig, 'Name'), numbers(order(k)));
    candidate = base;
    suffix = 2;
    while any(strcmp(used_names, candidate))
        candidate = sprintf('%s-%d', base, suffix);
        suffix = suffix + 1;
    end
    used_names{k} = candidate;
    png_file = fullfile(figure_folder, [candidate '.png']);
    svg_file = fullfile(figure_folder, [candidate '.svg']);
    fig_file = fullfile(figure_folder, [candidate '.fig']);
    print(fig, png_file, '-dpng', '-r200');
    try
        print(fig, '-vector', '-dsvg', svg_file);
    catch
        print(fig, svg_file, '-dsvg');
    end
    savefig(fig, fig_file);
    saved{k} = png_file;
end
end
