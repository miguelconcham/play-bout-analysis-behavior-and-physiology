%% Plot_playbout_behavior_variable_stacks
% Stack the peri-play-bout variables saved by
% Estimate_playbout_behavior_variables.m and plot one figure per variable.
%
% Variables plotted (each NP x T inside the saved matrices):
%   speed, acceleration, distance
% Vocalization is saved with those matrices and is left out of these figures.
% Window is the one stored with the matrices (onset at 0; default -20..+20 s).
%
% For each variable:
%   1) Z-score within an animal before stacking. The mean and SD are taken
%      over every finite peri-bout sample of that animal in that session,
%      then applied to all of that animal's play bouts. Animals therefore
%      share a color scale without one animal's raw units dominating.
%   2) Distance is stored twice (once per animal) and is the same trace.
%      It is kept once per play bout.
%   3) Figure, two rows:
%        top    every play bout of every animal (color = z-scored value)
%        bottom mean across those play bouts, with a 95% CI (mean +/- 1.96 SEM)
%
% The mean weights every play bout equally, so a recording with more bouts
% contributes more rows. Gray pixels in the image are missing samples.
% PNG files are written to a figures folder next to the saved .mat files.

%% Paths / options
experiment_folder = '\\experimentfs.home.bccn-berlin.de\experiment';
repo_root = [experiment_folder, ...
    '\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology'];
data_folder = fullfile(repo_root, 'Data', 'Analysis results', ...
    'playbout behavior regressors two animals');
figure_folder = fullfile(data_folder, 'figures');

color_limit = 3;   % symmetric z-score color range
ci_z = 1.96;       % normal 95% confidence interval
variables_to_plot = {'speed', 'acceleration', 'distance'};

%% Load every session saved by Estimate_playbout_behavior_variables
sessions = load_playbout_sessions(data_folder);
saved_names = sessions{1}.variable_names;
psth_time = sessions{1}.psth_time;
fprintf('Loaded %d sessions from\n%s\n', numel(sessions), data_folder);
fprintf('Plotting: %s\n', strjoin(variables_to_plot, ', '));
fprintf('Time: %.3g to %.3g s, %d bins\n', psth_time(1), psth_time(end), numel(psth_time));

if exist(figure_folder, 'dir') ~= 7
    mkdir(figure_folder);
end

for plot_index = 1:numel(variables_to_plot)
    variable_name = variables_to_plot{plot_index};
    variable_index = find(strcmp(saved_names, variable_name), 1);
    if isempty(variable_index)
        error('Saved matrices do not contain %s. Found: %s.', ...
            variable_name, strjoin(saved_names, ', '));
    end
    blocks = collect_variable_blocks(sessions, variable_index, variable_name);
    if isempty(blocks)
        warning('No finite play bouts for %s. Skipping.', variable_name);
        continue
    end
    plot_variable_stack(variable_name, psth_time, blocks, figure_folder, ...
        color_limit, ci_z);
end
fprintf('Done. Figures in:\n%s\n', figure_folder);

%% ========================= local helpers =========================

function sessions = load_playbout_sessions(data_folder)
files = dir(fullfile(data_folder, 'playbout_behavior_vars_*.mat'));
if ~isempty(files)
    names = {files.name};
    files = files(~strcmp(names, 'playbout_behavior_vars_session_list.mat'));
    [~, order] = sort({files.name});
    files = files(order);
end
if isempty(files)
    error(['No playbout_behavior_vars_*.mat files in:\n%s\n' ...
        'Run Estimate_playbout_behavior_variables.m first.'], data_folder);
end

sessions = cell(1, numel(files));
n_kept = 0;
reference_names = {};
reference_time = [];
for file_index = 1:numel(files)
    file_name = files(file_index).name;
    file_path = fullfile(files(file_index).folder, file_name);
    info = whos('-file', file_path);
    if ~ismember('playbout_behavior', {info.name})
        warning('Skipping %s (no playbout_behavior).', file_name);
        continue
    end
    loaded = load(file_path, 'playbout_behavior');
    pb = loaded.playbout_behavior;
    pb = standardize_session(pb, file_name);
    if isempty(reference_time)
        reference_names = pb.variable_names;
        reference_time = pb.psth_time;
    elseif ~isequal(pb.variable_names, reference_names)
        error('variable_names in %s do not match the first session.', file_name);
    elseif numel(pb.psth_time) ~= numel(reference_time) || ...
            max(abs(pb.psth_time - reference_time)) > 1e-6
        error('psth_time in %s does not match the first session.', file_name);
    end
    n_kept = n_kept + 1;
    sessions{n_kept} = pb;
end
sessions = sessions(1:n_kept);
if n_kept == 0
    error('No usable playbout_behavior files in %s.', data_folder);
end
end

function pb = standardize_session(pb, file_name)
if ~isfield(pb, 'session_name') || isempty(pb.session_name)
    pb.session_name = file_name;
end
if iscell(pb.session_name)
    pb.session_name = pb.session_name{1};
end
pb.session_name = char(pb.session_name);
pb.variable_names = cellstr(pb.variable_names);
pb.variable_names = pb.variable_names(:)';
pb.animal_names = cellstr(pb.animal_names);
pb.animal_names = pb.animal_names(:)';
for animal_index = 1:numel(pb.animal_names)
    pb.animal_names{animal_index} = char(pb.animal_names{animal_index});
end
pb.psth_time = double(pb.psth_time(:)');
if ~isfield(pb, 'vars_by_animal') || ~iscell(pb.vars_by_animal) || ...
        numel(pb.vars_by_animal) ~= numel(pb.animal_names)
    error('%s is missing vars_by_animal for each animal.', file_name);
end
if isempty(pb.variable_names) || numel(pb.psth_time) < 2 || any(diff(pb.psth_time) <= 0)
    error('%s has no variable names or a non-increasing time axis.', file_name);
end
end

function blocks = collect_variable_blocks(sessions, variable_index, variable_name)
blocks = struct('session_name', {}, 'animal_name', {}, 'shared', {}, 'values', {});
for session_index = 1:numel(sessions)
    pb = sessions{session_index};
    n_animals = numel(pb.animal_names);
    traces = cell(1, n_animals);
    for animal_index = 1:n_animals
        raw = pb.vars_by_animal{animal_index};
        if size(raw, 2) < variable_index || size(raw, 3) < 2
            error('vars for %s | %s do not contain %s.', ...
                pb.session_name, pb.animal_names{animal_index}, variable_name);
        end
        traces{animal_index} = reshape(raw(:, variable_index, :), ...
            size(raw, 1), size(raw, 3));
        if size(traces{animal_index}, 2) ~= numel(pb.psth_time)
            error('%s | %s has %d %s bins, but psth_time has %d.', ...
                pb.session_name, pb.animal_names{animal_index}, ...
                size(traces{animal_index}, 2), variable_name, numel(pb.psth_time));
        end
    end
    shared = traces_are_identical(traces);
    if shared
        animal_indices = 1;
    else
        animal_indices = 1:n_animals;
    end
    for animal_index = animal_indices
        values = finite_rows(zscore_matrix(traces{animal_index}));
        if isempty(values)
            warning('No finite %s samples for %s | %s.', ...
                variable_name, pb.session_name, pb.animal_names{animal_index});
            continue
        end
        block.session_name = pb.session_name;
        if shared
            block.animal_name = 'both';
        else
            block.animal_name = pb.animal_names{animal_index};
        end
        block.shared = shared;
        block.values = values;
        blocks(end + 1) = block; %#ok<AGROW>
    end
end
blocks = sort_blocks(blocks);
end

function shared = traces_are_identical(traces)
shared = false;
if numel(traces) < 2
    return
end
reference = traces{1};
for animal_index = 2:numel(traces)
    if ~isequaln(reference, traces{animal_index})
        return
    end
end
shared = true;
end

function Z = zscore_matrix(X)
% One mean and one SD over all finite samples, then applied in place.
Z = nan(size(X));
finite = isfinite(X);
if nnz(finite) < 2
    return
end
mu = mean(X(finite));
sigma = std(X(finite), 0);
if ~isfinite(sigma) || sigma == 0
    Z(finite) = 0;
    return
end
Z(finite) = (X(finite) - mu) ./ sigma;
end

function Y = finite_rows(X)
if isempty(X)
    Y = X;
    return
end
Y = X(any(isfinite(X), 2), :);
end

function blocks = sort_blocks(blocks)
n_blocks = numel(blocks);
if n_blocks < 2
    return
end
if all([blocks.shared])
    animal_key = {blocks.session_name};
else
    animal_key = {blocks.animal_name};
end
session_key = {blocks.session_name};
index_key = arrayfun(@(k) sprintf('%06d', k), 1:n_blocks, 'UniformOutput', false);
[~, order] = sortrows([animal_key(:), session_key(:), index_key(:)]);
blocks = blocks(order);
end

function plot_variable_stack(variable_name, time, blocks, figure_folder, ...
    color_limit, ci_z)
[stack, centers, labels, boundaries] = assemble_stack(blocks);
n_bouts = size(stack, 1);
n_blocks = numel(blocks);
shared = all([blocks.shared]);
if shared
    group_noun = 'sessions';
    image_title = sprintf(['Z-scored within each session, then stacked ' ...
        '(N = %d play bouts, %d sessions; same trace for both animals)'], ...
        n_bouts, n_blocks);
else
    group_noun = 'animal recordings';
    image_title = sprintf(['Z-scored within each animal, then stacked ' ...
        '(N = %d play bouts, %d animal recordings)'], n_bouts, n_blocks);
end
[mu, lo, hi] = mean_with_ci(stack, ci_z);
pretty = pretty_variable_name(variable_name);

fprintf('%s: %d play bouts from %d %s\n', pretty, n_bouts, n_blocks, group_noun);

fig = figure('Color', 'w', 'Position', [80 40 920 980], ...
    'NumberTitle', 'off', 'Name', ['Play-bout ' pretty]);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_image = nexttile(layout);
image_handle = imagesc(ax_image, time, 1:n_bouts, stack);
set(image_handle, 'AlphaData', double(isfinite(stack)));
hold(ax_image, 'on');
for boundary_index = 1:numel(boundaries)
    plot(ax_image, [time(1), time(end)], ...
        [boundaries(boundary_index), boundaries(boundary_index)], ...
        'Color', [0.15 0.15 0.15], 'LineWidth', 0.7, 'HandleVisibility', 'off');
end
xline(ax_image, 0, 'LineStyle', '--', 'Color', [0.15 0.15 0.15], ...
    'HandleVisibility', 'off');
colormap(ax_image, diverging_blue_red(256));
ax_image.CLim = [-color_limit, color_limit];
ax_image.Color = [0.86 0.86 0.86];
ax_image.YDir = 'reverse';
xlim(ax_image, [time(1), time(end)]);
ylim(ax_image, [0.5, n_bouts + 0.5]);
ax_image.YTick = centers;
ax_image.YTickLabel = labels;
ax_image.TickLabelInterpreter = 'none';
ax_image.FontSize = tick_font_size(n_blocks);
xlabel(ax_image, 'Time from play-bout onset (s)');
ylabel(ax_image, 'Play bout');
title(ax_image, image_title, 'Interpreter', 'none');
colorbar_handle = colorbar(ax_image);
colorbar_handle.Label.String = 'z-score';
colorbar_handle.Layout.Tile = 'east';
set(ax_image, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top');

ax_mean = nexttile(layout);
hold(ax_mean, 'on');
ci_handle = fill_confidence_interval(ax_mean, time, lo, hi);
mean_handle = plot(ax_mean, time, mu, 'k-', 'LineWidth', 1.6, 'DisplayName', 'Mean');
yline(ax_mean, 0, 'LineStyle', ':', 'Color', [0.45 0.45 0.45], ...
    'HandleVisibility', 'off');
xline(ax_mean, 0, 'LineStyle', '--', 'Color', [0.15 0.15 0.15], ...
    'HandleVisibility', 'off');
xlim(ax_mean, [time(1), time(end)]);
xlabel(ax_mean, 'Time from play-bout onset (s)');
ylabel(ax_mean, 'Mean z-score');
title(ax_mean, sprintf('Mean across play bouts with 95%% CI (mean +/- %.2f SEM)', ci_z), ...
    'Interpreter', 'none');
if isempty(ci_handle)
    legend(ax_mean, mean_handle, 'Location', 'best', 'Box', 'off');
else
    legend(ax_mean, [mean_handle, ci_handle], 'Location', 'best', 'Box', 'off');
end
set(ax_mean, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top');
linkaxes([ax_image, ax_mean], 'x');
title(layout, [pretty ' around play-bout onset'], 'Interpreter', 'none');

out_path = fullfile(figure_folder, ['playbout_stack_' file_token(variable_name) '.png']);
save_figure(fig, out_path);
fprintf('Saved %s\n', out_path);
end

function [stack, centers, labels, boundaries] = assemble_stack(blocks)
n_blocks = numel(blocks);
pieces = cell(1, n_blocks);
centers = zeros(1, n_blocks);
labels = cell(1, n_blocks);
boundaries = zeros(1, max(0, n_blocks - 1));
n_time = size(blocks(1).values, 2);
row = 0;
for block_index = 1:n_blocks
    values = blocks(block_index).values;
    if size(values, 2) ~= n_time
        error('Play-bout matrices do not share a time axis.');
    end
    n_rows = size(values, 1);
    pieces{block_index} = values;
    centers(block_index) = row + (n_rows + 1) / 2;
    row = row + n_rows;
    if block_index < n_blocks
        boundaries(block_index) = row + 0.5;
    end
    labels{block_index} = block_label(blocks(block_index));
end
stack = vertcat(pieces{:});
end

function label = block_label(block)
if block.shared
    label = block.session_name;
else
    label = sprintf('%s | %s', block.animal_name, block.session_name);
end
end

function [mu, lo, hi] = mean_with_ci(X, ci_z)
mu = mean(X, 1, 'omitnan');
n = sum(isfinite(X), 1);
sigma = std(X, 0, 1, 'omitnan');
sem = sigma ./ sqrt(n);
sem(n < 2) = NaN;
lo = mu - ci_z .* sem;
hi = mu + ci_z .* sem;
end

function ci_handle = fill_confidence_interval(ax, time, lo, hi)
ci_handle = gobjects(0);
valid = isfinite(time) & isfinite(lo) & isfinite(hi);
if ~any(valid)
    return
end
edges = diff([false, valid(:)', false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
color = [0.73 0.78 0.88];
first = true;
for segment_index = 1:numel(starts)
    idx = starts(segment_index):stops(segment_index);
    if numel(idx) < 2
        continue
    end
    tt = time(idx);
    patch_handle = fill(ax, [tt, fliplr(tt)], [lo(idx), fliplr(hi(idx))], color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9);
    if first
        patch_handle.DisplayName = '95% CI';
        ci_handle = patch_handle;
        first = false;
    else
        patch_handle.HandleVisibility = 'off';
    end
end
end

function cmap = diverging_blue_red(n_colors)
x = linspace(0, 1, n_colors)';
cmap = zeros(n_colors, 3);
left = x <= 0.5;
if any(left)
    fraction = x(left) / 0.5;
    cmap(left, :) = [fraction, fraction, ones(sum(left), 1)];
end
right = ~left;
if any(right)
    fraction = (x(right) - 0.5) / 0.5;
    cmap(right, :) = [ones(sum(right), 1), 1 - fraction, 1 - fraction];
end
end

function font_size = tick_font_size(n_blocks)
if n_blocks > 18
    font_size = 6;
elseif n_blocks > 10
    font_size = 7;
else
    font_size = 8;
end
end

function name = pretty_variable_name(variable_name)
switch lower(char(variable_name))
    case 'speed'
        name = 'Speed';
    case 'acceleration'
        name = 'Acceleration';
    case 'distance'
        name = 'Inter-animal distance';
    otherwise
        name = char(variable_name);
end
end

function name = file_token(variable_name)
name = regexprep(char(variable_name), '[^A-Za-z0-9_\-]', '_');
end

function save_figure(fig, file_path)
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, file_path, 'Resolution', 200);
else
    saveas(fig, file_path);
end
end
