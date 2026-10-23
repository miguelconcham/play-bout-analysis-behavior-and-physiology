%% Analyze_coherence_maps
% Delta phase difference versus distance from the LPAG-center channel.
%
% Rebuild with Estimate_coherence_maps.m (GENERATE_COHERNECE_MAPS_STRUCT).
% That Estimate stores pairwise delta (freq1, 1–5 Hz) and theta (freq2)
% PLI / wPLI / phase-lag matrices and saves coherence_pli_struct_2026.mat.
% This script does not re-estimate and does not restitch the old
% PAG/mPFC/NPX2 satellite files: it loads the newest merged / all-together
% / 2026 pairwise save. Do not load *lpag_ref* files here; those are
% plotted by Analyze_coherence_maps_lpag_ref.m.
%
% Plot 1 — every session is one line:
%   x = signed distance from that session's LPAG center (µm)
%   y = delta phase difference vs that LPAG center (radians)
%
% Plot 2 — one figure per animal, sessions overlaid, with vertical lines at
% region borders (same x/y as plot 1).
%
% Plot 3 — scatter of every stacked point, colored by area, with marker
% face blending. Only DR / PAG subregions / SupCol are colored. Points far
% from each session's phase-vs-distance trend are removed first.
%
% Plot 4 — mean and 95% CI of those colored areas after interpolating each
% session onto a shared distance axis.

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_coherence_maps.m';
end
repo_root = fileparts(fileparts(fileparts(this_file)));
addpath(fullfile(repo_root, 'Custom functions'));
addpath(fullfile(repo_root, 'Figure 2', 'Supporting codes v2'));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');
figure_folder = fullfile(saving_folder, 'delta phase vs lpag distance');
if exist(figure_folder, 'dir') ~= 7
    mkdir(figure_folder);
end

%% Load the newest merged / all-together / 2026 pairwise coherence structure
% Dual-animal PAG/mPFC slots were mismatched once and patched by stitching
% satellite mats. The combined result is saved as
% coherence_pli_struct_2026.mat (or an older *merged* / *all_together*
% file). Do not reload mpfc / B1D11007 / NPX2 pieces, and do not load
% *lpag_ref* vector saves here.
disp('loading coherence maps')
[structure_file, names_file] = find_latest_merged_coherence_files(saving_folder);
load(structure_file, 'pli_struct');
if ~isempty(names_file) && exist(names_file, 'file') == 2
    load(names_file, 'animal_names');
else
    animal_names = cell(numel(pli_struct), 1);
    for session_n = 1:numel(pli_struct)
        animal_names{session_n, 1} = sprintf('session %d', session_n);
    end
end
disp('ready')

%% One line per session / LPAG-containing probe
traces = extract_delta_phase_vs_lpag_distance(pli_struct, animal_names);
if isempty(traces)
    error('Analyze_coherence_maps:NoLPAG', ...
        'No LPAG-center traces found in the saved coherence maps.');
end
fprintf('Plotting %d LPAG-referenced delta phase traces\n', numel(traces));

%% Plot 1 — all sessions
session_fig = figure('Name', 'Delta phase vs distance from LPAG center', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.08 0.12 0.72 0.72]);
hold on
colors = generateDistinctColors(numel(traces));
for trace_n = 1:numel(traces)
    [distance_um, phase_rad] = local_sorted_trace(traces(trace_n));
    plot(distance_um, phase_rad, ...
        'Color', colors(trace_n, :), 'LineWidth', 1.2, ...
        'DisplayName', traces(trace_n).line_label)
end
xline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
yline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
xlabel('Distance from LPAG center (\mum)')
ylabel('Phase difference vs LPAG center (rad)')
title('Delta phase difference vs LPAG center')
legend('Location', 'eastoutside', 'Interpreter', 'none')
local_set_radian_ticks(gca)
box off
local_save_figure(session_fig, fullfile(figure_folder, 'all_sessions.png'));

%% Plot 2 — one figure per animal, sessions pooled, region borders
animal_ids = unique({traces.animal_id}, 'stable');
stack_table = table(zeros(0, 1), zeros(0, 1), cell(0, 1), cell(0, 1), cell(0, 1), zeros(0, 1), ...
    'VariableNames', {'distance_um', 'phase_rad', 'area', 'animal', 'session', 'probe'});
for animal_n = 1:numel(animal_ids)
    animal_id = animal_ids{animal_n};
    these = traces(strcmp({traces.animal_id}, animal_id));
    animal_fig = figure('Name', ['Delta phase vs distance — ', animal_id], ...
        'Color', 'w', 'Units', 'normalized', 'Position', [0.12 0.16 0.62 0.66]);
    hold on
    animal_colors = generateDistinctColors(numel(these));
    for trace_n = 1:numel(these)
        [distance_um, phase_rad, area] = local_sorted_trace(these(trace_n));
        plot(distance_um, phase_rad, ...
            'Color', animal_colors(trace_n, :), 'LineWidth', 1.2, ...
            'DisplayName', these(trace_n).line_label)
        n_row = numel(distance_um);
        animal = repmat({char(these(trace_n).animal_id)}, n_row, 1);
        session = repmat({char(these(trace_n).session_name)}, n_row, 1);
        if isfield(these(trace_n), 'probe_num') && ~isempty(these(trace_n).probe_num)
            probe = repmat(these(trace_n).probe_num, n_row, 1);
        else
            probe = ones(n_row, 1);
        end
        stack_table = [stack_table; table(distance_um, phase_rad, area, animal, session, probe, ...
            'VariableNames', {'distance_um', 'phase_rad', 'area', 'animal', 'session', 'probe'})]; %#ok<AGROW>
    end
    xline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
    yline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')

    [border_distance, border_label] = local_pooled_borders(these);
    for border_n = 1:numel(border_distance)
        xline(border_distance(border_n), '--', border_label{border_n}, ...
            'Color', [0.35 0.35 0.35], ...
            'LabelOrientation', 'aligned', ...
            'LabelHorizontalAlignment', 'center', ...
            'HandleVisibility', 'off')
    end

    xlabel('Distance from LPAG center (\mum)')
    ylabel('Phase difference vs LPAG center (rad)')
    title(['Delta phase difference vs LPAG center — ', animal_id])
    legend('Location', 'eastoutside', 'Interpreter', 'none')
    local_set_radian_ticks(gca)
    box off
    file_stem = regexprep(animal_id, '[^A-Za-z0-9]+', '_');
    local_save_figure(animal_fig, fullfile(figure_folder, [file_stem, '.png']));
end

%% Drop per-session points far from that recording's phase-vs-distance trend
stack_table = local_remove_session_trend_outliers(stack_table);

%% Plot 3 — all stacked points, colored by area
areas_in_color = {'DR', 'VLPAG', 'LPAG', 'DLPAG', 'SupCol'};
scatter_fig = figure('Name', 'Delta phase vs distance by area', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.12 0.68 0.7]);
hold on
area_labels = stack_table.area;
for k = 1:numel(area_labels)
    if isempty(area_labels{k})
        area_labels{k} = '(unlabeled)';
    end
end
area_names = unique(area_labels, 'stable');
area_colors = generateDistinctColors(numel(area_names));
index = ~ismember(area_names, areas_in_color);
area_colors(index, :) = 0.72;
for area_n = 1:numel(area_names)
    mask = strcmp(area_labels, area_names{area_n}) & ...
        isfinite(stack_table.distance_um) & isfinite(stack_table.phase_rad);
    if index(area_n)
        scatter(stack_table.distance_um(mask), stack_table.phase_rad(mask), ...
            16, area_colors(area_n, :), 'filled', ...
            'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.05, ...
            'HandleVisibility', 'off')
    else
        scatter(stack_table.distance_um(mask), stack_table.phase_rad(mask), ...
            22, area_colors(area_n, :), 'filled', ...
            'MarkerFaceAlpha', 0.28, 'MarkerEdgeAlpha', 0.12, ...
            'DisplayName', area_names{area_n})
    end
end
xline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
yline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
xlabel('Distance from LPAG center (\mum)')
ylabel('Phase difference vs LPAG center (rad)')
title('Delta phase vs distance, colored by area')
legend('Location', 'eastoutside', 'Interpreter', 'none')
local_set_radian_ticks(gca)
box off
local_save_figure(scatter_fig, fullfile(figure_folder, 'all_points_by_area.png'));

%% Plot 4 — session-interpolated mean and 95% CI for colored areas
colored_mask = ismember(stack_table.area, areas_in_color) & ...
    isfinite(stack_table.distance_um) & isfinite(stack_table.phase_rad);
if ~any(colored_mask)
    warning('Analyze_coherence_maps:NoColoredAreas', ...
        'No DR / PAG / SupCol points for the interpolated CI plot.');
else
    step_um = 5;
    dmin = min(stack_table.distance_um(colored_mask));
    dmax = max(stack_table.distance_um(colored_mask));
    distance_axis = (floor(dmin / step_um) * step_um : step_um : ...
        ceil(dmax / step_um) * step_um)';
    colored_names = area_names(~index);
    colored_colors = area_colors(~index, :);

    ci_fig = figure('Name', 'Delta phase vs distance, area CI', ...
        'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.12 0.68 0.7]);
    hold on
    y_scatter = [];
    y_ci = [];
    for area_n = 1:numel(colored_names)
        mask = strcmp(stack_table.area, colored_names{area_n}) & colored_mask;
        scatter(stack_table.distance_um(mask), stack_table.phase_rad(mask), ...
            14, colored_colors(area_n, :), 'filled', ...
            'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.05, ...
            'HandleVisibility', 'off')
        [mean_phase, ci_lo, ci_hi] = local_session_interpolated_ci( ...
            stack_table, colored_names{area_n}, distance_axis);
        local_fill_ci(distance_axis, ci_lo, ci_hi, colored_colors(area_n, :));
        plot(distance_axis, mean_phase, '-', ...
            'Color', colored_colors(area_n, :), 'LineWidth', 1.8, ...
            'DisplayName', colored_names{area_n})
        y_scatter = [y_scatter; stack_table.phase_rad(mask)]; %#ok<AGROW>
        y_ci = [y_ci; mean_phase(:); ci_lo(:); ci_hi(:)]; %#ok<AGROW>
    end
    xline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
    yline(0, '-', 'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off')
    xlabel('Distance from LPAG center (\mum)')
    ylabel('Phase difference vs LPAG center (rad)')
    title('Delta phase vs distance, 95% CI by area')
    legend('Location', 'eastoutside', 'Interpreter', 'none')
    local_set_observed_radian_ylim(gca, y_scatter, y_ci)
    box off
    local_save_figure(ci_fig, fullfile(figure_folder, 'area_interpolated_ci.png'));
end

%% -------------------------------------------------------------------------
function [distance_um, phase_rad, area] = local_sorted_trace(trace)
% Plot in channel-map depth order, not LFP-row order.
[distance_um, order] = sort(trace.distance_um);
distance_um = distance_um(:);
phase_rad = trace.phase_unwrapped(order);
phase_rad = phase_rad(:);
if isfield(trace, 'area') && numel(trace.area) >= numel(order)
    area = trace.area(order);
    area = area(:);
else
    area = repmat({''}, numel(order), 1);
end
if isstring(area)
    area = cellstr(area);
end
end

function stack_table = local_remove_session_trend_outliers(stack_table)
% Per session (and probe), fit a robust trend of phase vs distance and drop
% points whose residual is more than 3 scaled MAD from that trend, with a
% 0.5 rad floor so isolated jumps are still removed on a clean trace.
n_row = height(stack_table);
keep = true(n_row, 1);
session = stack_table.session;
if ismember('probe', stack_table.Properties.VariableNames)
    probe_str = arrayfun(@(p) num2str(p), stack_table.probe, 'UniformOutput', false);
    group = strcat(session, '::', probe_str);
else
    group = session;
end
group_ids = unique(group, 'stable');
for g = 1:numel(group_ids)
    idx = find(strcmp(group, group_ids{g}) & ...
        isfinite(stack_table.distance_um) & isfinite(stack_table.phase_rad));
    if numel(idx) < 8
        continue
    end
    distance_um = stack_table.distance_um(idx);
    phase_rad = stack_table.phase_rad(idx);
    [distance_um, order] = sort(distance_um);
    idx = idx(order);
    phase_rad = phase_rad(order);
    span = max(7, 2 * floor(numel(phase_rad) / 8) + 1);
    try
        trend = smooth(distance_um, phase_rad, span, 'rlowess');
    catch
        trend = movmedian(phase_rad, span);
    end
    residual = phase_rad - trend;
    med_res = median(residual);
    scaled_mad = 1.4826 * mad(residual, 1);
    threshold = max(3 * scaled_mad, 0.5);
    outlier = abs(residual - med_res) > threshold;
    keep(idx(outlier)) = false;
end
n_removed = sum(~keep);
fprintf('Removed %d / %d points far from the per-session phase-distance trend\n', ...
    n_removed, n_row);
stack_table = stack_table(keep, :);
end

function [mean_phase, ci_lo, ci_hi] = local_session_interpolated_ci( ...
    stack_table, area_name, distance_axis)
% Interpolate each session onto distance_axis, then mean ± 1.96 SEM.
sessions = unique(stack_table.session, 'stable');
phase_mat = nan(numel(distance_axis), numel(sessions));
for session_n = 1:numel(sessions)
    mask = strcmp(stack_table.session, sessions{session_n}) & ...
        strcmp(stack_table.area, area_name) & ...
        isfinite(stack_table.distance_um) & isfinite(stack_table.phase_rad);
    if sum(mask) < 2
        continue
    end
    distance_um = stack_table.distance_um(mask);
    phase_rad = stack_table.phase_rad(mask);
    [distance_um, order] = sort(distance_um);
    phase_rad = phase_rad(order);
    [distance_um, first] = unique(distance_um, 'first');
    phase_rad = phase_rad(first);
    if numel(distance_um) < 2
        continue
    end
    in_range = distance_axis >= distance_um(1) & distance_axis <= distance_um(end);
    phase_mat(in_range, session_n) = interp1(distance_um, phase_rad, ...
        distance_axis(in_range), 'linear');
end
n = sum(isfinite(phase_mat), 2);
mean_phase = nan(size(distance_axis));
ci_lo = nan(size(distance_axis));
ci_hi = nan(size(distance_axis));
enough = n >= 2;
if ~any(enough)
    return
end
mean_phase(enough) = mean(phase_mat(enough, :), 2, 'omitnan');
sem = std(phase_mat(enough, :), 0, 2, 'omitnan') ./ sqrt(n(enough));
half = 1.96 * sem;
ci_lo(enough) = mean_phase(enough) - half;
ci_hi(enough) = mean_phase(enough) + half;
end

function local_fill_ci(x, y_lo, y_hi, color)
ok = isfinite(x) & isfinite(y_lo) & isfinite(y_hi);
d = diff([false; ok(:); false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
for k = 1:numel(starts)
    sl = starts(k):ends(k);
    if numel(sl) < 2
        continue
    end
    fill([x(sl); flipud(x(sl))], [y_lo(sl); flipud(y_hi(sl))], color, ...
        'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off')
end
end

function [border_distance, border_label] = local_pooled_borders(these)
bin_um = 10;
buckets = containers.Map('KeyType', 'double', 'ValueType', 'any');
for trace_n = 1:numel(these)
    for border_n = 1:numel(these(trace_n).border_distance_um)
        distance = these(trace_n).border_distance_um(border_n);
        if ~isfinite(distance)
            continue
        end
        key = round(distance / bin_um) * bin_um;
        pair = {these(trace_n).border_from_area{border_n}, ...
            these(trace_n).border_to_area{border_n}};
        if isKey(buckets, key)
            bucket_pairs = buckets(key);
            bucket_pairs{end + 1} = pair; %#ok<AGROW>
            buckets(key) = bucket_pairs;
        else
            buckets(key) = {pair};
        end
    end
end
if buckets.Count == 0
    border_distance = [];
    border_label = {};
    return
end
keys_list = cell2mat(buckets.keys);
keys_list = sort(keys_list);
border_distance = zeros(numel(keys_list), 1);
border_label = cell(numel(keys_list), 1);
for k = 1:numel(keys_list)
    pairs = buckets(keys_list(k));
    counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for p = 1:numel(pairs)
        label = [pairs{p}{1}, ' | ', pairs{p}{2}];
        if isKey(counts, label)
            counts(label) = counts(label) + 1;
        else
            counts(label) = 1;
        end
    end
    labels = counts.keys;
    votes = cell2mat(counts.values);
    [~, winner] = max(votes);
    border_distance(k) = keys_list(k);
    border_label{k} = labels{winner};
end
end

function local_set_observed_radian_ylim(ax, y_scatter, y_ci)
% Fit y to the plotted CI and the bulk of the scatter (2–98th percentile).
y_scatter = y_scatter(isfinite(y_scatter));
y_ci = y_ci(isfinite(y_ci));
if isempty(y_scatter) && isempty(y_ci)
    return
end
if isempty(y_ci)
    low = min(y_scatter);
    high = max(y_scatter);
else
    low = min(y_ci);
    high = max(y_ci);
end
if ~isempty(y_scatter)
    if numel(y_scatter) >= 20
        edges = prctile(y_scatter, [2 98]);
    else
        edges = [min(y_scatter), max(y_scatter)];
    end
    low = min(low, edges(1));
    high = max(high, edges(2));
end
if ~isfinite(low) || ~isfinite(high) || high <= low
    return
end
pad = max(0.08 * (high - low), 0.05);
ylim(ax, [low - pad, high + pad])
local_set_radian_ticks(ax, 'tight')
end

function local_set_radian_ticks(ax, mode)
if nargin < 2 || isempty(mode)
    mode = 'atleast_pi';
end
ylim_now = ylim(ax);
if strcmp(mode, 'tight')
    low = ylim_now(1);
    high = ylim_now(2);
else
    low = min(-pi, ylim_now(1));
    high = max(pi, ylim_now(2));
    ylim(ax, [low, high])
end
tick = unique([-2*pi, -3*pi/2, -pi, -3*pi/4, -pi/2, -pi/4, 0, ...
    pi/4, pi/2, 3*pi/4, pi, 3*pi/2, 2*pi]);
tick = tick(tick >= low - 1e-9 & tick <= high + 1e-9);
if numel(tick) < 2
    return
end
yticks(ax, tick)
yticklabels(ax, arrayfun(@local_radian_label, tick, 'UniformOutput', false))
end

function label = local_radian_label(value)
if abs(value) < 1e-9
    label = '0';
elseif abs(value - pi) < 1e-9
    label = '\pi';
elseif abs(value + pi) < 1e-9
    label = '-\pi';
elseif abs(value - pi/2) < 1e-9
    label = '\pi/2';
elseif abs(value + pi/2) < 1e-9
    label = '-\pi/2';
elseif abs(value - pi/4) < 1e-9
    label = '\pi/4';
elseif abs(value + pi/4) < 1e-9
    label = '-\pi/4';
elseif abs(value - 3*pi/4) < 1e-9
    label = '3\pi/4';
elseif abs(value + 3*pi/4) < 1e-9
    label = '-3\pi/4';
elseif abs(value - 3*pi/2) < 1e-9
    label = '3\pi/2';
elseif abs(value + 3*pi/2) < 1e-9
    label = '-3\pi/2';
elseif abs(value - 2*pi) < 1e-9
    label = '2\pi';
elseif abs(value + 2*pi) < 1e-9
    label = '-2\pi';
else
    label = sprintf('%.2f', value);
end
end

function local_save_figure(fig, path)
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, path, 'Resolution', 200)
else
    saveas(fig, path)
end
disp(['saved ', path])
end
