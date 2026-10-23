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
stack_table = table(zeros(0, 1), zeros(0, 1), cell(0, 1), cell(0, 1), cell(0, 1), ...
    'VariableNames', {'distance_um', 'phase_rad', 'area', 'animal', 'session'});
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
        stack_table = [stack_table; table(distance_um, phase_rad, area, animal, session, ...
            'VariableNames', {'distance_um', 'phase_rad', 'area', 'animal', 'session'})]; %#ok<AGROW>
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

function local_set_radian_ticks(ax)
ylim_now = ylim(ax);
low = min(-pi, ylim_now(1));
high = max(pi, ylim_now(2));
ylim(ax, [low, high])
tick = unique([-2*pi, -pi, -pi/2, 0, pi/2, pi, 2*pi]);
tick = tick(tick >= low - 1e-9 & tick <= high + 1e-9);
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
