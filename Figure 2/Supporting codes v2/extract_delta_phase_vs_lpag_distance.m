function traces = extract_delta_phase_vs_lpag_distance(pli_struct, animal_names, phase_field)
%EXTRACT_DELTA_PHASE_VS_LPAG_DISTANCE  Delta phase vs LPAG-center distance.
%
% For every saved coherence-map session, take the pairwise phase-lag
% matrix and read one column: phase of each same-probe channel relative
% to that probe's LPAG-center site.
%
% phase_field selects the matrix. Default is 'mean_phaseLag_freq1'
% (delta, 1-5 Hz). Pass 'mean_phaseLag_freq2' for theta (6-12 Hz).
%
% LPAG center is the same geometry as Fig 2 PSTHs / GENERATE_PLAY_SESSION:
% mean (x, y) of LPAG sites on the probe, then the closest LPAG site.
% Sessions or probes without LPAG are skipped.
%
% Distance is signed probe y (µm) minus the LPAG-center y. Positive values
% are toward higher ycoords (typically dorsal).
%
% Channel maps are stored in LFP-row order (chanMap+1 in GENERATE), which
% is not anatomical. Same as the Estimate_coherence_maps heatmaps: sort
% each probe by ycoords (then x) before reading the phase-lag column,
% unwrapping, and drawing the line.
%
% traces(k) fields
%   session_name, animal_id, line_label, probe_num
%   lpag_center_index, lpag_center_y
%   channel_index, distance_um, phase_rad, phase_unwrapped, area
%   border_distance_um, border_from_area, border_to_area

if nargin < 3 || isempty(phase_field)
    phase_field = 'mean_phaseLag_freq1';
end
phase_field = char(phase_field);

if nargin < 2 || isempty(animal_names)
    animal_names = cell(numel(pli_struct), 1);
    for session_n = 1:numel(pli_struct)
        animal_names{session_n, 1} = sprintf('session %d', session_n);
    end
end

traces = struct( ...
    'session_name', {}, ...
    'animal_id', {}, ...
    'line_label', {}, ...
    'probe_num', {}, ...
    'lpag_center_index', {}, ...
    'lpag_center_y', {}, ...
    'channel_index', {}, ...
    'distance_um', {}, ...
    'phase_rad', {}, ...
    'phase_unwrapped', {}, ...
    'area', {}, ...
    'border_distance_um', {}, ...
    'border_from_area', {}, ...
    'border_to_area', {});

hard_coded_x_coords = [8 40; 258 290; 508 540; 758 790];

for session_n = 1:numel(pli_struct)
    entry = pli_struct(session_n);
    session_name = local_session_name(animal_names, session_n);
    animal_id = local_animal_id(session_name);

    if ~isfield(entry, phase_field)
        warning('extract_delta_phase_vs_lpag_distance:NoDeltaPhase', ...
            '%s: no %s — skipping', session_name, phase_field);
        continue
    end
    if ~isfield(entry, 'channel_map') || ~isfield(entry, 'areas_by_channel')
        warning('extract_delta_phase_vs_lpag_distance:NoMap', ...
            '%s: no channel map / area labels — skipping', session_name);
        continue
    end

    phase_matrix = entry.(phase_field);
    channel_map = entry.channel_map;
    areas = entry.areas_by_channel;
    n_chan = min([size(phase_matrix, 1), size(phase_matrix, 2), size(channel_map, 1), numel(areas)]);
    if n_chan < 2
        continue
    end

    x = channel_map(1:n_chan, 1);
    y = channel_map(1:n_chan, 2);
    valid = isfinite(x) & isfinite(y);
    if ~any(valid)
        continue
    end

    probe_num = local_probe_numbers(x, hard_coded_x_coords);
    area_names = cell(n_chan, 1);
    is_lpag = false(n_chan, 1);
    for ch = 1:n_chan
        area_names{ch} = local_channel_area(areas, ch);
        is_lpag(ch) = strcmp(area_names{ch}, 'LPAG') & valid(ch);
    end
    if ~any(is_lpag)
        continue
    end

    lpag_probes = unique(probe_num(is_lpag & isfinite(probe_num)));
    lpag_probes = lpag_probes(:)';
    n_lpag_probes = numel(lpag_probes);

    for probe_n = lpag_probes
        probe_mask = valid & (probe_num == probe_n);
        lpag_mask = probe_mask & is_lpag;
        if ~any(lpag_mask)
            continue
        end

        mean_x = mean(x(lpag_mask));
        mean_y = mean(y(lpag_mask));
        distance_to_mean = abs(x - mean_x) + abs(y - mean_y);
        distance_to_mean(~lpag_mask) = Inf;
        [~, ref] = min(distance_to_mean);
        if ~isfinite(distance_to_mean(ref))
            continue
        end

        probe_idx = find(probe_mask);
        probe_idx = local_sort_by_probe_depth(probe_idx, x, y);
        keep = isfinite(phase_matrix(probe_idx, ref));
        probe_idx = local_sort_by_probe_depth(probe_idx(keep), x, y);
        if numel(probe_idx) < 2
            continue
        end

        phase_rad = angle(exp(1i * phase_matrix(probe_idx, ref)));
        distance_um = y(probe_idx) - y(ref);
        area = area_names(probe_idx);

        [border_distance_um, border_from_area, border_to_area] = ...
            local_region_borders(distance_um, area);

        if n_lpag_probes > 1
            line_label = sprintf('%s probe %d', session_name, probe_n);
        else
            line_label = session_name;
        end

        trace = struct();
        trace.session_name = session_name;
        trace.animal_id = animal_id;
        trace.line_label = line_label;
        trace.probe_num = probe_n;
        trace.lpag_center_index = ref;
        trace.lpag_center_y = y(ref);
        trace.channel_index = probe_idx(:);
        trace.distance_um = distance_um(:);
        trace.phase_rad = phase_rad(:);
        trace.phase_unwrapped = unwrap(phase_rad(:));
        trace.area = area(:);
        trace.border_distance_um = border_distance_um(:);
        trace.border_from_area = border_from_area(:);
        trace.border_to_area = border_to_area(:);
        traces(end + 1) = trace; %#ok<AGROW>
    end
end
end

function session_name = local_session_name(animal_names, session_n)
if iscell(animal_names)
    if session_n > size(animal_names, 1)
        session_name = sprintf('session %d', session_n);
        return
    end
    session_name = animal_names{session_n, 1};
elseif isstring(animal_names) || ischar(animal_names)
    names = cellstr(animal_names);
    if session_n > numel(names)
        session_name = sprintf('session %d', session_n);
        return
    end
    session_name = names{session_n};
else
    session_name = sprintf('session %d', session_n);
    return
end
if isempty(session_name)
    session_name = sprintf('session %d', session_n);
else
    session_name = strtrim(char(session_name));
end
end

function animal_id = local_animal_id(session_name)
tokens = strsplit(strtrim(char(session_name)));
if isempty(tokens) || isempty(tokens{1})
    animal_id = session_name;
else
    animal_id = tokens{1};
end
end

function probe_idx = local_sort_by_probe_depth(probe_idx, x, y)
% Anatomical order along the shank: increasing y, then x at a given depth.
if isempty(probe_idx)
    return
end
coords = [y(probe_idx), x(probe_idx)];
[~, order] = sortrows(coords, [1 2]);
probe_idx = probe_idx(order);
end

function probe_num = local_probe_numbers(x, hard_coded_x_coords)
probe_num = nan(size(x));
is_npx2 = false;
for probe_n = 1:size(hard_coded_x_coords, 1)
    on_probe = ismember(x, hard_coded_x_coords(probe_n, :));
    if any(on_probe)
        is_npx2 = true;
        probe_num(on_probe) = probe_n;
    end
end
if ~is_npx2
    probe_num(isfinite(x)) = 1;
end
end

function label = local_channel_area(areas, ch)
label = '';
if ch > numel(areas)
    return
end
value = areas{ch};
if isempty(value)
    return
end
if iscell(value)
    if isempty(value)
        return
    end
    value = value{1};
end
if isstring(value) || ischar(value)
    label = strtrim(char(value));
end
end

function [border_distance_um, border_from_area, border_to_area] = ...
    local_region_borders(distance_um, area)
border_distance_um = [];
border_from_area = {};
border_to_area = {};
for ch = 2:numel(distance_um)
    from_area = area{ch - 1};
    to_area = area{ch};
    if isempty(from_area) || isempty(to_area) || strcmp(from_area, to_area)
        continue
    end
    border_distance_um(end + 1, 1) = ...
        0.5 * (distance_um(ch - 1) + distance_um(ch)); %#ok<AGROW>
    border_from_area{end + 1, 1} = from_area; %#ok<AGROW>
    border_to_area{end + 1, 1} = to_area; %#ok<AGROW>
end
end
