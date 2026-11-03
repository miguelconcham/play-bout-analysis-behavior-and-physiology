function traces = extract_delta_phase_vs_lpag_ref(pli_struct, animal_names)
%EXTRACT_DELTA_PHASE_VS_LPAG_REF  Delta phase vs mid-LPAG from vector saves.
%
% Companion to extract_delta_phase_vs_lpag_distance.m for the faster
% GENERATE_COHERENCE_MAPS_LPAG_REF output: mean_phaseLag_freq1 is already
% an nCh x 1 vector versus that probe's LPAG-center site (not a pairwise
% matrix). Do not pass full GENERATE_COHERNECE_MAPS_STRUCT maps here.
%
% Channel maps are stored in the shank / chanMap order used at Estimate
% time, which is not anatomical. Sort by y then x before unwrapping and
% locating region borders, same as the pairwise Analyze path.

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

for session_n = 1:numel(pli_struct)
    entry = pli_struct(session_n);
    session_name = local_session_name(animal_names, session_n);
    animal_id = local_animal_id(session_name);

    if ~isfield(entry, 'mean_phaseLag_freq1')
        warning('extract_delta_phase_vs_lpag_ref:NoDeltaPhase', ...
            '%s: no mean_phaseLag_freq1 — skipping', session_name);
        continue
    end
    if ~isfield(entry, 'channel_map') || ~isfield(entry, 'areas_by_channel')
        warning('extract_delta_phase_vs_lpag_ref:NoMap', ...
            '%s: no channel map / area labels — skipping', session_name);
        continue
    end
    if ~isfield(entry, 'lpag_center_index') || isempty(entry.lpag_center_index)
        warning('extract_delta_phase_vs_lpag_ref:NoRef', ...
            '%s: no lpag_center_index — skipping', session_name);
        continue
    end

    phase = entry.mean_phaseLag_freq1;
    if ~isvector(phase) || (size(phase, 1) > 1 && size(phase, 2) > 1)
        warning('extract_delta_phase_vs_lpag_ref:PairwiseMatrix', ...
            ['%s: mean_phaseLag_freq1 is a pairwise matrix. ', ...
            'Use extract_delta_phase_vs_lpag_distance / Analyze_coherence_maps.'], ...
            session_name);
        continue
    end
    phase = phase(:);
    channel_map = entry.channel_map;
    areas = entry.areas_by_channel;
    n_chan = min([numel(phase), size(channel_map, 1), numel(areas)]);
    if n_chan < 2
        continue
    end

    x = channel_map(1:n_chan, 1);
    y = channel_map(1:n_chan, 2);
    phase = phase(1:n_chan);
    valid = isfinite(x) & isfinite(y) & isfinite(phase);
    ref = entry.lpag_center_index;
    if ref < 1 || ref > n_chan || ~valid(ref)
        warning('extract_delta_phase_vs_lpag_ref:BadRef', ...
            '%s: LPAG-center index %g is not a valid channel — skipping', ...
            session_name, ref);
        continue
    end

    probe_idx = find(valid);
    probe_idx = local_sort_by_probe_depth(probe_idx, x, y);
    if numel(probe_idx) < 2
        continue
    end

    area_names = cell(n_chan, 1);
    for ch = 1:n_chan
        area_names{ch} = local_channel_area(areas, ch);
    end
    area = area_names(probe_idx);
    phase_rad = angle(exp(1i * phase(probe_idx)));
    distance_um = y(probe_idx) - y(ref);

    [border_distance_um, border_from_area, border_to_area] = ...
        local_region_borders(distance_um, area);

    if isfield(entry, 'probe_num') && ~isempty(entry.probe_num)
        probe_n = entry.probe_num;
    else
        probe_n = 1;
    end
    if isfield(entry, 'npx_type') && entry.npx_type == 2
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
if isempty(probe_idx)
    return
end
coords = [y(probe_idx), x(probe_idx)];
[~, order] = sortrows(coords, [1 2]);
probe_idx = probe_idx(order);
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
