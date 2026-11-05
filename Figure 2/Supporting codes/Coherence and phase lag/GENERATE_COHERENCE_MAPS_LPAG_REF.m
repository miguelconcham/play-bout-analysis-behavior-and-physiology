function pli_struct = GENERATE_COHERENCE_MAPS_LPAG_REF(current_dir, Hd_freq1, Hd_freq2, area2load)
%GENERATE_COHERENCE_MAPS_LPAG_REF  PLI / phase lag vs mid-LPAG only.
%
% Same geometry as Fig 2 PSTHs / GENERATE_PLAY_SESSION_LPAG_SPECTRA: one
% reference site at the center of LPAG on each probe that has LPAG. Only
% same-probe channels are compared to that reference, so this is much
% cheaper than GENERATE_COHERNECE_MAPS_STRUCT (no full channel-pair maps).
%
% Returns [] if LPAG is absent. Otherwise a struct array, one entry per
% LPAG-containing probe.

if nargin < 4 || isempty(area2load)
    area2load = 'PAG';
end

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
chan_map_folder = [data_root, '\NPX data\StarndarChannMap'];
area_limit_table = [data_root, '\Area_limits_GoodLooking.xlsx'];
hard_coded_x_coords = [8 40; 258 290; 508 540; 758 790];
sample_length_theta = 4200;
sample_length_delta = 15000;
n_samples = 500;
pctls = [1 5 95 99];

animal_code = strsplit(current_dir, {'\', '/'});
animal_code = animal_code(~cellfun(@isempty, animal_code));
animal_code = animal_code{end};
animal_code_params = strsplit(animal_code, ' ');
animal_batch = animal_code_params{1};
repeated_animal = animal_code_params{3};

%% Load LFP
disp(['LOADING LFP: ', animal_code])
if exist(fullfile(current_dir, ['LFP_', area2load, '.mat']), 'file') == 2
    NPX_Type = 2;
    load(fullfile(current_dir, ['LFP_', area2load, '.mat']), 'LFP');
elseif exist(fullfile(current_dir, ['LFP_', area2load, '.dat']), 'file') == 2
    NPX_Type = 1;
    file_pointer = fopen(fullfile(current_dir, ['LFP_', area2load, '.dat']), 'r');
    LFP = fread(file_pointer, 'int16');
    fclose(file_pointer);
    LFP = reshape(LFP, 384, numel(LFP) / 384);
else
    warning('GENERATE_COHERENCE_MAPS_LPAG_REF:NoLFP', ...
        '%s: no LFP_%s file — skipping', animal_code, area2load);
    pli_struct = [];
    return
end
LFP = double(LFP);
disp('LFP LOADED')

%% Channel map
if NPX_Type == 1
    load(fullfile(chan_map_folder, 'neuropixPhase3A_kilosortChanMap.mat'), ...
        'xcoords', 'ycoords', 'chanMap');
    if ~ismember(192, chanMap)
        pos_191 = find(chanMap == 191, 1);
        pos_193 = find(chanMap == 193, 1);
        if isempty(pos_191) || isempty(pos_193) || pos_193 ~= pos_191 + 1
            warning('GENERATE_COHERENCE_MAPS_LPAG_REF:ChanMap', ...
                '%s: inconsistent NPX1 channel map — skipping', animal_code);
            pli_struct = [];
            return
        end
        xcoords = [xcoords; NaN];
        xcoords(pos_193 + 1:end) = xcoords(pos_193:end - 1);
        xcoords(pos_193) = 43;
        ycoords = [ycoords; NaN];
        ycoords(pos_193 + 1:end) = ycoords(pos_193:end - 1);
        ycoords(pos_193) = 1900;
        chanMap = [chanMap; NaN];
        chanMap(pos_193 + 1:end) = chanMap(pos_193:end - 1);
        chanMap(pos_193) = 192;
    end
else
    load(fullfile(current_dir, ['chann_map_', area2load, '.mat']), ...
        'xcoords', 'ycoords', 'chanMap');
end

area_limit = readtable(area_limit_table);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);
lpag_rows = area_limit(ismember(area_limit.area, 'LPAG'), :);
if isempty(lpag_rows)
    warning('GENERATE_COHERENCE_MAPS_LPAG_REF:NoLPAG', ...
        '%s: no LPAG in area limits — skipping', animal_code);
    pli_struct = [];
    return
end

jobs = local_lpag_jobs(NPX_Type, lpag_rows, xcoords, ycoords, chanMap, ...
    hard_coded_x_coords, size(LFP, 1));
if isempty(jobs)
    warning('GENERATE_COHERENCE_MAPS_LPAG_REF:NoLPAGChannels', ...
        '%s: LPAG listed but no matching channels — skipping', animal_code);
    pli_struct = [];
    return
end

pli_struct = struct([]);
for job_n = 1:numel(jobs)
    job = jobs(job_n);
    shank_lfp = LFP(job.lfp_rows, :);
    disp(sprintf('  %s probe %d: %d channels, LPAG-center LFP row %d', ...
        animal_code, job.probe_num, numel(job.lfp_rows), job.detection_lfp_row));

    freq1 = Compute_pli_wpli_segments(shank_lfp, Hd_freq1, n_samples, ...
        sample_length_delta, pctls, 'ReferenceChannel', job.ref_in_shank, ...
        'Verbose', true);
    freq2 = Compute_pli_wpli_segments(shank_lfp, Hd_freq2, n_samples, ...
        sample_length_theta, pctls, 'ReferenceChannel', job.ref_in_shank, ...
        'Verbose', true);

    areas = local_areas_for_shank(job, area_limit, NPX_Type, area2load);
    entry = struct();
    entry.animal_code = animal_code;
    entry.npx_type = NPX_Type;
    entry.probe_num = job.probe_num;
    entry.lpag_center_index = job.ref_in_shank;
    entry.lpag_center_lfp_row = job.detection_lfp_row;
    entry.lpag_center_y = job.y_coords(job.ref_in_shank);
    entry.channel_map = [job.x_coords(:), job.y_coords(:)];
    entry.areas_by_channel = areas;
    entry.lfp_rows = job.lfp_rows(:);

    entry.PLI_vs_lpag_freq1 = freq1.PLI_vs_ref;
    entry.wPLI_vs_lpag_freq1 = freq1.wPLI_vs_ref;
    entry.mean_phaseLag_freq1 = freq1.mean_phaseLag_vs_ref;
    entry.pctls_freq1 = freq1.pctls;
    entry.freq_range_freq1 = freq1.freq_range;
    entry.sample_length_freq1 = sample_length_delta;
    entry.mean_overlap_freq1 = freq1.mean_overlap;

    entry.PLI_vs_lpag_freq2 = freq2.PLI_vs_ref;
    entry.wPLI_vs_lpag_freq2 = freq2.wPLI_vs_ref;
    entry.mean_phaseLag_freq2 = freq2.mean_phaseLag_vs_ref;
    entry.pctls_freq2 = freq2.pctls;
    entry.freq_range_freq2 = freq2.freq_range;
    entry.sample_length_freq2 = sample_length_theta;
    entry.mean_overlap_freq2 = freq2.mean_overlap;

    if isempty(pli_struct)
        pli_struct = entry;
    else
        pli_struct(end + 1) = entry; %#ok<AGROW>
    end
end
end

function jobs = local_lpag_jobs(NPX_Type, lpag_rows, xcoords, ycoords, chanMap, hard_x, n_lfp)
jobs = struct('probe_num', {}, 'detection_lfp_row', {}, 'ref_in_shank', {}, ...
    'lfp_rows', {}, 'x_coords', {}, 'y_coords', {});

if NPX_Type == 1
    depth_start = min(lpag_rows.depth_start);
    depth_end = max(lpag_rows.depth_end);
    in_lpag = ycoords >= depth_start & ycoords <= depth_end;
    if ~any(in_lpag)
        return
    end
    mean_loc = mean([xcoords(in_lpag), ycoords(in_lpag)], 1);
    distance = sum(abs([xcoords, ycoords] - mean_loc), 2);
    distance(~in_lpag) = Inf;
    [~, closest] = min(distance);
    valid = isfinite(xcoords) & isfinite(ycoords) & isfinite(chanMap);
    lfp_rows = chanMap(valid) + 1;
    keep = lfp_rows >= 1 & lfp_rows <= n_lfp;
    jobs(1).probe_num = 1;
    jobs(1).lfp_rows = lfp_rows(keep);
    jobs(1).x_coords = xcoords(valid);
    jobs(1).x_coords = jobs(1).x_coords(keep);
    jobs(1).y_coords = ycoords(valid);
    jobs(1).y_coords = jobs(1).y_coords(keep);
    jobs(1).detection_lfp_row = chanMap(closest) + 1;
    jobs(1).ref_in_shank = find(jobs(1).lfp_rows == jobs(1).detection_lfp_row, 1);
    if isempty(jobs(1).ref_in_shank)
        jobs = struct('probe_num', {}, 'detection_lfp_row', {}, 'ref_in_shank', {}, ...
            'lfp_rows', {}, 'x_coords', {}, 'y_coords', {});
    end
    return
end

probe_list = unique(lpag_rows.ProbeNum);
for p = 1:numel(probe_list)
    probe_n = probe_list(p);
    if probe_n < 1 || probe_n > size(hard_x, 1)
        continue
    end
    rows = lpag_rows(lpag_rows.ProbeNum == probe_n, :);
    depth_start = min(rows.depth_start);
    depth_end = max(rows.depth_end);
    shank_mask = ismember(xcoords, hard_x(probe_n, :));
    lpag_mask = shank_mask & ycoords >= depth_start & ycoords <= depth_end;
    if ~any(lpag_mask)
        continue
    end
    mean_loc = mean([xcoords(lpag_mask), ycoords(lpag_mask)], 1);
    distance = sum(abs([xcoords, ycoords] - mean_loc), 2);
    distance(~lpag_mask) = Inf;
    [~, closest] = min(distance);
    if ~shank_mask(closest)
        shank_idx = find(shank_mask);
        d = sum(abs([xcoords(shank_idx), ycoords(shank_idx)] - mean_loc), 2);
        [~, ii] = min(d);
        closest = shank_idx(ii);
    end
    lfp_rows = chanMap(shank_mask) + 1;
    x_sh = xcoords(shank_mask);
    y_sh = ycoords(shank_mask);
    keep = isfinite(lfp_rows) & lfp_rows >= 1 & lfp_rows <= n_lfp;
    lfp_rows = lfp_rows(keep);
    x_sh = x_sh(keep);
    y_sh = y_sh(keep);
    det_row = chanMap(closest) + 1;
    ref_in_shank = find(lfp_rows == det_row, 1);
    if isempty(ref_in_shank)
        continue
    end
    jobs(end + 1).probe_num = probe_n; %#ok<AGROW>
    jobs(end).detection_lfp_row = det_row;
    jobs(end).ref_in_shank = ref_in_shank;
    jobs(end).lfp_rows = lfp_rows(:);
    jobs(end).x_coords = x_sh(:);
    jobs(end).y_coords = y_sh(:);
end
end

function areas = local_areas_for_shank(job, area_limit, NPX_Type, area2load)
n_ch = numel(job.lfp_rows);
areas = cell(n_ch, 1);
for ch = 1:n_ch
    y = job.y_coords(ch);
    if NPX_Type == 1
        match = y >= area_limit.depth_start & y < area_limit.depth_end + 1 & ...
            ismember(area_limit.Probe_Area, area2load);
    else
        match = y >= area_limit.depth_start & y < area_limit.depth_end + 1 & ...
            area_limit.ProbeNum == job.probe_num & ...
            ismember(area_limit.Probe_Area, area2load);
    end
    if ~any(match)
        areas{ch} = '';
        continue
    end
    names = area_limit.area(match);
    if iscell(names)
        areas{ch} = names{1};
    else
        areas{ch} = char(names(1));
    end
end
end
