function results = GENERATE_DELTA_PEAK_TROUGH_CSD_LFP(npx_data_dir, Hd_delta, t1, t2)
%GENERATE_DELTA_PEAK_TROUGH_CSD_LFP  Raw LFP snippets around delta extrema.
%
% Detect delta peaks/troughs on the LPAG-center channel during Partners
% session intervals, cut raw (unreferenced) LFP on [t1, t2], and return
% one result per shank that has LPAG:
%   - NPX1: one shank (full channel count)
%   - NPX2: one result per ProbeNum with LPAG (reduced channel list)
%
% If LPAG is absent, results is empty.
%
% Inputs
%   npx_data_dir  session folder under Data/NPX data/NPX raw data
%   Hd_delta      FIR bandpass filter object (designfilt), e.g. 1–5 Hz
%   t1, t2        seconds around extremum (t1 < 0 < t2), e.g. -0.25, 0.25
%
% Output
%   results       struct array with fields for peak/trough stacks + coords

data_root       = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
synch_directory = fullfile(data_root, 'Synch data');
area_limit_table = fullfile(data_root, 'Area_limits_GoodLooking.xlsx');
behavior_data   = fullfile(data_root, 'Behavior backups');
area2analyse    = 'PAG';
sr_LFP          = 2500;
hard_coded_x_coords_NPX2 = [8 40; 258 290; 508 540; 758 790];

if t1 >= 0 || t2 <= 0
    error('GENERATE_DELTA_PEAK_TROUGH_CSD_LFP:window', 'Need t1 < 0 < t2');
end

animal_code        = strsplit(npx_data_dir, '\');
animal_code        = animal_code{end};
animal_code_params = strsplit(animal_code, ' ');
animal_batch       = animal_code_params{1};
repeated_animal    = animal_code_params{3};

%% Behavior + Partners session (play session) intervals
load(fullfile(synch_directory, animal_code, 'synch_model_video2NPX.mat'), 'synch_model_video2NPX');

Behavior_file = fullfile(behavior_data, [animal_code, '.txt']);
Behavior = readtable(Behavior_file);
Behavior(:, 2) = [];
Behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
Behavior.Type2 = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'}))   = {'Pounce'};
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
Behavior.Type2(ismember(Behavior.Type2, ''))                           = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'), :) = [];
Behavior.Start = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End   = predict(synch_model_video2NPX, Behavior.End);

play_session = Behavior(strcmp(Behavior.Type2, 'Partners session'), :);
play_session = play_session(:, {'Start', 'End'});
if isempty(play_session)
    warning('GENERATE_DELTA_PEAK_TROUGH_CSD_LFP:NoPlaySession', ...
        '%s: no Partners session intervals — skipping', animal_code);
    results = [];
    return
end

%% Channel map + NPX type
load(fullfile(npx_data_dir, ['chann_map_', area2analyse, '.mat']), 'chanMap', 'xcoords', 'ycoords');
if any(ismember(xcoords, hard_coded_x_coords_NPX2))
    NPX_Type = 2;
else
    NPX_Type = 1;
    if ~ismember(192, chanMap)
        pos_191 = find(chanMap == 191);
        pos_193 = find(chanMap == 193);
        if isempty(pos_191) || isempty(pos_193) || pos_193 ~= pos_191 + 1
            warning('GENERATE_DELTA_PEAK_TROUGH_CSD_LFP:ChanMap', ...
                '%s: inconsistent NPX1 channel map — skipping', animal_code);
            results = [];
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
end

area_limit = readtable(area_limit_table);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);
lpag_rows = area_limit(ismember(area_limit.area, {'LPAG'}), :);
if isempty(lpag_rows)
    warning('GENERATE_DELTA_PEAK_TROUGH_CSD_LFP:NoLPAG', ...
        '%s: no LPAG in area limits — skipping', animal_code);
    results = [];
    return
end

%% Load LFP (no re-referencing)
disp(['LOADING LFP: ', animal_code])
if NPX_Type == 2
    load(fullfile(npx_data_dir, ['LFP_', area2analyse, '.mat']), 'LFP');
else
    fid = fopen(fullfile(npx_data_dir, ['LFP_', area2analyse, '.dat']), 'r');
    LFP = fread(fid, 'int16');
    fclose(fid);
    LFP = reshape(LFP, 384, numel(LFP) / 384);
end
LFP = double(LFP);
n_time_total = size(LFP, 2);
lfp_time = (1:n_time_total) / sr_LFP;
disp('LFP LOADED')

pre_samp  = round(abs(t1) * sr_LFP);
post_samp = round(t2 * sr_LFP);
n_time_win = pre_samp + post_samp;
T0 = pre_samp; % 0-based sample index of the extremum inside each snippet
% MATLAB slice peak-pre_samp : peak+post_samp-1 → length n_time_win, peak at pre_samp+1 (1-based)

min_peak_dist = sr_LFP / Hd_delta.CutoffFrequency2;

%% Build per-shank jobs (only shanks with LPAG)
jobs = local_build_shank_jobs(NPX_Type, lpag_rows, xcoords, ycoords, chanMap, hard_coded_x_coords_NPX2);
if isempty(jobs)
    warning('GENERATE_DELTA_PEAK_TROUGH_CSD_LFP:NoLPAGChannels', ...
        '%s: LPAG listed but no matching channels — skipping', animal_code);
    results = [];
    return
end

results = struct([]);
for j = 1:numel(jobs)
    job = jobs(j);
    det_ch = job.detection_channel;
    shank_rows = job.lfp_rows;          % rows into LFP
    x_sh = job.x_coords(:);
    y_sh = job.y_coords(:);

    % Sort shank channels by depth for a stable CSD axis
    [y_sh, ord] = sort(y_sh);
    x_sh = x_sh(ord);
    shank_rows = shank_rows(ord);
    if ~ismember(det_ch, shank_rows)
        warning('GENERATE_DELTA_PEAK_TROUGH_CSD_LFP:DetNotOnShank', ...
            '%s probe %d: detection ch %d not on shank — skipping shank', ...
            animal_code, job.probe_num, det_ch);
        continue
    end

    disp(sprintf('  shank probe=%d  n_chan=%d  detect_ch=%d', ...
        job.probe_num, numel(shank_rows), det_ch));

    % Filter detection channel (full recording), then keep extrema in play sessions
    filtered = filtfilt(Hd_delta.Coefficients, 1, LFP(det_ch, :));
    amp = abs(hilbert(filtered));
    std_amp = std(amp);
    [~, peak_idx] = findpeaks(filtered, ...
        'MinPeakProminence', 0.5 * std_amp, ...
        'MinPeakDistance', min_peak_dist);
    [~, trough_idx] = findpeaks(-filtered, ...
        'MinPeakProminence', 0.5 * std_amp, ...
        'MinPeakDistance', min_peak_dist);

    peak_idx   = local_keep_in_sessions(peak_idx, lfp_time, play_session, pre_samp, post_samp, n_time_total);
    trough_idx = local_keep_in_sessions(trough_idx, lfp_time, play_session, pre_samp, post_samp, n_time_total);

    lfp_peak   = local_cut_snippets(LFP, shank_rows, peak_idx, pre_samp, post_samp, n_time_win);
    lfp_trough = local_cut_snippets(LFP, shank_rows, trough_idx, pre_samp, post_samp, n_time_win);

    entry = struct();
    entry.animal_code        = animal_code;
    entry.npx_type           = NPX_Type;
    entry.probe_num          = job.probe_num;
    entry.shank_label        = job.shank_label;
    entry.detection_channel  = det_ch;
    entry.channel_ids        = shank_rows(:);
    entry.x_coords           = x_sh;
    entry.y_coords           = y_sh;
    entry.lfp_peak           = lfp_peak;    % (n_events, n_chan, n_time)
    entry.lfp_trough         = lfp_trough;
    entry.peak_times_s       = lfp_time(peak_idx)';
    entry.trough_times_s     = lfp_time(trough_idx)';
    entry.T0                 = T0;
    entry.fs                 = sr_LFP;
    entry.t1                 = t1;
    entry.t2                 = t2;
    entry.play_session       = play_session;
    entry.n_events_peak      = size(lfp_peak, 1);
    entry.n_events_trough    = size(lfp_trough, 1);

    if isempty(results)
        results = entry;
    else
        results(end + 1) = entry; %#ok<AGROW>
    end
end
end

%% -------------------------------------------------------------------------
function jobs = local_build_shank_jobs(NPX_Type, lpag_rows, xcoords, ycoords, chanMap, hard_x)
jobs = struct('probe_num', {}, 'shank_label', {}, 'detection_channel', {}, ...
    'lfp_rows', {}, 'x_coords', {}, 'y_coords', {});

if NPX_Type == 1
    Y_Range = lpag_rows{1, {'ProbeNum', 'depth_start', 'depth_end'}};
    this_indexes = ycoords >= Y_Range(2) & ycoords <= Y_Range(3);
    if ~any(this_indexes)
        return
    end
    all_locs = [xcoords(this_indexes), ycoords(this_indexes)];
    mean_loc = mean(all_locs, 1);
    [~, closest] = min(sum(abs([xcoords, ycoords] - repmat(mean_loc, numel(ycoords), 1)), 2));
    det_ch = chanMap(closest);

    % Full probe (one shank)
    valid = ~isnan(xcoords) & ~isnan(ycoords) & ~isnan(chanMap);
    jobs(1).probe_num = 1;
    jobs(1).shank_label = 'NPX1';
    jobs(1).detection_channel = det_ch;
    jobs(1).lfp_rows = chanMap(valid);
    jobs(1).x_coords = xcoords(valid);
    jobs(1).y_coords = ycoords(valid);
    return
end

% NPX2: one job per ProbeNum that has LPAG
probe_list = unique(lpag_rows.ProbeNum);
for p = 1:numel(probe_list)
    probe_n = probe_list(p);
    rows = lpag_rows(lpag_rows.ProbeNum == probe_n, :);
    depth_start = min(rows.depth_start);
    depth_end   = max(rows.depth_end);
    shank_mask = ismember(xcoords, hard_x(probe_n, :));
    lpag_mask  = shank_mask & ycoords >= depth_start & ycoords <= depth_end;
    if ~any(lpag_mask)
        continue
    end
    all_locs = [xcoords(lpag_mask), ycoords(lpag_mask)];
    mean_loc = mean(all_locs, 1);
    [~, closest] = min(sum(abs([xcoords, ycoords] - repmat(mean_loc, numel(ycoords), 1)), 2));
    % Prefer closest site that is on this shank
    if ~shank_mask(closest)
        shank_idx = find(shank_mask);
        d = sum(abs([xcoords(shank_idx), ycoords(shank_idx)] - repmat(mean_loc, numel(shank_idx), 1)), 2);
        [~, ii] = min(d);
        closest = shank_idx(ii);
    end
    det_ch = chanMap(closest);

    jobs(end + 1).probe_num = probe_n; %#ok<AGROW>
    jobs(end).shank_label = sprintf('probe%d', probe_n);
    jobs(end).detection_channel = det_ch;
    jobs(end).lfp_rows = chanMap(shank_mask);
    jobs(end).x_coords = xcoords(shank_mask);
    jobs(end).y_coords = ycoords(shank_mask);
end
end

function idx = local_keep_in_sessions(idx, lfp_time, play_session, pre_samp, post_samp, n_time_total)
if isempty(idx)
    return
end
t = lfp_time(idx);
in_sess = false(size(t));
for s = 1:size(play_session, 1)
    in_sess = in_sess | (t >= play_session.Start(s) & t <= play_session.End(s));
end
ok_edge = idx > pre_samp & idx + post_samp - 1 <= n_time_total;
idx = idx(in_sess(:) & ok_edge(:));
end

function stack = local_cut_snippets(LFP, shank_rows, event_idx, pre_samp, post_samp, n_time_win)
n_ev = numel(event_idx);
n_ch = numel(shank_rows);
stack = zeros(n_ev, n_ch, n_time_win, 'single');
for e = 1:n_ev
    a = event_idx(e) - pre_samp;
    b = event_idx(e) + post_samp - 1;
    stack(e, :, :) = LFP(shank_rows, a:b);
end
end
