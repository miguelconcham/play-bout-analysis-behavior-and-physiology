function control_struct = GENERATE_PHASE_LOCKING_POWER_CONTROL(npx_data_dir, Hd_delta, options)
%GENERATE_PHASE_LOCKING_POWER_CONTROL Bout-resolved PPC and delta-power data.
%
% Each play bout is paired with an immediately preceding window of exactly
% the same duration. Pairs whose pre-window overlaps another play bout (or
% falls outside the LFP recording) are rejected. Unlike the original Figure
% 3 pipeline, spike phases are retained per bout so PPC can be recomputed
% after matching or resampling bouts by delta power.
%
% Delta power is mean(abs(hilbert(delta_lfp)).^2). Both linear power and dB
% power are stored. ADC calibration is unnecessary for within-channel dB
% differences, but linear values remain in ADC-count^2 unless the source
% LFP has already been converted.
%
% Required external functions are the same as the original pipeline:
% play_bout, readNPY, and circular_uniformize.

if nargin < 3
    options = struct;
end
if ~isfield(options, 'min_bout_duration')
    options.min_bout_duration = 0.25;
end
if ~isfield(options, 'phase_method')
    options.phase_method = 'uniformized';
end
if ~isfield(options, 'epoch_duration')
    options.epoch_duration = 1;
end

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
synch_directory = fullfile(data_root, 'Synch data');
behavior_data = fullfile(data_root, 'Behavior backups');
area_limit_file = fullfile(data_root, 'Area_limits_GoodLooking.xlsx');
play_behaviors = {'Pounce','CC','Boxing','Evasion','Pin','Escape','CB','CD'};
area_name = 'PAG';
sr_lfp = 2500;

[~, animal_code] = fileparts(npx_data_dir);
animal_parts = strsplit(animal_code, ' ');
animal_batch = animal_parts{1};
repeated_animal = animal_parts{3};

load(fullfile(synch_directory, animal_code, 'synch_model_video2NPX.mat'), ...
    'synch_model_video2NPX');
behavior = readtable(fullfile(behavior_data, [animal_code, '.txt']));
behavior(:,2) = [];
behavior.Properties.VariableNames = {'Animal','Start','End','Length','Type'};
behavior.Type2 = behavior.Type;
behavior.Type2(ismember(behavior.Type2, {'Pounce_A','Pounce_B'})) = {'Pounce'};
behavior.Type2(ismember(behavior.Type2, {'Pounce_Ai','Pounce_Bi'})) = {'PounceI'};
behavior.Type2(ismember(behavior.Type2, {''})) = {'Other'};
behavior(ismember(behavior.Animal, {'Reversal'}),:) = [];
behavior.Start = predict(synch_model_video2NPX, behavior.Start);
behavior.End = predict(synch_model_video2NPX, behavior.End);

animal_types = unique(behavior.Animal);
animal_types(ismember(animal_types, {'Session_structure'})) = [];
config.Behavior = behavior;
config.repeated_animal = repeated_animal;
config.animal_types = animal_types;
config.play_behaviors = play_behaviors;
config.beh_bin = 0.01;
config.conv_length = 1;
config.behavior_window = 0;
play_bouts = play_bout(config);
if istable(play_bouts)
    play_bouts = play_bouts{:,1:2};
else
    play_bouts = play_bouts(:,1:2);
end
play_bouts = sortrows(double(play_bouts), 1);

load(fullfile(npx_data_dir, ['chann_map_', area_name, '.mat']), ...
    'chanMap', 'xcoords', 'ycoords');
[areas_by_channel, channel_map, npx_type] = map_channels( ...
    chanMap, xcoords, ycoords, area_limit_file, animal_batch, repeated_animal);

spike_times = double(readNPY(fullfile(npx_data_dir, ...
    ['spike_times_', area_name, '.npy']))) / 30000;
spike_clusters = readNPY(fullfile(npx_data_dir, ...
    ['spike_clusters_', area_name, '.npy']));
cluster_info = readtable(fullfile(npx_data_dir, ...
    ['cluster_info_', area_name, '.tsv']), ...
    "FileType", "text", "Delimiter", '\t');
cluster_info = cluster_info(ismember(cluster_info.group, {'mua','good'}),:);
cluster_info.area = areas_by_channel(cluster_info.ch + 1);
cluster_info.ch = cluster_info.ch + 1;

mat_file = fullfile(npx_data_dir, ['LFP_', area_name, '.mat']);
dat_file = fullfile(npx_data_dir, ['LFP_', area_name, '.dat']);
if npx_type == 2
    if ~isfile(mat_file)
        error(['Session was identified as NPX2, but %s is missing. ' ...
            'The NPX1 binary fallback is intentionally disabled.'], mat_file);
    end
    loaded = load(mat_file, 'LFP');
    lfp = loaded.LFP;
    clear loaded
    if isempty(lfp) || size(lfp,1) ~= 384
        error('NPX2 LFP must be a nonempty 384-by-N array in %s.', mat_file);
    end
    n_lfp_samples = size(lfp,2);
else
    if ~isfile(dat_file)
        error('Session was identified as NPX1, but %s is missing.', dat_file);
    end
    file_info = dir(dat_file);
    bytes_per_sample = 2; % int16
    n_lfp_samples = file_info.bytes / (384*bytes_per_sample);
    if n_lfp_samples < 1 || n_lfp_samples ~= floor(n_lfp_samples)
        error('NPX1 LFP binary is empty or not divisible by 384 int16 channels.');
    end
    lfp = []; % NPX1 channels are streamed from disk inside the channel loop.
end

lfp_time = (1:n_lfp_samples) / sr_lfp;
[bout_table, rejected_bouts] = make_strict_pairs(play_bouts, ...
    lfp_time(1), lfp_time(end), options.min_bout_duration);

% Preserve the original Figure 3 duration-balancing logic for the candidate
% pools, then divide both conditions into complete fixed-duration epochs.
% Equal numbers of these epochs are selected downstream, so analyzed time
% is exactly equal even if clipping leaves different candidate-pool sizes.
pre_play_periods = build_balanced_periods(play_bouts);
pre_epoch_table = intervals_to_epochs(pre_play_periods, ...
    options.epoch_duration, lfp_time(1), lfp_time(end));
play_epoch_table = intervals_to_epochs(play_bouts, ...
    options.epoch_duration, lfp_time(1), lfp_time(end));
if isempty(pre_epoch_table) || isempty(play_epoch_table)
    error('No complete %.3g-second pre-play/play epochs for %s.', ...
        options.epoch_duration, animal_code);
end

channels_with_spikes = unique(cluster_info.ch)';
channel_data = repmat(struct( ...
    'channel', [], 'pre_power', [], 'play_power', [], ...
    'pre_power_db', [], 'play_power_db', [], ...
    'power_difference_db', [], 'amplitude_reference', [], ...
    'pre_epoch_power_db', [], 'play_epoch_power_db', [], ...
    'phase_method', ''), numel(channels_with_spikes), 1);
neurons = repmat(struct( ...
    'cluster_id', [], 'channel', [], 'area', '', ...
    'entire_preferred_angle', NaN, 'entire_ppc', NaN, ...
    'entire_ppc_p', NaN, 'pre_phase', {{}}, 'play_phase', {{}}, ...
    'pre_amplitude', {{}}, 'play_amplitude', {{}}, ...
    'pre_spike_count', [], 'play_spike_count', [], ...
    'pre_epoch_phase', {{}}, 'play_epoch_phase', {{}}, ...
    'pre_epoch_spike_count', [], 'play_epoch_spike_count', []), ...
    height(cluster_info), 1);

for cc = 1:numel(channels_with_spikes)
    channel = channels_with_spikes(cc);
    fprintf('Power-control analysis: channel %d (%d/%d)\n', ...
        channel, cc, numel(channels_with_spikes));

    if npx_type == 2
        raw_lfp_channel = double(lfp(channel,:));
    else
        raw_lfp_channel = double(read_npx1_lfp_channel( ...
            dat_file, channel, n_lfp_samples));
    end
    filtered_lfp = filtfilt(Hd_delta.Coefficients, 1, raw_lfp_channel);
    analytic_lfp = hilbert(filtered_lfp);
    amplitude = abs(analytic_lfp);
    phase = angle(analytic_lfp);

    [~, peak_index] = findpeaks(filtered_lfp, ...
        'MinPeakProminence', 0.5 * std(amplitude), ...
        'MinPeakDistance', sr_lfp / Hd_delta.CutoffFrequency2);
    if strcmpi(options.phase_method, 'uniformized')
        phase = circular_uniformize(phase);
    elseif ~strcmpi(options.phase_method, 'hilbert')
        error('phase_method must be ''uniformized'' or ''hilbert''.');
    end
    if ~isempty(peak_index)
        phase_origin = angle(mean(exp(1i * phase(peak_index))));
        phase = mod(phase - phase_origin + 5*pi, 2*pi) - pi;
    end

    n_bouts = height(bout_table);
    pre_power = nan(n_bouts,1);
    play_power = nan(n_bouts,1);
    for bb = 1:n_bouts
        pre_mask = lfp_time >= bout_table.PreStart(bb) & ...
            lfp_time < bout_table.PreEnd(bb);
        play_mask = lfp_time >= bout_table.PlayStart(bb) & ...
            lfp_time < bout_table.PlayEnd(bb);
        pre_power(bb) = mean(amplitude(pre_mask).^2, 'omitnan');
        play_power(bb) = mean(amplitude(play_mask).^2, 'omitnan');
    end
    channel_data(cc).channel = channel;
    channel_data(cc).pre_power = pre_power;
    channel_data(cc).play_power = play_power;
    channel_data(cc).pre_power_db = 10*log10(pre_power);
    channel_data(cc).play_power_db = 10*log10(play_power);
    channel_data(cc).power_difference_db = ...
        channel_data(cc).play_power_db - channel_data(cc).pre_power_db;
    channel_data(cc).pre_epoch_power_db = epoch_power_db( ...
        amplitude, lfp_time, pre_epoch_table);
    channel_data(cc).play_epoch_power_db = epoch_power_db( ...
        amplitude, lfp_time, play_epoch_table);
    channel_data(cc).amplitude_reference = single(prctile(amplitude, 0:100));
    channel_data(cc).phase_method = options.phase_method;

    neuron_index = find(cluster_info.ch == channel)';
    for nn = neuron_index
        cluster_id = cluster_info.cluster_id(nn);
        unit_spikes = spike_times(spike_clusters == cluster_id);
        finite_spikes = unit_spikes >= lfp_time(1) & unit_spikes <= lfp_time(end);
        unit_spikes = unit_spikes(finite_spikes);
        neurons(nn).cluster_id = cluster_id;
        neurons(nn).channel = channel;
        neurons(nn).area = cluster_info.area{nn};

        neurons(nn).pre_phase = cell(n_bouts,1);
        neurons(nn).play_phase = cell(n_bouts,1);
        neurons(nn).pre_amplitude = cell(n_bouts,1);
        neurons(nn).play_amplitude = cell(n_bouts,1);
        neurons(nn).pre_spike_count = zeros(n_bouts,1);
        neurons(nn).play_spike_count = zeros(n_bouts,1);
        for bb = 1:n_bouts
            pre_spikes = unit_spikes(unit_spikes >= bout_table.PreStart(bb) & ...
                unit_spikes < bout_table.PreEnd(bb));
            play_spikes = unit_spikes(unit_spikes >= bout_table.PlayStart(bb) & ...
                unit_spikes < bout_table.PlayEnd(bb));
            neurons(nn).pre_phase{bb} = interpolate_circular_phase( ...
                lfp_time, phase, pre_spikes);
            neurons(nn).play_phase{bb} = interpolate_circular_phase( ...
                lfp_time, phase, play_spikes);
            neurons(nn).pre_amplitude{bb} = interp1(lfp_time, amplitude, ...
                pre_spikes, 'linear');
            neurons(nn).play_amplitude{bb} = interp1(lfp_time, amplitude, ...
                play_spikes, 'linear');
            neurons(nn).pre_spike_count(bb) = numel(pre_spikes);
            neurons(nn).play_spike_count(bb) = numel(play_spikes);
        end

        [neurons(nn).pre_epoch_phase, ...
            neurons(nn).pre_epoch_spike_count] = phases_by_epoch( ...
            unit_spikes, lfp_time, phase, pre_epoch_table);
        [neurons(nn).play_epoch_phase, ...
            neurons(nn).play_epoch_spike_count] = phases_by_epoch( ...
            unit_spikes, lfp_time, phase, play_epoch_table);
    end
end

control_struct.animal = animal_batch;
control_struct.session = animal_code;
control_struct.cluster_info = cluster_info;
control_struct.bout_table = bout_table;
control_struct.rejected_bouts = rejected_bouts;
control_struct.pre_epoch_table = pre_epoch_table;
control_struct.play_epoch_table = play_epoch_table;
control_struct.channel_data = channel_data;
control_struct.neurons = neurons;
control_struct.channel_map = channel_map;
control_struct.areas_by_channel = areas_by_channel;
control_struct.npx_type = npx_type;
control_struct.sample_rate_lfp = sr_lfp;
control_struct.delta_band = [Hd_delta.CutoffFrequency1, Hd_delta.CutoffFrequency2];
control_struct.options = options;
end

function [pairs, rejected] = make_strict_pairs(play_bouts, first_time, last_time, min_duration)
duration = play_bouts(:,2) - play_bouts(:,1);
pre_start = play_bouts(:,1) - duration;
valid = isfinite(duration) & duration >= min_duration & ...
    pre_start >= first_time & play_bouts(:,2) <= last_time;
overlap = false(size(valid));
for bb = 1:size(play_bouts,1)
    if ~valid(bb), continue; end
    other = true(size(valid));
    other(bb) = false;
    overlap(bb) = any(play_bouts(other,1) < play_bouts(bb,1) & ...
        play_bouts(other,2) > pre_start(bb));
end
valid = valid & ~overlap;
bout_id = (1:size(play_bouts,1))';
pairs = table(bout_id(valid), pre_start(valid), play_bouts(valid,1), ...
    play_bouts(valid,1), play_bouts(valid,2), duration(valid), ...
    'VariableNames', {'BoutID','PreStart','PreEnd','PlayStart','PlayEnd','Duration'});
rejected = table(bout_id(~valid), play_bouts(~valid,1), play_bouts(~valid,2), ...
    overlap(~valid), 'VariableNames', {'BoutID','PlayStart','PlayEnd','PreOverlap'});
end

function [areas_by_channel, channel_map, npx_type] = map_channels( ...
    chan_map, xcoords, ycoords, area_limit_file, animal_batch, repeated_animal)
hard_coded_x = [8 40; 258 290; 508 540; 758 790];
chan_map = chan_map(:);
xcoords = xcoords(:);
ycoords = ycoords(:);
npx_type = 1 + any(ismember(xcoords, hard_coded_x));
if npx_type == 1 && ~ismember(192, chan_map)
    pos_191 = find(chan_map == 191, 1);
    pos_193 = find(chan_map == 193, 1);
    if ~isempty(pos_191) && pos_193 == pos_191 + 1
        xcoords = [xcoords; NaN];
        xcoords(pos_193+1:end) = xcoords(pos_193:end-1);
        xcoords(pos_193) = 43;
        ycoords = [ycoords; NaN];
        ycoords(pos_193+1:end) = ycoords(pos_193:end-1);
        ycoords(pos_193) = 1900;
        chan_map = [chan_map; NaN];
        chan_map(pos_193+1:end) = chan_map(pos_193:end-1);
        chan_map(pos_193) = 192;
    else
        error('Inconsistent NPX1 channel map: channel 192 is missing.');
    end
end
areas_by_channel = repmat({''}, 384, 1);
channel_map = nan(384,2);
area_limit = readtable(area_limit_file);
if strcmp(repeated_animal, 'Single2')
    animal_key = ['Batch', animal_batch(2), repeated_animal];
else
    animal_key = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, animal_key),:);

for ii = 1:numel(chan_map)
    channel = chan_map(ii);
    if channel < 1 || channel > 384, continue; end
    channel_map(channel,:) = [xcoords(ii), ycoords(ii)];
    if npx_type == 1
        match = ycoords(ii) >= area_limit.depth_start & ...
            ycoords(ii) < area_limit.depth_end + 1 & ...
            ismember(area_limit.Probe_Area, {'PAG'});
    else
        probe = find(any(ismember(hard_coded_x, xcoords(ii)),2), 1);
        match = ycoords(ii) >= area_limit.depth_start & ...
            ycoords(ii) < area_limit.depth_end + 1 & ...
            area_limit.ProbeNum == probe & ismember(area_limit.Probe_Area, {'PAG'});
    end
    if any(match)
        areas_by_channel{channel} = area_limit.area{find(match,1)};
    end
end
end

function phase_at_events = interpolate_circular_phase(time, phase, event_time)
if isempty(event_time)
    phase_at_events = zeros(0,1);
    return
end
phasor = interp1(time, exp(1i*phase), event_time, 'linear');
phase_at_events = angle(phasor);
end

function epoch_table = intervals_to_epochs(intervals, epoch_duration, first_time, last_time)
intervals = sortrows(double(intervals(:,1:2)), 1);
rows = zeros(0,3);
for ii = 1:size(intervals,1)
    interval_start = max(intervals(ii,1), first_time);
    interval_end = min(intervals(ii,2), last_time);
    n_epochs = floor((interval_end-interval_start) / epoch_duration);
    if n_epochs < 1, continue; end
    starts = interval_start + (0:n_epochs-1)' * epoch_duration;
    rows = [rows; [repmat(ii,n_epochs,1), starts, starts+epoch_duration]]; %#ok<AGROW>
end
epoch_table = array2table(rows, ...
    'VariableNames', {'SourceInterval','Start','End'});
end

function power_db = epoch_power_db(amplitude, lfp_time, epoch_table)
power_db = nan(height(epoch_table),1);
for ee = 1:height(epoch_table)
    mask = lfp_time >= epoch_table.Start(ee) & lfp_time < epoch_table.End(ee);
    power_db(ee) = 10*log10(mean(amplitude(mask).^2, 'omitnan'));
end
end

function [phase_cells, spike_count] = phases_by_epoch( ...
    spikes, lfp_time, phase, epoch_table)
phase_cells = cell(height(epoch_table),1);
spike_count = zeros(height(epoch_table),1);
for ee = 1:height(epoch_table)
    epoch_spikes = spikes(spikes >= epoch_table.Start(ee) & ...
        spikes < epoch_table.End(ee));
    phase_cells{ee} = interpolate_circular_phase( ...
        lfp_time, phase, epoch_spikes);
    spike_count(ee) = numel(epoch_spikes);
end
end

function channel_data = read_npx1_lfp_channel( ...
    dat_file, channel, n_samples)
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('Could not open NPX1 LFP binary: %s', dat_file);
end
cleanup = onCleanup(@() fclose(file_id));
status = fseek(file_id, 2*(channel-1), 'bof');
if status ~= 0
    error('Could not seek to channel %d in %s.', channel, dat_file);
end
[channel_data, count] = fread(file_id, n_samples, ...
    '*int16', 2*(384-1));
if count ~= n_samples
    error('Read %d/%d samples for channel %d in %s.', ...
        count, n_samples, channel, dat_file);
end
end
