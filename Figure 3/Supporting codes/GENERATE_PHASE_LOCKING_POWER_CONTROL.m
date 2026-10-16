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
% play_bout, readNPY, circular_uniformize, and CircStat.

if nargin < 3
    options = struct;
end
if ~isfield(options, 'n_rand'), options.n_rand = 1000; end
if ~isfield(options, 'min_bout_duration')
    options.min_bout_duration = 0.25;
end
if ~isfield(options, 'phase_method')
    options.phase_method = 'uniformized';
end

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
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
if isfile(mat_file)
    loaded = load(mat_file, 'LFP');
    lfp = double(loaded.LFP);
elseif isfile(dat_file)
    file_id = fopen(dat_file, 'r');
    cleanup = onCleanup(@() fclose(file_id));
    lfp = fread(file_id, 'int16');
    lfp = reshape(lfp, 384, []);
    clear cleanup
else
    error('No LFP_PAG.mat or LFP_PAG.dat found in %s.', npx_data_dir);
end

lfp_time = (1:size(lfp,2)) / sr_lfp;
[bout_table, rejected_bouts] = make_strict_pairs(play_bouts, ...
    lfp_time(1), lfp_time(end), options.min_bout_duration);
if isempty(bout_table)
    error('No valid one-to-one pre-play/play bout pairs for %s.', animal_code);
end

channels_with_spikes = unique(cluster_info.ch)';
channel_data = repmat(struct( ...
    'channel', [], 'pre_power', [], 'play_power', [], ...
    'pre_power_db', [], 'play_power_db', [], ...
    'power_difference_db', [], 'amplitude_reference', [], ...
    'phase_method', ''), numel(channels_with_spikes), 1);
neurons = repmat(struct( ...
    'cluster_id', [], 'channel', [], 'area', '', ...
    'entire_preferred_angle', NaN, 'entire_ppc', NaN, ...
    'entire_ppc_p', NaN, 'pre_phase', {{}}, 'play_phase', {{}}, ...
    'pre_amplitude', {{}}, 'play_amplitude', {{}}, ...
    'pre_spike_count', [], 'play_spike_count', []), height(cluster_info), 1);

for cc = 1:numel(channels_with_spikes)
    channel = channels_with_spikes(cc);
    fprintf('Power-control analysis: channel %d (%d/%d)\n', ...
        channel, cc, numel(channels_with_spikes));

    filtered_lfp = filtfilt(Hd_delta.Coefficients, 1, double(lfp(channel,:)));
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
    channel_data(cc).amplitude_reference = single(prctile(amplitude, 0:100));
    channel_data(cc).phase_method = options.phase_method;

    neuron_index = find(cluster_info.ch == channel)';
    for nn = neuron_index
        cluster_id = cluster_info.cluster_id(nn);
        unit_spikes = spike_times(spike_clusters == cluster_id);
        finite_spikes = unit_spikes >= lfp_time(1) & unit_spikes <= lfp_time(end);
        unit_spikes = unit_spikes(finite_spikes);
        all_phase = interpolate_circular_phase(lfp_time, phase, unit_spikes);

        neurons(nn).cluster_id = cluster_id;
        neurons(nn).channel = channel;
        neurons(nn).area = cluster_info.area{nn};
        neurons(nn).entire_preferred_angle = circ_mean(all_phase(:));
        neurons(nn).entire_ppc = compute_ppc_fast(all_phase);
        if options.n_rand > 0 && numel(all_phase) >= 2
            random_ppc = nan(options.n_rand,1);
            for rr = 1:options.n_rand
                sampled_phase = randsample(phase, numel(all_phase), false);
                random_ppc(rr) = compute_ppc_fast(sampled_phase);
            end
            neurons(nn).entire_ppc_p = ...
                (1 + sum(random_ppc >= neurons(nn).entire_ppc)) / ...
                (options.n_rand + 1);
        end

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
    end
end

control_struct.animal = animal_batch;
control_struct.session = animal_code;
control_struct.cluster_info = cluster_info;
control_struct.bout_table = bout_table;
control_struct.rejected_bouts = rejected_bouts;
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
