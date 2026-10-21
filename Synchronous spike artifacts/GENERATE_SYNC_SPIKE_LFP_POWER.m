function session_result = GENERATE_SYNC_SPIKE_LFP_POWER(npx_data_dir, options)
%GENERATE_SYNC_SPIKE_LFP_POWER  Raw-LFP power at synchronous spikes, one session.
%
% Called by Estimate_sync_spike_lfp_power. Do not start here.
%
% The pairs are the same play-bout coincidences as
% GENERATE_SYNC_SPIKE_WAVEFORMS. For each neuron the value at a spike is
% one sample of raw-LFP power on that neuron's channel (cluster ch + 1):
%
%   raw = LFP_PAG on that channel, unfiltered
%   power = abs(hilbert(double(raw))).^2
%   power_at_spike = interp1(lfp_time, power, spike_time, 'linear', NaN)
%   lfp_time = (1:numel(power))' / 2500
%
% This is not a cut of LFP samples around the spike, and it is not the
% 1-5 Hz delta power used for phase locking. The Hilbert transform uses
% the whole trace; the number kept for each spike is the power at that
% instant. The comparison pool is the spikes inside the play-bout windows,
% the same windows used for the waveform pipeline. The observed statistic
% is the mean of power at the synchronous spikes. Each of 500 nulls is the
% mean of the same count of this pair's non-synchronous spikes in those
% windows. A high percentile means the synchronous spikes sit on larger
% raw-LFP power. Each LFP channel is read once. mean_power_nonsync is the
% mean of every non-synchronous spike in the windows. mean_power_all is
% the mean of every spike in the windows, including the synchronous ones.

if nargin < 2 || isempty(options)
    options = struct();
end
options = sync_analysis_options(options);
[~, session_id] = fileparts(npx_data_dir);
session_result = blank_power_session(session_id, options);
if options.max_pairs < 1
    session_result.status = 'pair_budget_exhausted';
    return
end

selection = synchronous_play_pairs(npx_data_dir, options);
session_result.session_id = selection.session_id;
if isempty(selection.pairs)
    session_result.status = 'no_eligible_pair';
    fprintf('%s  no eligible synchronous pair\n', session_result.session_id);
    return
end
print_pairs(session_result.session_id, selection.pairs);
seed = session_rng_seed(options, session_result.session_id);
if ~isempty(seed)
    rng(seed, 'twister');
end
session_result.n_pairs = numel(selection.pairs);

cache = power_at_spikes(npx_data_dir, selection, options);
rows = empty_power_rows();
for pair_index = 1:numel(selection.pairs)
    pair = selection.pairs(pair_index);
    sides = [1 2];
    for side_index = 1:numel(sides)
        if sides(side_index) == 1
            cluster_id = pair.cluster_id_1;
            partner_id = pair.cluster_id_2;
            sync_index = pair.sync_index_1;
            area = pair.area_1;
            channel = pair.channel_1;
            n_spikes = pair.n_spikes_1;
        else
            cluster_id = pair.cluster_id_2;
            partner_id = pair.cluster_id_1;
            sync_index = pair.sync_index_2;
            area = pair.area_2;
            channel = pair.channel_2;
            n_spikes = pair.n_spikes_2;
        end
        neuron = neuron_by_id(selection.neurons, cluster_id);
        power = power_by_id(cache, cluster_id);
        if isempty(power) || isempty(neuron)
            note = sprintf('%s cluster %g pair %d: no LFP power', ...
                session_result.session_id, cluster_id, pair_index);
            session_result.notes{end + 1} = note; %#ok<AGROW>
            fprintf('  %s\n', note);
            continue
        end
        local_sync = local_sync_index(sync_index, neuron.interval_index);
        comparison = compare_sync_scalar(power, local_sync, options.n_repeats);
        if ~comparison.ok
            note = sprintf('%s cluster %g pair %d: %s', ...
                session_result.session_id, cluster_id, pair_index, comparison.reason);
            session_result.notes{end + 1} = note; %#ok<AGROW>
            fprintf('  %s\n', note);
            continue
        end
        rows = append_power_row(rows, power_row( ...
            session_result.session_id, pair_index, cluster_id, partner_id, ...
            area, channel, n_spikes, comparison, options.power_definition));
        fprintf(['%s  cluster %g  interval %d  sync %d  sync power %.3g  ' ...
            'nonsync power %.3g  percentile %.1f\n'], ...
            session_result.session_id, cluster_id, comparison.n_pool, comparison.n_sync, ...
            comparison.mean_sync, comparison.mean_nonsync, 100 * comparison.percentile);
    end
end
session_result.rows = rows;
if isempty(rows)
    session_result.status = 'no_saved_neuron';
else
    session_result.status = 'ok';
end
end

function cache = power_at_spikes(session_dir, selection, options)
cache = struct('cluster_id', {}, 'power', {});
channels = unique([selection.neurons.channel]);
for channel_index = 1:numel(channels)
    channel = channels(channel_index);
    started = tic;
    raw = read_lfp_channel(session_dir, channel, options.n_lfp_channels);
    envelope = abs(hilbert(double(raw(:))));
    clear raw
    power = envelope .^ 2;
    clear envelope
    lfp_time = (1:numel(power))' / options.lfp_sample_rate;
    fprintf('%s  LFP channel %g  %d samples  %.1f s\n', ...
        selection.session_id, channel, numel(power), toc(started));
    neuron_ids = neurons_on_channel(selection.neurons, channel);
    for neuron_index = 1:numel(neuron_ids)
        neuron = neuron_by_id(selection.neurons, neuron_ids(neuron_index));
        cache(end + 1).cluster_id = neuron.cluster_id; %#ok<AGROW>
        cache(end).power = interp1(lfp_time, power, ...
            neuron.spike_times(neuron.interval_index), 'linear', NaN);
    end
    clear power lfp_time
end
end

function ids = neurons_on_channel(neurons, channel)
ids = [];
for index = 1:numel(neurons)
    if neurons(index).channel == channel
        ids(end + 1) = neurons(index).cluster_id; %#ok<AGROW>
    end
end
ids = unique(ids);
end

function raw = read_lfp_channel(session_dir, channel, n_lfp_channels)
mat_file = fullfile(session_dir, 'LFP_PAG.mat');
dat_file = fullfile(session_dir, 'LFP_PAG.dat');
if exist(mat_file, 'file') == 2
    loaded = load(mat_file, 'LFP');
    if channel < 1 || channel > size(loaded.LFP, 1) || channel ~= round(channel)
        error('GENERATE_SYNC_SPIKE_LFP_POWER:lfp', ...
            'LFP channel %g is outside LFP_PAG.mat.', channel);
    end
    raw = loaded.LFP(channel, :);
    return
end
if exist(dat_file, 'file') ~= 2
    error('GENERATE_SYNC_SPIKE_LFP_POWER:lfp', ...
        'Missing LFP_PAG.mat and LFP_PAG.dat in %s.', session_dir);
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (n_lfp_channels * 2);
if n_samples ~= round(n_samples) || n_samples < 1
    error('GENERATE_SYNC_SPIKE_LFP_POWER:lfp', ...
        'LFP_PAG.dat is not a multiple of %d int16 channels.', n_lfp_channels);
end
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('GENERATE_SYNC_SPIKE_LFP_POWER:lfp', 'Could not open %s.', dat_file);
end
closer = onCleanup(@() fclose(file_id)); %#ok<NASGU>
if fseek(file_id, 2 * (channel - 1), 'bof') ~= 0
    error('GENERATE_SYNC_SPIKE_LFP_POWER:lfp', ...
        'Could not seek to LFP channel %d.', channel);
end
[raw, count] = fread(file_id, n_samples, 'int16=>double', 2 * (n_lfp_channels - 1));
if count ~= n_samples
    error('GENERATE_SYNC_SPIKE_LFP_POWER:lfp', ...
        'Read %d of %d LFP samples on channel %d.', count, n_samples, channel);
end
end

function local = local_sync_index(sync_index, interval_index)
[is_in, local] = ismember(sync_index(:), interval_index(:));
local = local(is_in);
end

function power = power_by_id(cache, cluster_id)
power = [];
for index = 1:numel(cache)
    if cache(index).cluster_id == cluster_id
        power = cache(index).power;
        return
    end
end
end

function neuron = neuron_by_id(neurons, cluster_id)
neuron = [];
for index = 1:numel(neurons)
    if neurons(index).cluster_id == cluster_id
        neuron = neurons(index);
        return
    end
end
end

function session_result = blank_power_session(session_id, options)
session_result.session_id = session_id;
session_result.status = '';
session_result.n_pairs = 0;
session_result.power_definition = options.power_definition;
session_result.notes = {};
session_result.options = options;
session_result.rows = empty_power_rows();
end

function rows = empty_power_rows()
rows = struct('session_id', {}, 'pair_index', {}, 'cluster_id', {}, ...
    'partner_id', {}, 'area', {}, 'channel', {}, ...
    'n_spikes', {}, 'n_interval', {}, 'n_sync', {}, 'n_nonsync', {}, ...
    'mean_power_all', {}, 'mean_power_nonsync', {}, 'mean_power_sync', {}, ...
    'null_mean_power', {}, 'percentile', {}, 'power_definition', {});
end

function row = power_row(session_id, pair_index, cluster_id, partner_id, ...
    area, channel, n_spikes, comparison, power_definition)
row.session_id = session_id;
row.pair_index = pair_index;
row.cluster_id = cluster_id;
row.partner_id = partner_id;
row.area = area;
row.channel = channel;
row.n_spikes = n_spikes;
row.n_interval = comparison.n_pool;
row.n_sync = comparison.n_sync;
row.n_nonsync = comparison.n_nonsync;
row.mean_power_all = comparison.mean_all;
row.mean_power_nonsync = comparison.mean_nonsync;
row.mean_power_sync = comparison.mean_sync;
row.null_mean_power = comparison.null_means;
row.percentile = comparison.percentile;
row.power_definition = power_definition;
end

function rows = append_power_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end + 1) = orderfields(row, rows); %#ok<AGROW>
end
end

function print_pairs(session_id, pairs)
for pair_index = 1:numel(pairs)
    pair = pairs(pair_index);
    fprintf(['%s  %s ch %g (%d spikes, %d in window, %d sync)  with  ' ...
        '%s ch %g (%d spikes, %d in window, %d sync)\n'], ...
        session_id, pair.area_1, pair.channel_1, pair.n_spikes_1, ...
        pair.n_interval_1, pair.n_sync_1, ...
        pair.area_2, pair.channel_2, pair.n_spikes_2, ...
        pair.n_interval_2, pair.n_sync_2);
end
end

function seed = session_rng_seed(options, session_id)
seed = [];
if ~isfield(options, 'random_seed') || isempty(options.random_seed) ...
        || ~isscalar(options.random_seed) || ~isfinite(options.random_seed)
    return
end
bytes = double(uint8(session_id(:)'));
seed = mod(abs(double(options.random_seed)) + sum(bytes .* (1:numel(bytes))), 2^31 - 1);
seed = floor(seed);
if seed < 1
    seed = 1;
end
end
