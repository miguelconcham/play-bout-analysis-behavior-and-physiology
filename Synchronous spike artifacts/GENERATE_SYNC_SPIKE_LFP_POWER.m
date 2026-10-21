function session_result = GENERATE_SYNC_SPIKE_LFP_POWER(npx_data_dir, options)
%GENERATE_SYNC_SPIKE_LFP_POWER  Raw-LFP power at synchronous spikes, one session.
%
% Called by Estimate_sync_spike_lfp_power. Do not start here.
%
% The pairs are the same play-bout coincidences as
% GENERATE_SYNC_SPIKE_WAVEFORMS. For each neuron the value at a synchronous
% spike is one sample of raw-LFP power on that neuron's channel
% (cluster ch + 1):
%
%   raw = LFP_PAG on that channel, unfiltered
%   power = abs(hilbert(double(raw))).^2
%   sample = round(spike_time * 2500)     nearest LFP sample, time k/2500
%
% This is not a cut of LFP samples around the spike, and it is not 1-5 Hz
% delta power. The 1-5 Hz band is used only to name the phase. Phase is
% angle(hilbert(filtfilt(b, 1, raw))) with the same FIR as the coincidence
% pipeline: bandpassfir, order 2000, Hamming, 1-5 Hz, 2500 Hz. The angle
% is not passed through circular_uniformize and is not rotated so that
% peaks share one phase. Angle 0 is the positive peak of that filtered
% oscillation, and +/- pi is the trough.
%
% The null is random LFP samples inside the same play-bout windows
% [bout_start-1.5, bout_start+1] s. Each synchronous event is matched to
% samples in the same 10-degree delta-phase bin. Samples within +/- 10 ms
% of that pair's synchronous spikes are left out, so the event is not
% compared with itself. A high percentile means the synchronous spikes sit
% on larger raw-LFP power than other moments at the same delta phase.
% Each LFP channel is read once.

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
session_result.coverage = selection.coverage;
print_coverage(session_result.session_id, selection.coverage, options.min_spikes);
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

coefficients = delta_fir_coefficients(options);
rows = empty_power_rows();
channels = unique([selection.neurons.channel]);
for channel_index = 1:numel(channels)
    channel = channels(channel_index);
    started = tic;
    raw = double(read_lfp_channel(npx_data_dir, channel, options.n_lfp_channels));
    raw = raw(:);
    power = abs(hilbert(raw)) .^ 2;
    phase = angle(hilbert(filtfilt(coefficients, 1, raw)));
    clear raw
    fprintf('%s  LFP channel %g  %d samples  %.1f s\n', ...
        session_result.session_id, channel, numel(power), toc(started));
    [rows, session_result] = record_channel_pairs(rows, session_result, selection, ...
        channel, power, phase, options);
    clear power phase
end
session_result.rows = rows;
if isempty(rows)
    session_result.status = 'no_saved_neuron';
else
    session_result.status = 'ok';
end
end

function [rows, session_result] = record_channel_pairs(rows, session_result, selection, ...
    channel, power, phase, options)
n_samples = numel(power);
for pair_index = 1:numel(selection.pairs)
    pair = selection.pairs(pair_index);
    for side_index = 1:2
        side = pair_side(pair, side_index);
        if side.channel ~= channel
            continue
        end
        neuron = neuron_by_id(selection.neurons, side.cluster_id);
        if isempty(neuron)
            note = sprintf('%s cluster %g pair %d: neuron missing from the selection', ...
                session_result.session_id, side.cluster_id, pair_index);
            session_result.notes{end + 1} = note; %#ok<AGROW>
            fprintf('  %s\n', note);
            continue
        end
        [pool_index, sync_sample] = lfp_phase_pool(n_samples, options.lfp_sample_rate, ...
            selection.bout_starts, options.hist_range, ...
            neuron.spike_times(side.sync_index), options.phase_guard_seconds);
        kept = sync_sample >= 1;
        comparison = compare_phase_matched_lfp( ...
            power(sync_sample(kept)), phase(sync_sample(kept)), ...
            power(pool_index), phase(pool_index), ...
            options.n_repeats, options.n_phase_bins);
        if ~comparison.ok
            note = sprintf('%s cluster %g pair %d: %s', ...
                session_result.session_id, side.cluster_id, pair_index, comparison.reason);
            session_result.notes{end + 1} = note; %#ok<AGROW>
            fprintf('  %s\n', note);
            continue
        end
        rows = append_power_row(rows, power_row( ...
            session_result.session_id, pair_index, side, neuron, ...
            comparison, options.power_definition));
        fprintf(['%s  cluster %g  interval spikes %d  sync %d  phase pool %d  ' ...
            'sync power %.3g  phase-matched %.3g  percentile %.1f\n'], ...
            session_result.session_id, side.cluster_id, numel(neuron.interval_index), ...
            comparison.n_sync, comparison.n_pool, comparison.mean_sync, ...
            comparison.mean_phase_matched, 100 * comparison.percentile);
    end
end
end

function side = pair_side(pair, side_index)
if side_index == 1
    side.cluster_id = pair.cluster_id_1;
    side.partner_id = pair.cluster_id_2;
    side.sync_index = pair.sync_index_1;
    side.area = pair.area_1;
    side.channel = pair.channel_1;
    side.n_spikes = pair.n_spikes_1;
else
    side.cluster_id = pair.cluster_id_2;
    side.partner_id = pair.cluster_id_1;
    side.sync_index = pair.sync_index_2;
    side.area = pair.area_2;
    side.channel = pair.channel_2;
    side.n_spikes = pair.n_spikes_2;
end
end

function coefficients = delta_fir_coefficients(options)
if exist('designfilt', 'file') ~= 2 && exist('designfilt', 'builtin') ~= 5
    error('GENERATE_SYNC_SPIKE_LFP_POWER:designfilt', ...
        'designfilt is required for the 1-5 Hz delta filter.');
end
    delta_filter = designfilt('bandpassfir', ...
    'FilterOrder', options.delta_filter_order, ...
    'CutoffFrequency1', options.delta_freq(1), ...
    'CutoffFrequency2', options.delta_freq(2), ...
    'SampleRate', options.lfp_sample_rate, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');
coefficients = delta_filter.Coefficients;
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
[~, ~, session_result.coverage] = select_synchronous_pairs(struct( ...
    'cluster_id', {}, 'channel', {}, 'area', {}, ...
    'spike_times', {}, 'spike_samples', {}, 'interval_index', {}), [], struct());
session_result.rows = empty_power_rows();
end

function rows = empty_power_rows()
rows = struct('session_id', {}, 'pair_index', {}, 'cluster_id', {}, ...
    'partner_id', {}, 'area', {}, 'channel', {}, ...
    'n_spikes', {}, 'n_interval', {}, 'n_sync', {}, 'n_phase_pool', {}, ...
    'min_phase_candidates', {}, 'phase_bin_width', {}, ...
    'mean_power_sync', {}, 'mean_power_phase_matched', {}, ...
    'null_mean_power', {}, 'percentile', {}, 'power_definition', {});
end

function row = power_row(session_id, pair_index, side, neuron, comparison, power_definition)
row.session_id = session_id;
row.pair_index = pair_index;
row.cluster_id = side.cluster_id;
row.partner_id = side.partner_id;
row.area = side.area;
row.channel = side.channel;
row.n_spikes = side.n_spikes;
row.n_interval = numel(neuron.interval_index);
row.n_sync = comparison.n_sync;
row.n_phase_pool = comparison.n_pool;
row.min_phase_candidates = comparison.min_candidates;
row.phase_bin_width = comparison.phase_bin_width;
row.mean_power_sync = comparison.mean_sync;
row.mean_power_phase_matched = comparison.mean_phase_matched;
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

function print_coverage(session_id, coverage, min_spikes)
fprintf(['%s  observed %d/%d neurons with >= %d spikes (%s), ' ...
    '%d/%d good units in SC, DLPAG, LPAG, VLPAG, and DR. Target %.0f%%.\n'], ...
    session_id, coverage.n_covered, coverage.n_with_spikes, min_spikes, ...
    percent_text(coverage.fraction_of_spiking), coverage.n_covered, ...
    coverage.n_area, 100 * coverage.target);
for index = 1:numel(coverage.areas)
    row = coverage.areas(index);
    fprintf('  %s  %d/%d with enough spikes (%s), %d good units\n', ...
        row.area, row.n_covered, row.n_with_spikes, ...
        percent_text(row.fraction_of_spiking), row.n_area);
end
if ~coverage.target_met
    fprintf(['  Short of %.0f%%: %d neurons have fewer than %d spikes, ' ...
        'and %d have enough spikes but no partner that can support a null.\n'], ...
        100 * coverage.target, coverage.n_below_min_spikes, min_spikes, ...
        coverage.n_spiking_without_partner);
end
end

function label = percent_text(fraction)
if ~isfinite(fraction)
    label = 'n/a';
else
    label = sprintf('%.0f%%', 100 * fraction);
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
