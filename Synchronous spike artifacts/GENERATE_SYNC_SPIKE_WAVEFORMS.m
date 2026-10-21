function session_result = GENERATE_SYNC_SPIKE_WAVEFORMS(npx_data_dir, options)
%GENERATE_SYNC_SPIKE_WAVEFORMS  Synchronous-spike waveforms for one session.
%
% Called by Estimate_sync_spike_waveforms. Do not start here.
%
% The comparison pool for each neuron is the spikes inside the play-bout
% windows, [bout_start - 1.5 s, bout_start + 1 s] for every play bout.
% Only those snippets are read from the continuous file. A second pair of
% the same neuron reuses that read. Synchronous spikes of this pair are
% excluded from the 500 draws; a spike that is synchronous only with a
% different partner can still be drawn.
%
% For each neuron in each pair the saved waveform detail is:
%   sync_mean            mean around the synchronous spikes
%   nonsync_mean         mean of the other spikes in the windows
%   percentile_time      fraction of the 500 null means strictly below
%                        sync_mean, one value per channel and time
%   percentile_abs_time  fraction of null means closer to nonsync_mean
%   null_p025, null_p975 pointwise null envelope
% The scalar distance of sync_mean from nonsync_mean, and the 500 scalar
% null distances, are saved as well. The 500 waveform tensors are not.
%
% session_result.rows is one struct per neuron in a pair.
% session_result.pairs keeps the two neurons of each pair together.
%
% The continuous path is column 2 of NPX_raw_data_folders.xlsx on the row
% whose SessionName is this NPX folder. A flat mean stops the session.

if nargin < 2 || isempty(options)
    options = struct();
end
options = sync_analysis_options(options);
[~, session_id] = fileparts(npx_data_dir);
session_result = blank_waveform_session(session_id, options);
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

first_neuron = neuron_by_id(selection.neurons, selection.pairs(1).cluster_id_1);
[continuous_file, checks] = choose_continuous_file( ...
    continuous_candidates(session_result.session_id, options), ...
    first_neuron.spike_samples(first_neuron.interval_index), ...
    first_neuron.channel, options.channel_range, file_choice(options));
for check_index = 1:numel(checks)
    if isempty(checks(check_index).message)
        fprintf('  %s  spike scale %.3g\n', checks(check_index).path, checks(check_index).scale);
    else
        fprintf('  %s  %s\n', checks(check_index).path, checks(check_index).message);
    end
end
fprintf('%s  continuous file %s\n', session_result.session_id, continuous_file);
session_result.continuous_file = continuous_file;
session_result.n_pairs = numel(selection.pairs);

cache = struct('cluster_id', {}, 'waveforms', {}, 'channels_kept', {}, 'time_samples', {});
timed_a_load = false;
rows = empty_waveform_rows();
pair_records = empty_waveform_pairs();
for pair_index = 1:numel(selection.pairs)
    pair = selection.pairs(pair_index);
    record = blank_pair_record(session_result.session_id, pair_index);
    sides = [1 2];
    for side_index = 1:numel(sides)
        if sides(side_index) == 1
            cluster_id = pair.cluster_id_1;
            partner_id = pair.cluster_id_2;
            sync_index = pair.sync_index_1;
            area = pair.area_1;
            channel = pair.channel_1;
            n_spikes = pair.n_spikes_1;
            n_interval = pair.n_interval_1;
        else
            cluster_id = pair.cluster_id_2;
            partner_id = pair.cluster_id_1;
            sync_index = pair.sync_index_2;
            area = pair.area_2;
            channel = pair.channel_2;
            n_spikes = pair.n_spikes_2;
            n_interval = pair.n_interval_2;
        end
        neuron = neuron_by_id(selection.neurons, cluster_id);
        local_sync = local_sync_index(sync_index, neuron.interval_index);
        [cache, waveforms, channels_kept, time_samples, elapsed] = cached_waveforms( ...
            cache, neuron, continuous_file, options);
        if ~timed_a_load
            timed_a_load = true;
            session_result.first_load_seconds = elapsed;
            fprintf('%s  cluster %g  read %d interval spikes in %.1f s\n', ...
                session_result.session_id, cluster_id, numel(neuron.interval_index), elapsed);
            if isfinite(options.max_first_load_seconds) ...
                    && options.max_first_load_seconds > 0 ...
                    && elapsed > options.max_first_load_seconds
                error('GENERATE_SYNC_SPIKE_WAVEFORMS:loadTime', ...
                    ['Cluster %g in %s took %.1f minutes to read %d interval spikes. ' ...
                     'Lower max_pairs, or raise max_first_load_seconds, before the full run.'], ...
                    cluster_id, session_result.session_id, elapsed / 60, ...
                    numel(neuron.interval_index));
            end
        end
        comparison = compare_sync_waveforms(waveforms, time_samples, local_sync, ...
            channel, channels_kept, options.n_repeats);
        side = waveform_side(cluster_id, partner_id, area, channel, ...
            n_spikes, n_interval, comparison, continuous_file);
        if sides(side_index) == 1
            record.neuron_1 = side;
        else
            record.neuron_2 = side;
        end
        if ~comparison.ok
            note = sprintf('%s cluster %g pair %d: %s', ...
                session_result.session_id, cluster_id, pair_index, comparison.reason);
            session_result.notes{end + 1} = note; %#ok<AGROW>
            fprintf('  %s\n', note);
            continue
        end
        rows = append_waveform_row(rows, waveform_row( ...
            session_result.session_id, pair_index, side));
        fprintf(['%s  cluster %g  interval %d  sync %d  nonsync %d  ' ...
            'distance %.3g  percentile %.1f  channel percentile %.1f  ' ...
            'extreme sample %.0f (abs percentile %.1f)\n'], ...
            session_result.session_id, cluster_id, n_interval, comparison.n_sync, ...
            comparison.n_nonsync, comparison.distance_sync, ...
            100 * comparison.percentile, 100 * comparison.percentile_channel, ...
            comparison.extreme_time_samples, 100 * comparison.extreme_abs_percentile);
    end
    pair_records = append_waveform_pair(pair_records, record);
end
session_result.rows = rows;
session_result.pairs = pair_records;
if isempty(rows)
    session_result.status = 'no_saved_neuron';
else
    session_result.status = 'ok';
end
end

function candidates = continuous_candidates(session_id, options)
if isfield(options, 'continuous_file') && ~isempty(options.continuous_file)
    candidates = {options.continuous_file};
    return
end
candidates = read_npx_continuous_locations(options.raw_folders_file, session_id);
end

function choice = file_choice(options)
choice.min_scale = options.min_spike_scale;
choice.n_channels = options.n_ap_channels;
choice.n_check = options.n_check_spikes;
choice.samples_before = options.samples_before;
choice.samples_after = options.samples_after;
end

function [cache, waveforms, channels_kept, time_samples, elapsed] = cached_waveforms( ...
    cache, neuron, continuous_file, options)
for index = 1:numel(cache)
    if cache(index).cluster_id == neuron.cluster_id
        waveforms = cache(index).waveforms;
        channels_kept = cache(index).channels_kept;
        time_samples = cache(index).time_samples;
        elapsed = NaN;
        return
    end
end
started = tic;
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    continuous_file, neuron.spike_samples(neuron.interval_index), neuron.channel, ...
    options.channel_range, ...
    'NChannels', options.n_ap_channels, ...
    'SamplesBefore', options.samples_before, ...
    'SamplesAfter', options.samples_after);
elapsed = toc(started);
cache(end + 1).cluster_id = neuron.cluster_id; %#ok<AGROW>
cache(end).waveforms = waveforms;
cache(end).channels_kept = channels_kept;
cache(end).time_samples = time_samples;
end

function local = local_sync_index(sync_index, interval_index)
[is_in, local] = ismember(sync_index(:), interval_index(:));
local = local(is_in);
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

function session_result = blank_waveform_session(session_id, options)
session_result.session_id = session_id;
session_result.status = '';
session_result.n_pairs = 0;
session_result.continuous_file = '';
session_result.first_load_seconds = NaN;
session_result.notes = {};
session_result.options = stored_options(options);
session_result.coverage = blank_session_coverage();
session_result.rows = empty_waveform_rows();
session_result.pairs = empty_waveform_pairs();
end

function coverage = blank_session_coverage()
[~, ~, coverage] = select_synchronous_pairs(struct( ...
    'cluster_id', {}, 'channel', {}, 'area', {}, ...
    'spike_times', {}, 'spike_samples', {}, 'interval_index', {}), [], struct());
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

function options = stored_options(options)
if isfield(options, 'continuous_locations')
    options = rmfield(options, 'continuous_locations');
end
end

function rows = empty_waveform_rows()
rows = struct('session_id', {}, 'pair_index', {}, 'cluster_id', {}, ...
    'partner_id', {}, 'area', {}, 'channel', {}, ...
    'n_spikes', {}, 'n_interval', {}, 'n_sync', {}, 'n_nonsync', {}, ...
    'channels_kept', {}, 'time_samples', {}, ...
    'nonsync_mean', {}, 'sync_mean', {}, ...
    'percentile_time', {}, 'percentile_abs_time', {}, ...
    'null_p025', {}, 'null_p975', {}, ...
    'distance_sync', {}, 'null_distances', {}, 'percentile', {}, ...
    'distance_sync_channel', {}, 'null_distances_channel', {}, ...
    'percentile_channel', {}, ...
    'extreme_time_samples', {}, 'extreme_percentile', {}, ...
    'extreme_abs_percentile', {}, 'neuron_channel_tail_fraction', {}, ...
    'continuous_file', {});
end

function row = waveform_row(session_id, pair_index, side)
row.session_id = session_id;
row.pair_index = pair_index;
row.cluster_id = side.cluster_id;
row.partner_id = side.partner_id;
row.area = side.area;
row.channel = side.channel;
row.n_spikes = side.n_spikes;
row.n_interval = side.n_interval;
row.n_sync = side.n_sync;
row.n_nonsync = side.n_nonsync;
row.channels_kept = side.channels_kept;
row.time_samples = side.time_samples;
row.nonsync_mean = side.nonsync_mean;
row.sync_mean = side.sync_mean;
row.percentile_time = side.percentile_time;
row.percentile_abs_time = side.percentile_abs_time;
row.null_p025 = side.null_p025;
row.null_p975 = side.null_p975;
row.distance_sync = side.distance_sync;
row.null_distances = side.null_distances;
row.percentile = side.percentile;
row.distance_sync_channel = side.distance_sync_channel;
row.null_distances_channel = side.null_distances_channel;
row.percentile_channel = side.percentile_channel;
row.extreme_time_samples = side.extreme_time_samples;
row.extreme_percentile = side.extreme_percentile;
row.extreme_abs_percentile = side.extreme_abs_percentile;
row.neuron_channel_tail_fraction = side.neuron_channel_tail_fraction;
row.continuous_file = side.continuous_file;
end

function rows = append_waveform_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end + 1) = orderfields(row, rows); %#ok<AGROW>
end
end

function pairs = empty_waveform_pairs()
pairs = struct('session_id', {}, 'pair_index', {}, 'neuron_1', {}, 'neuron_2', {});
end

function record = blank_pair_record(session_id, pair_index)
record.session_id = session_id;
record.pair_index = pair_index;
record.neuron_1 = empty_waveform_side();
record.neuron_2 = empty_waveform_side();
end

function side = empty_waveform_side()
side.ok = false;
side.reason = '';
side.cluster_id = NaN;
side.partner_id = NaN;
side.area = '';
side.channel = NaN;
side.n_spikes = 0;
side.n_interval = 0;
side.n_sync = 0;
side.n_nonsync = 0;
side.channels_kept = [];
side.time_samples = [];
side.nonsync_mean = [];
side.sync_mean = [];
side.percentile_time = [];
side.percentile_abs_time = [];
side.null_p025 = [];
side.null_p975 = [];
side.distance_sync = NaN;
side.null_distances = [];
side.percentile = NaN;
side.distance_sync_channel = NaN;
side.null_distances_channel = [];
side.percentile_channel = NaN;
side.extreme_time_samples = NaN;
side.extreme_percentile = NaN;
side.extreme_abs_percentile = NaN;
side.neuron_channel_tail_fraction = NaN;
side.continuous_file = '';
end

function side = waveform_side(cluster_id, partner_id, area, channel, ...
    n_spikes, n_interval, comparison, continuous_file)
side = empty_waveform_side();
side.ok = comparison.ok;
side.reason = comparison.reason;
side.cluster_id = cluster_id;
side.partner_id = partner_id;
side.area = area;
side.channel = channel;
side.n_spikes = n_spikes;
side.n_interval = n_interval;
side.n_sync = comparison.n_sync;
side.n_nonsync = comparison.n_nonsync;
side.channels_kept = comparison.channels_kept;
side.time_samples = comparison.time_samples;
side.nonsync_mean = comparison.nonsync_mean;
side.sync_mean = comparison.sync_mean;
side.percentile_time = comparison.percentile_time;
side.percentile_abs_time = comparison.percentile_abs_time;
side.null_p025 = comparison.null_p025;
side.null_p975 = comparison.null_p975;
side.distance_sync = comparison.distance_sync;
side.null_distances = comparison.null_distances;
side.percentile = comparison.percentile;
side.distance_sync_channel = comparison.distance_sync_channel;
side.null_distances_channel = comparison.null_distances_channel;
side.percentile_channel = comparison.percentile_channel;
side.extreme_time_samples = comparison.extreme_time_samples;
side.extreme_percentile = comparison.extreme_percentile;
side.extreme_abs_percentile = comparison.extreme_abs_percentile;
side.neuron_channel_tail_fraction = comparison.neuron_channel_tail_fraction;
side.continuous_file = continuous_file;
end

function pairs = append_waveform_pair(pairs, pair)
if isempty(pairs)
    pairs = pair;
else
    pairs(end + 1) = orderfields(pair, pairs); %#ok<AGROW>
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
