function session_result = GENERATE_SYNC_SPIKE_WAVEFORMS(npx_data_dir, options)
%GENERATE_SYNC_SPIKE_WAVEFORMS  Synchronous-spike waveforms for one session.
%
% Called by Estimate_sync_spike_waveforms. Do not start here.
%
% For each neuron in a play-bout coincident pair, every spike is read once
% from the continuous file (a second pair reuses that read). The template
% is the baseline-subtracted mean of all of those spikes. The synchronous
% mean uses every coincident spike. Five hundred null means use that same
% count of the neuron's other spikes. Distance is the mean absolute
% difference from the template, on the channel neighborhood and on the
% neuron channel alone. A high percentile means the synchronous waveform
% is more dissimilar than the null draws.
%
% The continuous path comes from NPX_raw_data_folders.xlsx, matched on the
% session folder name. Several paths are scored with a short read on the
% neuron channel; a flat mean stops the session. The null waveforms are
% not saved.

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
    first_neuron.spike_samples, first_neuron.channel, options.channel_range, ...
    file_choice(options));
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
        [cache, waveforms, channels_kept, time_samples, elapsed] = cached_waveforms( ...
            cache, selection.neurons, cluster_id, continuous_file, options);
        if ~timed_a_load
            timed_a_load = true;
            session_result.first_load_seconds = elapsed;
            fprintf('%s  cluster %g  read %d spikes in %.1f s\n', ...
                session_result.session_id, cluster_id, n_spikes, elapsed);
            if isfinite(options.max_first_load_seconds) ...
                    && options.max_first_load_seconds > 0 ...
                    && elapsed > options.max_first_load_seconds
                error('GENERATE_SYNC_SPIKE_WAVEFORMS:loadTime', ...
                    ['Cluster %g in %s took %.1f minutes to read %d spikes. ' ...
                     'Lower max_pairs, or raise max_first_load_seconds, before the full run.'], ...
                    cluster_id, session_result.session_id, elapsed / 60, n_spikes);
            end
        end
        comparison = compare_sync_waveforms(waveforms, time_samples, sync_index, ...
            channel, channels_kept, options.n_repeats);
        if ~comparison.ok
            note = sprintf('%s cluster %g pair %d: %s', ...
                session_result.session_id, cluster_id, pair_index, comparison.reason);
            session_result.notes{end + 1} = note; %#ok<AGROW>
            fprintf('  %s\n', note);
            continue
        end
        rows = append_waveform_row(rows, waveform_row( ...
            session_result.session_id, pair_index, cluster_id, partner_id, ...
            area, channel, n_spikes, comparison, continuous_file));
        fprintf(['%s  cluster %g  sync %d  distance %.3g  percentile %.1f' ...
            '  channel percentile %.1f\n'], ...
            session_result.session_id, cluster_id, comparison.n_sync, ...
            comparison.distance_sync, 100 * comparison.percentile, ...
            100 * comparison.percentile_channel);
    end
end
session_result.rows = rows;
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
if isfield(options, 'continuous_locations') && ~isempty(options.continuous_locations)
    locations = options.continuous_locations;
else
    locations = read_npx_continuous_locations(options.raw_folders_file);
end
candidates = {};
for index = 1:numel(locations)
    if ~strcmp(locations(index).session, session_id)
        continue
    end
    path = strtrim(locations(index).path);
    if isempty(path) || any(strcmp(candidates, path))
        continue
    end
    candidates{end + 1} = path; %#ok<AGROW>
end
if isempty(candidates)
    error('GENERATE_SYNC_SPIKE_WAVEFORMS:continuous', ...
        'No continuous path for %s in %s.', session_id, options.raw_folders_file);
end
end

function choice = file_choice(options)
choice.min_scale = options.min_spike_scale;
choice.n_channels = options.n_ap_channels;
choice.n_check = options.n_check_spikes;
choice.samples_before = options.samples_before;
choice.samples_after = options.samples_after;
end

function [cache, waveforms, channels_kept, time_samples, elapsed] = cached_waveforms( ...
    cache, neurons, cluster_id, continuous_file, options)
for index = 1:numel(cache)
    if cache(index).cluster_id == cluster_id
        waveforms = cache(index).waveforms;
        channels_kept = cache(index).channels_kept;
        time_samples = cache(index).time_samples;
        elapsed = NaN;
        return
    end
end
neuron = neuron_by_id(neurons, cluster_id);
if isempty(neuron)
    error('GENERATE_SYNC_SPIKE_WAVEFORMS:neuron', ...
        'Cluster %g was selected but has no spike list.', cluster_id);
end
started = tic;
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    continuous_file, neuron.spike_samples, neuron.channel, options.channel_range, ...
    'NChannels', options.n_ap_channels, ...
    'SamplesBefore', options.samples_before, ...
    'SamplesAfter', options.samples_after);
elapsed = toc(started);
cache(end + 1).cluster_id = cluster_id; %#ok<AGROW>
cache(end).waveforms = waveforms;
cache(end).channels_kept = channels_kept;
cache(end).time_samples = time_samples;
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
session_result.rows = empty_waveform_rows();
end

function options = stored_options(options)
if isfield(options, 'continuous_locations')
    options = rmfield(options, 'continuous_locations');
end
end

function rows = empty_waveform_rows()
rows = struct('session_id', {}, 'pair_index', {}, 'cluster_id', {}, ...
    'partner_id', {}, 'area', {}, 'channel', {}, 'n_spikes', {}, 'n_sync', {}, ...
    'channels_kept', {}, 'time_samples', {}, 'template_mean', {}, 'sync_mean', {}, ...
    'distance_sync', {}, 'null_distances', {}, 'percentile', {}, ...
    'distance_sync_channel', {}, 'null_distances_channel', {}, ...
    'percentile_channel', {}, 'continuous_file', {});
end

function row = waveform_row(session_id, pair_index, cluster_id, partner_id, ...
    area, channel, n_spikes, comparison, continuous_file)
row.session_id = session_id;
row.pair_index = pair_index;
row.cluster_id = cluster_id;
row.partner_id = partner_id;
row.area = area;
row.channel = channel;
row.n_spikes = n_spikes;
row.n_sync = comparison.n_sync;
row.channels_kept = comparison.channels_kept;
row.time_samples = comparison.time_samples;
row.template_mean = comparison.template_mean;
row.sync_mean = comparison.sync_mean;
row.distance_sync = comparison.distance_sync;
row.null_distances = comparison.null_distances;
row.percentile = comparison.percentile;
row.distance_sync_channel = comparison.distance_sync_channel;
row.null_distances_channel = comparison.null_distances_channel;
row.percentile_channel = comparison.percentile_channel;
row.continuous_file = continuous_file;
end

function rows = append_waveform_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end + 1) = orderfields(row, rows); %#ok<AGROW>
end
end

function print_pairs(session_id, pairs)
for pair_index = 1:numel(pairs)
    pair = pairs(pair_index);
    fprintf(['%s  %s ch %g (%d spikes, %d sync)  with  ' ...
        '%s ch %g (%d spikes, %d sync)\n'], ...
        session_id, pair.area_1, pair.channel_1, pair.n_spikes_1, pair.n_sync_1, ...
        pair.area_2, pair.channel_2, pair.n_spikes_2, pair.n_sync_2);
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
