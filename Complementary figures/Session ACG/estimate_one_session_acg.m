function summary = estimate_one_session_acg(session_name, spikes, clusters, neurons, lfp, lfp_time, options, save_dir)
%ESTIMATE_ONE_SESSION_ACG  Autocorrelograms and 500 nulls for one session.
%
% neurons is a struct array with cluster_id, channel, area, and lock_group.
% channel is the 1-based PAG LFP row. lfp_time uses the same clock as the
% spike times. When save_dir is set, each neuron is written to
% neuron_<id>.mat after its shuffles, and a finished neuron is skipped on
% a resumed run. Channel spectra are written once to channel_power.mat.
%
% The returned summary lists the neurons estimated in this call. It does
% not hold the null matrices.

if nargin < 8
    save_dir = '';
end
settings = session_acg_settings(options);
options.bin_edges = settings.bin_edges;
options.smooth_n = settings.smooth_n;
if ~isfield(options, 'resume') || isempty(options.resume)
    options.resume = true;
end
if ~isfield(options, 'log_file')
    options.log_file = '';
end

summary = blank_summary(session_name);
if isempty(neurons)
    summary.n_neurons = 0;
    return
end
t0 = lfp_time(1);
t1 = lfp_time(end);
[starts, ends] = session_time_windows(t0, t1, options.window_s);
summary.n_windows = numel(starts);
summary.n_neurons = numel(neurons);

if ~isempty(save_dir) && exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end
if ~isempty(save_dir)
    assert_saved_settings(save_dir, settings, options.resume);
end
if isfield(options, 'channel_power') && ~isempty(options.channel_power)
    power = options.channel_power;
    if size(power.fooof_delta, 2) ~= numel(starts) || size(power.hilbert_amp, 2) ~= numel(starts)
        error('session_acg:power', ...
            'Precomputed channel power does not have one column per window.');
    end
elseif ~isempty(save_dir)
    power = cached_channel_power(save_dir, neurons, lfp, lfp_time, starts, ends, options, settings);
else
    power = channel_power_for_neurons(neurons, lfp, starts, ends, options);
end

n_estimated = 0;
n_skipped = 0;
rows = cell(numel(neurons), 1);
for n = 1:numel(neurons)
    neuron_id = neurons(n).cluster_id;
    neuron_file = '';
    if ~isempty(save_dir)
        neuron_file = fullfile(save_dir, sprintf('neuron_%d.mat', round(neuron_id)));
        if options.resume && exist(neuron_file, 'file') == 2
            n_skipped = n_skipped + 1;
            session_acg_log(options.log_file, sprintf('%s id %g already saved', session_name, neuron_id));
            continue
        end
    end
    these_spikes = spikes(clusters == neuron_id);
    label = sprintf('%s id %g (%s)', session_name, neuron_id, neurons(n).lock_group);
    session_acg_log(options.log_file, sprintf('%s: %d spikes', label, numel(these_spikes)));
    options.seed = neuron_shuffle_seed(session_name, neuron_id);
    options.label = label;
    tail_index = tails_for_channel(power, neurons(n).channel, numel(starts), options.tail_percent);
    packed = neuron_acg_with_nulls(these_spikes, t0, t1, starts, ends, ...
        tail_index, options.bin_edges, options);
    neuron = pack_neuron(session_name, neurons(n), packed, power, tail_index);
    if ~isempty(neuron_file)
        temp_file = [neuron_file, '.writing'];
        if exist(temp_file, 'file') == 2
            delete(temp_file);
        end
        save(temp_file, 'neuron', 'settings', '-v7');
        movefile(temp_file, neuron_file);
    end
    n_estimated = n_estimated + 1;
    rows{n_estimated} = neuron_row(neuron);
    session_acg_log(options.log_file, sprintf( ...
        '%s done: session p %.3g oscillates %d; FOOOF high p %.3g; Hilbert high p %.3g', ...
        label, neuron.scores.session.first_p, neuron.scores.session.oscillates, ...
        neuron.scores.fooof_high.first_p, neuron.scores.hilbert_high.first_p));
end
summary.n_estimated = n_estimated;
summary.n_skipped = n_skipped;
if n_estimated == 0
    summary.rows = empty_rows();
else
    summary.rows = [rows{1:n_estimated}];
end
end


function summary = blank_summary(session_name)
summary.session = session_name;
summary.n_neurons = 0;
summary.n_windows = 0;
summary.n_estimated = 0;
summary.n_skipped = 0;
summary.rows = empty_rows();
end


function rows = empty_rows()
rows = struct('cluster_id', {}, 'lock_group', {}, 'session_oscillates', {}, ...
    'session_p', {}, 'session_lag', {}, 'session_prominence', {});
end


function row = neuron_row(neuron)
row.cluster_id = neuron.cluster_id;
row.lock_group = neuron.lock_group;
row.session_oscillates = neuron.scores.session.oscillates;
row.session_p = neuron.scores.session.first_p;
row.session_lag = neuron.scores.session.first_lag;
row.session_prominence = neuron.scores.session.first_prominence;
end


function assert_saved_settings(save_dir, settings, resume)
param_file = fullfile(save_dir, 'parameters.mat');
existing = dir(fullfile(save_dir, 'neuron_*.mat'));
if exist(param_file, 'file') == 2
    saved = load(param_file, 'settings');
    if ~session_acg_parameters_match(saved.settings, settings)
        error('session_acg:parameters', ...
            ['%s was estimated with different settings. Keep resume on the ' ...
            'same settings, or move this folder aside before starting again.'], save_dir);
    end
elseif ~isempty(existing)
    error('session_acg:parameters', ...
        '%s has neuron files but no parameters.mat.', save_dir);
else
    save(param_file, 'settings');
end
if ~resume && ~isempty(existing)
    error('session_acg:resume', ...
        ['%s already has neuron files. resume = true continues that run. ' ...
        'Move the folder aside to estimate it again.'], save_dir);
end
end


function power = cached_channel_power(save_dir, neurons, lfp, lfp_time, starts, ends, options, settings)
power_file = fullfile(save_dir, 'channel_power.mat');
if exist(power_file, 'file') == 2
    saved = load(power_file, 'power', 'settings');
    if isfield(saved, 'settings') && session_acg_parameters_match(saved.settings, settings) ...
            && isfield(saved, 'power') && size(saved.power.fooof_delta, 2) == numel(starts)
        power = saved.power;
        return
    end
end
if isempty(lfp)
    error('session_acg:lfp', ...
        'Session %s has no channel spectra and no LFP to estimate them.', save_dir);
end
if numel(lfp_time) ~= size(lfp, 2)
    error('session_acg:lfp', 'LFP time and LFP width do not match.');
end
power = channel_power_for_neurons(neurons, lfp, starts, ends, options);
save(power_file, 'power', 'settings', '-v7');
end


function power = channel_power_for_neurons(neurons, lfp, starts, ends, options)
channels = [neurons.channel];
channels = unique(channels(isfinite(channels)));
if ~isempty(lfp)
    channels = channels(channels >= 1 & channels <= size(lfp, 1));
end
n_windows = numel(starts);
power.channel = channels(:)';
power.fooof_delta = nan(numel(channels), n_windows);
power.hilbert_amp = nan(numel(channels), n_windows);
if isempty(channels) || isempty(lfp)
    return
end
for c = 1:numel(channels)
    session_acg_log(options.log_file, sprintf('channel %g (%d/%d), %d windows', ...
        channels(c), c, numel(channels), n_windows));
    measured = channel_window_delta(lfp(channels(c), :), options.lfp_sample_rate, ...
        starts, ends, options);
    power.fooof_delta(c, :) = measured.fooof_delta;
    power.hilbert_amp(c, :) = measured.hilbert_amp;
end
end


function tail_index = tails_for_channel(power, channel, n_windows, tail_percent)
tail_index.fooof_low = false(n_windows, 1);
tail_index.fooof_high = false(n_windows, 1);
tail_index.hilbert_low = false(n_windows, 1);
tail_index.hilbert_high = false(n_windows, 1);
if n_windows < 1 || ~isfinite(channel) || ~isfield(power, 'channel')
    return
end
hit = find(power.channel == channel, 1);
if isempty(hit) || size(power.fooof_delta, 2) ~= n_windows
    return
end
[tail_index.fooof_low, tail_index.fooof_high] = power_tail_index( ...
    power.fooof_delta(hit, :), tail_percent);
[tail_index.hilbert_low, tail_index.hilbert_high] = power_tail_index( ...
    power.hilbert_amp(hit, :), tail_percent);
end


function neuron = pack_neuron(session_name, identity, packed, power, tail_index)
neuron.session = session_name;
neuron.cluster_id = identity.cluster_id;
neuron.channel = identity.channel;
neuron.area = identity.area;
neuron.lock_group = identity.lock_group;
neuron.n_spikes = packed.n_spikes;
neuron.firing_rate = packed.firing_rate;
neuron.duration = packed.duration;
neuron.window_spike_count = packed.window_spike_count;
neuron.window_acg = packed.window_acg;
neuron.curves = packed.curves;
neuron.scores = packed.scores;
hit = [];
if isfinite(identity.channel) && isfield(power, 'channel')
    hit = find(power.channel == identity.channel, 1);
end
if isempty(hit)
    neuron.fooof_delta = nan(1, numel(packed.window_spike_count));
    neuron.hilbert_amp = nan(1, numel(packed.window_spike_count));
else
    neuron.fooof_delta = power.fooof_delta(hit, :);
    neuron.hilbert_amp = power.hilbert_amp(hit, :);
end
neuron.fooof_low_index = tail_index.fooof_low;
neuron.fooof_high_index = tail_index.fooof_high;
neuron.hilbert_low_index = tail_index.hilbert_low;
neuron.hilbert_high_index = tail_index.hilbert_high;
neuron.tail_used = packed.tail_used;
end
