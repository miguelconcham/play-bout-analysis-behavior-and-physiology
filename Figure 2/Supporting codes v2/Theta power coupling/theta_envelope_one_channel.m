function recording = theta_envelope_one_channel(lfp, selected, loaded, options)
%THETA_ENVELOPE_ONE_CHANNEL  Session spectra and behavior-aligned amplitude.
%
% The LFP is decimated and split at partner-session boundaries. Session
% estimates use the Hilbert envelope pooled across those intervals. Each
% annotated behavior is aligned separately, for the implanted animal and for
% the partner. Amplitude is the Hilbert envelope and the 6-12 Hz spectrogram.

if isfinite(options.random_seed)
    rng(options.random_seed + round(selected.lfp_row) + ...
        sum(double(loaded.animal_code)), 'twister');
end
rate_out = options.envelope_rate;
rate_in = options.sampling_rate;
decimated = theta_envelope_decimate(lfp, rate_in, rate_out);
[segments, full_envelope, full_power] = interval_signals( ...
    decimated, loaded.partner_ranges, options);
recording = [];
if isempty(segments)
    warning('theta_envelope_one_channel:interval', ...
        '%s channel %g: no interval long enough for the theta envelope.', ...
        loaded.animal_code, selected.channel);
    return
end

grids = struct('test', options.grid_test, 'fine', options.grid_fine);
welch = struct('window_seconds', options.welch_window_seconds, ...
    'overlap', options.welch_overlap, ...
    'resolution_hz', options.welch_resolution_hz);
envelopes = {segments.envelope};
lfps = {segments.lfp};
spectrum = theta_envelope_session_spectrum(envelopes, rate_out, grids, welch);
coherence = theta_envelope_session_coherence( ...
    envelopes, lfps, rate_out, grids, welch, false);
power_null = theta_envelope_spectrum_null(lfps, options.theta_coefficients, ...
    spectrum, rate_out, grids, welch, options.n_perm);
coherence_null = theta_envelope_coherence_null(envelopes, lfps, coherence, ...
    rate_out, grids, welch, options.n_perm);

catalog = theta_envelope_behavior_events(loaded.behavior, ...
    loaded.repeated_animal, options.play_behaviors, options.other_behaviors);
intervals = [[segments.t0].', [segments.t1].'];
behavior = align_behaviors(catalog, full_envelope, full_power, rate_out, ...
    intervals, rate_in, options, loaded, selected);
recording = pack_recording(loaded, selected, options, spectrum, coherence, ...
    power_null, coherence_null, behavior);
end

function [segments, full_envelope, full_power] = interval_signals(decimated, ranges, options)
full_envelope = nan(numel(decimated), 1);
full_power = nan(numel(decimated), 1);
segments = struct('envelope', {}, 'lfp', {}, 't0', {}, 't1', {});
rate_out = options.envelope_rate;
rate_in = options.sampling_rate;
if isempty(ranges)
    return
end
for interval_n = 1:size(ranges, 1)
    i1 = theta_envelope_time_index(ranges(interval_n, 1), rate_out, rate_in);
    i2 = theta_envelope_time_index(ranges(interval_n, 2), rate_out, rate_in);
    i1 = max(1, i1);
    i2 = min(numel(decimated), i2);
    if i2 - i1 < 10
        continue
    end
    lfp_segment = decimated(i1:i2);
    envelope = theta_envelope_signal(lfp_segment, options.theta_coefficients);
    if ~any(isfinite(envelope))
        continue
    end
    power = theta_envelope_band_power(lfp_segment, rate_out, ...
        options.theta_freq, options.power_window_seconds);
    full_envelope(i1:i2) = envelope;
    full_power(i1:i2) = power;
    segments(end + 1).envelope = envelope; %#ok<AGROW>
    segments(end).lfp = lfp_segment;
    segments(end).t0 = ranges(interval_n, 1);
    segments(end).t1 = ranges(interval_n, 2);
end
end

function behavior = align_behaviors(catalog, full_envelope, full_power, ...
        rate_out, intervals, rate_in, options, loaded, selected)
actors = {'self', 'partner'};
signals = {'envelope', 'power'};
traces = {full_envelope, full_power};
behavior = repmat(empty_behavior(), numel(catalog), 1);
for behavior_n = 1:numel(catalog)
    behavior(behavior_n).name = catalog(behavior_n).name;
    behavior(behavior_n).is_play = catalog(behavior_n).is_play;
    onsets = {catalog(behavior_n).self_onsets, catalog(behavior_n).partner_onsets};
    for actor_n = 1:numel(actors)
        actor = struct();
        actor.onsets = onsets{actor_n}(:);
        seed = behavior_seed(options, loaded, selected, ...
            catalog(behavior_n).name, actors{actor_n});
        for measure_n = 1:numel(signals)
            rng(seed, 'twister');
            actor.(signals{measure_n}) = one_measure(traces{measure_n}, ...
                rate_out, actor.onsets, intervals, rate_in, options);
        end
        behavior(behavior_n).(actors{actor_n}) = actor;
    end
end
end

function block = one_measure(signal, rate_out, onsets, intervals, rate_in, options)
n_freq = numel(options.grid_test.centers);
pre = repmat(empty_pre(n_freq, options.n_perm), size(options.pre_windows, 1), 1);
null_options = struct('min_shift_seconds', options.min_shift_seconds);
for window_n = 1:size(options.pre_windows, 1)
    null_options.n_perm = options.n_perm * (window_n == 1);
    null_options.compute_trial_summaries = window_n == 1;
    measured = theta_envelope_onset_null(signal, rate_out, onsets, ...
        options.pre_windows(window_n, :), options.grid_test.centers, ...
        options.grid_test.confirmatory, intervals, rate_in, null_options);
    pre(window_n) = pack_pre(measured);
end
[epochs, epoch_time] = theta_envelope_epochs(signal, rate_out, onsets, ...
    options.epoch_window, rate_in);
valid = all(isfinite(epochs), 2);
if any(valid)
    mean_epoch = mean(epochs(valid, :), 1);
else
    mean_epoch = nan(1, numel(epoch_time));
end
block = struct('epoch_time', epoch_time(:).', 'mean_epoch', mean_epoch, 'pre', pre);
end

function seed = behavior_seed(options, loaded, selected, behavior_name, actor)
base = 0;
if isfinite(options.random_seed)
    base = options.random_seed;
end
label = [loaded.animal_code, behavior_name, actor];
seed = mod(base + round(selected.lfp_row) + sum(double(label)), 2 ^ 31 - 1);
end

function item = empty_behavior()
item = struct('name', '', 'is_play', false, 'self', empty_actor(), 'partner', empty_actor());
end

function actor = empty_actor()
actor = struct('onsets', zeros(0, 1), ...
    'envelope', empty_measure(), 'power', empty_measure());
end

function block = empty_measure()
block = struct('epoch_time', zeros(1, 0), 'mean_epoch', zeros(1, 0), ...
    'pre', empty_pre(0, 0));
end

function recording = pack_recording(loaded, selected, options, spectrum, ...
        coherence, power_null, coherence_null, behavior)
recording.animal_name = loaded.animal_code;
recording.animal_area_table_name = loaded.table_name;
recording.repeated_animal = loaded.repeated_animal;
recording.lpag_channel = selected.channel;
recording.lpag_lfp_row = selected.lfp_row;
recording.lpag_probe_number = selected.probe_number;
recording.envelope_rate = options.envelope_rate;
recording.sampling_rate = options.sampling_rate;
recording.frequencies_test = options.grid_test.centers;
recording.frequencies_fine = options.grid_fine.centers;
recording.confirmatory_mask = options.grid_test.confirmatory;
recording.session = struct();
recording.session.frequencies = spectrum.frequencies;
recording.session.psd = spectrum.psd;
recording.session.coherence_frequencies = coherence.frequencies;
recording.session.coherence_spectrum = coherence.coherence;
recording.session.resolution_hz = spectrum.resolution_hz;
recording.session.n_segments = spectrum.n_segments;
recording.session.n_coherence_segments = coherence.n_segments;
recording.session.band_power_test = spectrum.band_power_test;
recording.session.band_power_fine = spectrum.band_power_fine;
recording.session.band_coherence_test = coherence.band_coherence_test;
recording.session.band_coherence_fine = coherence.band_coherence_fine;
recording.session.power_z_test = power_null.test.z;
recording.session.power_z_fine = power_null.fine.z;
recording.session.power_p_test = power_null.test.p;
recording.session.coherence_z_test = coherence_null.test.z;
recording.session.coherence_z_fine = coherence_null.fine.z;
recording.session.coherence_p_test = coherence_null.test.p;
recording.session.peak_power_hz = peak_hz(spectrum.band_power_fine, options.grid_fine.centers);
recording.session.peak_coherence_hz = peak_hz( ...
    coherence.band_coherence_fine, options.grid_fine.centers);
recording.session.peak_power_z_hz = peak_hz(power_null.fine.z, options.grid_fine.centers);
recording.session.peak_coherence_z_hz = peak_hz( ...
    coherence_null.fine.z, options.grid_fine.centers);
recording.behavior = behavior;
end

function pre = pack_pre(measured)
pre = struct( ...
    'window', measured.window, ...
    'coefficient', measured.coefficient, ...
    'n_valid', measured.n_valid, ...
    'mean_coefficient', measured.mean_coefficient, ...
    'ppc', measured.ppc, ...
    'loo_r', measured.loo_r, ...
    'loo_mean_r', measured.loo_mean_r, ...
    'phase_fraction', measured.phase_fraction, ...
    'null_mean_coefficient', measured.null_mean_coefficient, ...
    'null_ppc', measured.null_ppc, ...
    'null_loo_mean_r', measured.null_loo_mean_r, ...
    'null_phase_fraction', measured.null_phase_fraction);
end

function pre = empty_pre(n_freq, n_perm)
pre = struct( ...
    'window', [NaN NaN], ...
    'coefficient', zeros(0, n_freq), ...
    'n_valid', 0, ...
    'mean_coefficient', nan(1, n_freq), ...
    'ppc', nan(1, n_freq), ...
    'loo_r', zeros(0, 1), ...
    'loo_mean_r', NaN, ...
    'phase_fraction', NaN, ...
    'null_mean_coefficient', nan(n_perm, n_freq), ...
    'null_ppc', nan(n_perm, n_freq), ...
    'null_loo_mean_r', nan(n_perm, 1), ...
    'null_phase_fraction', nan(n_perm, 1));
end

function hz = peak_hz(z, centers)
z = double(z(:).');
if all(~isfinite(z))
    hz = NaN;
    return
end
z(~isfinite(z)) = -inf;
[~, idx] = max(z);
hz = centers(idx);
end
