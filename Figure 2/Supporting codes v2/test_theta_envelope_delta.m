function test_theta_envelope_delta()
%TEST_THETA_ENVELOPE_DELTA  Synthetic checks that do not need experiment files.
%
% Run this file by name after add_repo_paths. It checks the frequency grid,
% the surrogate spectrum, play-onset selection, the pre-onset Fourier
% measures, and a theta carrier whose amplitude oscillates at 2.5 Hz.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end
if exist('pkg', 'file') == 2
    try
        pkg load signal
    catch
    end
end
rng(1, 'twister');

test_options_and_grid();
test_phase_randomize();
test_play_events();
test_fourier_and_consistency();
test_signflip_and_cluster();
test_group_summaries();
test_onset_shift_null();
if signal_ready()
    test_sine_spectrum_and_coherence();
    test_modulated_theta();
    test_decimate_keeps_theta();
else
    disp('filtfilt or hilbert is missing; filter checks were skipped.');
end
disp('theta_envelope_delta tests passed');
end

function test_options_and_grid()
options = theta_envelope_delta_options(struct('n_perm', 0, 'min_events', 0));
assert(options.n_perm == 0, 'A permutation count of 0 was replaced.');
assert(options.min_events == 0, 'A minimum of 0 events was replaced.');
assert(isequal(options.modulation_band, [0.1 6]), 'The scan is not 0.1-6 Hz.');
assert(options.bandwidth_hz == 0.5, 'The window is not 0.5 Hz.');
assert(options.test_step_hz == 0.5, 'The test step is not 0.5 Hz.');
assert(options.fine_step_hz == 0.2, 'The fine step is not 0.2 Hz.');

test_grid = theta_envelope_frequency_grid([0.1 6], 0.5, 0.5, [2 3]);
fine_grid = theta_envelope_frequency_grid([0.1 6], 0.2, 0.5, [2 3]);
assert(test_grid.centers(1) == 0.1 && test_grid.centers(end) == 6, ...
    'The 0.5 Hz grid does not run from 0.1 to 6 Hz.');
assert(fine_grid.centers(1) == 0.1 && fine_grid.centers(end) == 6, ...
    'The 0.2 Hz grid does not run from 0.1 to 6 Hz.');
assert(numel(fine_grid.centers) > numel(test_grid.centers), ...
    'The fine grid is not finer than the test grid.');
assert(any(test_grid.centers >= 2 & test_grid.centers <= 3 & test_grid.confirmatory), ...
    '2-3 Hz is not marked on the test grid.');
flat = ones(100, 1);
frequencies = linspace(0.1, 6, 100).';
band = theta_envelope_band_mean(flat, frequencies, test_grid);
assert(max(abs(band - 1)) < 1e-12, 'A flat spectrum did not average to 1.');
end

function test_phase_randomize()
rng(2, 'twister');
x = randn(1000, 1);
y = theta_envelope_phase_randomize(x);
scale = max(abs(fft(x)));
assert(max(abs(abs(fft(x)) - abs(fft(y)))) / scale < 1e-6, ...
    'Phase randomization changed the Fourier magnitudes.');
assert(max(abs(imag(y))) < 1e-8, 'The surrogate is not real.');
end

function test_play_events()
behavior = struct();
behavior.Animal = {'Dual'; 'Dual'; 'Partner'; 'Dual'; 'Session_structure'; 'Dual'};
behavior.Start = [10; 12; 10.5; 30; 0; 11];
behavior.End = behavior.Start + 0.4;
behavior.Type = {'Pounce_A'; 'Sniffing'; 'Pounce'; 'Boxing'; 'Partners session'; 'Pin'};
events = theta_envelope_play_events(behavior, 'Dual', ...
    {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD'}, 2);
assert(numel(events.onsets) == 3, 'Expected three recorded-animal play onsets.');
assert(strcmp(events.behavior_type{1}, 'Pounce'), 'Pounce_A was not renamed.');
assert(events.clean_pre(1), 'The first onset should have a clean pre-window.');
assert(~events.clean_pre(3), 'An onset 1 s after play should not be clean.');
assert(events.clean_pre(2), 'The isolated boxing onset should be clean.');
end

function test_fourier_and_consistency()
rate = 50;
time = (-2:1 / rate:0).';
epoch = sin(2 * pi * 2.5 * time).';
epochs = repmat(epoch, 6, 1);
coefficient = theta_envelope_fourier(epochs, time, [2.5 4.5]);
assert(abs(coefficient(1, 1)) > 0.8, 'A 2.5 Hz sine was not recovered.');
assert(abs(coefficient(1, 2)) < 0.2, 'A 2.5 Hz sine leaked into 4.5 Hz.');
ppc = theta_envelope_phase_consistency(coefficient);
assert(ppc(1) > 0.99, 'Identical phases should have pairwise consistency near 1.');
flipped = coefficient;
flipped(4:6, 1) = -flipped(4:6, 1);
ppc_cancel = theta_envelope_phase_consistency(flipped);
assert(abs(ppc_cancel(1)) < 0.25, 'Opposite phases should cancel.');
correlation = theta_envelope_leave_one_out(epochs);
assert(all(correlation > 0.99), 'Identical windows should match the leave-one-out mean.');
epochs(6, :) = -epochs(6, :);
correlation = theta_envelope_leave_one_out(epochs);
assert(correlation(6) < -0.9, 'A flipped window should correlate negatively.');
end

function test_signflip_and_cluster()
rng(3, 'twister');
positive = theta_envelope_signflip([1; 1; 1; 1], 300);
negative = theta_envelope_signflip([-1; -1; -1; -1], 300);
assert(positive.p_greater < 0.2, 'A positive animal mean was not detected.');
assert(negative.p_greater > 0.8, 'A negative animal mean was treated as an increase.');
values = repmat([0 0 2 2 2 0], 5, 1);
cluster = theta_envelope_cluster_signflip(values, 200);
assert(numel(cluster.clusters) == 1, 'The positive bump was split or missed.');
assert(cluster.clusters(1).start == 3 && cluster.clusters(1).stop == 5, ...
    'The cluster is not the 2 Hz bump.');
assert(cluster.p < 0.2, 'The cluster sign-flip did not detect the bump.');
end

function test_group_summaries()
frequencies = [2.1 2.6 3.1];
mask = [true true true];
first = fake_recording('Batch1Dual1', [2 2 2], frequencies, mask, 1, 20);
second = fake_recording('Batch1Dual1', [4 0 0], frequencies, mask, 1, 20);
other = fake_recording('Batch2Single', [0 0 0], frequencies, mask, 0, 3);
sessions = theta_envelope_group_sessions([first, second, other]);
dual = sessions.power_z_test(strcmp(sessions.animals, 'Batch1Dual1'), :);
assert(max(abs(dual - [3 1 1])) < 1e-12, 'Within-animal z-scores were not averaged.');
onsets = theta_envelope_group_onsets([first, second, other], 15);
assert(onsets.included(strcmp(onsets.animals, 'Batch1Dual1')), ...
    'An animal with enough onsets was excluded.');
assert(~onsets.included(strcmp(onsets.animals, 'Batch2Single')), ...
    'An animal below the onset minimum was included.');
assert(onsets.locked_p == 0, 'A null of zeros should give p = 0 for a locked animal.');
end

function test_onset_shift_null()
rng(4, 'twister');
rate = 50;
n = rate * 50;
time = (1:n).' / rate;
onsets = (8:6:38).';
envelope = zeros(n, 1);
relative = (-2:1 / rate:0).';
for onset = onsets.'
    idx = theta_envelope_time_index(onset + relative, rate, rate);
    envelope(idx) = sin(2 * pi * 2.5 * relative);
end
measured = theta_envelope_onset_null(envelope, rate, onsets, [-2 0], ...
    [2.5 4.5], [true false], [time(1), time(end)], rate, ...
    struct('n_perm', 30, 'min_shift_seconds', 3, 'compute_trial_summaries', true));
assert(measured.n_valid == numel(onsets), 'Phase-locked windows were dropped.');
assert(abs(measured.mean_coefficient(1)) > 0.8, ...
    'The 2.5 Hz coefficient of identical windows is too small.');
assert(abs(measured.mean_coefficient(1)) > 4 * abs(measured.mean_coefficient(2)), ...
    'The onset waveform is not selective for 2.5 Hz.');
assert(measured.ppc(1) > 0.9, 'Phase consistency of identical windows is too low.');
assert(measured.loo_mean_r > 0.9, 'Leave-one-out correlation should be near 1.');
assert(measured.phase_fraction > 0.9, 'Every window should share the same phase.');
null_size = abs(measured.null_mean_coefficient(:, 1));
assert(mean(null_size < abs(measured.mean_coefficient(1))) > 0.9, ...
    'Onset shifts did not break the 2.5 Hz alignment.');
end

function test_sine_spectrum_and_coherence()
rng(5, 'twister');
rate = 50;
time = (1:rate * 40).' / rate;
envelope = sin(2 * pi * 2.5 * time);
lfp = sin(2 * pi * 2.5 * time + 0.4);
grids = test_grids();
welch = struct('window_seconds', 2, 'overlap', 0.5, 'resolution_hz', 0.1);
spectrum = theta_envelope_session_spectrum(envelope, rate, grids, welch);
[~, idx] = max(spectrum.band_power_fine);
assert(abs(grids.fine.centers(idx) - 2.5) <= 0.2, ...
    'The sine spectrum did not peak at 2.5 Hz.');
coherence = theta_envelope_session_coherence(envelope, lfp, rate, grids, welch, false);
shuffled = theta_envelope_session_coherence(envelope, lfp, rate, grids, welch, true);
band = @(grid, values) values(abs(grid.centers - 2.5) < 0.15);
assert(band(grids.fine, coherence.band_coherence_fine) > 0.8, ...
    'Two 2.5 Hz sines should be coherent.');
assert(band(grids.fine, shuffled.band_coherence_fine) < 0.5, ...
    'Segment shifts did not reduce coherence.');
end

function test_modulated_theta()
rng(6, 'twister');
rate = 50;
time = (1:rate * 40).' / rate;
options = theta_envelope_delta_options(struct('envelope_rate', rate));
coefficients = theta_envelope_filter_coefficients(options);
carrier = filtfilt(coefficients.', 1, randn(size(time)));
carrier = carrier / std(carrier);
modulator = sin(2 * pi * 2.5 * time);
lfp = (1 + 0.9 * modulator) .* carrier + 0.8 * modulator;
envelope = theta_envelope_signal(lfp, coefficients);
grids = test_grids();
welch = struct('window_seconds', 2, 'overlap', 0.5, 'resolution_hz', 0.1);
spectrum = theta_envelope_session_spectrum(envelope, rate, grids, welch);
[~, idx] = max(spectrum.band_power_fine);
assert(abs(grids.fine.centers(idx) - 2.5) <= 0.3, ...
    'The theta envelope did not oscillate at the modulation frequency.');
power_null = theta_envelope_spectrum_null(lfp, coefficients, spectrum, rate, grids, welch, 12);
coherence = theta_envelope_session_coherence(envelope, lfp, rate, grids, welch, false);
coherence_null = theta_envelope_coherence_null(envelope, lfp, coherence, rate, grids, welch, 12);
power_z = power_null.test.z(grids.test.confirmatory);
coherence_z = coherence_null.test.z(grids.test.confirmatory);
assert(max(power_z) > 2, 'Envelope power at 2-3 Hz did not beat the phase-randomized null.');
assert(max(coherence_z) > 2, 'Envelope-LFP coherence did not beat the segment-shift null.');
end

function test_decimate_keeps_theta()
rate_in = 2500;
rate_out = 50;
time = (1:rate_in * 8).' / rate_in;
kept = theta_envelope_decimate(sin(2 * pi * 9 * time), rate_in, rate_out);
removed = theta_envelope_decimate(sin(2 * pi * 40 * time), rate_in, rate_out);
center = 80:(numel(kept) - 80);
assert(std(kept(center)) > 5 * std(removed(center)), ...
    'Decimation did not preserve 9 Hz relative to 40 Hz.');
index = theta_envelope_time_index(8, rate_out, rate_in);
assert(abs(index - round((round(8 * rate_in) - 1) / (rate_in / rate_out)) - 1) < 1e-9, ...
    'Decimated time index does not follow the LFP clock.');
end

function recording = fake_recording(animal, z_values, frequencies, mask, coefficient, n_valid)
recording.animal_area_table_name = animal;
recording.frequencies_test = frequencies;
recording.frequencies_fine = frequencies;
recording.confirmatory_mask = mask;
recording.session.power_z_test = z_values;
recording.session.coherence_z_test = z_values;
recording.session.power_z_fine = z_values;
recording.session.coherence_z_fine = z_values;
recording.session.band_power_fine = z_values;
recording.session.band_coherence_fine = z_values;
recording.onset.epoch_time = -1:0.5:1;
recording.onset.epochs = coefficient * sin(2 * pi * 2.5 * recording.onset.epoch_time);
recording.onset.epochs = repmat(recording.onset.epochs, n_valid, 1);
recording.onset.clean_pre = true(n_valid, 1);
pre = struct();
pre.n_valid = n_valid;
pre.mean_coefficient = coefficient * ones(1, numel(frequencies));
pre.ppc = ones(1, numel(frequencies));
pre.loo_mean_r = coefficient;
pre.phase_fraction = double(coefficient > 0);
pre.null_mean_coefficient = zeros(5, numel(frequencies));
pre.null_ppc = zeros(5, numel(frequencies));
pre.null_loo_mean_r = zeros(5, 1);
pre.null_phase_fraction = zeros(5, 1);
recording.onset.pre = pre;
end

function grids = test_grids()
grids.test = theta_envelope_frequency_grid([0.1 6], 0.5, 0.5, [2 3]);
grids.fine = theta_envelope_frequency_grid([0.1 6], 0.2, 0.5, [2 3]);
end

function ready = signal_ready()
ready = exist('filtfilt', 'file') == 2 && exist('hilbert', 'file') == 2;
end
