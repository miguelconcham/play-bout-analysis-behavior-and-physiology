function test_session_acg_shuffle()
%TEST_SESSION_ACG_SHUFFLE  Checks that do not need the experiment files.
%
% Pair counts, windowing, power tails, the delta-peak call, FOOOF delta,
% one resumed session, and the grand-average null. Run this file by name
% after add_repo_paths, or run it directly.

this_dir = fileparts(mfilename('fullpath'));
complementary_dir = fileparts(this_dir);
repo_root = fileparts(complementary_dir);
addpath(this_dir);
addpath(complementary_dir);
addpath(genpath(fullfile(repo_root, 'Custom functions')));
addpath(genpath(fullfile(repo_root, 'Figure 2', 'Supporting codes')));

assert(strcmp(canonical_kept_area('SupCol', {'SC', 'LPAG'}), 'SC'), ...
    'SupCol was not read as superior colliculus.');
assert(strcmp(canonical_kept_area('VLPAG', {'SC', 'LPAG', 'VLPAG'}), 'VLPAG'), ...
    'VLPAG was folded into LPAG.');
assert(isempty(canonical_kept_area('DMPAG', {'SC', 'LPAG'})), ...
    'DMPAG was kept.');

session_column = {'A'; 'B'; {'A'}};
assert(isequal(session_name_matches(session_column, 'A'), [true; false; true]), ...
    'Nested session names were not matched.');

edges = -1:0.01:1;
spikes = [0; 0.25; 0.5; 1; 1];
counts = pair_acg_counts(spikes, edges);
reference = histcounts(autocorrelogram(spikes, [edges(1), edges(end)]), edges);
assert(isequal(counts, reference), ...
    'Pair counts did not match autocorrelogram at the outer lag.');
rate = interval_spike_acg(spikes, 0, 2, edges, true);
reference_rate = run_window_acgs(spikes, 0, 2, edges);
assert(max(abs(rate - reference_rate)) < 1e-9, ...
    'The session autocorrelogram did not match the firing-run helper.');

one_window = pool_interval_acgs([0.1; 0.4], 0, 1, edges);
two_windows = pool_interval_acgs([0.1; 0.4; 10.1; 10.4], [0; 10], [1; 11], edges);
assert(max(abs(one_window - two_windows)) < 1e-9, ...
    'Two copies of a window did not pool to the same autocorrelogram.');

[starts, ends] = session_time_windows(0, 25, 10);
assert(isequal(starts, [0; 10]) && isequal(ends, [10; 20]), ...
    'A 5 s leftover was kept.');
window_counts = window_spike_counts([0; 9.9; 10; 19.999; 20], 0, 10, 2);
assert(isequal(window_counts, [2; 2]), ...
    'A spike on a window boundary was counted twice or dropped.');

[low_index, high_index] = power_tail_index((1:8)', 25);
assert(isequal(find(low_index), [1; 2]) && isequal(find(high_index), [7; 8]), ...
    'The 25% tails were not the two lowest and two highest scores.');
[empty_low, empty_high] = power_tail_index([1; 2; NaN], 25);
assert(~any(empty_low) && ~any(empty_high), ...
    'Three scores were split into tails.');

[scaled, mass] = sum_to_one_acg([1 2 3 4], [-0.4 -0.2 0.2 0.4], 0);
assert(abs(sum(scaled) - 1) < 1e-12 && abs(mass - 10) < 1e-12, ...
    'The autocorrelogram was not scaled by its own mass.');

centers = -1:0.01:1;
t = abs(centers);
peak = 0.15 * exp(-t / 0.35) + 0.55 * exp(-((t - 0.30) / 0.035) .^ 2);
decay = 0.15 * exp(-t / 0.35);
null_decay = repmat(decay, 20, 1);
scored = score_acg_oscillation(peak, null_decay, centers, 1, 0.05);
assert(scored.oscillates && scored.first_p <= 0.05 && scored.first_lag > 0.25 && scored.first_lag < 0.36, ...
    'A 300 ms peak did not clear the null.');
assert(max(scored.pctl) > 90, ...
    'The delta peak was not in the top of the null.');
flat_score = score_acg_oscillation(decay, null_decay, centers, 1, 0.05);
assert(~flat_score.oscillates, ...
    'A smooth decay was called a delta oscillation.');

options = estimate_options();
options.fooof_max_peaks = 6;
options.fooof_fit_range = [0.5 40];
frequency = (0.5:0.25:40)';
aperiodic = 10 .^ (1.5 - 1.2 * log10(frequency));
with_peak = aperiodic .* 10 .^ (0.8 * exp(-((frequency - 3) / 0.5) .^ 2));
delta_background = fooof_band_power(frequency, aperiodic, [1 6], options);
delta_bump = fooof_band_power(frequency, with_peak, [1 6], options);
assert(isfinite(delta_background) && isfinite(delta_bump) && delta_bump > delta_background + 0.05, ...
    'Removing the aperiodic fit did not leave the 3 Hz bump.');

sample_rate = 250;
sample_index = (1:20 * sample_rate)';
sample_time = sample_index / sample_rate;
rng(1, 'twister');
signal = 0.05 * randn(size(sample_time));
signal(sample_time >= 10) = 4 * sin(2 * pi * 3 * sample_time(sample_time >= 10));
window_options = options;
window_options.delta_band = [1 6];
window_options.filter_order = 64;
window_options.lfp_sample_rate = sample_rate;
measured = channel_window_delta(signal, sample_rate, [0; 10], [10; 20], window_options);
assert(measured.hilbert_amp(2) > 5 * measured.hilbert_amp(1), ...
    'The Hilbert amplitude did not follow the larger delta sine.');
assert(isfinite(measured.fooof_delta(1)) && isfinite(measured.fooof_delta(2)) ...
    && measured.fooof_delta(2) > measured.fooof_delta(1), ...
    'FOOOF delta did not rise in the window that contains the sine.');

options.n_shuffles = 25;
options.window_s = 10;
options.min_spikes = 20;
options.smooth_half = 0.01;
options.resume = false;
options.areas_to_keep = {'SC'};
settings = session_acg_settings(options);
options.bin_edges = settings.bin_edges;
options.smooth_n = settings.smooth_n;
root = tempname;
session_dir = fullfile(root, 'SessionA');
mkdir(root);
cleanup_root = onCleanup(@() rmdir(root, 's')); %#ok<NASGU>
save(fullfile(root, 'parameters.mat'), 'settings');
options.log_file = fullfile(root, 'estimate_log.txt');
channel_power.channel = 1;
channel_power.fooof_delta = [1 2 3 10];
channel_power.hilbert_amp = [1 2 3 10];
options.channel_power = channel_power;
spikes = (0.15:0.3:39.9)';
neurons = struct('cluster_id', 7, 'channel', 1, 'area', 'SC', 'lock_group', 'trough');
summary = estimate_one_session_acg('SessionA', spikes, ones(size(spikes)), ...
    neurons, [], [0 40], options, session_dir);
assert(summary.n_estimated == 1 && summary.n_skipped == 0, ...
    'The first pass did not estimate the one neuron.');
assert(summary.rows(1).session_oscillates && summary.rows(1).session_p <= 0.05, ...
    'The regular train was not called a delta oscillation.');
saved_neuron = load(fullfile(session_dir, 'neuron_7.mat'), 'neuron');
assert(size(saved_neuron.neuron.curves.session.null, 1) == options.n_shuffles, ...
    'The saved null did not keep every shuffle.');
assert(saved_neuron.neuron.scores.fooof_high.oscillates, ...
    'The high FOOOF-delta tail of a regular train did not oscillate.');
assert(saved_neuron.neuron.tail_used.fooof_low && saved_neuron.neuron.tail_used.fooof_high, ...
    'The power tails were dropped although each had enough spikes.');

options.resume = true;
resumed = estimate_one_session_acg('SessionA', spikes, ones(size(spikes)), ...
    neurons, [], [0 40], options, session_dir);
assert(resumed.n_estimated == 0 && resumed.n_skipped == 1, ...
    'The resumed run estimated the neuron again.');
grand = session_acg_grand_average(root);
assert(grand.groups.trough.session.n_neurons == 1, ...
    'The trough grand mean did not include the saved neuron.');
assert(grand.groups.delta.session.n_neurons == 1, ...
    'Delta-locked did not include the trough neuron.');
assert(grand.counts.nonentrained.session.n == 0, ...
    'A missing group was given neurons.');
assert(grand.groups.trough.session.result.oscillates, ...
    'The one-neuron grand mean did not keep that neuron''s delta peak.');

alignment_root = fullfile(root, 'alignment');
align_dir = fullfile(alignment_root, 'SessionB');
mkdir(alignment_root);
mkdir(align_dir);
align_options = options;
align_options.n_shuffles = 2;
align_options.bin_w = 0.5;
align_options.acg_lag = 1;
align_options.bin_edges = -1:0.5:1;
align_settings = session_acg_settings(align_options);
settings = align_settings; %#ok<NASGU>
save(fullfile(alignment_root, 'parameters.mat'), 'settings');
save(fullfile(align_dir, 'parameters.mat'), 'settings');
save_alignment_neuron(align_dir, 'SessionB', 1, 'trough', ...
    [1 0 0 0; 0 1 0 0], true);
save_alignment_neuron(align_dir, 'SessionB', 2, 'peak', ...
    [1 0 0 0; 0 0 1 0], false);
aligned = session_acg_grand_average(alignment_root);
delta_block = aligned.groups.delta.session;
assert(delta_block.n_neurons == 2, ...
    'Delta-locked did not pool trough and peak.');
assert(abs(delta_block.sum_null(1, 1) - 2) < 1e-9, ...
    'Shuffle 1 was not averaged with shuffle 1.');
assert(abs(delta_block.sum_null(2, 2) - 1) < 1e-9 && abs(delta_block.sum_null(2, 3) - 1) < 1e-9, ...
    'Shuffle 2 mixed the two neurons.');
assert(aligned.counts.delta.session.n == 2 && aligned.counts.delta.session.n_oscillates == 1, ...
    'The oscillation count mixed the two neurons.');
assert(aligned.counts.trough.session.n_oscillates == 1 && aligned.counts.peak.session.n_oscillates == 0, ...
    'Group counts did not stay with their own neurons.');

outside = crosscorr_acg_counts([0.2; 9.8], 0, 10, -1:0.1:1, true);
assert(all(isnan(outside)), ...
    'Spikes outside [1, 9] still formed a reference autocorrelogram.');
interior = crosscorr_acg_counts([0.5; 1.2], 0, 10, -1:0.1:1, true);
assert(isequal(interior, histcounts(-0.7, -1:0.1:1)), ...
    'A spike before 1 s was used as a reference, or the partner in the window was missed.');
both_sides = crosscorr_acg_counts([2; 2.4], 0, 10, -1:0.1:1, true);
assert(isequal(both_sides, histcounts([-0.4 0.4], -1:0.1:1)), ...
    'An interior pair did not appear on both sides of the crosscorrelogram.');
one_window = crosscorr_acg_counts([2; 2.4], 0, 10, -1:0.1:1, false);
pooled_windows = pool_crosscorr_acgs([2; 2.4; 12; 12.4], [0; 10], [10; 20], -1:0.1:1);
assert(isequal(pooled_windows, 2 * one_window), ...
    'Two interior windows did not add their crosscorrelogram counts.');

missing = load_pag_recording(tempname);
assert(~missing.ok, 'A missing session folder was treated as a recording.');

changed = options;
changed.min_spikes = options.min_spikes + 1;
assert(~session_acg_parameters_match(options, changed), ...
    'A changed minimum spike count still matched the saved run.');
changed = options;
changed.note = 'not part of the estimate';
assert(session_acg_parameters_match(options, changed), ...
    'An extra note changed the estimate identity.');
changed = options;
changed.acg_method = 'crosscorr';
assert(~session_acg_parameters_match(options, changed), ...
    'The crosscorrelogram method matched the edge-corrected run.');

disp('test_session_acg_shuffle passed')
end


function options = estimate_options()
options = struct( ...
    'version', 1, ...
    'n_shuffles', 25, ...
    'window_s', 10, ...
    'min_spikes', 20, ...
    'tail_percent', 25, ...
    'alpha', 0.01, ...
    'non_entrained_p', 0.1, ...
    'oscillation_alpha', 0.05, ...
    'acg_lag', 1, ...
    'bin_w', 0.01, ...
    'delta_band', [1 6], ...
    'smooth_half', 0.01, ...
    'acg_refractory_s', 0.005, ...
    'areas_to_keep', {{'SC'}}, ...
    'fooof_max_peaks', 6, ...
    'fooof_fit_range', [0.5 40], ...
    'filter_order', 64, ...
    'lfp_sample_rate', 250, ...
    'acg_method', 'edge', ...
    'resume', false);
end


function save_alignment_neuron(session_dir, session_name, cluster_id, lock_group, null_rows, oscillates)
names = session_acg_names('curves');
n_bins = size(null_rows, 2);
neuron.session = session_name;
neuron.cluster_id = cluster_id;
neuron.channel = 1;
neuron.area = 'SC';
neuron.lock_group = lock_group;
neuron.curves = struct();
neuron.scores = struct();
for i = 1:numel(names)
    neuron.curves.(names{i}).observed = nan(1, n_bins, 'single');
    neuron.curves.(names{i}).null = nan(size(null_rows, 1), n_bins, 'single');
    neuron.scores.(names{i}).oscillates = false;
end
neuron.curves.session.observed = single(ones(1, n_bins) / n_bins);
neuron.curves.session.null = single(null_rows);
neuron.scores.session.oscillates = oscillates;
save(fullfile(session_dir, sprintf('neuron_%d.mat', cluster_id)), 'neuron');
end
