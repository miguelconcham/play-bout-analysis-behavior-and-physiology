function test_sync_event_statistics()
%TEST_SYNC_EVENT_STATISTICS  Checks that do not need the experiment file system.
%
% Coincidence mask, pair caps, waveform distance, phase-matched LFP power,
% and continuous-file choice. Run this file by name after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(this_dir), 'Custom functions'));
end

mask = sync_spike_mask([0; 1; 1.004; 5], 1.003, 1, [-1.5 1], 0.005);
assert(isequal(mask, [false; true; true; false]), ...
    'Play-bout coincidence mask did not mark the two spikes near 1 s.');

options = sync_analysis_options(struct('max_pairs', 0, 'min_spikes', 10));
assert(options.max_pairs == 0, 'A pair budget of 0 was replaced by the default.');
assert(options.min_spikes == 10, 'min_spikes override was lost.');

[pairs, kept] = select_synchronous_pairs(example_neurons(), 10, ...
    sync_analysis_options(struct('max_pairs_per_neuron', 1, 'max_pairs', 2)));
assert(numel(pairs) == 1, 'Neuron cap kept more than one pair.');
assert(pairs(1).cluster_id_1 == 11 && pairs(1).cluster_id_2 == 22, ...
    'The pair with more spikes was not taken first.');
assert(pairs(1).n_sync_1 == 20 && numel(pairs(1).sync_index_1) == 20, ...
    'Synchronous spikes were capped or dropped.');
assert(numel(kept) == 2, 'Kept neurons should be the two members of the pair.');
neuron_b = kept([kept.cluster_id] == 22);
assert(numel(neuron_b.interval_index) == 50, ...
    'The bout-window pool should be the synchronous spikes plus the in-window controls.');
assert(all(neuron_b.spike_times(neuron_b.interval_index) >= 8.5 ...
    & neuron_b.spike_times(neuron_b.interval_index) <= 11), ...
    'Spikes outside the bout window entered the comparison pool.');
assert(all(ismember(pairs(1).sync_index_2, neuron_b.interval_index)), ...
    'A synchronous spike fell outside the bout-window pool.');
neuron_a = kept([kept.cluster_id] == 11);
assert(numel(setdiff(neuron_a.interval_index, pairs(1).sync_index_1)) >= pairs(1).n_sync_1, ...
    'Spikes synchronous only with the other partner were removed from this pair''s pool.');

[pairs, ~, coverage] = select_synchronous_pairs(example_neurons(), 10, ...
    sync_analysis_options(struct('max_pairs_per_neuron', 2, 'max_pairs', 2)));
assert(numel(pairs) == 2, 'Two uses of the best-sampled neuron were expected.');
assert(pairs(2).cluster_id_1 == 11 && pairs(2).cluster_id_2 == 33, ...
    'The second pair was not the neuron still missing from the sample.');
assert(coverage.target_met && coverage.n_covered == 3, ...
    'Three areas with one neuron each were not all observed.');

same_area = example_neurons();
for neuron_index = 1:numel(same_area)
    same_area(neuron_index).area = 'LPAG';
end
[pairs, ~, coverage] = select_synchronous_pairs(same_area, 10, ...
    sync_analysis_options(struct('min_neuron_fraction', 0.5)));
assert(numel(pairs) == 1 && coverage.n_covered == 2 && coverage.target_met, ...
    'Pairs were added after the area had already reached the target.');

[pairs, kept, coverage] = select_synchronous_pairs(hub_neurons(), 10, ...
    sync_analysis_options(struct('min_spikes', 50, 'min_sync', 10)));
assert(numel(kept) == 5 && coverage.n_covered == 5 && coverage.target_met, ...
    'Reusing the hub did not bring every partner neuron into the sample.');
assert(numel(pairs) == 4, ...
    'The hub sample kept combinations that did not add a neuron.');

[pairs, ~] = select_synchronous_pairs(example_neurons(), 10, ...
    sync_analysis_options(struct('min_sync', 100)));
assert(isempty(pairs), 'A pair below min_sync was kept.');

short_options = sync_analysis_options(struct( ...
    'min_spikes', 30, 'min_sync', 10, 'max_pairs', 5));
[pairs, ~] = select_synchronous_pairs(too_few_other_spikes(), 10, short_options);
assert(isempty(pairs), 'A pair without enough non-synchronous spikes was kept.');

high_power = compare_sync_scalar([10 * ones(30, 1); ones(200, 1)], 1:30, 40);
assert(high_power.ok && high_power.percentile >= 0.99, ...
    'High synchronous values were not above the null.');
assert(all(abs(high_power.null_means - 1) < 1e-12), ...
    'Null draws included the synchronous spikes.');
assert(abs(high_power.mean_nonsync - 1) < 1e-12, ...
    'The non-synchronous mean included the synchronous spikes.');
low_power = compare_sync_scalar([zeros(30, 1); ones(200, 1)], 1:30, 40);
assert(low_power.ok && low_power.percentile <= 0.05, ...
    'Low synchronous values were not below the null.');
blocked = compare_sync_scalar(ones(20, 1), 1:20, 10);
assert(~blocked.ok, 'A neuron with no non-synchronous spikes should be skipped.');

phase_options = sync_analysis_options();
assert(isequal(phase_options.delta_freq, [1 5]), 'Delta band is not 1-5 Hz.');
assert(phase_options.delta_filter_order == 2000, 'Delta FIR order is not 2000.');
assert(phase_options.n_phase_bins == 36, 'Delta phase is not split into 36 bins.');
assert(phase_options.phase_guard_seconds == 0.01, 'The phase guard is not 10 ms.');
matched_phase = 0.05 * ones(80, 1);
other_phase = (pi - 0.05) * ones(80, 1);
pool_phase = [matched_phase; other_phase];
pool_power = [ones(80, 1); 50 * ones(80, 1)];
high_phase = compare_phase_matched_lfp(10 * ones(20, 1), 0.05 * ones(20, 1), ...
    pool_power, pool_phase, 40, 36);
assert(high_phase.ok, high_phase.reason);
assert(all(abs(high_phase.null_means - 1) < 1e-12), ...
    'The phase-matched null used LFP from a different delta phase.');
assert(high_phase.percentile == 1, ...
    'Synchronous power above the matched phase was not percentile 1.');
assert(abs(high_phase.mean_phase_matched - 1) < 1e-12, ...
    'The phase-matched mean included the other phase bin.');
assert(high_phase.min_candidates == 80 && high_phase.n_pool == 160, ...
    'The phase pool count is not the full window sample.');
low_phase = compare_phase_matched_lfp(zeros(20, 1), 0.05 * ones(20, 1), ...
    pool_power, pool_phase, 40, 36);
assert(low_phase.ok && low_phase.percentile == 0, ...
    'Synchronous power below the matched phase was not percentile 0.');
empty_phase = compare_phase_matched_lfp(10, pi / 2, pool_power, pool_phase, 10, 36);
assert(~empty_phase.ok, 'A synchronous phase with no pool samples was accepted.');
trough = compare_phase_matched_lfp(10, pi, 1, -pi, 5, 36);
assert(trough.ok && abs(trough.mean_phase_matched - 1) < 1e-12, ...
    '+pi and -pi were treated as different delta phases.');
same_peak = compare_phase_matched_lfp(4, 0, 2, 0.05, 5, 36);
assert(same_peak.ok && abs(same_peak.mean_phase_matched - 2) < 1e-12, ...
    'Phase 0 and a sample 0.05 rad later fell in different bins.');
[pool_index, sync_sample] = lfp_phase_pool(100, 10, 5, [-1 1], [5; 0; 20], 0.1);
assert(sync_sample(1) == 50 && sync_sample(2) == 0 && sync_sample(3) == 0, ...
    'Synchronous spikes were not mapped to the nearest in-trace LFP sample.');
assert(isequal(pool_index(:), [40:48, 52:60]'), ...
    'The phase pool is not the bout window with 10 ms around the event removed.');

time_samples = -20:10;
[~, peak_at] = ismember(0, time_samples);
matched = zeros(80, 1, numel(time_samples));
matched(:, 1, peak_at) = 10;
matched(1:20, 1, :) = matched(1:20, 1, :) + 50;
matched_result = compare_sync_waveforms(matched, time_samples, 1:20, 1, 1, 25);
assert(matched_result.ok, matched_result.reason);
assert(matched_result.distance_sync < 1e-9, ...
    'Baseline subtraction left a DC offset in the distance.');
assert(matched_result.percentile == 0, ...
    'Identical shapes were treated as more dissimilar than the null.');
assert(all(matched_result.percentile_time(:) == 0), ...
    'Identical samples were ranked above the null means.');
assert(matched_result.extreme_abs_percentile == 0, ...
    'An identical neuron channel was marked as an extreme departure.');

artifact = zeros(80, 3, numel(time_samples));
artifact(:, 2, peak_at) = 20;
artifact(1:20, 1, peak_at) = 50;
artifact_result = compare_sync_waveforms(artifact, time_samples, 1:20, 5, [4 5 6], 25);
assert(artifact_result.ok, artifact_result.reason);
assert(artifact_result.percentile == 1, ...
    'A side-channel difference was not farther than the null.');
assert(artifact_result.percentile_channel == 0, ...
    'The neuron channel was marked dissimilar when its shape matched.');
assert(isequal(size(artifact_result.nonsync_mean), [3, numel(time_samples)]), ...
    'Non-synchronous mean is not channels by time.');
assert(artifact_result.percentile_time(1, peak_at) == 1, ...
    'The side-channel peak was not above every null mean.');
assert(artifact_result.percentile_time(2, peak_at) == 0, ...
    'The matching neuron-channel sample was ranked above the null.');
assert(artifact_result.percentile_abs_time(1, peak_at) == 1, ...
    'The absolute side-channel departure was not above the null.');
assert(artifact_result.n_pool == 80 && artifact_result.n_nonsync == 60, ...
    'The waveform pool did not keep the non-synchronous spikes.');
assert(artifact_result.null_p975(1, peak_at) < 1, ...
    'The null envelope on the unchanged side channel was not near zero.');

[chosen, ~] = choose_continuous_file( ...
    {spike_file(), zero_file()}, 50, 3, 2, check_choice());
assert(strcmp(chosen, spike_file()), ...
    'The flat continuous file was chosen over the file with a spike.');
flat_failed = false;
try
    choose_continuous_file(zero_file(), 50, 3, 2, check_choice());
catch
    flat_failed = true;
end
assert(flat_failed, 'A flat continuous file was accepted.');

waveform_session = GENERATE_SYNC_SPIKE_WAVEFORMS( ...
    fullfile(tempdir, 'missing-sync-session'), struct('max_pairs', 0));
assert(strcmp(waveform_session.status, 'pair_budget_exhausted'), ...
    'Waveform GENERATE did not stop when the pair budget was already filled.');
power_session = GENERATE_SYNC_SPIKE_LFP_POWER( ...
    fullfile(tempdir, 'missing-sync-session'), struct('max_pairs', 0));
assert(strcmp(power_session.status, 'pair_budget_exhausted'), ...
    'LFP-power GENERATE did not stop when the pair budget was already filled.');

delete_if_present(spike_file());
delete_if_present(zero_file());
older = struct('session_id', 'old', 'rows', struct('percentile', 0.2));
newer = struct('session_id', 'new', 'pairs', struct('pair_index', 1), 'rows', ...
    struct('percentile', 0.9));
combined = append_struct_array(older, newer);
assert(numel(combined) == 2, 'Sessions with different fields were not combined.');
assert(isempty(combined(1).pairs) && combined(2).pairs.pair_index == 1, ...
    'The older session did not keep an empty field for the new one.');
rows = append_struct_array(struct('percentile', 0.2), struct('percentile', 0.9, 'nonsync_mean', 1));
assert(isempty(rows(1).nonsync_mean) && rows(2).nonsync_mean == 1, ...
    'Neuron rows with different fields were not combined.');

[centers, counts] = sync_event_acg([0; 0.002; 0.004], [-0.05 0.05], 0.001);
assert(sum(counts(centers > 0 & centers < 0.003)) == 2, ...
    'The autocorrelogram missed a pair of synchronous spikes 2 ms apart.');
assert(sum(counts(centers > 0.003 & centers < 0.005)) == 1, ...
    'The autocorrelogram missed the pair 4 ms apart.');
assert(sum(counts(abs(centers) <= 0.0005)) == 0, ...
    'The zero lag was left in the autocorrelogram.');
summary = normalized_acg_stack([0 2 0; 0 0 4; 0 0 0]);
assert(summary.n == 2, 'A neuron with no spike pairs stayed in the normalized stack.');
assert(abs(sum(summary.normalized(1, :)) - 1) < 1e-12, ...
    'A normalized autocorrelogram does not sum to 1.');
assert(all(abs(summary.mean_acg - [0 0.5 0.5]) < 1e-12), ...
    'The mean normalized autocorrelogram is wrong.');
assert(summary.normalized(1, 2) == 1 && summary.normalized(2, 3) == 1, ...
    'Rows were not ordered by the lag of the peak.');
half_width = 1.96 * std([1 0], 0) / sqrt(2);
assert(abs(summary.ci_high(2) - (0.5 + half_width)) < 1e-12, ...
    'The autocorrelogram confidence interval is not 1.96 standard errors.');

check_location_lookup();
check_spatial_range();
check_spatial_stack();
fprintf('Synchronous-event statistic checks passed.\n');
end

function check_spatial_range()
chan_map = [10; 11; 12];
xcoords = [0; 32; 0];
ycoords = [0; 20; -20];
waveform = [0 10 0; 0 4 0; 0 2 -1];
profile = sync_waveform_spatial_range([12 10 11], waveform([3 1 2], :), 10, ...
    chan_map, xcoords, ycoords);
assert(profile.neuron_found, 'The neuron channel was missing from the probe map.');
assert(issorted(profile.distance_um), 'Spatial samples were not ordered along the probe.');
assert(profile.center_channel == 10, ...
    'The largest-range channel was not selected as the origin.');
assert(profile.distance_um(profile.channel == 10) == 0, ...
    'The largest-range channel is not at 0 um.');
assert(abs(profile.distance_um(profile.channel == 11) - sqrt(32 ^ 2 + 20 ^ 2)) < 1e-9, ...
    'Euclidean distance missed the lateral offset.');
assert(profile.distance_um(profile.channel == 12) == -20, ...
    'A site above the neuron was not given a negative distance.');
assert(profile.range(profile.channel == 10) == 10, 'Peak-to-peak range is wrong.');
assert(profile.range(profile.channel == 12) == 3, 'Peak-to-peak range dropped the negative sample.');
assert(abs(pair_channel_distance(10, 11, chan_map, xcoords, ycoords) - sqrt(32 ^ 2 + 20 ^ 2)) < 1e-9, ...
    'Pair separation missed the Euclidean distance.');
assert(pair_channel_distance(10, 12, chan_map, xcoords, ycoords) == 20, ...
    'Pair separation along depth was not 20 um.');
assert(pair_channel_distance(10, 10, chan_map, xcoords, ycoords) == 0, ...
    'A neuron paired with its own channel was not at 0 um.');
assert(isnan(pair_channel_distance(10, 99, chan_map, xcoords, ycoords)), ...
    'A missing partner channel produced a distance.');
shifted = sync_waveform_spatial_range([10 11 12], [0 4 0; 0 10 0; 0 2 -1], 10, ...
    chan_map, xcoords, ycoords);
assert(shifted.center_channel == 11, ...
    'The range plot stayed on the assigned channel when another channel was larger.');
assert(shifted.distance_um(shifted.channel == 11) == 0, ...
    'The largest-range channel is not at 0 um after recentering.');
assert(abs(shifted.distance_um(shifted.channel == 10) + sqrt(32 ^ 2 + 20 ^ 2)) < 1e-9, ...
    'Distance from the peak channel missed the assigned channel.');
assert(abs(shifted.distance_um(shifted.channel == 12) + sqrt(32 ^ 2 + 40 ^ 2)) < 1e-9, ...
    'Distance from the peak channel missed the far site.');
lateral = sync_waveform_spatial_range([10 13], [0 1; 0 1], 10, [10; 13], [0; 16], [100; 100]);
assert(abs(lateral.distance_um(lateral.channel == 13) - 16) < 1e-9, ...
    'A site at the same depth did not take its sign from x.');
dropped = sync_waveform_spatial_range([10 99], [0 5; 0 9], 10, chan_map, xcoords, ycoords);
assert(isequal(dropped.channel, 10), 'A channel missing from the map was plotted.');

session_dir = fullfile(tempdir, 'sync_geom_npx1');
if exist(session_dir, 'dir') ~= 7
    mkdir(session_dir);
end
chanMap = [(1:191).'; (193:384).']; %#ok<NASGU>
xcoords = repmat([16; 48], ceil(numel(chanMap) / 2), 1); %#ok<NASGU>
xcoords = xcoords(1:numel(chanMap)); %#ok<NASGU>
ycoords = (1:numel(chanMap)).' * 20; %#ok<NASGU>
save(fullfile(session_dir, 'chann_map_PAG.mat'), 'chanMap', 'xcoords', 'ycoords');
[restored, restored_x, restored_y] = read_npx_channel_geometry(session_dir);
assert(restored(192) == 192 && restored_x(192) == 43 && restored_y(192) == 1900, ...
    'NPX1 channel 192 was not inserted at its probe position.');
assert(restored(193) == 193, 'Channels after 192 shifted onto the wrong site.');
delete(fullfile(session_dir, 'chann_map_PAG.mat'));

npx2_dir = fullfile(tempdir, 'sync_geom_npx2');
if exist(npx2_dir, 'dir') ~= 7
    mkdir(npx2_dir);
end
xcoords(1:2) = [8; 40]; %#ok<NASGU>
save(fullfile(npx2_dir, 'chann_map_PAG.mat'), 'chanMap', 'xcoords', 'ycoords');
[npx2_map, ~, ~] = read_npx_channel_geometry(npx2_dir);
assert(~any(npx2_map == 192) && numel(npx2_map) == numel(chanMap), ...
    'An NPX2 map was given the NPX1 channel-192 site.');
delete(fullfile(npx2_dir, 'chann_map_PAG.mat'));
end

function check_location_lookup()
if exist('readtable', 'file') ~= 2 && exist('readtable', 'builtin') ~= 5
    return
end
location_file = fullfile(tempdir, 'npx_raw_data_folders.csv');
write_text(location_file, sprintf('%s\n%s\n%s\n', ...
    'SessionName,Location', ...
    'B2S2 1110 Single,server/B2S2/continuous.dat', ...
    'B4D4 0826 Dual,server/B4D4/continuous.dat'));
paths = read_npx_continuous_locations(location_file, 'B2S2 1110 Single');
delete_if_present(location_file);
assert(isequal(paths, {'server/B2S2/continuous.dat'}), ...
    'The SessionName row did not return its Location.');
end

function write_text(path, text)
file_id = fopen(path, 'w');
assert(file_id >= 0, 'Could not write %s.', path);
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
fwrite(file_id, text);
end

function neurons = hub_neurons()
blocks = zeros(20, 4);
for block_index = 1:4
    blocks(:, block_index) = 8.70 + (block_index - 1) * 0.30 + (0:19)' * 0.01;
end
controls = zeros(30, 5);
for control_index = 1:5
    controls(:, control_index) = 10.00 + (control_index - 1) * 0.12 + (0:29)' * 0.002;
end
areas = {'DLPAG', 'VLPAG', 'DR', 'SupCol'};
neurons = empty_neurons(5);
neurons(1) = make_neuron(1, 1, 'LPAG', [blocks(:); controls(:, 1); 100 + (1:40)']);
for partner_index = 1:4
    neurons(partner_index + 1) = make_neuron(10 + partner_index, 10 + partner_index, ...
        areas{partner_index}, [blocks(:, partner_index); controls(:, partner_index + 1); ...
        1000 * partner_index + (1:40)']);
end
end

function neurons = example_neurons()
ab = 9.00 + (0:19)' * 0.02;
ac = 9.50 + (0:19)' * 0.02;
control_a = 10.40 + (0:29)' * 0.01;
control_b = 10.00 + (0:29)' * 0.01;
control_c = 8.60 + (0:29)' * 0.01;
neurons = empty_neurons(3);
neurons(1) = make_neuron(11, 10, 'LPAG', [ab; ac; control_a; 100 + (1:4960)']);
neurons(2) = make_neuron(22, 20, 'DR', [ab; control_b; 20000 + (1:3980)']);
neurons(3) = make_neuron(33, 30, 'SupCol', [ac; control_c; 40000 + (1:2980)']);
end

function neurons = too_few_other_spikes()
sync = 9 + (0:24)' * 0.01;
neurons = empty_neurons(2);
neurons(1) = make_neuron(1, 1, 'LPAG', [sync; 100 + (1:15)']);
neurons(2) = make_neuron(2, 2, 'DR', [sync; 200 + (1:15)']);
end

function neurons = empty_neurons(n)
neurons = repmat(struct('cluster_id', [], 'channel', [], 'area', '', ...
    'spike_times', [], 'spike_samples', []), 1, n);
end

function neuron = make_neuron(cluster_id, channel, area, spike_times)
neuron.cluster_id = cluster_id;
neuron.channel = channel;
neuron.area = area;
neuron.spike_times = spike_times(:);
neuron.spike_samples = round(spike_times(:) * 30000);
end

function path = spike_file()
path = fullfile(tempdir, 'sync_stat_spike_continuous.dat');
raw = zeros(384, 100);
raw(3, 50) = 80;
write_raw(path, raw);
end

function path = zero_file()
path = fullfile(tempdir, 'sync_stat_zero_continuous.dat');
write_raw(path, zeros(384, 100));
end

function write_raw(path, raw)
file_id = fopen(path, 'w', 'ieee-le');
assert(file_id >= 0, 'Could not write %s.', path);
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
fwrite(file_id, raw, 'int16');
end

function choice = check_choice()
choice.min_scale = 1;
choice.n_channels = 384;
choice.n_check = 5;
choice.samples_before = 30;
choice.samples_after = 30;
end

function delete_if_present(path)
if exist(path, 'file') == 2
    delete(path);
end
end

function check_spatial_stack()
local = {[-20; 0; 40], [1; 4; 1]};
flat = {[-20; 0; 40], [2; 2; 2]};
summary = spatial_range_stack({local{1}, flat{1}}, {local{2}, flat{2}}, -40:20:40);
assert(summary.n == 2, 'Both spatial profiles should enter the stack.');
center = summary.centers == 0;
assert(abs(summary.mean_range(center) - 1) < 1e-12, ...
    'The mean range at the neuron channel is not 1.');
side = summary.centers == -20;
assert(abs(summary.mean_range(side) - 0.625) < 1e-12, ...
    'The mean of the normalized side channels is wrong.');
local_row = summary.normalized(:, summary.centers == 40) == 0.25;
assert(any(local_row) && abs(summary.normalized(local_row, center) - 1) < 1e-12, ...
    'The local profile was not divided by the range at its own channel.');
dropped = spatial_range_stack({0}, {0}, -40:20:40);
assert(dropped.n == 0, 'A flat neuron channel stayed in the spatial stack.');
times = [9 + (0:24)' * 0.05; 100 + (0:9)'];
samples = (1:numel(times))';
picked = template_spike_samples(times, samples, 10, [-1.5 1], 5, 20);
assert(numel(picked) == 5 && all(picked <= 25), ...
    'The template did not keep evenly spaced spikes from the play-bout windows.');
few = template_spike_samples(times(1:10), samples(1:10), 10, [-1.5 1], 100, 20);
assert(numel(few) == 10, ...
    'A short train was replaced when it had fewer than the window minimum.');
end
