function test_sync_event_statistics()
%TEST_SYNC_EVENT_STATISTICS  Checks that do not need the experiment file system.
%
% Coincidence mask, pair caps, waveform distance, raw-power percentile,
% and continuous-file choice. Run this file by name after add_repo_paths.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
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

[pairs, ~] = select_synchronous_pairs(example_neurons(), 10, ...
    sync_analysis_options(struct('max_pairs_per_neuron', 2, 'max_pairs', 2)));
assert(numel(pairs) == 2, 'Two disjoint uses of the best-sampled neuron were expected.');
assert(pairs(2).cluster_id_1 == 11 && pairs(2).cluster_id_2 == 33, ...
    'The second pair was not the next-largest spike count.');

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
low_power = compare_sync_scalar([zeros(30, 1); ones(200, 1)], 1:30, 40);
assert(low_power.ok && low_power.percentile <= 0.05, ...
    'Low synchronous values were not below the null.');
blocked = compare_sync_scalar(ones(20, 1), 1:20, 10);
assert(~blocked.ok, 'A neuron with no non-synchronous spikes should be skipped.');

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

artifact = zeros(80, 3, numel(time_samples));
artifact(:, 2, peak_at) = 20;
artifact(1:20, 1, peak_at) = 50;
artifact_result = compare_sync_waveforms(artifact, time_samples, 1:20, 5, [4 5 6], 25);
assert(artifact_result.ok, artifact_result.reason);
assert(artifact_result.percentile == 1, ...
    'A side-channel difference was not farther than the null.');
assert(artifact_result.percentile_channel == 0, ...
    'The neuron channel was marked dissimilar when its shape matched.');
assert(isequal(size(artifact_result.template_mean), [3, numel(time_samples)]), ...
    'Template mean is not channels by time.');

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
check_location_lookup();
fprintf('Synchronous-event statistic checks passed.\n');
end

function check_location_lookup()
if exist('readtable', 'file') ~= 2 && exist('readtable', 'builtin') ~= 5
    return
end
split_file = fullfile(tempdir, 'npx_locations_split.csv');
full_file = fullfile(tempdir, 'npx_locations_full.csv');
write_text(split_file, sprintf('%s\n%s\n%s\n', ...
    'Animal,Session,Location', ...
    'B2S2,1110 Single2,server/probeA/continuous.dat', ...
    'B4D4,0826 Dual,server/other/continuous.dat'));
write_text(full_file, sprintf('%s\n%s\n', ...
    'Session,Path', ...
    'B2S2 1110 Single2,D:/data/continuous.dat'));
split_paths = read_npx_continuous_locations(split_file, 'B2S2 1110 Single2');
full_paths = read_npx_continuous_locations(full_file, 'B2S2 1110 Single2');
delete_if_present(split_file);
delete_if_present(full_file);
assert(numel(split_paths) == 1 && ~isempty(strfind(split_paths{1}, 'probeA')), ...
    'Animal plus session did not return that row''s location.');
assert(numel(full_paths) == 1 && ~isempty(strfind(full_paths{1}, 'continuous.dat')), ...
    'The session-folder column did not return its path.');
end

function write_text(path, text)
file_id = fopen(path, 'w');
assert(file_id >= 0, 'Could not write %s.', path);
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
fwrite(file_id, text);
end

function neurons = example_neurons()
ab = 9.00 + (0:19)' * 0.02;
ac = 9.50 + (0:19)' * 0.02;
neurons = empty_neurons(3);
neurons(1) = make_neuron(11, 10, 'LPAG', [ab; ac; 100 + (1:4960)']);
neurons(2) = make_neuron(22, 20, 'DR', [ab; 20000 + (1:3980)']);
neurons(3) = make_neuron(33, 30, 'SupCol', [ac; 40000 + (1:2980)']);
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
