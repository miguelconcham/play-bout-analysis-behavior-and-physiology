function test_plot_random_good_cluster_checks()
%TEST_PLOT_RANDOM_GOOD_CLUSTER_CHECKS  Selection, waveform channels, and ACG.

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end
addpath(fullfile(fileparts(this_dir), 'Custom functions'));

cluster_info.id = [10; 11; 12];
cluster_info.ch = [4; 5; 6];
cluster_info.group = {'good'; ' mua '; 'good'};
spike_times = [0 10 20 30 40];
spike_clusters = [10 10 11 12 12];
clusters = good_clusters_from_sorting(cluster_info, spike_times, spike_clusters, 'PAG');
assert(numel(clusters) == 2, 'Only good clusters should be kept.');
assert(clusters(1).id == 10, 'First good cluster id was not kept.');
assert(clusters(1).channel == 5, 'Phy channel 4 should be continuous channel 5.');
assert(isequal(clusters(1).spike_samples, [1; 11]), ...
    'A spike_times file that starts at 0 should be shifted to sample 1.');
assert(isequal(clusters(2).spike_samples, [31; 41]), ...
    'Spikes of the second good cluster were not collected.');
assert(strcmp(clusters(1).source, 'PAG'), 'Cluster source was not kept.');

one_based.id = 3;
one_based.ch = 0;
one_based.group = {'good'};
kept = good_clusters_from_sorting(one_based, [1 5 9], [3 3 3]);
assert(isequal(kept.spike_samples, [1; 5; 9]), ...
    'Spike times that already start at sample 1 should stay unchanged.');
assert(kept.channel == 1, 'Phy channel 0 should be continuous channel 1.');

n_channels = 32;
n_samples = 100;
dat_file = write_coded_continuous(n_channels, n_samples);
file_cleaner = onCleanup(@() delete_if_present(dat_file)); %#ok<NASGU>
spike_samples = [40; 70];
neuron_channel = 8;
result = plot_cluster_waveform_and_acg(dat_file, spike_samples, neuron_channel, ...
    'ChannelRange', 5, ...
    'NWaveforms', 10, ...
    'NChannels', n_channels, ...
    'SampleRate', 1000, ...
    'AcgWindow', [-0.1 0.1], ...
    'AcgBin', 0.01, ...
    'SessionName', 'B2S2 1110 Single2', ...
    'ClusterID', 10, ...
    'Plot', false);

assert(isequal(result.channels_kept, 3:13), ...
    'Channel range 5 around channel 8 should keep channels 3:13.');
center = find(result.time_samples == 0);
assert(numel(center) == 1, 'The spike sample should sit at time 0.');
for spike_index = 1:numel(result.spike_samples)
    expected = (result.spike_samples(spike_index) - 1) * n_channels + result.channels_kept(:);
    got = reshape(result.waveforms(spike_index, :, center), [], 1);
    assert(isequal(got, expected), 'Loaded waveform does not match the continuous file.');
end

expected_lags = autocorrelogram(spike_samples / 1000, [-0.1 0.1]);
assert(numel(result.lags) == numel(expected_lags), 'ACG lag count does not match autocorrelogram.m.');
assert(max(abs(sort(result.lags) - sort(expected_lags(:)))) < 1e-12, ...
    'ACG lags do not match autocorrelogram.m.');

report.ok = false;
report.rows = [];
summary = plot_random_good_cluster_checks('Report', report);
assert(isempty(summary), 'A failed folder match should not plot clusters.');

fprintf('plot_random_good_cluster_checks tests passed.\n');
end

function dat_file = write_coded_continuous(n_channels, n_samples)
codes = int16(bsxfun(@plus, (1:n_channels)', (0:n_samples - 1) * n_channels));
dat_file = [tempname '.dat'];
file_id = fopen(dat_file, 'wb', 'ieee-le');
assert(file_id > 0, 'Could not create the synthetic continuous file.');
fwrite(file_id, codes, 'int16');
fclose(file_id);
end

function delete_if_present(path)
if exist(path, 'file') == 2
    delete(path);
end
end
