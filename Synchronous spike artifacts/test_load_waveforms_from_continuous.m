function test_load_waveforms_from_continuous()
%TEST_LOAD_WAVEFORMS_FROM_CONTINUOUS  Synthetic continuous-file checks.
%
% Builds a small little-endian int16 file whose value at each
% (channel, sample) is unique, then checks channel capping, spike
% alignment, and windows that run off either end of the recording.

n_channels = 384;
n_samples = 40;
dat_file = write_coded_continuous(n_channels, n_samples);
file_cleaner = onCleanup(@() delete_if_present(dat_file)); %#ok<NASGU>

samples_before = 10;
samples_after = 8;

% Interior channel: full neighborhood, both sides.
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    dat_file, [30 12 30], 100, 5, ...
    'NChannels', n_channels, ...
    'SamplesBefore', samples_before, ...
    'SamplesAfter', samples_after);
assert(isequal(channels_kept, 95:105), 'Interior channel window was not 95:105.');
assert(isequal(time_samples, -samples_before:samples_after), 'Time axis is not centered on the spike.');
assert(isequal(size(waveforms), [3, numel(channels_kept), numel(time_samples)]), ...
    'Waveform array is not n_spikes x n_channels x n_times.');
assert_snippet(waveforms, 1, 30, channels_kept, samples_before, samples_after, n_channels, n_samples);
assert_snippet(waveforms, 2, 12, channels_kept, samples_before, samples_after, n_channels, n_samples);
assert_snippet(waveforms, 3, 30, channels_kept, samples_before, samples_after, n_channels, n_samples);

% Low edge: drop channels below 1, do not extend the high side to compensate.
[waveforms, channels_kept] = load_waveforms_from_continuous( ...
    dat_file, 20, 3, 5, ...
    'NChannels', n_channels, ...
    'SamplesBefore', samples_before, ...
    'SamplesAfter', samples_after);
assert(isequal(channels_kept, 1:8), 'Channel 3 range 5 should cap at 1:8.');
assert(numel(channels_kept) < 11, 'Low-edge window was recentered to a full neighborhood.');
assert_snippet(waveforms, 1, 20, channels_kept, samples_before, samples_after, n_channels, n_samples);
assert(channels_kept(3) == 3, 'Neuron channel 3 should be the third kept channel.');

% High edge: drop channels above 384, do not extend the low side.
[waveforms, channels_kept] = load_waveforms_from_continuous( ...
    dat_file, int64(22), 382, 5, ...
    'NChannels', n_channels, ...
    'SamplesBefore', 4, ...
    'SamplesAfter', 4);
assert(isequal(channels_kept, 377:384), 'Channel 382 range 5 should cap at 377:384.');
assert(numel(channels_kept) < 11, 'High-edge window was recentered to a full neighborhood.');
assert_snippet(waveforms, 1, 22, channels_kept, 4, 4, n_channels, n_samples);
assert(channels_kept(end) == 384, 'Capped window should still include channel 384.');

% First and last channel, and a range of zero.
[~, only_first] = load_waveforms_from_continuous( ...
    dat_file, 5, 1, 0, 'NChannels', n_channels, 'SamplesBefore', 0, 'SamplesAfter', 0);
[waveforms, only_last, time_samples] = load_waveforms_from_continuous( ...
    dat_file, 5, 384, 0, 'NChannels', n_channels, 'SamplesBefore', 0, 'SamplesAfter', 0);
assert(isequal(only_first, 1), 'Range 0 on channel 1 should keep channel 1.');
assert(isequal(only_last, 384), 'Range 0 on channel 384 should keep channel 384.');
assert(isequal(time_samples, 0), 'A zero-length margin should keep the spike sample only.');
assert(waveforms(1, 1, 1) == code_at(5, 384, n_channels), ...
    'Spike sample on channel 384 does not match the file.');

% Range wider than the probe.
[~, all_channels] = load_waveforms_from_continuous( ...
    dat_file, 8, 10, 1000, ...
    'NChannels', n_channels, 'SamplesBefore', 1, 'SamplesAfter', 1);
assert(isequal(all_channels, 1:n_channels), 'A huge range should keep every channel.');

% Window hanging off the start and the end of the recording.
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    dat_file, [3, n_samples - 1], 20, 1, ...
    'NChannels', n_channels, 'SamplesBefore', 6, 'SamplesAfter', 5);
assert_snippet(waveforms, 1, 3, channels_kept, 6, 5, n_channels, n_samples);
assert_snippet(waveforms, 2, n_samples - 1, channels_kept, 6, 5, n_channels, n_samples);
assert(all(isnan(reshape(waveforms(1, :, 1:3), 1, []))), ...
    'Samples before the file should be NaN.');
assert(all(isnan(reshape(waveforms(2, :, end-3:end), 1, []))), ...
    'Samples after the file should be NaN.');
assert(~any(isnan(waveforms(1, :, time_samples == 0))), ...
    'The spike sample itself should be present at the start edge.');
assert(~any(isnan(waveforms(2, :, time_samples == 0))), ...
    'The spike sample itself should be present at the end edge.');

% No spikes: still report the channel window, do not read a waveform.
[waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    dat_file, [], 15, 2, ...
    'NChannels', n_channels, 'SamplesBefore', 7, 'SamplesAfter', 9);
assert(isequal(size(waveforms), [0, 5, 17]), 'Empty spike list should keep the array shape.');
assert(isequal(channels_kept, 13:17), 'Empty spike list should still cap the channel window.');
assert(isequal(time_samples, -7:9), 'Empty spike list should still return the time axis.');

% A file that is not 384 channels.
small_channels = 16;
small_samples = 25;
small_file = write_coded_continuous(small_channels, small_samples);
small_cleaner = onCleanup(@() delete_if_present(small_file)); %#ok<NASGU>
[waveforms, channels_kept] = load_waveforms_from_continuous( ...
    small_file, [2, 25, 10], 2, 4, ...
    'NChannels', small_channels, 'SamplesBefore', 3, 'SamplesAfter', 3);
assert(isequal(channels_kept, 1:6), 'Small file should cap at its own channel count.');
assert_snippet(waveforms, 1, 2, channels_kept, 3, 3, small_channels, small_samples);
assert_snippet(waveforms, 2, 25, channels_kept, 3, 3, small_channels, small_samples);
assert_snippet(waveforms, 3, 10, channels_kept, 3, 3, small_channels, small_samples);

assert_error_id(@() load_waveforms_from_continuous( ...
    [tempname, '.dat'], 1, 1, 1), 'load_waveforms_from_continuous:file');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 0, 10, 1, 'NChannels', n_channels), ...
    'load_waveforms_from_continuous:spikeSample');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, n_samples + 1, 10, 1, 'NChannels', n_channels), ...
    'load_waveforms_from_continuous:spikeSample');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 1.5, 10, 1, 'NChannels', n_channels), ...
    'load_waveforms_from_continuous:spikeSample');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 4, 0, 1, 'NChannels', n_channels), ...
    'load_waveforms_from_continuous:channel');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 4, n_channels + 1, 1, 'NChannels', n_channels), ...
    'load_waveforms_from_continuous:channel');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 4, [10 11], 1, 'NChannels', n_channels), ...
    'load_waveforms_from_continuous:channel');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 4, 10, -1, 'NChannels', n_channels), ...
    'load_waveforms_from_continuous:channelRange');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 4, 10, 1, 'SamplesBefore', -2), ...
    'load_waveforms_from_continuous:window');
assert_error_id(@() load_waveforms_from_continuous( ...
    dat_file, 4, 10, 1, 'NotAParameter', 1), ...
    'load_waveforms_from_continuous:nameValue');

bad_file = [tempname, '.dat'];
bad_id = fopen(bad_file, 'wb');
fwrite(bad_id, uint8(1:10));
fclose(bad_id);
bad_cleaner = onCleanup(@() delete_if_present(bad_file)); %#ok<NASGU>
assert_error_id(@() load_waveforms_from_continuous( ...
    bad_file, 1, 1, 0), 'load_waveforms_from_continuous:size');

fprintf('load_waveforms_from_continuous tests passed.\n');
end

function assert_snippet(waveforms, spike_index, spike_sample, channels_kept, ...
    samples_before, samples_after, n_channels, n_samples)
got = reshape(waveforms(spike_index, :, :), numel(channels_kept), []);
expected = expected_snippet(spike_sample, channels_kept, ...
    samples_before, samples_after, n_channels, n_samples);
assert(isequal(size(got), size(expected)), 'Snippet size does not match the requested window.');
missing = isnan(expected);
assert(isequal(isnan(got), missing), 'NaN padding does not match samples outside the file.');
assert(all(got(~missing) == expected(~missing)), ...
    'Snippet values do not match the continuous file.');
end

function snippet = expected_snippet(spike_sample, channels_kept, ...
    samples_before, samples_after, n_channels, n_samples)
file_samples = spike_sample + (-samples_before:samples_after);
snippet = nan(numel(channels_kept), numel(file_samples));
for time_index = 1:numel(file_samples)
    sample = file_samples(time_index);
    if sample < 1 || sample > n_samples
        continue
    end
    snippet(:, time_index) = code_at(sample, channels_kept(:), n_channels);
end
end

function value = code_at(sample, channel, n_channels)
value = (sample - 1) * n_channels + channel;
end

function dat_file = write_coded_continuous(n_channels, n_samples)
codes = int16(bsxfun(@plus, (1:n_channels)', (0:n_samples - 1) * n_channels));
assert(max(codes(:)) <= intmax('int16'), 'Test code does not fit in int16.');
dat_file = [tempname, '.dat'];
file_id = fopen(dat_file, 'wb', 'ieee-le');
assert(file_id > 0, 'Could not create the synthetic continuous file.');
fwrite(file_id, codes, 'int16');
fclose(file_id);
end

function assert_error_id(call, identifier)
caught = '';
threw = false;
try
    call();
catch ME
    threw = true;
    caught = ME.identifier;
end
assert(threw, ['Expected error ', identifier, '.']);
assert(strcmp(caught, identifier), ...
    ['Expected ', identifier, ', got ', caught, '.']);
end

function delete_if_present(path)
if exist(path, 'file') == 2
    delete(path);
end
end
