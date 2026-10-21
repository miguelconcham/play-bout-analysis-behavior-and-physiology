function [waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
    continuous_file, spike_samples, neuron_channel, channel_range, varargin)
%LOAD_WAVEFORMS_FROM_CONTINUOUS  Spike snippets from a multiplexed int16 file.
%
%   [waveforms, channels_kept, time_samples] = ...
%       load_waveforms_from_continuous(continuous_file, spike_samples, ...
%           neuron_channel, channel_range)
%
% Cuts one waveform per spike from a continuous Neuropixels file without
% loading the recording. The file is little-endian int16, multiplexed as
% channel-then-sample, the same layout as
%
%   reshape(fread(fopen(continuous_file, 'r', 'ieee-le'), 'int16'), 384, [])
%
% in the GENERATE scripts. Each spike is centered on the sample you pass.
%
% Inputs
%   continuous_file   Path to the continuous file (for example continuous.dat).
%   spike_samples     Spike times in samples, one per spike. Sample 1 is the
%                     first frame in the file (the first column of the
%                     reshaped array above). Pass kilosort spike_times as
%                     stored; if a file is 0-based, add 1 before calling.
%   neuron_channel    Channel of this neuron, 1-based, matching rows of the
%                     reshaped continuous array. cluster_info.ch in this
%                     repository is 0-based, so use that value plus 1.
%   channel_range     Channels kept on each side of neuron_channel. The
%                     neighborhood is
%                     (neuron_channel-channel_range):(neuron_channel+channel_range),
%                     the same span the waveform checker built with a fixed
%                     range of 5. The window stops at channel 1 and at the
%                     last channel. Channel 3 with range 5 on a 384-channel
%                     file keeps 1:8. Channel 382 with range 5 keeps 377:384.
%
% Name-value pairs
%   'NChannels'       Channels per sample in the file. Default 384.
%   'SamplesBefore'   Samples kept before the spike sample. Default 30
%                     (1 ms at 30 kHz), the half-width the waveform checker
%                     used around the spike.
%   'SamplesAfter'    Samples kept after the spike sample. Default 30.
%
% Outputs
%   waveforms         n_spikes x n_channels_kept x n_times, double, in the
%                     raw int16 counts of the file. Index it as
%                     waveforms(spike, channel, time). channel indexes
%                     channels_kept, not the probe. A window that runs off
%                     the start or end of the file is filled with NaN.
%   channels_kept     1 x n_channels_kept probe channels that were read,
%                     after the end-of-probe cap.
%   time_samples      1 x n_times offsets from the spike sample. 0 is the
%                     spike. Divide by the sampling rate and multiply by
%                     1000 for milliseconds.
%
% Call once per neuron. For a synchronous pair, use the same spike_samples
% and each neuron's own channel, then compare the two waveform stacks.
%
% Example
%   neuron_channel = double(cluster_info.ch(cluster_info.id == id)) + 1;
%   spike_samples = double(spike_times(spike_clusters == id));
%   [waveforms, channels_kept, time_samples] = ...
%       load_waveforms_from_continuous( ...
%           fullfile(npx_data_dir, 'continuous.dat'), ...
%           spike_samples, neuron_channel, 5);

if is_string(continuous_file)
    if numel(continuous_file) ~= 1
        error('load_waveforms_from_continuous:file', ...
            'continuous_file must be one path.');
    end
    continuous_file = char(continuous_file);
end
if ~ischar(continuous_file) || isempty(strtrim(continuous_file))
    error('load_waveforms_from_continuous:file', ...
        'continuous_file must be a path to a continuous int16 file.');
end

options = parse_options(varargin);
n_channels = options.n_channels;
samples_before = options.samples_before;
samples_after = options.samples_after;

neuron_channel = require_scalar_integer(neuron_channel, ...
    'load_waveforms_from_continuous:channel', 'neuron_channel');
channel_range = require_scalar_integer(channel_range, ...
    'load_waveforms_from_continuous:channelRange', 'channel_range');
if neuron_channel < 1 || neuron_channel > n_channels
    error('load_waveforms_from_continuous:channel', ...
        'neuron_channel must be from 1 to %d.', n_channels);
end
if channel_range < 0
    error('load_waveforms_from_continuous:channelRange', ...
        'channel_range must be a non-negative integer.');
end

spike_samples = require_sample_index(spike_samples);
channels_kept = channel_window(neuron_channel, channel_range, n_channels);
time_samples = -samples_before:samples_after;
n_times = numel(time_samples);
n_spikes = numel(spike_samples);

file_id = fopen(continuous_file, 'r', 'ieee-le');
if file_id < 0
    error('load_waveforms_from_continuous:file', ...
        'Could not open continuous file: %s', continuous_file);
end
file_closer = onCleanup(@() fclose(file_id)); %#ok<NASGU>

n_samples = continuous_sample_count(file_id, n_channels, continuous_file);
if any(spike_samples < 1) || any(spike_samples > n_samples)
    bad = find(spike_samples < 1 | spike_samples > n_samples, 1, 'first');
    error('load_waveforms_from_continuous:spikeSample', ...
        'Spike sample %g is outside 1..%d.', spike_samples(bad), n_samples);
end

waveforms = nan(n_spikes, numel(channels_kept), n_times);
if n_spikes == 0
    return
end

% Read in time order. Results stay in the order of spike_samples.
[sorted_samples, sort_order] = sort(spike_samples);
for spike_index = 1:n_spikes
    snippet = read_spike_window( ...
        file_id, sorted_samples(spike_index), channels_kept, ...
        samples_before, samples_after, n_channels, n_samples);
    waveforms(sort_order(spike_index), :, :) = reshape( ...
        snippet, 1, size(snippet, 1), size(snippet, 2));
end

% file_closer releases the file when this function returns.
end

function yes = is_string(value)
yes = false;
if exist('isstring', 'builtin')
    yes = isstring(value);
end
end

function options = parse_options(name_value)
options.n_channels = 384;
options.samples_before = 30;
options.samples_after = 30;

if mod(numel(name_value), 2) ~= 0
    error('load_waveforms_from_continuous:nameValue', ...
        'Optional arguments must be name-value pairs.');
end

for pair_index = 1:2:numel(name_value)
    name = name_value{pair_index};
    value = name_value{pair_index + 1};
    if is_string(name)
        name = char(name);
    end
    if ~ischar(name)
        error('load_waveforms_from_continuous:nameValue', ...
            'Optional argument names must be text.');
    end
    switch lower(strtrim(name))
        case 'nchannels'
            options.n_channels = require_scalar_integer(value, ...
                'load_waveforms_from_continuous:nChannels', 'NChannels');
            if options.n_channels < 1
                error('load_waveforms_from_continuous:nChannels', ...
                    'NChannels must be a positive integer.');
            end
        case 'samplesbefore'
            options.samples_before = require_scalar_integer(value, ...
                'load_waveforms_from_continuous:window', 'SamplesBefore');
            if options.samples_before < 0
                error('load_waveforms_from_continuous:window', ...
                    'SamplesBefore must be a non-negative integer.');
            end
        case 'samplesafter'
            options.samples_after = require_scalar_integer(value, ...
                'load_waveforms_from_continuous:window', 'SamplesAfter');
            if options.samples_after < 0
                error('load_waveforms_from_continuous:window', ...
                    'SamplesAfter must be a non-negative integer.');
            end
        otherwise
            error('load_waveforms_from_continuous:nameValue', ...
                'Unknown parameter "%s".', name);
    end
end
end

function channels_kept = channel_window(neuron_channel, channel_range, n_channels)
% Stop at channel 1 and at the last channel.
first_channel = max(1, neuron_channel - channel_range);
last_channel = min(n_channels, neuron_channel + channel_range);
channels_kept = first_channel:last_channel;
end

function spike_samples = require_sample_index(spike_samples)
if ~isnumeric(spike_samples) || ~isreal(spike_samples)
    error('load_waveforms_from_continuous:spikeSample', ...
        'spike_samples must be real sample indices.');
end
spike_samples = double(spike_samples(:));
if any(~isfinite(spike_samples)) || any(spike_samples ~= round(spike_samples))
    error('load_waveforms_from_continuous:spikeSample', ...
        'spike_samples must be finite integer sample indices, not seconds.');
end
spike_samples = round(spike_samples);
end

function value = require_scalar_integer(value, identifier, name)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) ...
        || value ~= round(value)
    error(identifier, '%s must be a finite integer scalar.', name);
end
value = round(double(value));
end

function n_samples = continuous_sample_count(file_id, n_channels, continuous_file)
start_position = ftell(file_id);
fseek(file_id, 0, 'eof');
n_bytes = ftell(file_id);
fseek(file_id, start_position, 'bof');

bytes_per_frame = n_channels * 2;
if n_bytes <= 0 || mod(n_bytes, bytes_per_frame) ~= 0
    error('load_waveforms_from_continuous:size', ...
        ['%s is %d bytes, which is not a positive multiple of %d int16 ' ...
         'channels. Set NChannels if this file is not 384 channels.'], ...
        continuous_file, n_bytes, n_channels);
end
n_samples = n_bytes / bytes_per_frame;
end

function snippet = read_spike_window(file_id, spike_sample, channels_kept, ...
    samples_before, samples_after, n_channels, n_samples)
n_kept = numel(channels_kept);
n_times = samples_before + samples_after + 1;
snippet = nan(n_kept, n_times);

requested_start = spike_sample - samples_before;
requested_end = spike_sample + samples_after;
read_start = max(1, requested_start);
read_end = min(n_samples, requested_end);
n_read = read_end - read_start + 1;

byte_offset = (read_start - 1) * n_channels * 2;
if fseek(file_id, byte_offset, 'bof') ~= 0
    error('load_waveforms_from_continuous:read', ...
        'Could not seek to sample %d.', read_start);
end

[raw, count] = fread(file_id, [n_channels, n_read], 'int16=>double');
if count ~= n_channels * n_read
    error('load_waveforms_from_continuous:read', ...
        'Read %d of %d values at sample %d.', count, n_channels * n_read, read_start);
end

destination = (read_start - requested_start) + 1;
snippet(:, destination:destination + n_read - 1) = raw(channels_kept, :);
end
