function profile = sync_waveform_spatial_range(channels, waveform, neuron_channel, ...
    chan_map, xcoords, ycoords)
%SYNC_WAVEFORM_SPATIAL_RANGE  Peak-to-peak range versus distance on the probe.
%
% channels is the saved channel list. waveform has one row per channel and
% one column per time sample of the averaged synchronous waveform. The
% range on a channel is the maximum of that row minus the minimum.
%
% distance_um is the Euclidean distance in the probe plane from the channel
% with the largest peak-to-peak range, in micrometers. Channels with a
% larger y than that channel are positive. Channels with a smaller y are
% negative. A channel at the same depth takes its sign from x. The
% largest-range channel is at 0. When several channels share that range,
% the assigned neuron channel is kept if it is one of them; otherwise the
% tied channel closest to the assigned channel is used. Rows are sorted
% from negative distance to positive distance.

profile = empty_profile();
channels = double(channels(:));
waveform = double(waveform);
if isempty(channels) || isempty(waveform)
    return
end
if size(waveform, 1) ~= numel(channels) && size(waveform, 2) == numel(channels)
    waveform = waveform.';
end
if size(waveform, 1) ~= numel(channels)
    profile.reason = 'waveform rows do not match the channel list';
    return
end

n_channels = numel(channels);
range_value = nan(n_channels, 1);
x_site = nan(n_channels, 1);
y_site = nan(n_channels, 1);
on_map = false(n_channels, 1);
for index = 1:n_channels
    range_value(index) = peak_to_peak(waveform(index, :));
    [x_site(index), y_site(index), on_map(index)] = site_xy( ...
        channels(index), chan_map, xcoords, ycoords);
end
[x_neuron, y_neuron, neuron_found] = site_xy(neuron_channel, chan_map, xcoords, ycoords);
profile.neuron_found = neuron_found;
origin = peak_channel_index(range_value, on_map, channels, x_site, y_site, ...
    neuron_channel, neuron_found, x_neuron, y_neuron);
if origin < 1
    return
end
profile.center_channel = channels(origin);
distance = nan(n_channels, 1);
for index = 1:n_channels
    if on_map(index)
        distance(index) = signed_distance(x_site(index), y_site(index), ...
            x_site(origin), y_site(origin));
    end
end
keep = isfinite(distance) & isfinite(range_value);
profile.channel = channels(keep);
profile.distance_um = distance(keep);
profile.range = range_value(keep);
[profile.distance_um, order] = sort(profile.distance_um);
profile.channel = profile.channel(order);
profile.range = profile.range(order);
end

function origin = peak_channel_index(range_value, on_map, channels, x_site, y_site, ...
    neuron_channel, neuron_found, x_neuron, y_neuron)
origin = 0;
candidates = find(on_map & isfinite(range_value));
if isempty(candidates)
    return
end
peak = max(range_value(candidates));
tied = candidates(range_value(candidates) == peak);
neuron_index = find(channels == neuron_channel, 1);
if ~isempty(neuron_index) && any(tied == neuron_index)
    origin = neuron_index;
    return
end
origin = tied(1);
if ~neuron_found
    [~, pick] = min(channels(tied));
    origin = tied(pick);
    return
end
best_distance = inf;
for index = 1:numel(tied)
    separation = hypot(x_site(tied(index)) - x_neuron, y_site(tied(index)) - y_neuron);
    if separation < best_distance ...
            || (separation == best_distance && channels(tied(index)) < channels(origin))
        origin = tied(index);
        best_distance = separation;
    end
end
end

function distance = signed_distance(x_site, y_site, x0, y0)
dx = x_site - x0;
dy = y_site - y0;
distance = sqrt(dx .^ 2 + dy .^ 2);
direction = sign(dy);
if dy == 0
    direction = sign(dx);
end
distance = direction * distance;
end

function value = peak_to_peak(samples)
samples = samples(isfinite(samples));
if isempty(samples)
    value = NaN;
    return
end
value = max(samples) - min(samples);
end

function [x_site, y_site, found] = site_xy(channel, chan_map, xcoords, ycoords)
x_site = NaN;
y_site = NaN;
found = false;
if isempty(chan_map) || ~isfinite(channel)
    return
end
match = find(chan_map(:) == channel, 1);
if isempty(match)
    return
end
x_site = xcoords(match);
y_site = ycoords(match);
found = isfinite(x_site) && isfinite(y_site);
end

function profile = empty_profile()
profile.channel = zeros(0, 1);
profile.distance_um = zeros(0, 1);
profile.range = zeros(0, 1);
profile.center_channel = NaN;
profile.neuron_found = false;
profile.reason = '';
end
