function profile = sync_waveform_spatial_range(channels, waveform, neuron_channel, ...
    chan_map, xcoords, ycoords)
%SYNC_WAVEFORM_SPATIAL_RANGE  Peak-to-peak range versus distance on the probe.
%
% channels is the saved channel list. waveform has one row per channel and
% one column per time sample of the averaged synchronous waveform. The
% range on a channel is the maximum of that row minus the minimum.
%
% distance_um is the Euclidean distance in the probe plane from the neuron
% channel, in micrometers. Channels with a larger y than the neuron are
% positive. Channels with a smaller y are negative. A channel at the same
% depth takes its sign from x. The neuron channel itself is 0. Rows are
% sorted from negative distance to positive distance so a line through
% them peaks at the neuron when the spike is local.

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

[x0, y0, neuron_found] = site_xy(neuron_channel, chan_map, xcoords, ycoords);
profile.neuron_found = neuron_found;
n_channels = numel(channels);
distance = nan(n_channels, 1);
range_value = nan(n_channels, 1);
for index = 1:n_channels
    range_value(index) = peak_to_peak(waveform(index, :));
    if ~neuron_found
        continue
    end
    [x_site, y_site, found] = site_xy(channels(index), chan_map, xcoords, ycoords);
    if found
        distance(index) = signed_distance(x_site, y_site, x0, y0);
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
profile.neuron_found = false;
profile.reason = '';
end
