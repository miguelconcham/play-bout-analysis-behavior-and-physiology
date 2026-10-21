function distance_um = pair_channel_distance(channel_a, channel_b, chan_map, xcoords, ycoords)
%PAIR_CHANNEL_DISTANCE  Euclidean distance between two probe channels.
%
% Distance is in micrometers, in the probe plane of xcoords and ycoords.
% A channel missing from chan_map returns NaN. The same channel returns 0.

distance_um = NaN;
[x_a, y_a, found_a] = site_xy(channel_a, chan_map, xcoords, ycoords);
[x_b, y_b, found_b] = site_xy(channel_b, chan_map, xcoords, ycoords);
if ~found_a || ~found_b
    return
end
distance_um = hypot(x_a - x_b, y_a - y_b);
end

function [x_site, y_site, found] = site_xy(channel, chan_map, xcoords, ycoords)
x_site = NaN;
y_site = NaN;
found = false;
if isempty(chan_map) || isempty(channel) || ~isfinite(channel)
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
