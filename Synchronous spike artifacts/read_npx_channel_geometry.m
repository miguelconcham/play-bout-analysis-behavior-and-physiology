function [chan_map, xcoords, ycoords] = read_npx_channel_geometry(session_dir)
%READ_NPX_CHANNEL_GEOMETRY  Probe x,y in micrometers for one NPX session.
%
% Loads chann_map_PAG.mat. chan_map(site) is the channel number used in the
% waveform mats (cluster ch + 1). xcoords and ycoords are the site position
% in micrometers. On NPX1, a map that skips channel 192 is repaired the
% same way as the area assignment: channel 192 is inserted between 191 and
% 193 at x = 43 um, y = 1900 um.

map_file = fullfile(session_dir, 'chann_map_PAG.mat');
if exist(map_file, 'file') ~= 2
    error('read_npx_channel_geometry:map', ...
        'Missing chann_map_PAG.mat in %s.', session_dir);
end
loaded = load(map_file, 'chanMap', 'xcoords', 'ycoords');
chan_map = double(loaded.chanMap(:));
xcoords = double(loaded.xcoords(:));
ycoords = double(loaded.ycoords(:));
if numel(chan_map) ~= numel(xcoords) || numel(chan_map) ~= numel(ycoords)
    error('read_npx_channel_geometry:map', ...
        'chanMap, xcoords, and ycoords differ in length in %s.', map_file);
end
if is_npx2(xcoords) || any(chan_map == 192)
    return
end
pos_191 = find(chan_map == 191, 1);
pos_193 = find(chan_map == 193, 1);
if isempty(pos_191) || isempty(pos_193) || pos_193 ~= pos_191 + 1
    error('read_npx_channel_geometry:map', ...
        'Inconsistent channel map in %s.', map_file);
end
xcoords = [xcoords; NaN];
xcoords(pos_193 + 1:end) = xcoords(pos_193:end - 1);
xcoords(pos_193) = 43;
ycoords = [ycoords; NaN];
ycoords(pos_193 + 1:end) = ycoords(pos_193:end - 1);
ycoords(pos_193) = 1900;
chan_map = [chan_map; NaN];
chan_map(pos_193 + 1:end) = chan_map(pos_193:end - 1);
chan_map(pos_193) = 192;
end

function npx2 = is_npx2(xcoords)
npx2_x = [8 40; 258 290; 508 540; 758 790];
npx2 = any(ismember(xcoords, npx2_x(:)));
end
