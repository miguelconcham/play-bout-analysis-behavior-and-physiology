function selected = select_mid_lpag_channels(xcoords, ycoords, chanMap, lpag_limits, npx_type)
%SELECT_MID_LPAG_CHANNELS  Site nearest the middle of each LPAG depth range.
%
% This is the geometry rule in GENERATE_PSTH_PLAY_BOUT and
% GENERATE_PLAY_SESSION_LPAG_SPECTRA. lpag_limits is n-by-3:
% [ProbeNum, depth_start, depth_end], in the same units as ycoords.
%
% The returned channel is chanMap at that site. LFP rows in this project
% are indexed by that value (GENERATE_PSTH_PLAY_BOUT and
% GENERATE_SPIKE_TRAIN_PARSING). read_npx_channel_geometry documents
% chanMap as the 1-based channel, equal to cluster ch + 1.
%
% NPX1 maps that skip channel 192 are repaired the same way as the
% existing PAG code before the site is chosen.

xcoords = double(xcoords(:));
ycoords = double(ycoords(:));
chanMap = double(chanMap(:));
lpag_limits = double(lpag_limits);
if numel(xcoords) ~= numel(ycoords) || numel(xcoords) ~= numel(chanMap)
    error('select_mid_lpag_channels:map', ...
        'xcoords, ycoords, and chanMap must have the same length.');
end
if isempty(lpag_limits) || size(lpag_limits, 2) ~= 3
    error('select_mid_lpag_channels:limits', ...
        'LPAG limits must be [ProbeNum depth_start depth_end].');
end
if ~ismember(npx_type, [1 2])
    error('select_mid_lpag_channels:npx', 'npx_type must be 1 or 2.');
end

if npx_type == 1
    [xcoords, ycoords, chanMap] = repair_npx1_missing_channel( ...
        xcoords, ycoords, chanMap);
end

hard_coded_x = [8 40; 258 290; 508 540; 758 790];
selected = repmat(struct( ...
    'channel', NaN, ...
    'lfp_row', NaN, ...
    'map_index', NaN, ...
    'probe_number', NaN, ...
    'x_um', NaN, ...
    'y_um', NaN, ...
    'depth_range', [NaN NaN]), size(lpag_limits, 1), 1);

for row_n = 1:size(lpag_limits, 1)
    probe_number = lpag_limits(row_n, 1);
    depth_range = lpag_limits(row_n, 2:3);
    in_area = ycoords >= depth_range(1) & ycoords <= depth_range(2);
    if npx_type == 2
        if probe_number < 1 || probe_number > size(hard_coded_x, 1)
            error('select_mid_lpag_channels:probe', ...
                'Invalid probe number %g in LPAG area limits.', probe_number);
        end
        in_area = in_area & ismember(xcoords, hard_coded_x(probe_number, :));
    end
    in_area = in_area & isfinite(xcoords) & isfinite(ycoords) & isfinite(chanMap);
    if ~any(in_area)
        error('select_mid_lpag_channels:empty', ...
            'No mapped channels fall in LPAG limits for row %d.', row_n);
    end

    mean_location = mean([xcoords(in_area), ycoords(in_area)], 1);
    distance = abs(xcoords - mean_location(1)) + abs(ycoords - mean_location(2));
    distance(~in_area) = Inf;
    [~, map_index] = min(distance);
    channel = chanMap(map_index);
    if ~isfinite(channel) || channel ~= round(channel) || channel < 1
        error('select_mid_lpag_channels:channel', ...
            'Selected LPAG channel %g is not a positive LFP row.', channel);
    end

    selected(row_n).channel = channel;
    selected(row_n).lfp_row = channel;
    selected(row_n).map_index = map_index;
    selected(row_n).probe_number = probe_number;
    selected(row_n).x_um = xcoords(map_index);
    selected(row_n).y_um = ycoords(map_index);
    selected(row_n).depth_range = depth_range;
end
end

function [xcoords, ycoords, chanMap] = repair_npx1_missing_channel(xcoords, ycoords, chanMap)
if any(chanMap == 192)
    return
end
pos_191 = find(chanMap == 191, 1);
pos_193 = find(chanMap == 193, 1);
if isempty(pos_191) && isempty(pos_193)
    return
end
if isempty(pos_191) || isempty(pos_193) || pos_193 ~= pos_191 + 1
    error('select_mid_lpag_channels:chanMap', ...
        'Inconsistent NPX1 channel map around missing channel 192.');
end
xcoords = [xcoords; NaN];
xcoords(pos_193 + 1:end) = xcoords(pos_193:end - 1);
xcoords(pos_193) = 43;
ycoords = [ycoords; NaN];
ycoords(pos_193 + 1:end) = ycoords(pos_193:end - 1);
ycoords(pos_193) = 1900;
chanMap = [chanMap; NaN];
chanMap(pos_193 + 1:end) = chanMap(pos_193:end - 1);
chanMap(pos_193) = 192;
end
