function options = distance_coverage_options(overrides)
%DISTANCE_COVERAGE_OPTIONS  Cuts used to match play and non-play on distance.
%
% The numbers are the same cuts as the high-speed section of
% Analyze_speed_freq_relation.m. high_z = 1 means one within-session
% standard deviation from that session's own mean distance.
% Far samples have distance z >= high_z. Close samples have z <= -high_z.
%
% Durations in the saved mats are spans of overlapping spectrogram frames.
% min_nonplay_s and min_play_s apply to a whole tail (far or close).
% The ceiling bin uses min_ceiling_play_s and min_ceiling_nonplay_s.
% The count curve uses min_bin_play_s and min_bin_nonplay_s.

if nargin < 1
    overrides = [];
end

options = struct( ...
    'high_z', 1, ...
    'ceiling_bin_width', 0.5, ...
    'min_nonplay_s', 2, ...
    'min_play_s', 1, ...
    'min_ceiling_nonplay_s', 1, ...
    'min_ceiling_play_s', 0.5, ...
    'curve_bin_width', 0.25, ...
    'min_bin_nonplay_s', 0.5, ...
    'min_bin_play_s', 0.5, ...
    'curve_z_limit', 6);

if isempty(overrides)
    return
end

fields = fieldnames(overrides);
for field_n = 1:numel(fields)
    name = fields{field_n};
    if ~isfield(options, name)
        error('distance_coverage_options:UnknownOption', ...
            'Unknown distance-coverage option %s.', name);
    end
    options.(name) = overrides.(name);
end
end
