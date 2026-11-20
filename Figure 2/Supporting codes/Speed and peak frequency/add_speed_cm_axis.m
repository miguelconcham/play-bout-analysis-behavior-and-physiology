function add_speed_cm_axis(pixel_axis, scale)
%ADD_SPEED_CM_AXIS  Top ruler in cm/s for a pixels-per-second speed axis.
%
% The plotted points are not moved. linear multiplies the limits by
% 0.1334. log shifts the natural-log limits by log(0.1334).

if isempty(pixel_axis) || ~isgraphics(pixel_axis)
    return
end
cm_limits = speed_axis_in_cm(xlim(pixel_axis), scale);
overlay = axes('Parent', pixel_axis.Parent, ...
    'Position', pixel_axis.Position, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'Color', 'none', ...
    'YTick', [], ...
    'Box', 'off', ...
    'XLim', cm_limits, ...
    'HitTest', 'off');
overlay.YAxis.Visible = 'off';
if strcmp(scale, 'log')
    xlabel(overlay, 'log speed (cm/s)')
else
    xlabel(overlay, 'Speed (cm/s)')
end
linkprop([pixel_axis, overlay], 'Position');
end
