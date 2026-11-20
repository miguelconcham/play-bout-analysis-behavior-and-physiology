function cm_limits = speed_axis_in_cm(speed_limits, scale)
%SPEED_AXIS_IN_CM  Same speed axis, labelled in cm/s.
%
% linear: the plotted values are pixels per second, and cm/s is that
% value times 0.1334. log: the plotted values are natural log of pixels
% per second, and log(cm/s) is that value plus log(0.1334). Either way
% the data stay on the original coordinates.

speed_limits = double(speed_limits(:)).';
if numel(speed_limits) ~= 2 || any(~isfinite(speed_limits))
    error('speed_axis_in_cm:limits', 'speed_limits must be two finite values.');
end
if strcmp(scale, 'linear')
    cm_limits = speed_limits * pixel_to_cm();
elseif strcmp(scale, 'log')
    cm_limits = speed_limits + log(pixel_to_cm());
else
    error('speed_axis_in_cm:scale', 'scale must be linear or log.');
end
end
