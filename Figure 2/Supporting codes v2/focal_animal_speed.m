function [speed, speed_time] = focal_animal_speed(position, time)
%FOCAL_ANIMAL_SPEED  Speed from successive focal-animal positions.
%
% position is n x 2 in image pixels. time is the NPX time of each position,
% in seconds. speed(k) is the pixel displacement from sample k to k+1
% divided by the time step, so the unit is pixels per second. A non-positive
% time step is returned as NaN rather than an infinite speed.

position = double(position);
time = double(time(:));
if size(position, 2) ~= 2 || size(position, 1) ~= numel(time)
    error('focal_animal_speed:size', ...
        'position must be n x 2 and time must have one entry per position.');
end
if any(~isfinite(position(:))) || any(~isfinite(time))
    error('focal_animal_speed:finite', ...
        'position and time must be finite. Drop missing clicks before this call.');
end
if numel(time) < 2
    speed = zeros(0, 1);
    speed_time = zeros(0, 1);
    return
end

dt = diff(time);
displacement = sqrt(sum(diff(position, 1, 1).^2, 2));
speed = displacement ./ dt;
speed(~isfinite(speed) | dt <= 0) = NaN;
speed_time = time(2:end);
end
