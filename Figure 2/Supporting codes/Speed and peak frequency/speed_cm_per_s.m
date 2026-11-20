function cm_per_s = speed_cm_per_s(pixel_displacement, frame_dt)
%SPEED_CM_PER_S  Centimetres per second from a pixel step and its duration.
%
% pixel_displacement is the distance between two tracked positions, in
% image pixels. frame_dt is the NPX time between those positions, in
% seconds. cm/s = (pixels / seconds) * 0.1334. A non-positive time step
% is NaN, matching focal_animal_speed.

displacement = double(pixel_displacement);
frame_dt = double(frame_dt);
cm_per_s = displacement ./ frame_dt * pixel_to_cm();
cm_per_s(~isfinite(frame_dt) | frame_dt <= 0) = NaN;
end
