function cm_per_s = frame_speed_to_cm_s(pixels_per_frame)
%FRAME_SPEED_TO_CM_S  Convert frame-to-frame displacement to cm/s.
%
% pixels_per_frame is the displacement between consecutive video frames.
% The tracking time base is tracking_video_fps frames per second.

cm_per_s = double(pixels_per_frame) * tracking_video_fps() * pixel_to_cm();
end
