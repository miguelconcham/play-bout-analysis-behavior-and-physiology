function frames_per_second = tracking_video_fps()
%TRACKING_VIDEO_FPS  Video frame rate used to build the tracking time base.
%
% GENERATE_SPEED_FREQ_RELATION places a sample on every video frame at
% frames/30 before the video-to-NPX sync. Speed stored beside that sample
% is the pixel displacement from one frame to the next.

frames_per_second = 30;
end
