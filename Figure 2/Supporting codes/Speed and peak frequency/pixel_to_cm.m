function cm_per_pixel = pixel_to_cm()
%PIXEL_TO_CM  Image-pixel scale for the play-tracking videos.
%
% Positions in the tracking files are image pixels. Multiply by this
% value to get centimetres. Frame-to-frame speed uses tracking_video_fps
% as well: cm/s = pixels per frame * frames per second * cm per pixel.

cm_per_pixel = 0.1334;
end
