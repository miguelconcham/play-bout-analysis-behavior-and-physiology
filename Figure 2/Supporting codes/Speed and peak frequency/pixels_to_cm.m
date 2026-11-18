function centimetres = pixels_to_cm(pixels)
%PIXELS_TO_CM  Convert an image distance from pixels to centimetres.

centimetres = double(pixels) * pixel_to_cm();
end
