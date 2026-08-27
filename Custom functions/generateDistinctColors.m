function rgbColors = generateDistinctColors(n)
%GENERATEDISTINCTCOLORS Generate n distinct RGB colors avoiding grayscale and black
% Output: rgbColors is an n×3 matrix with RGB rows in [0,1]

    if n < 1
        rgbColors = [];
        return;
    end

    % Generate colors by spacing hues evenly
    hues = linspace(0, 1, n + 1);  % add 1 to avoid wrapping around to the first hue
    hues(end) = [];               % remove the extra hue

    saturation = 0.8;             % high saturation to avoid gray
    value = 0.9;                  % high value to avoid black/dark colors

    rgbColors = hsv2rgb([hues', repmat(saturation, n, 1), repmat(value, n, 1)]);
end