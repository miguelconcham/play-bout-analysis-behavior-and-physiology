function band_power = theta_envelope_band_mean(power, frequencies, grid)
%THETA_ENVELOPE_BAND_MEAN  Mean of a spectrum inside each grid window.

power = power(:);
frequencies = frequencies(:);
n_centers = numel(grid.centers);
band_power = nan(1, n_centers);
if isempty(power) || isempty(frequencies) || numel(power) ~= numel(frequencies)
    return
end
for k = 1:n_centers
    in_window = frequencies >= grid.windows(k, 1) & ...
        frequencies <= grid.windows(k, 2);
    if any(in_window)
        band_power(k) = mean(power(in_window));
    end
end
end
