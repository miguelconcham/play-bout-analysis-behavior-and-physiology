function peak = select_band_peak(peak_params, band)
%SELECT_BAND_PEAK  Strongest FOOOF peak whose center lies inside a band.
%
% peak_params is the FOOOF matrix [center frequency, peak power, bandwidth].
% The selected peak is the one with the largest peak power. Equal powers
% keep the lower center frequency. An empty band returns present = false.

if nargin < 2 || numel(band) ~= 2 || ~all(isfinite(band)) || band(2) < band(1)
    error('select_band_peak:band', 'band must be [low high] in Hz.');
end

peak = struct( ...
    'present', false, ...
    'center_frequency', NaN, ...
    'power', NaN, ...
    'bandwidth', NaN, ...
    'n_peaks_in_band', 0);

if isempty(peak_params)
    return
end
if size(peak_params, 2) ~= 3
    error('select_band_peak:params', ...
        'peak_params must be [center frequency, power, bandwidth].');
end

in_band = peak_params(:, 1) >= band(1) & peak_params(:, 1) <= band(2) & ...
    isfinite(peak_params(:, 1)) & isfinite(peak_params(:, 2));
peak.n_peaks_in_band = nnz(in_band);
if peak.n_peaks_in_band == 0
    return
end

candidates = peak_params(in_band, :);
[~, order] = sortrows(candidates, [-2 1]);
chosen = candidates(order(1), :);
peak.present = true;
peak.center_frequency = chosen(1);
peak.power = chosen(2);
peak.bandwidth = chosen(3);
end
