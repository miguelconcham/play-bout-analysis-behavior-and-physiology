function grid = theta_envelope_frequency_grid(band, step, bandwidth, confirmatory_band)
%THETA_ENVELOPE_FREQUENCY_GRID  Centers and integration windows.
%
% Centers run from band(1) to band(2) in steps of `step`, and the endpoint
% is kept when the step does not land on it. The integration window is
% bandwidth wide and is clipped to `band`. confirmatory_band marks windows
% that overlap that range; pass [] to leave the mask empty.

if numel(band) ~= 2 || band(1) >= band(2)
    error('theta_envelope_frequency_grid:band', 'band must be [low high].');
end
if step <= 0 || bandwidth <= 0
    error('theta_envelope_frequency_grid:step', ...
        'step and bandwidth must be positive.');
end

start_milli = round(band(1) * 1000);
stop_milli = round(band(2) * 1000);
step_milli = round(step * 1000);
if step_milli < 1
    error('theta_envelope_frequency_grid:step', ...
        'step is finer than the 0.001 Hz grid.');
end
centers_milli = start_milli:step_milli:stop_milli;
if isempty(centers_milli) || centers_milli(end) ~= stop_milli
    centers_milli = unique([centers_milli, stop_milli]);
end
centers = centers_milli / 1000;
half = bandwidth / 2;
windows = zeros(numel(centers), 2);
for k = 1:numel(centers)
    windows(k, 1) = max(band(1), centers(k) - half);
    windows(k, 2) = min(band(2), centers(k) + half);
end

grid = struct();
grid.centers = centers(:).';
grid.windows = windows;
grid.bandwidth_hz = bandwidth;
grid.step_hz = step;
grid.band_hz = band(:).';
grid.confirmatory = false(1, numel(centers));
if nargin >= 4 && ~isempty(confirmatory_band)
    overlaps = windows(:, 1) <= confirmatory_band(2) & ...
        windows(:, 2) >= confirmatory_band(1);
    grid.confirmatory = overlaps(:).';
end
end
