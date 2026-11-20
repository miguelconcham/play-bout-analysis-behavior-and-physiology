function settings = stored_band_settings(entry, fallback)
%STORED_BAND_SETTINGS  Frequency grid already saved with a speed mat.
%
% Play minus non-play is computed from that power. Rebuilding SpeedCm
% must keep f, freq_range, and the spectrogram window, or the points
% above and below zero change for a reason other than the centimetre scale.

if nargin < 2 || isempty(fallback)
    fallback = struct( ...
        'freq_range', [1 5], ...
        'f', 1:.1:16, ...
        'wind_length', 0.250, ...
        'wind_overlap', 0.249);
end
settings = fallback;
if isempty(entry)
    settings.reused = false;
    return
end
settings.reused = false;
if isfield(entry, 'freq_range') && numel(entry.freq_range) == 2 ...
        && all(isfinite(entry.freq_range))
    settings.freq_range = double(entry.freq_range(:)).';
    settings.reused = true;
end
if isfield(entry, 'f') && ~isempty(entry.f) && all(isfinite(entry.f(:)))
    settings.f = double(entry.f(:)).';
    settings.reused = true;
end
if isfield(entry, 'wind_length') && isscalar(entry.wind_length) ...
        && isfinite(entry.wind_length) && entry.wind_length > 0
    settings.wind_length = double(entry.wind_length);
end
if isfield(entry, 'wind_overlap') && isscalar(entry.wind_overlap) ...
        && isfinite(entry.wind_overlap) && entry.wind_overlap >= 0
    settings.wind_overlap = double(entry.wind_overlap);
end
end
