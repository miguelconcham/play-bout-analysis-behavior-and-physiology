function delta = fooof_band_power(frequency, power, band, options)
%FOOOF_BAND_POWER  Mean aperiodic-removed log power in one frequency band.
%
% frequency is in Hz and power is linear. FOOOF is fit on options.fooof_fit_range
% (default 0.5-40 Hz). delta is the mean of the flattened spectrum, log10
% power minus the aperiodic fit, inside band. A fit that does not run
% returns NaN. The periodic peaks are left in the flattened spectrum, so a
% delta bump raises this value above the 1/f background.

if nargin < 4 || isempty(options)
    options = struct();
end
if ~isfield(options, 'fooof_max_peaks') || isempty(options.fooof_max_peaks)
    options.fooof_max_peaks = 6;
end
if ~isfield(options, 'fooof_fit_range') || isempty(options.fooof_fit_range)
    options.fooof_fit_range = [0.5 40];
end
delta = NaN;
[frequency, power] = even_positive_spectrum(frequency, power, options.fooof_fit_range);
if numel(frequency) < 10
    return
end
warn_state = warning;
cleanup = onCleanup(@() warning(warn_state)); %#ok<NASGU>
warning('off', 'all');
try
    fitted = fooof_matlab(frequency, power, ...
        'MaxNPeaks', options.fooof_max_peaks, ...
        'PeakWidthLimits', [0.5 12], ...
        'AperiodicMode', 'fixed');
catch
    return
end
if ~isfield(fitted, 'flat_spectrum') || numel(fitted.flat_spectrum) ~= numel(fitted.f)
    return
end
in_band = fitted.f >= band(1) & fitted.f <= band(2);
if ~any(in_band)
    return
end
delta = mean(fitted.flat_spectrum(in_band));
end


function [frequency, power] = even_positive_spectrum(frequency, power, fit_range)
frequency = double(frequency(:));
power = double(power(:));
ok = isfinite(frequency) & isfinite(power) & frequency > 0 & power > 0;
frequency = frequency(ok);
power = power(ok);
if numel(frequency) < 2
    frequency = [];
    power = [];
    return
end
[frequency, order] = sort(frequency);
power = power(order);
[frequency, unique_index] = unique(frequency);
power = power(unique_index);
keep = frequency >= fit_range(1) & frequency <= fit_range(2);
frequency = frequency(keep);
power = power(keep);
if numel(frequency) < 10
    frequency = [];
    power = [];
    return
end
spacing = diff(frequency);
step = median(spacing);
if any(abs(spacing - step) > 1e-6 * max(1, abs(step)))
    even_frequency = (frequency(1):step:frequency(end))';
    if numel(even_frequency) < 10
        frequency = [];
        power = [];
        return
    end
    power = interp1(frequency, power, even_frequency, 'linear');
    frequency = even_frequency;
end
end
