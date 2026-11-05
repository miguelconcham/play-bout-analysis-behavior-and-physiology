function envelope = theta_envelope_signal(signal, coefficients)
%THETA_ENVELOPE_SIGNAL  Z-scored Hilbert amplitude of a theta-filtered trace.
%
% `signal` is already at the rate the coefficients were designed for.
% Samples within one FIR order of each end are left NaN. The finite center
% is z-scored so sessions share a scale. The mean of a z-scored envelope is
% zero; slow 0.1 Hz power is kept because the modulation scan starts there.

signal = double(signal(:));
coefficients = coefficients(:);
if any(~isfinite(signal))
    error('theta_envelope_signal:input', ...
        'The theta filter input contains non-finite samples.');
end
order = numel(coefficients) - 1;
envelope = nan(size(signal));
if numel(signal) <= 3 * order
    return
end
filtered = filtfilt(coefficients.', 1, signal);
envelope = abs(hilbert(filtered));
envelope(1:order) = NaN;
envelope(end - order + 1:end) = NaN;
finite = isfinite(envelope);
center = envelope(finite);
scale = std(center);
if isempty(center) || ~isfinite(scale) || scale == 0
    envelope(:) = NaN;
    return
end
envelope(finite) = (center - mean(center)) / scale;
end
