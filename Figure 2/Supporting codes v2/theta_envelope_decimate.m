function y = theta_envelope_decimate(signal, rate_in, rate_out)
%THETA_ENVELOPE_DECIMATE  Anti-aliased integer downsampling.
%
% Output sample k is the filtered input sample (k-1)*factor+1. The lowpass
% sits at 80 percent of the new Nyquist frequency.

signal = double(signal(:));
if rate_out == rate_in
    y = signal;
    return
end
factor = rate_in / rate_out;
if factor ~= round(factor) || factor < 2
    error('theta_envelope_decimate:rate', ...
        'rate_in must be an integer multiple of rate_out.');
end
factor = round(factor);
new_nyquist = rate_out / 2;
cutoff = 0.8 * new_nyquist;
transition = 0.2 * new_nyquist;
order = 2 * ceil((3.3 * rate_in / transition) / 2);
coefficients = windowed_lowpass(order, cutoff, rate_in);
if numel(signal) <= 3 * (numel(coefficients) - 1)
    error('theta_envelope_decimate:length', ...
        'The LFP is shorter than the anti-alias filter.');
end
filtered = filtfilt(coefficients(:).', 1, signal);
y = filtered(1:factor:end);
end

function coefficients = windowed_lowpass(order, cutoff, fs)
half = order / 2;
n = (-half:half).';
fc = cutoff / fs;
h = 2 * fc * normalized_sinc(2 * fc * n);
window = 0.54 + 0.46 * cos(2 * pi * n / order);
coefficients = h .* window;
coefficients = coefficients / sum(coefficients);
end

function y = normalized_sinc(x)
y = ones(size(x));
nonzero = x ~= 0;
y(nonzero) = sin(pi * x(nonzero)) ./ (pi * x(nonzero));
end
