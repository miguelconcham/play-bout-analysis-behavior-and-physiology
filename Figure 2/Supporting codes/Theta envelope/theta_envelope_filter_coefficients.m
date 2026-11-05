function coefficients = theta_envelope_filter_coefficients(options)
%THETA_ENVELOPE_FILTER_COEFFICIENTS  Even-order Hamming FIR for the theta band.
%
% The order follows the 1.5 Hz transition at the envelope sampling rate, so
% the passband stays clear of the confirmatory 2-3 Hz modulation band.
% designfilt matches the other Hamming filters in this project. fir1, then a
% windowed-sinc difference of lowpasses, are used when designfilt is absent.

if nargin < 1 || isempty(options)
    options = struct();
end
options = theta_envelope_delta_options(options);
fs = options.envelope_rate;
f1 = options.theta_freq(1);
f2 = options.theta_freq(2);
order = options.theta_filter_order;
if exist('designfilt', 'file') == 2 || exist('designfilt', 'builtin') == 5
    designed = designfilt('bandpassfir', ...
        'FilterOrder', order, ...
        'CutoffFrequency1', f1, ...
        'CutoffFrequency2', f2, ...
        'SampleRate', fs, ...
        'DesignMethod', 'window', ...
        'Window', 'hamming');
    coefficients = designed.Coefficients(:);
    return
end
if exist('fir1', 'file') == 2
    coefficients = fir1(order, [f1, f2] / (fs / 2), ...
        theta_envelope_hamming_window(order + 1));
    coefficients = coefficients(:);
    return
end
coefficients = windowed_bandpass(order, f1, f2, fs);
end

function coefficients = windowed_bandpass(order, f1, f2, fs)
low = windowed_lowpass(order, f1, fs);
high = windowed_lowpass(order, f2, fs);
coefficients = high - low;
center = (f1 + f2) / 2;
omega = 2 * pi * center / fs;
samples = (0:numel(coefficients) - 1).';
response = abs(sum(coefficients .* exp(-1i * omega * samples)));
if response > 0
    coefficients = coefficients / response;
end
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
