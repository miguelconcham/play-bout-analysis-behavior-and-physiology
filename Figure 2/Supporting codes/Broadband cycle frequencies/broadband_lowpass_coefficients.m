function coefficients = broadband_lowpass_coefficients(options)
%BROADBAND_LOWPASS_COEFFICIENTS  Hamming FIR for cycle detection.

if nargin < 1 || isempty(options)
    options = struct();
end
options = broadband_cycle_options(options);
order = options.lowpass_order;
if exist('designfilt', 'file') == 2 || exist('designfilt', 'builtin') == 5
    lowpass = designfilt('lowpassfir', ...
        'FilterOrder', order, ...
        'CutoffFrequency', options.lowpass_hz, ...
        'SampleRate', options.sampling_rate, ...
        'DesignMethod', 'window', ...
        'Window', 'hamming');
    coefficients = lowpass.Coefficients;
    return
end
if exist('fir1', 'file') == 2
    coefficients = fir1(order, ...
        options.lowpass_hz / (options.sampling_rate / 2), ...
        'low', hamming(order + 1));
    return
end
error('broadband_lowpass_coefficients:filter', ...
    'designfilt or fir1 is required for the cycle low-pass.');
end
