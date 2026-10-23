function coefficients = delta_lowpass_coefficients(options)
%DELTA_LOWPASS_COEFFICIENTS  Low-pass used to locate cycle extrema.
%
% Cole & Voytek (2019) low-pass at about four times the oscillation
% frequency so the waveform survives, and they give the simulation
% low-pass a duration of 100 ms. cycle_lowpass_hz and cycle_lowpass_order
% carry that choice (20 Hz and order 250 at 2500 Hz for a 1-5 Hz band).
% MATLAB uses designfilt. Without it, fir1 builds the same Hamming window.

if nargin < 1 || isempty(options)
    options = struct();
end
options = delta_episode_options(options);
order = options.cycle_lowpass_order;
if exist('designfilt', 'file') == 2 || exist('designfilt', 'builtin') == 5
    lowpass = designfilt('lowpassfir', ...
        'FilterOrder', order, ...
        'CutoffFrequency', options.cycle_lowpass_hz, ...
        'SampleRate', options.sampling_rate, ...
        'DesignMethod', 'window', ...
        'Window', 'hamming');
    coefficients = lowpass.Coefficients;
    return
end
if exist('fir1', 'file') == 2
    coefficients = fir1(order, ...
        options.cycle_lowpass_hz / (options.sampling_rate / 2), ...
        'low', hamming(order + 1));
    return
end
error('delta_lowpass_coefficients:filter', ...
    'designfilt or fir1 is required for the cycle low-pass.');
end
