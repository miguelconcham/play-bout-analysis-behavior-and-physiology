function coefficients = delta_bandpass_coefficients(options)
%DELTA_BANDPASS_COEFFICIENTS  Hamming FIR used for the 1-5 Hz delta band.
%
% Same design as Estimate_delta_peak_trough_csd, Figure 2, and the
% synchronous-spike LFP filter: bandpassfir, window method, Hamming.

if nargin < 1 || isempty(options)
    options = struct();
end
options = delta_episode_options(options);
if exist('designfilt', 'file') ~= 2 && exist('designfilt', 'builtin') ~= 5
    error('delta_bandpass_coefficients:designfilt', ...
        'designfilt is required for the 1-5 Hz delta filter.');
end

delta_filter = designfilt('bandpassfir', ...
    'FilterOrder', options.delta_filter_order, ...
    'CutoffFrequency1', options.delta_freq(1), ...
    'CutoffFrequency2', options.delta_freq(2), ...
    'SampleRate', options.sampling_rate, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');
coefficients = delta_filter.Coefficients;
end
