function detected = detect_delta_episodes(lfp, sampling_rate, options)
%DETECT_DELTA_EPISODES  Delta-filter one continuous LFP interval and detect.
%
% lfp is the minimally processed interval (double, constant removed by
% the caller). The filter is the existing 1-5 Hz Hamming FIR, applied
% with filtfilt. Instantaneous amplitude is abs(hilbert(filtered)).
% Episodes are then built by characterize_delta_episodes.

if nargin < 3 || isempty(options)
    options = struct();
end
options = delta_episode_options(options);
if isempty(options.filter_coefficients)
    options.filter_coefficients = delta_bandpass_coefficients(options);
end

lfp = double(lfp(:).');
if any(~isfinite(lfp))
    error('detect_delta_episodes:lfp', 'LFP contains non-finite samples.');
end
filter_order = numel(options.filter_coefficients) - 1;
if numel(lfp) <= 3 * filter_order
    error('detect_delta_episodes:length', ...
        ['Interval has %d samples. filtfilt needs more than three times ' ...
         'the %d-point delta filter.'], numel(lfp), filter_order);
end

filtered = filtfilt(options.filter_coefficients, 1, lfp);
amplitude = abs(hilbert(filtered));
detected = characterize_delta_episodes(filtered, amplitude, sampling_rate, options);
detected.lfp = lfp;
detected.filtered = filtered;
detected.amplitude = amplitude;
end
