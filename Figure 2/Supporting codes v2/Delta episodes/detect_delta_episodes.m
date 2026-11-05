function detected = detect_delta_episodes(lfp, sampling_rate, options)
%DETECT_DELTA_EPISODES  Delta-filter one continuous LFP interval and detect.
%
% lfp is the minimally processed interval (double, constant removed by
% the caller). The filter is the existing 1-5 Hz Hamming FIR, applied
% with filtfilt. Instantaneous amplitude is abs(hilbert(filtered)).
% Episodes are then built by characterize_delta_episodes.
%
% Cycle shape is measured by characterize_cycle_features on a low-pass
% of this same LFP. The 1-5 Hz signal is used only to place zero-crossings.
% Peak-to-peak episode fields are left as the amplitude detector defined them.

if nargin < 3 || isempty(options)
    options = struct();
end
options = delta_episode_options(options);
if isempty(options.filter_coefficients)
    options.filter_coefficients = delta_bandpass_coefficients(options);
end
if isempty(options.cycle_lowpass_coefficients)
    options.cycle_lowpass_coefficients = delta_lowpass_coefficients(options);
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

lowpass_order = numel(options.cycle_lowpass_coefficients) - 1;
if numel(lfp) <= 3 * lowpass_order
    error('detect_delta_episodes:length', ...
        ['Interval has %d samples. filtfilt needs more than three times ' ...
         'the %d-point cycle low-pass.'], numel(lfp), lowpass_order);
end
broad = filtfilt(options.cycle_lowpass_coefficients, 1, lfp);
[detected.episodes, detected.cycles] = characterize_cycle_features( ...
    broad, filtered, sampling_rate, detected.episodes, options);
detected.lfp = lfp;
detected.filtered = filtered;
detected.amplitude = amplitude;
end
