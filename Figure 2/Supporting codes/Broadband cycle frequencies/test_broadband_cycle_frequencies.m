function test_broadband_cycle_frequencies()
%TEST_BROADBAND_CYCLE_FREQUENCIES  Synthetic low-pass cycle checks.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
    addpath(fullfile(fileparts(this_dir), 'Delta episodes'));
end
if exist('pkg', 'file') == 2
    try
        pkg load signal
    catch
    end
end

options = broadband_cycle_options(struct('lowpass_order', 100));
sampling_rate = options.sampling_rate;
time_5 = (0:6 * sampling_rate - 1) / sampling_rate;
detected_5 = detect_broadband_cycles(sin(2 * pi * 5 * time_5), ...
    sampling_rate, options);
accepted_5 = detected_5.frequency_hz(detected_5.accepted);
assert(~isempty(accepted_5), 'A 5 Hz sine produced no accepted cycles.');
assert(all(abs(accepted_5 - 5) < 0.25), ...
    'Accepted 5 Hz cycles are farther than 0.25 Hz from 5 Hz.');

time_10 = (0:4 * sampling_rate - 1) / sampling_rate;
combined = [sin(2 * pi * 5 * time_5), sin(2 * pi * 10 * time_10)];
detected_both = detect_broadband_cycles(combined, sampling_rate, options);
accepted_both = detected_both.frequency_hz(detected_both.accepted);
assert(nnz(accepted_both > 4 & accepted_both < 6) >= 10, ...
    'The concatenated signal lost the 5 Hz cycles.');
assert(nnz(accepted_both > 9 & accepted_both < 11) >= 10, ...
    'The concatenated signal lost the 10 Hz cycles.');

disp('broadband_cycle_frequencies tests passed');
end
