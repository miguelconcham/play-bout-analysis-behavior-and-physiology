function test_playbout_erp()
%TEST_PLAYBOUT_ERP  Epoch cut, saturated-trial rejection, and session ERPs.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

test_defaults_cover_the_playbout_window();
test_epoch_uses_the_entire_window_and_the_recording_clock();
test_ties_use_the_earlier_sample();
test_epochs_outside_the_recording_are_missing();
test_saturation_flags_both_int16_rails();
test_raw_mean_keeps_every_finite_trial();
test_clean_mean_drops_saturated_trials_and_subtracts_baseline();
test_session_pool_keeps_one_electrode_and_counts_rejections();
test_grand_average_matches_the_balanced_mean();
disp('test_playbout_erp passed');
end

function test_defaults_cover_the_playbout_window()
options = playbout_erp_defaults();
assert(isequal(options.window_s, [-2 2]), 'The default epoch is not -2 s to +2 s.');
assert(options.sampling_rate == 2500, 'The LFP sampling rate is not 2500 Hz.');
assert(options.saturation_limit == 32767, 'The int16 rail is not the saturation limit.');
assert(isequal(options.baseline_window, [-2 0]), 'The baseline is not [-2, 0).');
assert(strcmp(options.formula, 'Voltage ~ 1 + (1|Animal)'), ...
    'The grand-average formula changed.');
end

function test_epoch_uses_the_entire_window_and_the_recording_clock()
lfp = 1:20000;
onset_s = 4;
[traces, time_s, complete] = cut_event_lfp(lfp, onset_s, [-2 2], 2500);

assert(complete, 'A window inside the recording was marked incomplete.');
assert(numel(time_s) == 10001, 'The -2 s to +2 s epoch is not the full sample grid.');
assert(time_s(1) == -2 && time_s(end) == 2, 'The epoch does not span -2 s to +2 s.');
assert(any(time_s == 0), 'The onset sample is missing.');
assert(traces(time_s == 0) == 10000, ...
    'The onset sample is not the recording sample at time onset.');
assert(traces(1) == 10000 - 5000, 'The first sample is not 2 s before onset.');
assert(traces(end) == 10000 + 5000, 'The last sample is not 2 s after onset.');
end

function test_ties_use_the_earlier_sample()
lfp = 1:30;
[earlier, time_s] = cut_event_lfp(lfp, 1.05, [-0.2 0.2], 10);
[later, ~] = cut_event_lfp(lfp, 1.06, [-0.2 0.2], 10);
assert(earlier(time_s == 0) == 10, 'An exact sample tie did not use the earlier sample.');
assert(later(time_s == 0) == 11, 'A time just after a tie did not move to the next sample.');
end

function test_epochs_outside_the_recording_are_missing()
lfp = 1:20000;
[traces, ~, complete] = cut_event_lfp(lfp, 0.5, [-2 2], 2500);
assert(~complete, 'An epoch hanging off the start was kept.');
assert(all(isnan(traces)), 'An incomplete epoch was not NaN.');
end

function test_saturation_flags_both_int16_rails()
traces = [
    1 2 3
    10 32767 4
    -32768 0 1
    NaN NaN NaN
    32766 -32766 0
    ];
mask = saturated_lfp_trials(traces, 32767);
assert(isequal(mask, [false; true; true; false; false]), ...
    'Saturated trials were not limited to the int16 rails.');
end

function test_raw_mean_keeps_every_finite_trial()
erp = mean_event_erp([1 1; 3 5; NaN NaN]);
assert(isequal(erp, [2 3]), 'The raw session ERP dropped a finite trial or kept a NaN trial.');
empty_erp = mean_event_erp(nan(2, 4));
assert(isequal(size(empty_erp), [1 4]) && all(isnan(empty_erp)), ...
    'A recording with no finite trial did not return a NaN ERP.');
end

function test_clean_mean_drops_saturated_trials_and_subtracts_baseline()
time_s = [-1 -0.5 0 0.5 1];
traces = [
    2 4 8 8 8
    1 1 32767 1 1
    NaN NaN NaN NaN NaN
    ];
[erp, info] = clean_session_erp(traces, time_s, 32767, [-1 0]);
assert(isequal(info.n_complete, 2), 'The incomplete trial was counted as complete.');
assert(isequal(info.n_saturated, 1), 'The rail trial was kept.');
assert(isequal(info.n_kept, 1), 'The clean trial was not the only trial in the ERP.');
assert_close(erp, [-1 1 5 5 5]);
end

function test_session_pool_keeps_one_electrode_and_counts_rejections()
time_s = [-1 -0.5 0 0.5 1];
good = [2 4 8 8 8];
saturated = [1 1 32767 1 1];
other_electrode = good + 100;

first = new_playbout_erp_session();
first.animal_id = 'B1D1';
first.session_id = 'B1D1 1013 Dual';
first.electrode_index = 1;
first.time_s = time_s;
first.lfp_trials = [good; saturated];

second = first;
second.animal_id = 'B2S2';
second.session_id = 'B2S2 1014 Dual';
second.lfp_trials = good + 2;

third = first;
third.session_id = 'B1D1 1013 Dual electrode 2';
third.electrode_index = 2;
third.lfp_trials = other_electrode;

options = playbout_erp_defaults();
options.baseline_window = [-1 0];
assembled = session_erps_for_lme([first second third], options);

assert(assembled.n_recordings_used == 2, 'The second electrode entered the grand average.');
assert(isequal(assembled.animal_id, {'B1D1'; 'B2S2'}), 'Recording order or identity changed.');
assert(assembled.n_trials_complete == 3, 'Complete trials were miscounted.');
assert(assembled.n_trials_saturated == 1, 'The saturated trial was not counted.');
assert(assembled.n_trials_kept == 2, 'A rejected or extra-electrode trial stayed in.');
assert_close(assembled.session_erp(1, :), [-1 1 5 5 5]);
assert_close(assembled.session_erp(2, :), [-1 1 5 5 5]);
assert(assembled.lfp_abs_max >= 32767, 'The rail sample was ignored while scanning amplitude.');
end

function test_grand_average_matches_the_balanced_mean()
if exist('fitlme', 'file') ~= 2
    disp('Skipping the mixed-model grand average because fitlme is not on the path.');
    return
end
session_erp = [
    1 1 1
    3 3 3
    1 1 1
    3 3 3
    ];
animal_id = {'B1D1'; 'B1D1'; 'B2S2'; 'B2S2'};
results = grand_average_erp_lme(session_erp, animal_id, [-1 0 1], 0.05);
assert(isequal(results.formula, 'Voltage ~ 1 + (1|Animal)'), ...
    'The fitted formula is not the grand-mean model.');
assert_close(results.mean.', [2 2 2]);
assert(all(results.ci(:, 1) < 2 & results.ci(:, 2) > 2), ...
    'The confidence interval does not contain the balanced grand mean.');
assert(results.n_sessions == 4 && results.n_animals == 2, ...
    'Session or animal counts were not kept.');
end

function assert_close(actual, expected)
if numel(actual) ~= numel(expected) || any(abs(actual(:) - expected(:)) > 1e-10)
    error('test_playbout_erp:assert', ...
        'Expected [%s], got [%s].', num2str(expected(:).'), num2str(actual(:).'));
end
end
