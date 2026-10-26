function test_session_partner_mean_psth()
%TEST_SESSION_PARTNER_MEAN_PSTH  Day-by-partner means without experiment files.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

test_means_ignore_outliers_partner3_and_other_days();
test_length_limit_drops_short_bouts();
test_limit_keeps_values_equal_to_the_cutoff();
test_bout_end_is_excluded();
test_empty_input();
disp('test_session_partner_mean_psth passed');
end

function test_means_ignore_outliers_partner3_and_other_days()
time = [-1 0 1];
power = [
    1 1 1
    3 3 9
    10 0 0
    2 2 2
    4 4 4
    7 7 7
    ];
partner = [1; 1; 1; 2; 1; 3];
session_id = {'dayA'; 'dayA'; 'dayA'; 'dayA'; 'dayB'; 'dayA'};
bout_length = [2; 0.5; 2; 2; 2; 2];

obs = session_partner_mean_psth(power, partner, session_id, bout_length, time, 4, 0);

assert(isequal(obs.session_id, {'dayA'; 'dayA'; 'dayB'}), 'Sessions were pooled or reordered.');
assert(isequal(obs.partner, [1; 2; 1]), 'Partner 3 entered the mean, or partner order changed.');
assert_close(obs.power(1, :), [2, 4/3, 0.5]);
assert_close(obs.n_bouts(1, :), [2 3 2]);
assert_close(obs.power(2, :), [2 2 2]);
assert_close(obs.power(3, :), [4 4 4]);
assert(isequal(obs.time, time), 'The time axis was not kept.');
end

function test_length_limit_drops_short_bouts()
time = [-1 0 1];
power = [
    1 1 1
    3 5 9
    10 0 0
    ];
partner = [1; 1; 1];
session_id = {'dayA'; 'dayA'; 'dayA'};
bout_length = [2; 0.5; 2];

obs = session_partner_mean_psth(power, partner, session_id, bout_length, time, 4, 1);

assert(isequal(obs.session_id, {'dayA'}), 'The short-bout session was split.');
assert_close(obs.power, [1 0.5 0.5]);
assert_close(obs.n_bouts, [1 2 2]);
end

function test_limit_keeps_values_equal_to_the_cutoff()
obs = session_partner_mean_psth([4 0 0], 1, {'day'}, 2, [-1 0 1], 4, 0);
assert(obs.power(1) == 4, 'A sample at the |z| limit was dropped.');
assert(obs.n_bouts(1) == 1, 'The sample at the |z| limit was not counted.');

obs = session_partner_mean_psth([4.1 0 0; 1 0 0], [1; 1], {'day'; 'day'}, [2; 2], [-1 0 1], 4, 0);
assert(abs(obs.power(1) - 1) < 1e-12, 'A sample above the |z| limit stayed in the mean.');
assert(obs.n_bouts(1) == 1, 'The outlier was still counted.');
end

function test_bout_end_is_excluded()
obs = session_partner_mean_psth([2 8], 1, {'day'}, 1, [0 1], 4, 0);
assert(obs.power(1) == 2, 'A sample before the bout end was dropped.');
assert(isnan(obs.power(2)), 'The sample at the bout end was kept.');
assert(obs.n_bouts(2) == 0, 'The sample at the bout end was counted.');
end

function test_empty_input()
obs = session_partner_mean_psth(zeros(0, 3), zeros(0, 1), {}, zeros(0, 1), [1 2 3], 4, 0);
assert(isempty(obs.power), 'An empty bout list produced a mean.');
assert(isempty(obs.session_id), 'An empty bout list produced a session.');
assert(isequal(obs.time, [1 2 3]), 'The time axis was dropped for an empty bout list.');
end

function assert_close(actual, expected)
if ~isequal(size(actual), size(expected))
    error('test_session_partner_mean_psth:Size', 'Sizes differ.');
end
if any(abs(actual - expected) > 1e-12)
    error('test_session_partner_mean_psth:Value', 'Values differ.');
end
end
