function test_assign_partner_index()
%TEST_ASSIGN_PARTNER_INDEX  Partner windows and bout assignment.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

test_partners_session_rows_ignore_tickling();
test_onset_inside_a_window_keeps_a_bout_that_runs_past_the_end();
test_non_partner_type_falls_back_to_session_blocks();
test_original_window_rule();
disp('test_assign_partner_index passed');
end

function test_original_window_rule()
windows = [0 10; 10 30];
partner = partner_index_in_windows([1; 12; 1], [5; 20; 15], windows(:, 1), windows(:, 2));
assert(isequal(partner, [1; 2; 1]), ...
    'The original whole-bout partner index changed.');
end

function test_partners_session_rows_ignore_tickling()
behavior.Animal = {'Session_structure'; 'Session_structure'; 'Session_structure'; 'B1'};
behavior.Start = [100; 0; 300; 5];
behavior.End = [200; 90; 400; 10];
behavior.Type = {'Partners session'; 'Tickling'; 'Partners session'; 'Pounce'};

windows = partner_windows_from_behavior(behavior);

assert(isequal(size(windows), [2 2]), 'Tickling or the animal row was kept.');
assert(isequal(windows, [100 200; 300 400]), 'Partner windows were not ordered by start time.');
end

function test_onset_inside_a_window_keeps_a_bout_that_runs_past_the_end()
windows = [0 10; 10 20];
partner = assign_partner_index([1; 12; 9.5; 25], windows(:, 1), windows(:, 2));
assert(isequaln(partner, [1; 2; 1; NaN]), 'Bout onsets were assigned to the wrong partner.');
end

function test_non_partner_type_falls_back_to_session_blocks()
behavior.Animal = {'Session_structure'; 'Session_structure'; 'Session_structure'};
behavior.Start = [20; 0; 40];
behavior.End = [30; 10; 50];
behavior.Type = {'Play'; 'Intro'; 'Tickling'};

windows = partner_windows_from_behavior(behavior);

assert(isequal(windows, [0 10; 20 30]), 'The non-tickling session blocks were not used.');
end
