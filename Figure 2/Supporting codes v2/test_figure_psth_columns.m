function test_figure_psth_columns()
%TEST_FIGURE_PSTH_COLUMNS  Paper range (4001) and the shorter 801-bin grid.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

test_paper_range_keeps_all_4001_bins();
test_columns_past_the_paper_axis_are_dropped();
test_preallocated_4000_matches_the_paper_axis();
test_shorter_grid_is_not_indexed_at_4000();
disp('test_figure_psth_columns passed');
end

function test_paper_range_keeps_all_4001_bins()
psth = repmat(1:4001, 2, 1);
aligned = figure_psth_columns(psth, 4001);
assert(isequal(aligned, psth), 'The 4001-bin paper range was changed.');
% Samples the old 1:4000 slice kept are still the leading columns.
assert(isequal(aligned(:, 1:4000), psth(:, 1:4000)), ...
    'The leading 4000 paper bins changed.');
end

function test_columns_past_the_paper_axis_are_dropped()
psth = repmat(1:4002, 2, 1);
aligned = figure_psth_columns(psth, 4001);
assert(isequal(size(aligned), [2 4001]), 'The paper axis was not 4001 bins.');
assert(isequal(aligned, psth(:, 1:4001)), 'A bin inside the paper range was dropped.');
end

function test_preallocated_4000_matches_the_paper_axis()
psth = repmat(1:4000, 2, 1);
aligned = figure_psth_columns(psth, 4001);
assert(isequal(size(aligned), [2 4001]), 'The 4000-bin matrix was not placed on the 4001-bin axis.');
assert(isequal(aligned(:, 1:4000), psth), 'The original 4000 bins changed.');
assert(all(isnan(aligned(:, 4001))), 'The missing paper endpoint was not empty.');
end

function test_shorter_grid_is_not_indexed_at_4000()
psth = reshape(1:(3 * 801), 3, 801);
on_its_own_axis = figure_psth_columns(psth, 801);
against_paper_axis = figure_psth_columns(psth, 4001);
assert(isequal(on_its_own_axis, psth), 'An 801-bin PSTH was changed on its own axis.');
assert(isequal(against_paper_axis, psth), 'An 801-bin PSTH was indexed on the 4001-bin axis.');
end
