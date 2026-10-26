function test_match_psth_columns()
%TEST_MATCH_PSTH_COLUMNS  Widths used when pooling behavior-onset PSTHs.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

test_equal_width_is_kept();
test_extra_endpoint_is_dropped();
test_one_missing_endpoint_is_nan();
test_fifty_millisecond_grid_is_not_cut_to_4000();
test_different_grids_are_reported();
disp('test_match_psth_columns passed');
end

function test_equal_width_is_kept()
psth = reshape(1:10, 2, 5);
assert(isequal(match_psth_columns(psth, 5), psth), 'An equal width was changed.');
end

function test_extra_endpoint_is_dropped()
psth = [(1:4000); (1:4000)];
psth(:, 4001) = 9;
aligned = match_psth_columns(psth, 4000);
assert(isequal(size(aligned), [2 4000]), 'The 10 ms grid was not trimmed to 4000.');
assert(isequal(aligned(1, :), 1:4000), 'Columns other than the endpoint were changed.');
end

function test_one_missing_endpoint_is_nan()
psth = [1 2 3; 4 5 6];
aligned = match_psth_columns(psth, 4);
assert(isequal(size(aligned), [2 4]), 'The missing endpoint was not appended.');
assert(isequal(aligned(:, 1:3), psth), 'Existing samples were changed.');
assert(all(isnan(aligned(:, 4))), 'The appended endpoint was not NaN.');
end

function test_fifty_millisecond_grid_is_not_cut_to_4000()
% -20:0.05:20 has 801 samples. Column 4000 does not exist on that grid.
psth = ones(3, 801);
aligned = match_psth_columns(psth, 801);
assert(isequal(size(aligned), [3 801]), 'The 50 ms grid was not kept intact.');
end

function test_different_grids_are_reported()
threw = false;
try
    match_psth_columns(ones(2, 801), 4000);
catch err
    threw = strcmp(err.identifier, 'match_psth_columns:Width') ...
        || ~isempty(strfind(err.message, '801'));
end
assert(threw, 'An 801-bin PSTH was forced onto the 4000-bin grid.');
end
