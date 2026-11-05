function test_psth_plot_time()
%TEST_PSTH_PLOT_TIME  Paper 4001-bin axis and the 801-bin grid.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

test_paper_axis_is_unchanged();
test_preallocated_4000_uses_the_leading_paper_bins();
test_801_bins_span_the_histogram_range();
disp('test_psth_plot_time passed');
end

function test_paper_axis_is_unchanged()
hist_range = [-20 20];
bin_size = 1 - 0.99;
time = hist_range(1):bin_size:hist_range(2) + bin_size;
assert(numel(time) == 4001, 'The paper time axis was not 4001 bins.');
plot_time = psth_plot_time(time, 4001, hist_range, bin_size);
assert(isequal(plot_time, time), 'The 4001-bin paper axis was replaced.');
end

function test_preallocated_4000_uses_the_leading_paper_bins()
hist_range = [-20 20];
bin_size = 1 - 0.99;
time = hist_range(1):bin_size:hist_range(2) + bin_size;
plot_time = psth_plot_time(time, 4000, hist_range, bin_size);
assert(isequal(plot_time, time(1:4000)), 'The 4000 stored bins left the paper clock.');
end

function test_801_bins_span_the_histogram_range()
hist_range = [-20 20];
bin_size = 1 - 0.99;
time = hist_range(1):bin_size:hist_range(2) + bin_size;
plot_time = psth_plot_time(time, 801, hist_range, bin_size);
expected = hist_range(1):0.05:hist_range(2);
assert(numel(plot_time) == 801, 'The 801-bin axis changed length.');
assert(max(abs(plot_time - expected)) < 1e-9, '801 bins were not -20:0.05:20.');
end
