function plot_time = psth_plot_time(time, n_columns, hist_range, bin_size)
%PSTH_PLOT_TIME  Time axis with one sample per PSTH column.
%
%   When numel(time) equals n_columns, plot_time is that vector. On the
%   paper range this is the 4001-bin axis from hist_range, unchanged.
%   A matrix one sample shorter uses the leading bins of that axis (the
%   4000-bin preallocation). A different grid, such as 801 bins across the
%   same histogram range, is spaced from hist_range(1) to hist_range(2).

if nargin < 4 || isempty(bin_size) || ~isscalar(bin_size) || bin_size <= 0
    error('psth_plot_time:Bin', 'bin_size must be a positive step in seconds.');
end
if nargin < 3 || numel(hist_range) ~= 2
    error('psth_plot_time:Range', 'hist_range must be [start end] in seconds.');
end
if ~isscalar(n_columns) || n_columns < 1 || n_columns ~= floor(n_columns)
    error('psth_plot_time:Width', 'n_columns must be a positive integer.');
end

time = double(time(:))';
if numel(time) == n_columns
    plot_time = time;
    return
end

if n_columns > 1
    inclusive_step = (hist_range(2) - hist_range(1)) / (n_columns - 1);
    if abs(inclusive_step - bin_size) > max(1e-6, 0.2 * bin_size)
        plot_time = hist_range(1) + (0:n_columns - 1) * inclusive_step;
        return
    end
end

if numel(time) > n_columns
    plot_time = time(1:n_columns);
    return
end

plot_time = hist_range(1) + (0:n_columns - 1) * bin_size;
end
