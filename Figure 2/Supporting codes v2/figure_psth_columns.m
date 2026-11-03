function aligned = figure_psth_columns(psth, n_axis)
%FIGURE_PSTH_COLUMNS  Columns of a PSTH on the Figure 2 time axis.
%
%   n_axis is numel(time) in Figure 2. For the paper range, hist_range
%   [-20 20] and a 10 ms step from wind_length - wind_overlap, that vector
%   has 4001 bins. A PSTH at least that wide keeps exactly those bins, so
%   column 4001 stays. The old pool took columns 1:4000 and dropped it.
%   A matrix one bin shorter (the 4000-bin preallocation) is padded with one
%   NaN so it still sits on that axis. A shorter grid, such as 801 bins, is
%   returned unchanged.

if nargin < 2 || ~isscalar(n_axis) || n_axis < 1 || n_axis ~= floor(n_axis)
    error('figure_psth_columns:Width', ...
        'n_axis must be a positive integer time-bin count.');
end

n_have = size(psth, 2);
if n_have >= n_axis
    aligned = psth(:, 1:n_axis);
    return
end
if n_have + 1 == n_axis
    aligned = psth;
    aligned(:, n_axis) = NaN;
    return
end

aligned = psth;
end
