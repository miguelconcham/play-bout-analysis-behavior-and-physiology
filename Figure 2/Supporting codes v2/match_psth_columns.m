function aligned = match_psth_columns(psth, n_columns)
%MATCH_PSTH_COLUMNS  Align a behavior-onset PSTH to the play-bout time grid.
%
%   psth is nBout x nTime. n_columns is the width of play_bout_onset for the
%   same session. Equal widths are kept. A one-column difference is the
%   inclusive endpoint of hist_range(1):bin:hist_range(2) (4001 vs 4000 at
%   10 ms, 801 vs 800 at 50 ms): the extra column is dropped, or one NaN is
%   appended. Any larger gap is left untouched and reported, because those
%   columns are not the same time axis.

if nargin < 2 || ~isscalar(n_columns) || n_columns < 1 || n_columns ~= floor(n_columns)
    error('match_psth_columns:Width', ...
        'n_columns must be a positive integer time-bin count.');
end

n_have = size(psth, 2);
if n_have == n_columns
    aligned = psth;
    return
end
if n_have == n_columns + 1
    aligned = psth(:, 1:n_columns);
    return
end
if n_have + 1 == n_columns
    aligned = psth;
    aligned(:, n_columns) = NaN;
    return
end

error('match_psth_columns:Width', ...
    ['Behavior onset has %d time bins and play-bout onset has %d. ' ...
     'The pooling line used to take columns 1:4000 (the ±20 s grid at 10 ms). ' ...
     '801 bins is -20:0.05:20 s. Those widths cannot be stacked.'], ...
    n_have, n_columns);
end
