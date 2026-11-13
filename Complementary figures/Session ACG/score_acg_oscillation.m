function score = score_acg_oscillation(observed_row, null_rows, centers, smooth_n, alpha)
%SCORE_ACG_OSCILLATION  Delta-peak test of one autocorrelogram against its null.
%
% The curve and each null row are scored with acg_delta_peaks. The
% one-sided p-value is (1 + number of null prominences at least as large)
% / (number of finite nulls + 1). A missing null peak counts as prominence
% 0 when the lag window was observed. oscillates is true only when the
% observed curve has a delta peak and that prominence clears alpha.
% pctl is the smoothed folded curve as a percentile of the null curves.
% 100 means that lag is above every shuffle.

score = blank_oscillation_score();
observed_row = double(observed_row(:))';
if isempty(observed_row) || ~any(isfinite(observed_row))
    return
end
observed = acg_delta_peaks(observed_row, centers, smooth_n);
score.first_present = observed.first_present;
score.first_lag = observed.first_lag;
score.first_prominence = observed.first_prominence;
score.second_present = observed.second_present;
score.second_lag = observed.second_lag;
score.second_prominence = observed.second_prominence;
score.second_measurable = observed.second_measurable;
score.pctl_lag = observed.fold_t;
if isempty(null_rows)
    return
end
null_rows = double(null_rows);
if size(null_rows, 2) ~= numel(observed_row)
    error('score_acg_oscillation:width', ...
        'Null autocorrelograms must have one column per lag bin.');
end
n_shuffles = size(null_rows, 1);
null_fold = nan(n_shuffles, numel(observed.fold_t));
null_first = nan(n_shuffles, 1);
null_second = nan(n_shuffles, 1);
for s = 1:n_shuffles
    if ~any(isfinite(null_rows(s, :)))
        continue
    end
    peaks = acg_delta_peaks(null_rows(s, :), centers, smooth_n);
    null_first(s) = null_prominence(peaks, 'first');
    null_second(s) = null_prominence(peaks, 'second');
    if numel(peaks.fold_y) == numel(observed.fold_t)
        null_fold(s, :) = peaks.fold_y;
    end
end
score.first_p = exception_p(null_prominence(observed, 'first'), null_first);
if observed.second_measurable
    score.second_p = exception_p(null_prominence(observed, 'second'), null_second);
end
score.n_shuffles = sum(isfinite(null_first));
score.pctl = amplitude_percentile(observed.fold_y, null_fold);
score.oscillates = observed.first_present && isfinite(score.first_p) && score.first_p <= alpha;
end


function value = null_prominence(peaks, which)
if ~peaks.first_measurable
    value = NaN;
    return
end
if strcmp(which, 'first')
    value = peaks.first_prominence;
    if ~isfinite(value)
        value = 0;
    end
    return
end
if peaks.second_measurable
    value = peaks.second_prominence;
    if ~isfinite(value)
        value = 0;
    end
else
    value = 0;
end
end


function p = exception_p(observed, nulls)
nulls = nulls(isfinite(nulls));
if ~isfinite(observed) || isempty(nulls)
    p = NaN;
    return
end
p = (1 + sum(nulls >= observed)) / (numel(nulls) + 1);
end


function percentile = amplitude_percentile(observed, null_rows)
observed = double(observed(:))';
percentile = nan(size(observed));
if isempty(null_rows) || isempty(observed)
    return
end
finite = isfinite(null_rows);
n = sum(finite, 1);
below = sum(finite & null_rows < observed, 1);
equal = sum(finite & null_rows == observed, 1);
percentile = 100 * (below + 0.5 * equal) ./ n;
percentile(n == 0 | ~isfinite(observed)) = NaN;
end


function score = blank_oscillation_score()
score.first_present = false;
score.first_lag = NaN;
score.first_prominence = NaN;
score.first_p = NaN;
score.second_present = false;
score.second_lag = NaN;
score.second_prominence = NaN;
score.second_measurable = false;
score.second_p = NaN;
score.oscillates = false;
score.pctl_lag = [];
score.pctl = [];
score.n_shuffles = 0;
end
