function scores = score_acg_curves(acg, centers, smooth_n, delta_lag)
%SCORE_ACG_CURVES  Burst, delta-cycle peaks, and folded delta power.
%
% Each row is one neuron. Peak fields come from acg_delta_peaks. Folded
% delta power is the 2-5 Hz power of the half-mirrored autocorrelogram
% after 100 ms. Both are the measures used for the locking-period
% paired tests.

if nargin < 3 || isempty(smooth_n)
    smooth_n = 1;
end
if nargin < 4 || isempty(delta_lag)
    delta_lag = [0.2 0.5];
end
acg = double(acg);
if isvector(acg)
    acg = acg(:)';
end
n = size(acg, 1);
scores = empty_acg_scores(n);
if n < 1 || isempty(centers)
    return
end
scores.acg_delta = folded_acg_delta_power(acg, centers, delta_lag);
if ~delta_peak_finder_available()
    return
end
for row = 1:n
    peaks = acg_delta_peaks(acg(row, :), centers, smooth_n);
    scores.burst(row) = peaks.burst;
    scores.first_measurable(row) = peaks.first_measurable;
    scores.first_present(row) = peaks.first_present;
    scores.first_lag(row) = peaks.first_lag;
    scores.first_prominence(row) = peaks.first_prominence;
    scores.second_measurable(row) = peaks.second_measurable;
    scores.second_present(row) = peaks.second_present;
    scores.second_lag(row) = peaks.second_lag;
    scores.second_prominence(row) = peaks.second_prominence;
    scores.n_cycles(row) = peaks.n_cycles;
end
end


function scores = empty_acg_scores(n)
scores.burst = nan(n, 1);
scores.first_measurable = false(n, 1);
scores.first_present = false(n, 1);
scores.first_lag = nan(n, 1);
scores.first_prominence = nan(n, 1);
scores.second_measurable = false(n, 1);
scores.second_present = false(n, 1);
scores.second_lag = nan(n, 1);
scores.second_prominence = nan(n, 1);
scores.n_cycles = nan(n, 1);
scores.acg_delta = nan(n, 1);
end


function ok = delta_peak_finder_available()
persistent found
if ~isempty(found)
    ok = found;
    return
end
found = false;
if exist('findpeaks', 'file') ~= 2
    ok = false;
    return
end
try
    [~, locs] = findpeaks([0 1 0], [0 1 2], ...
        'MinPeakDistance', 1, 'MinPeakWidth', 0.5);
    found = isnumeric(locs);
catch
    found = false;
end
ok = found;
end
