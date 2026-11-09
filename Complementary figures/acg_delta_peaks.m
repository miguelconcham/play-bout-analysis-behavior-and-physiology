function peaks = acg_delta_peaks(acg_row, centers, smooth_n)
%ACG_DELTA_PEAKS  Burst level and delta-cycle prominences of one autocorrelogram.
%
% The positive and negative lags are averaged. The burst is the mean level
% from the first positive bin through 50 ms, measured before smoothing.
% Oscillatory peaks are local maxima of the Gaussian-smoothed curve from
% 100 ms onward. smooth_n is the moving-average width in bins. The
% Gaussian standard deviation is one quarter of that width.
%
% A local maximum is a peak only when all three of these hold:
%   the width at half-prominence is at least 40 ms, so a one-bin
%   deflection is ignored;
%   the prominence is at least a fifth of the 10-90% range of the
%   smoothed curve in the search window;
%   the curve does not rise above the maximum on both sides within
%   200 ms. A bump inside a trough fails that last check, because both
%   shoulders of the trough are higher. A peak on a decaying tail does
%   not: only the side toward zero is higher.
% The first peak is the earliest maximum in the delta period, 200-500 ms.
% A maximum before 200 ms is a faster event. It is kept as the first peak
% only when no maximum lies in that delta period; otherwise the delta
% maximum is the first peak and the early one is not its second cycle.
% A second peak counts only when it lies near twice the chosen lag, passes
% the same checks, and the whole search window was observed. Peaks within
% 50 ms of the right edge are left out, because findpeaks treats the
% endpoint as a valley.

if nargin < 3 || isempty(smooth_n)
    smooth_n = 1;
end
burst_end = 0.05;
search_start = 0.10;
edge_margin = 0.05;
min_distance = 0.15;
min_width = 0.04;
shoulder_window = 0.20;
prominence_fraction = 0.20;
delta_period = [0.20 0.50];
second_scale = [1.7 2.3];

peaks = empty_peaks();
[t_pos, y_raw] = fold_acg_row(acg_row, centers);
if isempty(t_pos)
    return
end
peaks.fold_t = t_pos;
peaks.fold_y = y_raw;
burst_idx = t_pos <= burst_end & isfinite(y_raw);
if any(burst_idx)
    peaks.burst = mean(y_raw(burst_idx));
end

y = smooth_positive_acg(y_raw, t_pos, smooth_n);
peaks.fold_y = y;
search = t_pos >= search_start & t_pos <= (t_pos(end) - edge_margin);
if ~any(search)
    return
end
t_seg = t_pos(search);
y_seg = y(search);
stop = find(~isfinite(y_seg), 1, 'first');
if ~isempty(stop)
    t_seg = t_seg(1:stop - 1);
    y_seg = y_seg(1:stop - 1);
end
if numel(y_seg) < 5
    return
end
peaks.first_measurable = true;
level = y_seg(isfinite(y_seg));
scale = prctile(level, 90) - prctile(level, 10);

[pk, locs, ~, prom] = findpeaks(y_seg, t_seg, ...
    'MinPeakDistance', min_distance, 'MinPeakWidth', min_width);
if isempty(locs)
    peaks.first_prominence = 0;
    peaks.n_cycles = 0;
    return
end
[locs, order] = sort(locs(:));
pk = pk(order);
prom = prom(order);
keep = false(numel(locs), 1);
for i = 1:numel(locs)
    keep(i) = peak_is_reliable(t_seg, y_seg, locs(i), pk(i), prom(i), ...
        scale, shoulder_window, prominence_fraction);
end
locs = locs(keep);
pk = pk(keep);
prom = prom(keep);
if isempty(locs)
    peaks.first_prominence = 0;
    peaks.n_cycles = 0;
    return
end
in_delta = locs >= delta_period(1) & locs <= delta_period(2);
if any(in_delta)
    first_pick = find(in_delta, 1, 'first');
else
    first_pick = 1;
end
peaks.first_present = true;
peaks.first_lag = locs(first_pick);
peaks.first_height = pk(first_pick);
peaks.first_prominence = prom(first_pick);
tau = locs(first_pick);
second_lo = second_scale(1) * tau;
second_hi = second_scale(2) * tau;
peaks.second_measurable = second_hi <= t_seg(end);
if ~peaks.second_measurable
    peaks.n_cycles = 1;
    peaks = add_third_peak(peaks, locs, pk, prom, t_seg(end));
    return
end
in_second = locs >= second_lo & locs <= second_hi;
if ~any(in_second)
    peaks.second_prominence = 0;
    peaks.n_cycles = 1;
    peaks = add_third_peak(peaks, locs, pk, prom, t_seg(end));
    return
end
second_locs = locs(in_second);
second_prom = prom(in_second);
second_pk = pk(in_second);
[~, pick] = min(abs(second_locs - 2 * tau));
peaks.second_present = true;
peaks.second_lag = second_locs(pick);
peaks.second_height = second_pk(pick);
peaks.second_prominence = second_prom(pick);
peaks.n_cycles = 2;
peaks = add_third_peak(peaks, locs, pk, prom, t_seg(end));
end


function peaks = add_third_peak(peaks, locs, pk, prom, t_end)
peaks.third_measurable = false;
peaks.third_present = false;
peaks.third_lag = NaN;
peaks.third_height = NaN;
peaks.third_prominence = NaN;
if ~peaks.first_present || ~isfinite(peaks.first_lag) || isempty(locs)
    return
end
tau = peaks.first_lag;
third_lo = 2.6 * tau;
third_hi = 3.4 * tau;
peaks.third_measurable = third_hi <= t_end;
if ~peaks.third_measurable
    return
end
in_third = locs >= third_lo & locs <= third_hi;
in_third = in_third & abs(locs - peaks.first_lag) > 1e-6;
if peaks.second_present && isfinite(peaks.second_lag)
    in_third = in_third & abs(locs - peaks.second_lag) > 1e-6;
end
if ~any(in_third)
    return
end
third_locs = locs(in_third);
third_prom = prom(in_third);
third_pk = pk(in_third);
[~, pick] = min(abs(third_locs - 3 * tau));
peaks.third_present = true;
peaks.third_lag = third_locs(pick);
peaks.third_height = third_pk(pick);
peaks.third_prominence = third_prom(pick);
end


function keep = peak_is_reliable(t, y, loc, height, prom, scale, shoulder_window, prominence_fraction)
if ~isfinite(prom) || ~isfinite(height) || ~isfinite(scale) || scale <= 0
    keep = false;
    return
end
if prom < prominence_fraction * scale
    keep = false;
    return
end
left = t < loc & t >= loc - shoulder_window;
right = t > loc & t <= loc + shoulder_window;
keep = ~(any(y(left) > height) && any(y(right) > height));
end


function peaks = empty_peaks()
peaks.burst = NaN;
peaks.first_measurable = false;
peaks.first_present = false;
peaks.first_lag = NaN;
peaks.first_height = NaN;
peaks.first_prominence = NaN;
peaks.second_measurable = false;
peaks.second_present = false;
peaks.second_lag = NaN;
peaks.second_height = NaN;
peaks.second_prominence = NaN;
peaks.third_measurable = false;
peaks.third_present = false;
peaks.third_lag = NaN;
peaks.third_height = NaN;
peaks.third_prominence = NaN;
peaks.n_cycles = NaN;
peaks.fold_t = [];
peaks.fold_y = [];
end


function [t_pos, y_fold] = fold_acg_row(y, centers)
y = double(y(:))';
centers = double(centers(:))';
t_pos = centers(centers > 0);
y_fold = nan(size(t_pos));
for k = 1:numel(t_pos)
    [~, ip] = min(abs(centers - t_pos(k)));
    [~, ineg] = min(abs(centers + t_pos(k)));
    y_fold(k) = mean([y(ip), y(ineg)], 'omitnan');
end
end


function y = smooth_positive_acg(y_raw, t_pos, smooth_n)
y = y_raw;
if smooth_n <= 1 || isempty(y_raw)
    return
end
keep = t_pos > 0.02;
if ~any(keep)
    return
end
% Sigma is a quarter of the moving-average width, in bins.
y(keep) = gaussian_smooth_acg(y_raw(keep), smooth_n / 4);
end
