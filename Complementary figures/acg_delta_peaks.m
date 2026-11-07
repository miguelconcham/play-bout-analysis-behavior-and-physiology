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
% The first peak is the earliest maximum that passes. A second peak
% counts only when it lies near twice that lag, passes the same checks,
% and the whole search window was observed. Peaks within 50 ms of the
% right edge are left out, because findpeaks treats the endpoint as a valley.

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
peaks.first_present = true;
peaks.first_lag = locs(1);
peaks.first_height = pk(1);
peaks.first_prominence = prom(1);
tau = locs(1);
second_lo = second_scale(1) * tau;
second_hi = second_scale(2) * tau;
peaks.second_measurable = second_hi <= t_seg(end);
if ~peaks.second_measurable
    peaks.n_cycles = 1;
    return
end
in_second = locs >= second_lo & locs <= second_hi;
if ~any(in_second)
    peaks.second_prominence = 0;
    peaks.n_cycles = 1;
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
