function stats = speed_peak_correlation(speed, delta_cf, theta_cf, min_n)
%SPEED_PEAK_CORRELATION  Within-session correlations of peak frequency.
%
% speed, delta_cf, and theta_cf are one entry per speed bin. Pearson and
% Spearman correlations are estimated separately for delta versus speed,
% theta versus speed, and delta versus theta. A pair is omitted when either
% value is missing. Fewer than min_n finite pairs, or a constant series,
% returns NaN rather than a correlation.

speed = double(speed(:));
delta_cf = double(delta_cf(:));
theta_cf = double(theta_cf(:));
if numel(speed) ~= numel(delta_cf) || numel(speed) ~= numel(theta_cf)
    error('speed_peak_correlation:size', ...
        'speed, delta peak frequency, and theta peak frequency must match.');
end
if nargin < 4 || isempty(min_n)
    min_n = 3;
end
if ~isscalar(min_n) || min_n < 3 || min_n ~= round(min_n)
    error('speed_peak_correlation:min', ...
        'min_n must be an integer of at least 3.');
end

[n_delta, pearson_delta, pearson_delta_p, spearman_delta, spearman_delta_p] = ...
    pair_correlation(speed, delta_cf, min_n);
[n_theta, pearson_theta, pearson_theta_p, spearman_theta, spearman_theta_p] = ...
    pair_correlation(speed, theta_cf, min_n);
both = isfinite(speed) & isfinite(delta_cf) & isfinite(theta_cf);
[n_delta_theta, pearson_delta_theta, pearson_delta_theta_p, ...
    spearman_delta_theta, spearman_delta_theta_p] = pair_correlation( ...
    delta_cf(both), theta_cf(both), min_n);

stats = struct( ...
    'min_n', min_n, ...
    'n_delta', n_delta, ...
    'pearson_delta', pearson_delta, ...
    'pearson_delta_p', pearson_delta_p, ...
    'spearman_delta', spearman_delta, ...
    'spearman_delta_p', spearman_delta_p, ...
    'n_theta', n_theta, ...
    'pearson_theta', pearson_theta, ...
    'pearson_theta_p', pearson_theta_p, ...
    'spearman_theta', spearman_theta, ...
    'spearman_theta_p', spearman_theta_p, ...
    'n_delta_theta', n_delta_theta, ...
    'pearson_delta_theta', pearson_delta_theta, ...
    'pearson_delta_theta_p', pearson_delta_theta_p, ...
    'spearman_delta_theta', spearman_delta_theta, ...
    'spearman_delta_theta_p', spearman_delta_theta_p);
end

function [n, pearson_r, pearson_p, spearman_r, spearman_p] = ...
        pair_correlation(x, y, min_n)

x = x(:);
y = y(:);
if numel(x) ~= numel(y)
    error('speed_peak_correlation:size', 'The two series must match.');
end
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
n = numel(x);
pearson_r = NaN;
pearson_p = NaN;
spearman_r = NaN;
spearman_p = NaN;
if n < min_n || std(x, 0) == 0 || std(y, 0) == 0
    return
end
[pearson_r, pearson_p] = correlation_and_p(x, y);
[spearman_r, spearman_p] = correlation_and_p(average_rank(x), average_rank(y));
end

function [r, p] = correlation_and_p(x, y)
% Student-t p-value used by MATLAB corr for Pearson and Spearman.
x = x(:) - mean(x);
y = y(:) - mean(y);
r = (x.' * y) / sqrt((x.' * x) * (y.' * y));
n = numel(x);
if ~isfinite(r)
    r = NaN;
    p = NaN;
    return
end
if abs(r) >= 1
    r = sign(r);
    p = 0;
    return
end
degrees = n - 2;
t_stat = r * sqrt(degrees / (1 - r * r));
p = betainc(degrees / (degrees + t_stat * t_stat), degrees / 2, 0.5);
end

function ranks = average_rank(x)
% Average ranks for ties, matching MATLAB tiedrank.
[sorted, order] = sort(x(:));
ranks = zeros(numel(x), 1);
n = numel(x);
start = 1;
while start <= n
    stop = start;
    while stop < n && sorted(stop + 1) == sorted(start)
        stop = stop + 1;
    end
    ranks(order(start:stop)) = (start + stop) / 2;
    start = stop + 1;
end
end
