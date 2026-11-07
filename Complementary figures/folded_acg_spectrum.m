function [freq, spec, band_power] = folded_acg_spectrum(acg, centers, delta_lag)
%FOLDED_ACG_SPECTRUM  One-sided spectrum of the half-mirrored autocorrelogram.
%
% Positive and negative lags are averaged. Only lags after 100 ms are
% kept, so the burst does not enter the spectrum. freq and spec are cell
% columns, one row of the autocorrelogram each. spec holds the power of
% each positive Fourier component, 2*|coef|^2/n^2. band_power sums the
% components inside the delta band, which is the value returned by
% folded_acg_delta_power. delta_lag is the delta period window in
% seconds, [short long]. The band is 1/long to 1/short Hz.

if nargin < 3 || isempty(delta_lag)
    delta_lag = [0.2 0.5];
end
acg = double(acg);
centers = double(centers(:))';
if isvector(acg)
    acg = acg(:)';
end
n_rows = size(acg, 1);
freq = cell(n_rows, 1);
spec = cell(n_rows, 1);
band_power = nan(n_rows, 1);
t_pos = centers(centers > 0);
if numel(t_pos) < 8 || numel(delta_lag) ~= 2 || delta_lag(2) <= delta_lag(1)
    return
end
folded = nan(n_rows, numel(t_pos));
for k = 1:numel(t_pos)
    [~, ip] = min(abs(centers - t_pos(k)));
    [~, ineg] = min(abs(centers + t_pos(k)));
    folded(:, k) = mean(acg(:, [ip, ineg]), 2, 'omitnan');
end
keep = t_pos > 0.1;
t_pos = t_pos(keep);
folded = folded(:, keep);
f_lo = 1 / delta_lag(2);
f_hi = 1 / delta_lag(1);
for row = 1:n_rows
    y = folded(row, :);
    stop = find(~isfinite(y), 1, 'first');
    if ~isempty(stop)
        y = y(1:stop - 1);
    end
    if numel(y) < 8
        continue
    end
    t = t_pos(1:numel(y));
    dt = median(diff(t));
    if ~isfinite(dt) || dt <= 0
        continue
    end
    y = y - mean(y);
    n = numel(y);
    coef = fft(y);
    row_freq = (0:n - 1) / (n * dt);
    df = 1 / (n * dt);
    nyquist = 1 / (2 * dt);
    positive = row_freq > 0 & row_freq < nyquist;
    freq{row} = row_freq(positive);
    spec{row} = (2 * abs(coef(positive)) .^ 2) / (n ^ 2);
    % Half a bin past each edge keeps a rhythm at 2 Hz or 5 Hz when that
    % frequency falls between Fourier bins of this short segment.
    band = row_freq >= (f_lo - df / 2) & row_freq <= (f_hi + df / 2) ...
        & positive;
    if ~any(band)
        continue
    end
    band_power(row) = sum((2 * abs(coef(band)) .^ 2) / (n ^ 2));
end
end
