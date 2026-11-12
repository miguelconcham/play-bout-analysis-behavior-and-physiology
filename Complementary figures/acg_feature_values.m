function [values, names] = acg_feature_values(logit, lag)
%ACG_FEATURE_VALUES  Autocorrelogram features for one neuron per row.
%
% logit is neurons by lags. lag is in seconds. It is not assumed to start
% at 0, so a first bin at 0.05 s does not shift an interval by one bin.
% Only positive lags are used.
%
% The first trough is the minimum on (0, 0.25] s. Each later extremum is
% the max or min on the following 0.40 s. Modulation1 is the first peak
% minus the first trough. Modulation2 is the second peak minus the second
% trough. PeakInterval and TroughInterval are those time gaps, in seconds.
% The frequencies are their reciprocals. Persistence is Modulation2 divided
% by Modulation1. OscillationScore is the mean of the two modulations.
% Baseline is the mean logit at lags of at least 0.60 s.

names = {'Modulation1', 'Modulation2', 'Persistence', ...
    'PeakInterval', 'PeakFrequency', 'TroughInterval', 'TroughFrequency', ...
    'OscillationScore', 'Baseline'};
logit = double(logit);
lag = double(lag(:).');
if isempty(logit)
    values = zeros(0, numel(names));
    return
end
if size(logit, 2) ~= numel(lag)
    error('acg_feature_values:axis', ...
        'The logit matrix has %d lags and the axis has %d.', size(logit, 2), numel(lag));
end
n = size(logit, 1);
values = nan(n, numel(names));
for i = 1:n
    values(i, :) = one_neuron_features(logit(i, :), lag);
end
end


function row = one_neuron_features(y, lag)
row = nan(1, 9);
[trough1_t, trough1_y] = lag_extremum(lag, y, 0, 0.25, false);
[peak1_t, peak1_y] = lag_extremum(lag, y, trough1_t, trough1_t + 0.40, true);
[trough2_t, trough2_y] = lag_extremum(lag, y, peak1_t, peak1_t + 0.40, false);
[peak2_t, peak2_y] = lag_extremum(lag, y, trough2_t, trough2_t + 0.40, true);
row(1) = peak1_y - trough1_y;
row(2) = peak2_y - trough2_y;
if isfinite(row(1)) && row(1) ~= 0 && isfinite(row(2))
    row(3) = row(2) / row(1);
end
row(4) = peak2_t - peak1_t;
if isfinite(row(4)) && row(4) > 0
    row(5) = 1 / row(4);
end
row(6) = trough2_t - trough1_t;
if isfinite(row(6)) && row(6) > 0
    row(7) = 1 / row(6);
end
if isfinite(row(1)) && isfinite(row(2))
    row(8) = mean([row(1), row(2)]);
end
late = lag >= 0.60 & isfinite(y);
if any(late)
    row(9) = mean(y(late));
end
end


function [t_ext, y_ext] = lag_extremum(lag, y, t_start, t_stop, want_max)
t_ext = NaN;
y_ext = NaN;
if ~isfinite(t_start) || ~isfinite(t_stop)
    return
end
idx = find(lag > t_start & lag <= t_stop & isfinite(y));
if isempty(idx)
    return
end
if want_max
    [y_ext, local] = max(y(idx));
else
    [y_ext, local] = min(y(idx));
end
t_ext = lag(idx(local));
end
