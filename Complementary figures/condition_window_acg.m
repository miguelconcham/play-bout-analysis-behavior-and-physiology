function block = condition_window_acg(spike_times, epochs, bin_edges, power, method, count_method)
%CONDITION_WINDOW_ACG  Autocorrelogram of one play-power condition.
%
% 'overlap' is the default. count_method chooses that overlap count:
% 'pairs' divides pair counts by T-|lag|, 'fft' uses the same correction
% on a binned FFT, and 'crosscorr' insets the reference spikes by the lag
% limit. 'eligible' replaces the curve with the mean partner count of
% spikes that can see each whole bin. The firing rate is still spikes
% divided by epoch duration, including silent epochs. The eligible curve
% is not divided by T-|lag|.

if nargin < 4
    power = [];
end
if nargin < 5 || isempty(method)
    method = 'overlap';
end
if nargin < 6 || isempty(count_method)
    count_method = 'pairs';
else
    count_method = canonical_count_method(count_method);
end
method = canonical_acg_method(method);
if strcmp(method, 'eligible') && ~strcmp(count_method, 'pairs')
    error('condition_window_acg:count', ...
        ['Eligible-spike autocorrelograms keep their own count. ' ...
        'acg_count applies to overlap estimates.']);
end
block = window_condition_acg(spike_times, epochs, bin_edges, power, count_method);
if strcmp(method, 'overlap') || block.n_epochs < 1
    return
end
epochs = double(epochs(:, 1:2));
keep = all(isfinite(epochs), 2) & epochs(:, 2) > epochs(:, 1);
epochs = epochs(keep, :);
[partner_counts, n_eligible] = eligible_window_counts( ...
    spike_times, epochs(:, 1), epochs(:, 2), bin_edges, false);
n_eligible(~isfinite(n_eligible)) = 0;
partner_counts(~isfinite(partner_counts)) = 0;
den = sum(n_eligible, 1);
num = sum(partner_counts, 1);
acg = nan(1, numel(den));
ok = den > 0;
acg(ok) = num(ok) ./ den(ok);
block.acg = acg;
end
