function block = condition_window_acg(spike_times, epochs, bin_edges, power, method)
%CONDITION_WINDOW_ACG  Autocorrelogram of one play-power condition.
%
% 'overlap' is pair counts / (T-|lag|), pooled across half-open epochs.
% 'eligible' replaces that curve with the mean partner count of spikes
% that can see each whole bin. The firing rate is still spikes divided by
% epoch duration, including silent epochs. The eligible curve is not
% divided by T-|lag|.

if nargin < 4
    power = [];
end
if nargin < 5 || isempty(method)
    method = 'eligible';
end
method = canonical_acg_method(method);
block = window_condition_acg(spike_times, epochs, bin_edges, power);
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
