function plan = power_run_strata(power, spikes, duration, min_spikes, min_lag_s)
%POWER_RUN_STRATA  Contiguous power strata with about equal spike counts.
%
% Runs are ordered from low to high delta power. The cut points divide the
% cumulative spike count into equal shares. A run is never split. The
% result uses 4 strata, then 2, then 1, and keeps the first partition in
% which every stratum has at least min_spikes spikes and at least one run
% longer than min_lag_s. plan.k is 0 when no partition qualifies.
% plan.index is 0 for a run left out, and 1..k for a run in a stratum.

if nargin < 4 || isempty(min_spikes)
    min_spikes = 50;
end
if nargin < 5 || isempty(min_lag_s)
    min_lag_s = 0.5;
end
power = double(power(:));
spikes = double(spikes(:));
duration = double(duration(:));
n = numel(power);
plan.k = 0;
plan.index = zeros(n, 1);
if n < 1 || numel(spikes) ~= n || numel(duration) ~= n
    return
end
usable = isfinite(power) & isfinite(spikes) & spikes > 0 ...
    & isfinite(duration) & duration > 0;
if ~any(usable)
    return
end
for k = [4 2 1]
    index = zeros(n, 1);
    index(usable) = bins_by_spike_mass(power(usable), spikes(usable), k);
    if partition_is_usable(index, spikes, duration, k, min_spikes, min_lag_s)
        plan.k = k;
        plan.index = index;
        return
    end
end
end


function index = bins_by_spike_mass(power, spikes, k)
n = numel(power);
[~, order] = sortrows([power, (1:n)']);
ordered_spikes = spikes(order);
cumulative = cumsum(ordered_spikes);
edges = (1:(k - 1)) .* (cumulative(end) / k);
bin_sorted = ones(n, 1);
for i = 1:n
    bin_sorted(i) = 1 + sum(cumulative(i) > edges);
end
index = zeros(n, 1);
index(order) = bin_sorted;
end


function ok = partition_is_usable(index, spikes, duration, k, min_spikes, min_lag_s)
ok = false;
if k < 1
    return
end
for b = 1:k
    in = index == b;
    if ~any(in) || sum(spikes(in)) < min_spikes || max(duration(in)) <= min_lag_s
        return
    end
end
ok = true;
end
