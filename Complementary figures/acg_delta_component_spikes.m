function component = acg_delta_component_spikes(spike_times, run_starts, run_ends, lag_range)
%ACG_DELTA_COMPONENT_SPIKES  Spikes that make one autocorrelogram side peak.
%
% The autocorrelogram counts pairs. With a symmetric window each pair is
% entered twice, once at +lag and once at -lag. In a train spaced about one
% delta cycle apart, every interior spike is the earlier member of one pair
% and the later member of another, so the histogram lists that spike twice.
% Phase belongs to the spike, not to the histogram entry.
%
% component.mask is true once for each spike that belongs to at least one
% pair inside a single firing run whose lag is inside lag_range. Pairs that
% cross a gap between runs are left out, because the run autocorrelograms
% never saw them. component.listings counts how many times the symmetric
% autocorrelogram would list each spike as a reference in that lag window.
% component.n_acg_entries is the sum of those listings and equals twice the
% number of pairs. The phase sample is mask, not listings.

spike_times = double(spike_times(:));
n_raw = numel(spike_times);
component = empty_component(n_raw);
if nargin < 4 || numel(lag_range) ~= 2 || ~all(isfinite(lag_range)) ...
        || lag_range(1) < 0 || ~(lag_range(2) > lag_range(1))
    error('acg_delta_component_spikes:lag', ...
        'lag_range must be [min max] seconds, with 0 <= min < max.');
end
if n_raw < 1
    return
end

finite = isfinite(spike_times);
if ~any(finite)
    return
end
raw_index = find(finite);
[t, order] = sort(spike_times(finite));
sorted_raw = raw_index(order);
n = numel(t);
run_id = assign_runs(t, run_starts, run_ends);

lo = lag_range(1);
hi = lag_range(2);
listings = zeros(n, 1);
n_pairs = 0;
right = 1;
for i = 1:n
    if run_id(i) == 0
        continue
    end
    while right <= n && t(right) < t(i) + lo
        right = right + 1;
    end
    k = right;
    while k <= n && t(k) <= t(i) + hi
        % Forward only, so the -lag copy of this pair is not a second pair.
        if run_id(k) == run_id(i) && t(k) > t(i)
            n_pairs = n_pairs + 1;
            listings(i) = listings(i) + 1;
            listings(k) = listings(k) + 1;
        end
        k = k + 1;
    end
end

component.listings(sorted_raw) = listings;
component.in_run(sorted_raw) = run_id > 0;
component.mask(sorted_raw) = listings > 0;
component.n_pairs = n_pairs;
component.n_acg_entries = sum(listings);
component.n_spikes = sum(component.mask);
component.n_run_spikes = sum(component.in_run);
end


function component = empty_component(n_raw)
component.mask = false(n_raw, 1);
component.in_run = false(n_raw, 1);
component.listings = zeros(n_raw, 1);
component.n_pairs = 0;
component.n_acg_entries = 0;
component.n_spikes = 0;
component.n_run_spikes = 0;
end


function run_id = assign_runs(t, run_starts, run_ends)
n = numel(t);
run_id = zeros(n, 1);
if nargin < 2 || isempty(run_starts) || isempty(run_ends)
    return
end
run_starts = double(run_starts(:));
run_ends = double(run_ends(:));
n_runs = min(numel(run_starts), numel(run_ends));
for r = 1:n_runs
    if ~isfinite(run_starts(r)) || ~isfinite(run_ends(r)) || run_ends(r) <= run_starts(r)
        continue
    end
    in = t >= run_starts(r) & t <= run_ends(r) & run_id == 0;
    run_id(in) = r;
end
end
