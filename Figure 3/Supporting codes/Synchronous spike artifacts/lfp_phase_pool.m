function [pool_index, sync_sample] = lfp_phase_pool(n_samples, sample_rate, ...
    bout_starts, hist_range, sync_times, guard_seconds)
%LFP_PHASE_POOL  LFP samples eligible for a phase-matched null.
%
% Sample k has time k/sample_rate, matching lfp_time = (1:N)'/sample_rate.
% The pool is every sample inside any play-bout window
% [bout_start + hist_range(1), bout_start + hist_range(2)]. Samples within
% guard_seconds of a synchronous spike are removed. sync_sample is the
% nearest LFP index of each synchronous spike, or 0 when that spike falls
% outside the trace.

window_mask = bout_window_mask(bout_starts, hist_range, n_samples, sample_rate);
sync_sample = nearest_lfp_sample(sync_times, n_samples, sample_rate);
pool_index = phase_pool_index(window_mask, sync_sample, ...
    round(guard_seconds * sample_rate));
end

function mask = bout_window_mask(bout_starts, hist_range, n_samples, sample_rate)
mask = false(n_samples, 1);
if isempty(bout_starts) || n_samples < 1
    return
end
starts = bout_starts(:) + hist_range(1);
stops = bout_starts(:) + hist_range(2);
first = max(1, ceil(starts * sample_rate));
last = min(n_samples, floor(stops * sample_rate));
for index = 1:numel(first)
    if first(index) <= last(index)
        mask(first(index):last(index)) = true;
    end
end
end

function sample = nearest_lfp_sample(spike_time, n_samples, sample_rate)
sample = round(double(spike_time(:)) * sample_rate);
outside = ~isfinite(sample) | sample < 1 | sample > n_samples;
sample(outside) = 0;
end

function pool_index = phase_pool_index(window_mask, sync_sample, guard_samples)
exclude = false(size(window_mask));
sync_sample = sync_sample(sync_sample >= 1);
n_samples = numel(window_mask);
guard_samples = max(0, round(guard_samples));
for index = 1:numel(sync_sample)
    first = max(1, sync_sample(index) - guard_samples);
    last = min(n_samples, sync_sample(index) + guard_samples);
    exclude(first:last) = true;
end
pool_index = find(window_mask & ~exclude);
end
