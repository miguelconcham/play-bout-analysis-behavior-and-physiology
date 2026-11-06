function shuffled = shuffle_spikes_within_runs(spike_times, run_starts, run_ends)
%SHUFFLE_SPIKES_WITHIN_RUNS  Uniform spike times that destroy the autocorrelogram.
%
% Each firing run keeps its spike count. Those spikes are replaced by
% independent uniform times inside that run. Interspike intervals are not
% preserved, so a rhythmic autocorrelogram is flattened. A shuffle of the
% intervals would keep a peak that comes from a preferred interval.

spike_times = double(spike_times(:));
run_starts = double(run_starts(:));
run_ends = double(run_ends(:));
n_runs = numel(run_starts);
counts = zeros(n_runs, 1);
if numel(run_ends) ~= n_runs
    error('shuffle_spikes_within_runs:runs', ...
        'Run starts and run ends must have the same length.');
end
for j = 1:n_runs
    if ~isfinite(run_starts(j)) || ~isfinite(run_ends(j)) || run_ends(j) <= run_starts(j)
        continue
    end
    counts(j) = sum(spike_times >= run_starts(j) & spike_times <= run_ends(j));
end
shuffled = zeros(sum(counts), 1);
pos = 1;
for j = 1:n_runs
    n = counts(j);
    if n < 1
        continue
    end
    width = run_ends(j) - run_starts(j);
    shuffled(pos:pos + n - 1) = run_starts(j) + width * rand(n, 1);
    pos = pos + n;
end
end
