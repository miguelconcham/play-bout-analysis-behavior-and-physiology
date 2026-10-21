function comparison = compare_sync_scalar(values, sync_index, n_repeats)
%COMPARE_SYNC_SCALAR  Mean at synchronous spikes versus matched random means.
%
% values is one number per spike in the comparison pool. For the LFP
% pipeline that pool is the spikes inside the play-bout windows, and
% sync_index indexes that pool. The observed value is the mean of the
% finite synchronous samples. Each of the n_repeats nulls is the mean of
% that same count of the other finite values, drawn without replacement
% inside the draw. Synchronous spikes of this pair are excluded. percentile
% is the fraction of null means strictly smaller than the observed mean, so
% a high percentile means the synchronous spikes sit higher than the
% non-synchronous spikes in the pool. mean_all is the mean of every finite
% value in the pool, including the synchronous spikes. mean_nonsync is the
% mean of the other finite values.

if nargin < 3 || isempty(n_repeats)
    n_repeats = 500;
end
n_repeats = round(double(n_repeats));
comparison.ok = false;
comparison.reason = '';
comparison.n_pool = 0;
comparison.n_sync = 0;
comparison.n_nonsync = 0;
comparison.mean_all = NaN;
comparison.mean_nonsync = NaN;
comparison.mean_sync = NaN;
comparison.null_means = nan(max(n_repeats, 0), 1);
comparison.percentile = NaN;

if ~isscalar(n_repeats) || ~isfinite(n_repeats) || n_repeats < 1
    comparison.reason = 'n_repeats must be a positive integer';
    return
end

values = double(values(:));
n_spikes = numel(values);
is_sync = false(n_spikes, 1);
sync_index = unique(round(double(sync_index(:)')));
sync_index = sync_index(sync_index >= 1 & sync_index <= n_spikes);
if ~isempty(sync_index)
    is_sync(sync_index) = true;
end

sync_values = values(is_sync);
sync_values = sync_values(isfinite(sync_values));
other = values(~is_sync);
other = other(isfinite(other));
comparison.n_pool = n_spikes;
comparison.n_sync = numel(sync_values);
comparison.n_nonsync = numel(other);
comparison.mean_all = mean_finite(values);
comparison.mean_nonsync = mean_finite(other);
if isempty(sync_values)
    comparison.reason = 'no finite values at the synchronous spikes';
    return
end
if numel(other) < numel(sync_values)
    comparison.reason = 'fewer other spikes than synchronous spikes';
    return
end

comparison.mean_sync = mean(sync_values);
n_sync = numel(sync_values);
for repeat_index = 1:n_repeats
    draw = other(randperm(numel(other), n_sync));
    comparison.null_means(repeat_index) = mean(draw);
end
comparison.percentile = mean(comparison.null_means < comparison.mean_sync);
comparison.ok = true;
end

function value = mean_finite(values)
values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = mean(values);
end
end
