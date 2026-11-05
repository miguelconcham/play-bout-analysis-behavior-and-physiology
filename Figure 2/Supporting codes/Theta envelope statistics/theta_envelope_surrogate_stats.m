function stats = theta_envelope_surrogate_stats(observed, null_draws)
%THETA_ENVELOPE_SURROGATE_STATS  Z-score and right-tail p against a null.
%
% observed is 1-by-frequency. null_draws is n_perm-by-frequency.
% p is the fraction of finite null draws at least as large as observed.

observed = double(observed(:).');
n_freq = numel(observed);
stats.null_mean = nan(1, n_freq);
stats.null_std = nan(1, n_freq);
stats.z = nan(1, n_freq);
stats.p = nan(1, n_freq);
if isempty(null_draws)
    return
end
if size(null_draws, 2) ~= n_freq
    error('theta_envelope_surrogate_stats:size', ...
        'Null draws must have one column per frequency.');
end
for f = 1:n_freq
    draws = null_draws(:, f);
    draws = draws(isfinite(draws));
    if isempty(draws) || ~isfinite(observed(f))
        continue
    end
    stats.null_mean(f) = mean(draws);
    if numel(draws) > 1
        stats.null_std(f) = std(draws);
    end
    if isfinite(stats.null_std(f)) && stats.null_std(f) > 0
        stats.z(f) = (observed(f) - stats.null_mean(f)) / stats.null_std(f);
    end
    stats.p(f) = mean(draws >= observed(f));
end
end
