function pooled = pool_period_acgs(acg, n_eligible, bin_edges, run_length, score, already_corrected, low_percent, high_percent, method)
%POOL_PERIOD_ACGS  Pool firing-run rows with the selected autocorrelogram.
%
% 'eligible' reads partner counts in acg and eligible-spike counts in
% n_eligible. 'overlap' reads pair counts, or overlap-corrected rates when
% already_corrected is true, and ignores n_eligible.

method = canonical_acg_method(method);
if strcmp(method, 'eligible')
    pooled = pool_eligible_period_acgs(acg, n_eligible, run_length, score, ...
        low_percent, high_percent);
else
    pooled = pool_locking_period_acgs(acg, bin_edges, run_length, score, ...
        already_corrected, low_percent, high_percent);
end
end
