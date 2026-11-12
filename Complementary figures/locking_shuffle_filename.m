function file_name = locking_shuffle_filename(n_shuffles, method)
%LOCKING_SHUFFLE_FILENAME  Checkpoint name for one autocorrelogram method.
%
% Overlap keeps locking_period_acg_shuffle_nN.mat. Eligible uses a
% different name so a new null does not replace the overlap checkpoint.

method = canonical_acg_method(method);
if strcmp(method, 'eligible')
    file_name = sprintf('locking_period_acg_shuffle_eligible_n%d.mat', n_shuffles);
else
    file_name = sprintf('locking_period_acg_shuffle_n%d.mat', n_shuffles);
end
end
