function file_name = locking_shuffle_filename(n_shuffles, method, count_method)
%LOCKING_SHUFFLE_FILENAME  Checkpoint name for one autocorrelogram method.
%
% Overlap pair counts keep locking_period_acg_shuffle_nN.mat. Eligible,
% FFT, and crosscorr use a different name so a new null does not replace
% another count.

if nargin < 3
    count_method = [];
end
token = acg_count_token(method, count_method);
if isempty(token)
    file_name = sprintf('locking_period_acg_shuffle_n%d.mat', n_shuffles);
else
    file_name = sprintf('locking_period_acg_shuffle%s_n%d.mat', token, n_shuffles);
end
end
