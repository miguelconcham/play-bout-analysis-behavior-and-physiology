function file_name = shuffle_checkpoint_filename(n_shuffles, method, count_method)
%SHUFFLE_CHECKPOINT_FILENAME  Overlap pair counts are the standard autocorrelogram.
%
% Overlap pair counts keep locking_period_acg_shuffle_nN.mat, the name
% used before a method selector existed. Eligible, FFT, and crosscorr
% nulls use a different name so they do not replace that file. A missing
% method is overlap. A missing count is pairs.

if nargin < 2 || isempty(method)
    method = 'overlap';
end
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
