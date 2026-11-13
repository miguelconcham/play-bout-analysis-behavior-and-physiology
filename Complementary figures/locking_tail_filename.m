function file_name = locking_tail_filename(tail_percent, method, count_method)
%LOCKING_TAIL_FILENAME  Tail-comparison file for one autocorrelogram method.
%
% Overlap pair counts keep locking_period_acg_tails_N.mat. Eligible
% inserts _eligible. FFT and crosscorr insert their own token so a new
% tail estimate does not replace another count.

if nargin < 3
    count_method = [];
end
token = acg_count_token(method, count_method);
file_name = sprintf('locking_period_acg_tails%s_%g.mat', token, tail_percent);
end
