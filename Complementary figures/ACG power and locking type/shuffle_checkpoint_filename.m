function file_name = shuffle_checkpoint_filename(n_shuffles, method)
%SHUFFLE_CHECKPOINT_FILENAME  Overlap is the standard autocorrelogram.
%
% Overlap keeps locking_period_acg_shuffle_nN.mat, the name used before a
% method selector existed. Eligible-spike nulls use a different name so
% they do not replace that file. A missing method is overlap.

if nargin < 2 || isempty(method)
    method = 'overlap';
end
method = lower(strtrim(char(method)));
if strcmp(method, 'eligible')
    file_name = sprintf('locking_period_acg_shuffle_eligible_n%d.mat', n_shuffles);
elseif strcmp(method, 'overlap')
    file_name = sprintf('locking_period_acg_shuffle_n%d.mat', n_shuffles);
else
    error('shuffle_checkpoint_filename:method', ...
        'acg_method must be ''overlap'' or ''eligible''.');
end
end
