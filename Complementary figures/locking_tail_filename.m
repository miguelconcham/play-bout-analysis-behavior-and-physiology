function file_name = locking_tail_filename(tail_percent, method)
%LOCKING_TAIL_FILENAME  Tail-comparison file for one autocorrelogram method.
%
% Overlap keeps locking_period_acg_tails_N.mat. Eligible inserts
% _eligible so a new tail estimate does not replace an overlap file.

method = canonical_acg_method(method);
if strcmp(method, 'eligible')
    file_name = sprintf('locking_period_acg_tails_eligible_%g.mat', tail_percent);
else
    file_name = sprintf('locking_period_acg_tails_%g.mat', tail_percent);
end
end
