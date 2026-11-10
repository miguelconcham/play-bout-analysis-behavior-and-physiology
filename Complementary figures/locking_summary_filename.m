function file_name = locking_summary_filename(maximize_strata_number, method)
%LOCKING_SUMMARY_FILENAME  Summary file for the strata rule and ACG method.
%
% Overlap keeps the original names. Eligible inserts _eligible so a new
% estimate does not replace a summary written by the overlap method.

method = canonical_acg_method(method);
if strcmp(method, 'eligible')
    method_token = '_eligible';
else
    method_token = '';
end
if maximize_strata_number
    file_name = ['locking_period_acg_summary', method_token, '_max_strata.mat'];
else
    file_name = ['locking_period_acg_summary', method_token, '.mat'];
end
end
