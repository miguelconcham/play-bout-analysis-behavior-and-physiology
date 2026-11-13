function file_name = locking_summary_filename(maximize_strata_number, method, count_method)
%LOCKING_SUMMARY_FILENAME  Summary file for the strata rule and ACG method.
%
% Overlap pair counts keep the original names. Eligible inserts
% _eligible. FFT and crosscorr insert _fft or _crosscorr, so a new count
% does not replace a summary written by another count.

if nargin < 3
    count_method = [];
end
token = acg_count_token(method, count_method);
if maximize_strata_number
    file_name = ['locking_period_acg_summary', token, '_max_strata.mat'];
else
    file_name = ['locking_period_acg_summary', token, '.mat'];
end
end
