function file_name = play_power_retry_filename(method)
%PLAY_POWER_RETRY_FILENAME  Failed-session log for one autocorrelogram method.
%
% Overlap keeps the original log name. Eligible writes its own log.

method = canonical_acg_method(method);
if strcmp(method, 'eligible')
    file_name = 'play_power_acg_eligible_retry_sessions.mat';
else
    file_name = 'play_power_acg_retry_sessions.mat';
end
end
