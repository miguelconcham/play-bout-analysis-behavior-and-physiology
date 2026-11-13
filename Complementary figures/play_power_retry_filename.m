function file_name = play_power_retry_filename(method, count_method)
%PLAY_POWER_RETRY_FILENAME  Failed-session log for one autocorrelogram method.
%
% Overlap pair counts keep the original log name. Eligible, FFT, and
% crosscorr each write their own log.

if nargin < 2
    count_method = [];
end
token = acg_count_token(method, count_method);
file_name = ['play_power_acg', token, '_retry_sessions.mat'];
end
