function session_acg_log(log_file, message)
%SESSION_ACG_LOG  Print one line and append it to the overnight log.
%
% The stamp is Europe/Berlin. The file is opened and closed on each line
% so a stopped run still leaves the lines that were written.

stamp = char(datetime('now', 'TimeZone', 'Europe/Berlin', 'Format', 'yyyy-MM-dd HH:mm:ss'));
line = sprintf('[%s Europe/Berlin] %s', stamp, message);
fprintf('%s\n', line);
if nargin < 1 || isempty(log_file)
    return
end
fid = fopen(log_file, 'a');
if fid < 0
    return
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', line);
end
