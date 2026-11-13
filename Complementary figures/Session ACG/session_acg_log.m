function session_acg_log(log_file, message)
%SESSION_ACG_LOG  Print one line and append it to the overnight log.
%
% The stamp is Europe/Berlin when datetime is available. The file is
% opened and closed on each line so a stopped run still leaves the lines
% that were written.

use_datetime = exist('datetime', 'class') == 8 || exist('datetime', 'builtin') == 5 ...
    || exist('datetime', 'file') == 2;
if use_datetime
    stamp = char(datetime('now', 'TimeZone', 'Europe/Berlin', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    zone_label = 'Europe/Berlin';
else
    stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    zone_label = 'local';
end
line = sprintf('[%s %s] %s', stamp, zone_label, message);
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
