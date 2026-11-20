function output_path = path_without_overwrite(requested_path)
%PATH_WITHOUT_OVERWRITE  A save path that leaves an existing file untouched.
%
% If requested_path is free, it is returned. If a file is already there,
% the path becomes name_rebuilt.ext, then name_rebuilt_2.ext, and so on.
% The existing file is not opened for writing.

if nargin < 1 || isempty(requested_path)
    error('path_without_overwrite:Path', 'A file path is required.');
end
if exist('isstring', 'builtin') && isstring(requested_path)
    requested_path = char(requested_path);
end
if ~ischar(requested_path)
    error('path_without_overwrite:Path', 'A file path is required.');
end
if exist(requested_path, 'file') ~= 2
    output_path = requested_path;
    return
end

[folder, name, ext] = fileparts(requested_path);
output_path = fullfile(folder, [name '_rebuilt' ext]);
copy_n = 2;
while exist(output_path, 'file') == 2
    output_path = fullfile(folder, sprintf('%s_rebuilt_%d%s', name, copy_n, ext));
    copy_n = copy_n + 1;
end
fprintf(['Left untouched: %s\n' ...
    'New result written beside it: %s\n'], requested_path, output_path);
end
