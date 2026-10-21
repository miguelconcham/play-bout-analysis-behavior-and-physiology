function paths = read_npx_continuous_locations(xlsx_path, session_id)
%READ_NPX_CONTINUOUS_LOCATIONS  Location for one NPX raw-data folder.
%
% NPX_raw_data_folders.xlsx has two columns:
%   SessionName   the folder under NPX raw data, e.g. 'B2S2 1110 Single'
%   Location      the continuous-file directory for that folder
% session_id is that folder name. paths is the Location cell on the
% matching row.

if nargin < 2 || isempty(session_id)
    error('read_npx_continuous_locations:session', ...
        'Pass the NPX folder name, for example ''B2S2 1110 Single''.');
end
if exist(xlsx_path, 'file') ~= 2
    error('read_npx_continuous_locations:file', 'Could not find %s', xlsx_path);
end
session_id = strtrim(char(session_id));

sheet = readtable(xlsx_path);
if size(sheet, 2) < 2
    error('read_npx_continuous_locations:columns', ...
        'NPX_raw_data_folders.xlsx needs SessionName in column 1 and Location in column 2.');
end
session_names = column_text(sheet{:, 1});
locations = column_text(sheet{:, 2});
paths = locations(strcmp(session_names, session_id));
paths = paths(~cellfun('isempty', paths));
if isempty(paths)
    error('read_npx_continuous_locations:missing', ...
        'No row for "%s" in %s.', session_id, xlsx_path);
end
end

function texts = column_text(column)
n = numel(column);
texts = cell(n, 1);
for index = 1:n
    if iscell(column)
        value = column{index};
    else
        value = column(index);
    end
    if exist('isstring', 'builtin') && isstring(value)
        value = char(value);
    end
    if ischar(value)
        texts{index} = strtrim(value);
    else
        texts{index} = '';
    end
end
end
