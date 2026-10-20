function locations = read_npx_continuous_locations(path)
%READ_NPX_CONTINUOUS_LOCATIONS  Session names and continuous paths.
%
% Reads NPX_raw_data_folders.xlsx. locations(k).session is the session
% folder name and locations(k).path is the continuous file or its folder.

sheet = read_location_sheet(path);
headers = sheet(1, :);
if size(sheet, 1) < 2
    rows = cell(0, size(sheet, 2));
else
    rows = sheet(2:end, :);
end
session_col = 1;
if ~isempty(rows)
    sample = canonical_text(rows{1, 1});
    if ~isempty(regexp(sample, '[\\/]', 'once'))
        session_col = find_header(headers, 'session');
        if isempty(session_col)
            error('read_npx_continuous_locations:columns', ...
                ['The first raw-folder column looks like a path, ' ...
                 'and no session column was found.']);
        end
    end
end
path_col = [];
for header_index = 1:numel(headers)
    if header_index == session_col
        continue
    end
    if strcmp(header_role(headers{header_index}), 'path')
        path_col = header_index;
        break
    end
end
if isempty(path_col)
    path_col = first_path_column(rows, session_col);
end
if isempty(path_col)
    error('read_npx_continuous_locations:columns', ...
        ['NPX_raw_data_folders needs a session name and a file path. ' ...
         'Headers found: %s'], strjoin(headers, ', '));
end

locations = struct('session', {}, 'path', {});
for row_index = 1:size(rows, 1)
    session = canonical_text(rows{row_index, session_col});
    file_path = canonical_text(rows{row_index, path_col});
    if isempty(session) && isempty(file_path)
        continue
    end
    locations(end + 1).session = session; %#ok<AGROW>
    locations(end).path = file_path;
end
end

function column = find_header(headers, role_name)
column = [];
for index = 1:numel(headers)
    if strcmp(header_role(headers{index}), role_name)
        column = index;
        return
    end
end
end

function role = header_role(header)
key = lower(regexprep(canonical_text(header), '[^A-Za-z0-9]', ''));
role = '';
if isempty(key)
    return
end
if any(strcmp(key, {'session', 'sessionname', 'sessionfolder', 'animalfolder'})) ...
        || ~isempty(strfind(key, 'session'))
    role = 'session';
elseif ~isempty(strfind(key, 'path')) || ~isempty(strfind(key, 'location')) ...
        || ~isempty(strfind(key, 'folder')) || ~isempty(strfind(key, 'file')) ...
        || ~isempty(strfind(key, 'continuous')) || ~isempty(strfind(key, 'directory'))
    role = 'path';
end
end

function column = first_path_column(rows, session_col)
column = [];
if isempty(rows)
    return
end
for column_index = 1:size(rows, 2)
    if column_index == session_col
        continue
    end
    n_path = 0;
    n_filled = 0;
    for row_index = 1:size(rows, 1)
        text = canonical_text(rows{row_index, column_index});
        if isempty(text)
            continue
        end
        n_filled = n_filled + 1;
        if ~isempty(regexp(text, '[\\/]', 'once')) || ~isempty(strfind(lower(text), '.dat'))
            n_path = n_path + 1;
        end
    end
    if n_filled > 0 && n_path == n_filled
        column = column_index;
        return
    end
end
end

function sheet = read_location_sheet(path)
if exist(path, 'file') ~= 2
    error('read_npx_continuous_locations:file', 'Could not find %s', path);
end
try
    table_data = readtable(path, 'VariableNamingRule', 'preserve');
catch
    try
        table_data = readtable(path, 'PreserveVariableNames', true);
    catch
        table_data = readtable(path);
    end
end
headers = cellstr(table_data.Properties.VariableNames);
if exist('table2cell', 'file') == 2 || exist('table2cell', 'builtin') == 5
    rows = table2cell(table_data);
else
    rows = cell(height(table_data), width(table_data));
    for column_index = 1:numel(headers)
        rows(:, column_index) = column_to_cell(table_data.(headers{column_index}));
    end
end
if isempty(rows)
    rows = cell(0, numel(headers));
end
sheet = [headers; rows];
end

function cells = column_to_cell(column)
n = numel(column);
cells = cell(n, 1);
if iscell(column)
    for index = 1:n
        cells{index} = column{index};
    end
elseif isnumeric(column) || islogical(column)
    for index = 1:n
        cells{index} = column(index);
    end
elseif exist('isstring', 'builtin') && isstring(column)
    cells = cellstr(column);
elseif ischar(column)
    cells = cellstr(column);
else
    for index = 1:n
        cells{index} = column(index);
    end
end
end

function text = canonical_text(value)
text = '';
if iscell(value)
    if numel(value) ~= 1
        return
    end
    text = canonical_text(value{1});
    return
end
if isnumeric(value) && isscalar(value)
    if ~isfinite(value)
        return
    end
    if value == round(value)
        text = sprintf('%.0f', value);
    else
        text = sprintf('%g', value);
    end
    return
end
if exist('isstring', 'builtin') && isstring(value)
    if numel(value) ~= 1
        return
    end
    text = strtrim(char(value));
    return
end
if ischar(value)
    text = strtrim(value);
    if size(text, 1) > 1
        text = '';
    end
end
end
