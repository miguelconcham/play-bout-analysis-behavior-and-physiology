function paths = read_npx_continuous_locations(xlsx_path, session_id)
%READ_NPX_CONTINUOUS_LOCATIONS  Continuous-file location for one session.
%
% session_id is the NPX raw-data folder, for example 'B2S2 1110 Single2'.
% The animal name is the first word. The session name is that full folder
% name. NPX_raw_data_folders.xlsx has those names in a column and the
% continuous-file location in another. paths is a cell array of the
% matching locations. One session can have more than one row.

if nargin < 2 || isempty(session_id)
    error('read_npx_continuous_locations:session', ...
        'Pass the session folder name, for example ''B2S2 1110 Single2''.');
end
if exist(xlsx_path, 'file') ~= 2
    error('read_npx_continuous_locations:file', 'Could not find %s', xlsx_path);
end
session_id = strtrim(char(session_id));
name_parts = regexp(session_id, '\s+', 'split');
name_parts = name_parts(~cellfun('isempty', name_parts));
animal = name_parts{1};
if numel(name_parts) > 1
    session_rest = strjoin(name_parts(2:end), ' ');
else
    session_rest = '';
end

sheet = readtable(xlsx_path);
headers = sheet.Properties.VariableNames;
animal_col = header_index(headers, {'animal'}, false);
session_col = header_index(headers, {'session'}, false);
location_col = header_index(headers, ...
    {'location', 'path', 'folder', 'file', 'continuous', 'directory'}, true);
if isempty(session_col)
    session_col = 1;
end
if isempty(location_col)
    if size(sheet, 2) < 2
        error('read_npx_continuous_locations:columns', ...
            'NPX_raw_data_folders needs a session column and a location column.');
    end
    location_col = 2;
    if location_col == session_col
        location_col = setdiff(1:size(sheet, 2), session_col);
        location_col = location_col(1);
    end
end

sessions = column_text(sheet{:, session_col});
if isempty(animal_col)
    keep = strcmp(sessions, session_id);
else
    animals = column_text(sheet{:, animal_col});
    keep = strcmp(animals, animal) & ...
        (strcmp(sessions, session_id) | strcmp(sessions, session_rest));
end
locations = column_text(sheet{:, location_col});
paths = locations(keep);
paths = paths(~cellfun('isempty', paths));
paths = unique_paths(paths);
if isempty(paths)
    error('read_npx_continuous_locations:missing', ...
        'No location for session "%s" in %s.', session_id, xlsx_path);
end
end

function index = header_index(headers, tokens, skip_name_columns)
index = [];
for column = 1:numel(headers)
    key = lower(regexprep(headers{column}, '[^A-Za-z0-9]', ''));
    if isempty(key)
        continue
    end
    if skip_name_columns && (~isempty(strfind(key, 'session')) || ~isempty(strfind(key, 'animal')))
        continue
    end
    for token_index = 1:numel(tokens)
        if ~isempty(strfind(key, tokens{token_index}))
            index = column;
            return
        end
    end
end
end

function texts = column_text(column)
n = numel(column);
texts = cell(n, 1);
for index = 1:n
    if iscell(column)
        value = column{index};
    elseif isnumeric(column) || islogical(column)
        value = column(index);
    else
        value = column(index);
    end
    if iscell(value)
        if numel(value) == 1
            value = value{1};
        else
            texts{index} = '';
            continue
        end
    end
    if isnumeric(value) || islogical(value)
        if isempty(value) || any(~isfinite(value))
            texts{index} = '';
        else
            texts{index} = strtrim(num2str(value));
        end
        continue
    end
    if exist('isstring', 'builtin') && isstring(value)
        value = char(value);
    end
    if isempty(value)
        texts{index} = '';
    else
        texts{index} = strtrim(char(value));
    end
end
end

function paths = unique_paths(paths)
kept = {};
for index = 1:numel(paths)
    if ~any(strcmp(kept, paths{index}))
        kept{end + 1} = paths{index}; %#ok<AGROW>
    end
end
paths = kept;
end
