function report = compare_npx_folder_tables(experiment_info, raw_folders)
%COMPARE_NPX_FOLDER_TABLES  Check raw-file paths against experiment info.
%
%   report = compare_npx_folder_tables(experiment_info, raw_folders)
%
% experiment_info is Experiment_Info_GoodLooking.xlsx, as a table or as a
% cell array whose first row is the header. Columns, by header or by order:
%   1 animal code          e.g. B2S2
%   2 date folder          e.g. 20231110, 2023-11-10, or 1110
%   3 NPX recording folder
%   4 probe                folder name, or a probe number
%   5 recording number     folder name, or a number
%
% raw_folders is NPX_raw_data_folders.xlsx. Column 1 is the session name
% used under Data\NPX data\NPX raw data, e.g. 'B2S2 1110 Single2'. Another
% column is the full path of the continuous file or its folder.
%
% A row agrees when:
%   the session name is {animal} {MMDD} {Single|Dual|Single2}
%   MMDD is the date folder (0826 and 20230826 are the same day)
%   the label follows the animal code (D -> Dual, S -> Single or Single2)
%   the path contains, in order, the date folder, then the NPX recording
%   folder, then both the probe and the recording number
%
% Probe number 1 matches a path token ProbeA (Open Ephys, A = 1) or imec0
% (SpikeGLX). Recording number 2 matches a path token recording2.
% Probe and recording number may occur in either order after the recording
% folder. Pass the folder text itself in those columns when the name does
% not follow that pattern.
%
% report.ok is true when every experiment row has a matching path and every
% raw-folder row is listed in the experiment info. report.rows lists each
% experiment row and each raw-folder row that failed.
%
% This is not Area_limits_GoodLooking.xlsx. That workbook is the depth table.

[info_headers, info_rows] = sheet_cells(experiment_info);
[folder_headers, folder_rows] = sheet_cells(raw_folders);

if looks_like_area_limits(info_headers) || looks_like_area_limits(folder_headers)
    error('compare_npx_folder_tables:areaLimits', ...
        ['This sheet is Area_limits_GoodLooking (depth limits). ' ...
         'The experiment sheet is Experiment_Info_GoodLooking.xlsx.']);
end

[info_col, info_note] = experiment_columns(info_headers);
[folder_col, folder_note] = raw_folder_columns(folder_headers, folder_rows);

experiments = parse_experiments(info_rows, info_col);
raw = parse_raw_folders(folder_rows, folder_col);
[exp_results, raw_results] = match_folders(experiments, raw);

report.notes = notes_from(info_note, folder_note);
report.rows = [exp_results, raw_results];
report.n_experiment = numel(experiments);
report.n_raw = numel(raw);
if isempty(report.rows)
    report.n_ok = 0;
    report.n_bad = 0;
else
    report.n_ok = sum(strcmp({report.rows.status}, 'ok'));
    report.n_bad = numel(report.rows) - report.n_ok;
end
report.ok = report.n_bad == 0 && report.n_experiment > 0;
end

function [headers, rows] = sheet_cells(sheet)
if iscell(sheet)
    if isempty(sheet) || size(sheet, 1) < 1
        error('compare_npx_folder_tables:sheet', ...
            'The sheet is empty.');
    end
    headers = cell(1, size(sheet, 2));
    for c = 1:size(sheet, 2)
        headers{c} = canonical_text(sheet{1, c});
    end
    if size(sheet, 1) == 1
        rows = cell(0, size(sheet, 2));
    else
        rows = sheet(2:end, :);
    end
    return
end
if exist('istable', 'builtin') && istable(sheet)
    headers = cellstr(sheet.Properties.VariableNames);
    n = height(sheet);
    p = width(sheet);
    rows = cell(n, p);
    for c = 1:p
        rows(:, c) = column_to_cell(sheet.(headers{c}));
    end
    return
end
error('compare_npx_folder_tables:sheet', ...
    'Pass a table or a cell array with a header row.');
end

function cells = column_to_cell(col)
n = numel(col);
cells = cell(n, 1);
if iscell(col)
    for i = 1:n
        cells{i} = col{i};
    end
elseif isnumeric(col) || islogical(col)
    for i = 1:n
        cells{i} = col(i);
    end
elseif exist('isdatetime', 'builtin') && isdatetime(col)
    for i = 1:n
        cells{i} = col(i);
    end
elseif exist('isstring', 'builtin') && isstring(col)
    cells = cellstr(col);
elseif ischar(col)
    cells = cellstr(col);
else
    for i = 1:n
        cells{i} = col(i);
    end
end
end

function tf = looks_like_area_limits(headers)
keys = cellfun(@header_key, headers, 'UniformOutput', false);
tf = any(strcmp(keys, 'depthstart')) && any(strcmp(keys, 'animalname'));
end

function [col, note] = experiment_columns(headers)
roles = cellfun(@experiment_role, headers, 'UniformOutput', false);
col.animal = find(strcmp(roles, 'animal'), 1);
col.date = find(strcmp(roles, 'date'), 1);
col.recording_folder = find(strcmp(roles, 'recording_folder'), 1);
col.probe = find(strcmp(roles, 'probe'), 1);
col.recording_number = find(strcmp(roles, 'recording_number'), 1);
note = '';
if any(structfun(@isempty, col))
    if numel(headers) < 5
        error('compare_npx_folder_tables:columns', ...
            ['Experiment info needs animal, date, NPX recording folder, ' ...
             'probe, and recording number. Headers found: %s'], ...
            strjoin(headers, ', '));
    end
    col.animal = 1;
    col.date = 2;
    col.recording_folder = 3;
    col.probe = 4;
    col.recording_number = 5;
    note = ['Experiment info headers were not recognized (' ...
        strjoin(headers, ', ') '); used columns 1-5 in that order.'];
end
end

function role = experiment_role(header)
key = header_key(header);
role = '';
if isempty(key)
    return
end
if any(strcmp(key, {'recordingnumber', 'recnumber', 'recnum', ...
        'recordingnum', 'recordingn', 'recn'})) ...
        || ~isempty(strfind(key, 'recordingnumber')) ...
        || ~isempty(strfind(key, 'recnum'))
    role = 'recording_number';
elseif any(strcmp(key, {'animalcode', 'animal', 'animalid', 'mouse', 'subject', 'rat'})) ...
        || ~isempty(strfind(key, 'animalcode'))
    role = 'animal';
elseif ~isempty(strfind(key, 'probe'))
    role = 'probe';
elseif ~isempty(strfind(key, 'date'))
    role = 'date';
elseif ~isempty(strfind(key, 'folder')) || strcmp(key, 'npxrecording') ...
        || ~isempty(strfind(key, 'recordingfolder'))
    role = 'recording_folder';
end
end

function [col, note] = raw_folder_columns(headers, rows)
col.session = 1;
note = '';
if ~isempty(rows)
    sample = rows{1, 1};
    if ischar(sample) && ~isempty(regexp(sample, '[\\/]', 'once'))
        col.session = find(cellfun(@(h) strcmp(raw_role(h), 'session'), headers), 1);
        if isempty(col.session)
            error('compare_npx_folder_tables:columns', ...
                'The first raw-folder column looks like a path, and no session column was found.');
        end
    end
end
col.path = find(cellfun(@(h) strcmp(raw_role(h), 'path'), headers), 1);
if isempty(col.path)
    col.path = first_path_column(rows, col.session);
end
if isempty(col.path)
    error('compare_npx_folder_tables:columns', ...
        ['NPX_raw_data_folders needs a session name and a file path. ' ...
         'Headers found: %s'], strjoin(headers, ', '));
end
if col.session ~= 1
    note = 'Session name was not column 1 of the raw-folder sheet.';
end
end

function role = raw_role(header)
key = header_key(header);
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
for c = 1:size(rows, 2)
    if c == session_col
        continue
    end
    n_path = 0;
    n_filled = 0;
    for r = 1:size(rows, 1)
        text = canonical_text(rows{r, c});
        if isempty(text)
            continue
        end
        n_filled = n_filled + 1;
        if ~isempty(regexp(text, '[\\/]', 'once')) || ~isempty(strfind(lower(text), '.dat'))
            n_path = n_path + 1;
        end
    end
    if n_filled > 0 && n_path == n_filled
        column = c;
        return
    end
end
end

function experiments = parse_experiments(rows, col)
experiments = empty_experiment();
keep = 0;
for r = 1:size(rows, 1)
    animal = canonical_text(rows{r, col.animal});
    date_info = interpret_date(rows{r, col.date});
    recording_folder = canonical_text(rows{r, col.recording_folder});
    probe = canonical_text(rows{r, col.probe});
    recording_number = canonical_text(rows{r, col.recording_number});
    if isempty(animal) && isempty(date_info.mmdd) && isempty(recording_folder) ...
            && isempty(probe) && isempty(recording_number)
        continue
    end
    keep = keep + 1;
    experiments(keep).row = r + 1; %#ok<AGROW>
    experiments(keep).animal = animal;
    experiments(keep).date = date_info;
    experiments(keep).recording_folder = recording_folder;
    experiments(keep).probe = probe;
    experiments(keep).probe_tokens = probe_tokens(probe);
    experiments(keep).recording_number = recording_number;
    experiments(keep).recording_tokens = recording_tokens(recording_number);
end
if keep == 0
    experiments = empty_experiment();
end
end

function raw = parse_raw_folders(rows, col)
raw = empty_raw();
keep = 0;
for r = 1:size(rows, 1)
    session = canonical_text(rows{r, col.session});
    path = canonical_text(rows{r, col.path});
    if isempty(session) && isempty(path)
        continue
    end
    keep = keep + 1;
    parsed = parse_session_name(session);
    raw(keep).row = r + 1; %#ok<AGROW>
    raw(keep).session = session;
    raw(keep).path = path;
    raw(keep).segments = split_path(path);
    raw(keep).parse_ok = parsed.ok;
    raw(keep).animal = parsed.animal;
    raw(keep).mmdd = parsed.mmdd;
    raw(keep).label = parsed.label;
    raw(keep).label_ok = parsed.label_ok;
    raw(keep).parse_detail = parsed.detail;
end
if keep == 0
    raw = empty_raw();
end
end

function [exp_results, raw_results] = match_folders(experiments, raw)
exp_results = empty_result();
raw_ok = false(1, numel(raw));

for i = 1:numel(experiments)
    item = experiments(i);
    result = blank_result('experiment', item.row, session_hint(item), item);
    matches = matching_raw(item, raw);
    if isempty(item.animal) || isempty(item.recording_folder) ...
            || isempty(item.probe) || isempty(item.recording_number)
        result.status = 'bad path';
        result.detail = 'Animal, date, NPX recording folder, probe, and recording number all have to be filled in.';
    elseif isempty(item.date.mmdd)
        result.status = 'bad path';
        result.detail = 'Could not read the date folder.';
    elseif isempty(matches)
        result.status = 'no session';
        result.detail = sprintf('No raw-folder session starts with "%s %s".', ...
            item.animal, item.date.mmdd);
    else
        hit = [];
        label_hit = [];
        first_reason = '';
        for k = 1:numel(matches)
            j = matches(k);
            [ok, reason] = path_matches(raw(j).segments, item);
            if ok && raw(j).label_ok
                raw_ok(j) = true;
                if isempty(hit)
                    hit = j;
                end
            elseif ok && isempty(label_hit)
                label_hit = j;
            elseif isempty(first_reason)
                first_reason = reason;
            end
        end
        if ~isempty(hit)
            result.status = 'ok';
            result.session_name = raw(hit).session;
            result.actual_path = raw(hit).path;
            result.detail = '';
        elseif ~isempty(label_hit)
            result.status = 'bad session name';
            result.session_name = raw(label_hit).session;
            result.actual_path = raw(label_hit).path;
            result.detail = raw(label_hit).parse_detail;
        else
            result.status = 'bad path';
            result.session_name = raw(matches(1)).session;
            result.actual_path = raw(matches(1)).path;
            result.detail = first_reason;
        end
    end
    exp_results(end + 1) = result; %#ok<AGROW>
end

raw_results = empty_result();
for j = 1:numel(raw)
    if raw_ok(j)
        continue
    end
    item_matches = experiment_matches_raw(raw(j), experiments);
    result = blank_result('raw', raw(j).row, raw(j).session, []);
    result.actual_path = raw(j).path;
    if ~raw(j).parse_ok
        result.status = 'bad session name';
        result.detail = raw(j).parse_detail;
    elseif isempty(item_matches)
        result.status = 'unlisted';
        result.detail = sprintf('Session "%s" is not in the experiment info.', raw(j).session);
    elseif ~raw(j).label_ok
        result.status = 'bad session name';
        result.detail = raw(j).parse_detail;
    else
        result.status = 'bad path';
        [~, reason] = path_matches(raw(j).segments, experiments(item_matches(1)));
        result.detail = reason;
    end
    raw_results(end + 1) = result; %#ok<AGROW>
end
end

function matches = matching_raw(item, raw)
matches = [];
for j = 1:numel(raw)
    if raw(j).parse_ok && strcmpi(raw(j).animal, item.animal) ...
            && strcmp(raw(j).mmdd, item.date.mmdd)
        matches(end + 1) = j; %#ok<AGROW>
    end
end
end

function matches = experiment_matches_raw(raw_row, experiments)
matches = [];
if ~raw_row.parse_ok
    return
end
for i = 1:numel(experiments)
    if strcmpi(experiments(i).animal, raw_row.animal) ...
            && strcmp(experiments(i).date.mmdd, raw_row.mmdd)
        matches(end + 1) = i; %#ok<AGROW>
    end
end
end

function [ok, reason] = path_matches(segments, item)
ok = false;
reason = '';
date_at = 0;
for s = 1:numel(segments)
    if segment_matches_date(segments{s}, item.date)
        date_at = s;
        break
    end
end
if date_at == 0
    reason = sprintf('Path has no date folder for %s.', date_label(item.date));
    return
end
if date_at < numel(segments) && any(strcmpi(segments(date_at + 1:end), item.animal))
    reason = sprintf('Animal code %s appears after the date folder.', item.animal);
    return
end
folder_parts = split_path(item.recording_folder);
if isempty(folder_parts)
    folder_parts = {item.recording_folder};
end
folder_at = 0;
last_needed = numel(segments) - numel(folder_parts) + 1;
for s = date_at + 1:max(date_at, last_needed)
    if s < 1 || s + numel(folder_parts) - 1 > numel(segments)
        continue
    end
    if all(strcmpi(segments(s:s + numel(folder_parts) - 1), folder_parts))
        folder_at = s + numel(folder_parts) - 1;
        break
    end
end
if folder_at == 0
    reason = sprintf('Path has no NPX recording folder "%s" after the date.', ...
        item.recording_folder);
    return
end
probe_at = find_token_after(segments, folder_at, item.probe_tokens);
if probe_at == 0
    reason = sprintf('Path has no probe folder for "%s" (%s) after "%s".', ...
        item.probe, strjoin(item.probe_tokens, ' or '), item.recording_folder);
    return
end
number_at = find_token_after(segments, folder_at, item.recording_tokens);
if number_at == 0
    reason = sprintf('Path has no recording-number folder for "%s" (%s) after "%s".', ...
        item.recording_number, strjoin(item.recording_tokens, ' or '), item.recording_folder);
    return
end
ok = true;
end

function index = find_token_after(segments, folder_at, tokens)
index = 0;
for s = folder_at + 1:numel(segments)
    for t = 1:numel(tokens)
        if token_in_segment(segments{s}, tokens{t})
            index = s;
            return
        end
    end
end
end

function parsed = parse_session_name(session)
parsed.ok = false;
parsed.animal = '';
parsed.mmdd = '';
parsed.label = '';
parsed.label_ok = false;
parsed.detail = '';
tokens = regexp(strtrim(session), '\s+', 'split');
tokens = tokens(~cellfun('isempty', tokens));
if numel(tokens) ~= 3
    parsed.detail = sprintf(['Session name "%s" is not ' ...
        '"{animal} {MMDD} {Single|Dual|Single2}".'], session);
    return
end
parsed.animal = tokens{1};
parsed.label = tokens{3};
date_info = interpret_date(tokens{2});
parsed.mmdd = date_info.mmdd;
if isempty(parsed.mmdd)
    parsed.detail = sprintf('Session name "%s" has no MMDD date.', session);
    return
end
parsed.label_ok = label_follows_animal(parsed.animal, parsed.label);
parsed.ok = true;
if ~parsed.label_ok
    parsed.detail = sprintf('Session "%s" uses %s for animal %s.', ...
        session, parsed.label, parsed.animal);
end
end

function tf = label_follows_animal(animal, label)
tf = false;
token = upper(canonical_text(animal));
letter = regexp(token, '^B\d+([DS])\d+$', 'tokens', 'once');
label_key = lower(canonical_text(label));
if isempty(letter)
    tf = any(strcmp(label_key, {'single', 'single2', 'dual'}));
    return
end
if strcmp(letter{1}, 'D')
    tf = strcmp(label_key, 'dual');
else
    tf = any(strcmp(label_key, {'single', 'single2'}));
end
end

function info = interpret_date(value)
info.mmdd = '';
info.yyyymmdd = '';
info.exact = {};
if is_datetime_value(value)
    info.yyyymmdd = datestr(value, 'yyyymmdd');
    info.mmdd = datestr(value, 'mmdd');
    info.exact = {info.yyyymmdd, datestr(value, 'yyyy-mm-dd')};
    return
end
number = numeric_value(value);
if ~isempty(number) && number == round(number)
    if number >= 1e7 && number < 3e7
        info = date_from_yyyymmdd(sprintf('%.0f', number));
        return
    end
    if number >= 20000 && number <= 80000
        info.yyyymmdd = datestr(number + datenum(1899, 12, 30), 'yyyymmdd');
        info.mmdd = info.yyyymmdd(5:8);
        info.exact = {info.yyyymmdd, datestr(number + datenum(1899, 12, 30), 'yyyy-mm-dd')};
        return
    end
    if number >= 101 && number <= 1231
        info.mmdd = sprintf('%04.0f', number);
        info.exact = {info.mmdd};
        return
    end
end
text = canonical_text(value);
if isempty(text)
    return
end
digits = regexprep(text, '\D', '');
info.exact = {text};
if numel(digits) >= 8
    info = date_from_yyyymmdd(digits(1:8));
    if ~any(strcmpi(info.exact, text))
        info.exact = [{text}, info.exact];
    end
elseif numel(digits) == 3 || numel(digits) == 4
    info.mmdd = sprintf('%04d', str2double(digits));
    info.exact = unique_text([{text}, {info.mmdd}]);
end
end

function info = date_from_yyyymmdd(digits)
info.yyyymmdd = digits(1:8);
info.mmdd = digits(5:8);
info.exact = {info.yyyymmdd, [digits(1:4) '-' digits(5:6) '-' digits(7:8)]};
end

function tf = segment_matches_date(segment, info)
tf = false;
for k = 1:numel(info.exact)
    if strcmpi(segment, info.exact{k})
        tf = true;
        return
    end
end
if ~isempty(info.yyyymmdd)
    tf = ~isempty(regexpi(segment, ['(^|[^0-9])' info.yyyymmdd '([^0-9]|$)'], 'once'));
    if tf
        return
    end
end
if numel(info.mmdd) == 4
    tf = ~isempty(regexpi(segment, ...
        ['(^|[^0-9])20[0-9]{2}' info.mmdd '([^0-9]|$)'], 'once'));
    if tf
        return
    end
    tf = ~isempty(regexpi(segment, ...
        ['(^|[^0-9])20[0-9]{2}-' info.mmdd(1:2) '-' info.mmdd(3:4) '([^0-9]|$)'], 'once'));
end
end

function tokens = probe_tokens(probe)
tokens = {};
if isempty(probe)
    return
end
tokens = {probe};
number = str2double(probe);
if isnan(number) || number ~= round(number) || number < 0 || number > 16
    return
end
letters = 'ABCDEFGHIJKLMNOP';
if number >= 1 && number <= numel(letters)
    tokens{end + 1} = ['Probe' letters(number)]; %#ok<AGROW>
    tokens{end + 1} = ['imec' num2str(number - 1)]; %#ok<AGROW>
elseif number == 0
    tokens{end + 1} = 'imec0'; %#ok<AGROW>
end
end

function tokens = recording_tokens(recording_number)
tokens = {};
if isempty(recording_number)
    return
end
tokens = {recording_number};
number = str2double(recording_number);
if isnan(number) || number ~= round(number) || number < 0
    return
end
tokens{end + 1} = ['recording' num2str(number)]; %#ok<AGROW>
end

function tf = token_in_segment(segment, token)
if isempty(token)
    tf = false;
    return
end
pattern = ['(^|[^A-Za-z0-9])' escape_regex(token) '([^A-Za-z0-9]|$)'];
tf = ~isempty(regexpi(segment, pattern, 'once'));
end

function text = escape_regex(text)
text = regexprep(text, '([\\^$.|?*+(){}\[\]])', '\\$1');
end

function segments = split_path(path_str)
segments = regexp(strtrim(path_str), '[\\/]+', 'split');
segments = segments(~cellfun('isempty', segments));
end

function text = date_label(info)
forms = info.exact;
if ~isempty(info.mmdd)
    forms{end + 1} = info.mmdd;
end
if ~isempty(info.yyyymmdd)
    forms{end + 1} = info.yyyymmdd;
end
forms = unique_text(forms);
text = strjoin(forms, ' / ');
end

function text = session_hint(item)
if isempty(item.date.mmdd)
    text = item.animal;
else
    text = [item.animal ' ' item.date.mmdd];
end
end

function key = header_key(header)
key = lower(regexprep(canonical_text(header), '[^A-Za-z0-9]', ''));
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
if is_datetime_value(value)
    text = datestr(value, 'yyyy-mm-dd');
    return
end
number = numeric_value(value);
if ~isempty(number)
    if ~isfinite(number)
        return
    end
    if number == round(number)
        text = sprintf('%.0f', number);
    else
        text = sprintf('%g', number);
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

function number = numeric_value(value)
number = [];
if isnumeric(value) && isscalar(value)
    number = double(value);
end
end

function tf = is_datetime_value(value)
tf = false;
if exist('isdatetime', 'builtin') && isdatetime(value) && numel(value) == 1
    tf = true;
end
end

function values = unique_text(values)
keep = {};
for k = 1:numel(values)
    if isempty(values{k})
        continue
    end
    if ~any(strcmpi(keep, values{k}))
        keep{end + 1} = values{k}; %#ok<AGROW>
    end
end
values = keep;
end

function notes = notes_from(info_note, folder_note)
notes = {};
if ~isempty(info_note)
    notes{end + 1} = info_note; %#ok<AGROW>
end
if ~isempty(folder_note)
    notes{end + 1} = folder_note; %#ok<AGROW>
end
end

function item = empty_experiment()
item = struct('row', {}, 'animal', {}, 'date', {}, 'recording_folder', {}, ...
    'probe', {}, 'probe_tokens', {}, 'recording_number', {}, 'recording_tokens', {});
end

function item = empty_raw()
item = struct('row', {}, 'session', {}, 'path', {}, 'segments', {}, ...
    'parse_ok', {}, 'animal', {}, 'mmdd', {}, 'label', {}, 'label_ok', {}, ...
    'parse_detail', {});
end

function item = empty_result()
item = struct('source', {}, 'sheet_row', {}, 'session_name', {}, ...
    'animal_code', {}, 'date_mmdd', {}, 'recording_folder', {}, ...
    'probe', {}, 'recording_number', {}, 'actual_path', {}, ...
    'status', {}, 'detail', {});
end

function result = blank_result(source, sheet_row, session_name, item)
result.source = source;
result.sheet_row = sheet_row;
result.session_name = session_name;
result.animal_code = '';
result.date_mmdd = '';
result.recording_folder = '';
result.probe = '';
result.recording_number = '';
result.actual_path = '';
result.status = '';
result.detail = '';
if ~isempty(item)
    result.animal_code = item.animal;
    result.date_mmdd = item.date.mmdd;
    result.recording_folder = item.recording_folder;
    result.probe = item.probe;
    result.recording_number = item.recording_number;
end
end
