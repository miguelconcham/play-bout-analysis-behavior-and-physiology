% Check_NPX_folders_and_load_waveforms
%
% One script, three sections. Run this file from the Editor, or type its name
% in the command window. Later sections use the variables created above them.
% Supporting functions are at the bottom of this file.
%
% 1. Read Experiment_Info_GoodLooking.xlsx and NPX_raw_data_folders.xlsx.
%    Print every place the session name and the continuous-file path disagree.
% 2. Open each session folder under NPX raw data, pick one random good unit,
%    and store it in selected_units:
%       cluster_id, session_id, animal_name, spike_times (samples), channel
%    channel is cluster_info.ch + 1 (that table is 0-based).
% 3. For each selected unit, load a subset of spike_times from the continuous
%    file listed for that session in NPX_raw_data_folders.xlsx. The load uses
%    load_waveforms_from_continuous with five channels on each side of the
%    neuron. A figure shows that stack and the autocorrelogram.
%
% Section 1 only prints. It does not stop the later sections, so read the
% dissimilarities before treating a waveform as proof of the folder assignment.
%
% Area_limits_GoodLooking.xlsx is the depth table. This script reads
% Experiment_Info_GoodLooking.xlsx.

%% Paths
experiment_info_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\Experiment_Info_GoodLooking.xlsx'];
raw_folders_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\NPX data\' ...
    'NPX_raw_data_folders.xlsx'];
npx_raw_data = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\NPX data\NPX raw data'];

repo_root = fileparts(fileparts(mfilename('fullpath')));
if isempty(repo_root)
    repo_root = fileparts(experiment_info_file);
end
run(fullfile(repo_root, 'add_repo_paths.m'));

%% 1. Read both workbooks and report dissimilarities
experiment_sheet = read_location_sheet(experiment_info_file);
raw_sheet = read_location_sheet(raw_folders_file);
[dissimilarities, raw_locations, notes, n_experiment] = compare_location_sheets( ...
    experiment_sheet, raw_sheet);

for note_index = 1:numel(notes)
    fprintf('Note: %s\n', notes{note_index});
end
fprintf('Experiment info: %d rows. Raw-folder sheet: %d rows.\n', ...
    n_experiment, numel(raw_locations));
if isempty(dissimilarities)
    fprintf('No dissimilarities.\n');
else
    fprintf('%d dissimilarities.\n', numel(dissimilarities));
    for diff_index = 1:numel(dissimilarities)
        fprintf('%s\n  %s\n', dissimilarities(diff_index).where, ...
            dissimilarities(diff_index).detail);
    end
end

%% 2. One random good unit per NPX session folder
rng('shuffle');
seed_state = rng;
fprintf('Random seed: %d\n', seed_state.Seed);

selected_units = struct('cluster_id', {}, 'session_id', {}, ...
    'animal_name', {}, 'spike_times', {}, 'channel', {});

if exist(npx_raw_data, 'dir') ~= 7
    fprintf('NPX raw data folder was not found:\n  %s\n', npx_raw_data);
else
    listing = dir(npx_raw_data);
    for folder_index = 1:numel(listing)
        if ~listing(folder_index).isdir
            continue
        end
        session_id = listing(folder_index).name;
        if strcmp(session_id, '.') || strcmp(session_id, '..')
            continue
        end
        session_dir = fullfile(npx_raw_data, session_id);
        name_tokens = regexp(strtrim(session_id), '\s+', 'split');
        name_tokens = name_tokens(~cellfun('isempty', name_tokens));
        if isempty(name_tokens)
            animal_name = session_id;
        else
            animal_name = name_tokens{1};
        end
        units = good_units_in_session(session_dir);
        if isempty(units)
            fprintf('%s  no good unit\n', session_id);
            continue
        end
        pick = units(randi(numel(units)));
        selected_units(end + 1).cluster_id = pick.cluster_id; %#ok<SAGROW>
        selected_units(end).session_id = session_id;
        selected_units(end).animal_name = animal_name;
        selected_units(end).spike_times = pick.spike_times;
        selected_units(end).channel = pick.channel;
        fprintf('%s  cluster %g  channel %g  %d spikes\n', ...
            session_id, pick.cluster_id, pick.channel, numel(pick.spike_times));
    end
end
fprintf('Selected %d units.\n', numel(selected_units));

%% 3. Load those spikes from the continuous files in NPX_raw_data_folders.xlsx
channel_range = 5;
n_waveforms = 40;
fs = 30000;
acg_window = [-0.05 0.05];
acg_bin = 0.001;
n_acg_spikes = 3000;

for unit_index = 1:numel(selected_units)
    session_id = selected_units(unit_index).session_id;
    paths = paths_for_session(raw_locations, session_id);
    if isempty(paths)
        fprintf('%s  no continuous path in NPX_raw_data_folders.xlsx\n', session_id);
        continue
    end
    n_draw = min(n_waveforms, numel(selected_units(unit_index).spike_times));
    draw_index = randperm(numel(selected_units(unit_index).spike_times), n_draw);
    waveform_samples = selected_units(unit_index).spike_times(draw_index);
    stored = false;
    for path_index = 1:numel(paths)
        continuous_file = continuous_file_from_location(paths{path_index});
        fprintf('%s  cluster %g  %s\n', session_id, ...
            selected_units(unit_index).cluster_id, continuous_file);
        try
            [waveforms, channels_kept, time_samples] = load_waveforms_from_continuous( ...
                continuous_file, waveform_samples, ...
                selected_units(unit_index).channel, channel_range);
        catch load_error
            fprintf('  %s\n', load_error.message);
            continue
        end
        if ~stored
            selected_units(unit_index).continuous_file = continuous_file;
            selected_units(unit_index).waveform_samples = waveform_samples;
            selected_units(unit_index).waveforms = waveforms;
            selected_units(unit_index).channels_kept = channels_kept;
            selected_units(unit_index).time_samples = time_samples;
            stored = true;
        end
        try
            show_waveform_and_acg(selected_units(unit_index), continuous_file, ...
                waveforms, channels_kept, time_samples, fs, acg_window, ...
                acg_bin, n_acg_spikes);
        catch plot_error
            fprintf('  waveforms were loaded; the figure failed: %s\n', ...
                plot_error.message);
        end
    end
end

%% Local functions

function [headers, rows] = split_sheet(sheet)
headers = sheet(1, :);
if size(sheet, 1) < 2
    rows = cell(0, size(sheet, 2));
else
    rows = sheet(2:end, :);
end
end

function sheet = read_location_sheet(path)
if exist(path, 'file') ~= 2
    error('Check_NPX_folders_and_load_waveforms:file', ...
        'Could not find %s', path);
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
elseif exist('isdatetime', 'builtin') && isdatetime(column)
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

function [dissimilarities, raw_locations, notes, n_experiment] = compare_location_sheets(experiment_sheet, raw_sheet)
[info_headers, info_rows] = split_sheet(experiment_sheet);
[folder_headers, folder_rows] = split_sheet(raw_sheet);
if looks_like_area_limits(info_headers) || looks_like_area_limits(folder_headers)
    error('Check_NPX_folders_and_load_waveforms:areaLimits', ...
        ['This sheet is Area_limits_GoodLooking (depth limits). ' ...
         'The experiment sheet is Experiment_Info_GoodLooking.xlsx.']);
end

[info_col, info_note] = experiment_columns(info_headers);
[folder_col, folder_note] = raw_folder_columns(folder_headers, folder_rows);
notes = {};
if ~isempty(info_note)
    notes{end + 1} = info_note; %#ok<AGROW>
end
if ~isempty(folder_note)
    notes{end + 1} = folder_note; %#ok<AGROW>
end

experiments = struct('row', {}, 'animal', {}, 'mmdd', {}, 'yyyymmdd', {}, ...
    'exact', {}, 'recording_folder', {}, 'probe', {}, 'probe_tokens', {}, ...
    'recording_number', {}, 'recording_tokens', {});
for row_index = 1:size(info_rows, 1)
    animal = canonical_text(info_rows{row_index, info_col.animal});
    date_info = interpret_date(info_rows{row_index, info_col.date});
    recording_folder = canonical_text(info_rows{row_index, info_col.recording_folder});
    probe = canonical_text(info_rows{row_index, info_col.probe});
    recording_number = canonical_text(info_rows{row_index, info_col.recording_number});
    if isempty(animal) && isempty(date_info.mmdd) && isempty(recording_folder) ...
            && isempty(probe) && isempty(recording_number)
        continue
    end
    experiments(end + 1).row = row_index + 1; %#ok<AGROW>
    experiments(end).animal = animal;
    experiments(end).mmdd = date_info.mmdd;
    experiments(end).yyyymmdd = date_info.yyyymmdd;
    experiments(end).exact = date_info.exact;
    experiments(end).recording_folder = recording_folder;
    experiments(end).probe = probe;
    experiments(end).probe_tokens = probe_tokens(probe);
    experiments(end).recording_number = recording_number;
    experiments(end).recording_tokens = recording_tokens(recording_number);
end

raw_locations = struct('row', {}, 'session', {}, 'path', {}, 'segments', {}, ...
    'parse_ok', {}, 'animal', {}, 'mmdd', {}, 'label_ok', {}, 'detail', {});
for row_index = 1:size(folder_rows, 1)
    session = canonical_text(folder_rows{row_index, folder_col.session});
    path = canonical_text(folder_rows{row_index, folder_col.path});
    if isempty(session) && isempty(path)
        continue
    end
    parsed = parse_session_name(session);
    raw_locations(end + 1).row = row_index + 1; %#ok<AGROW>
    raw_locations(end).session = session;
    raw_locations(end).path = path;
    raw_locations(end).segments = split_path(path);
    raw_locations(end).parse_ok = parsed.ok;
    raw_locations(end).animal = parsed.animal;
    raw_locations(end).mmdd = parsed.mmdd;
    raw_locations(end).label_ok = parsed.label_ok;
    raw_locations(end).detail = parsed.detail;
end

dissimilarities = struct('where', {}, 'detail', {});
n_experiment = numel(experiments);
if isempty(experiments) || isempty(raw_locations)
    dissimilarities = push_difference(dissimilarities, 'sheets', ...
        'The experiment sheet or the raw-folder sheet has no data rows.');
end

raw_ok = false(1, numel(raw_locations));
for item_index = 1:numel(experiments)
    item = experiments(item_index);
    where = sprintf('Experiment_Info row %d (%s)', item.row, session_hint(item));
    matches = matching_raw(item, raw_locations);
    if isempty(item.animal) || isempty(item.recording_folder) ...
            || isempty(item.probe) || isempty(item.recording_number)
        dissimilarities = push_difference(dissimilarities, where, ...
            'Animal, date, NPX recording folder, probe, and recording number all have to be filled in.');
    elseif isempty(item.mmdd)
        dissimilarities = push_difference(dissimilarities, where, ...
            'Could not read the date folder.');
    elseif isempty(matches)
        dissimilarities = push_difference(dissimilarities, where, ...
            sprintf('No raw-folder session starts with "%s %s".', item.animal, item.mmdd));
    else
        hit = [];
        label_hit = [];
        first_reason = '';
        for match_index = 1:numel(matches)
            raw_index = matches(match_index);
            [ok, reason] = path_matches(raw_locations(raw_index).segments, item);
            if ok && raw_locations(raw_index).label_ok
                raw_ok(raw_index) = true;
                if isempty(hit)
                    hit = raw_index;
                end
            elseif ok && isempty(label_hit)
                label_hit = raw_index;
            elseif isempty(first_reason)
                first_reason = reason;
            end
        end
        if isempty(hit) && ~isempty(label_hit)
            dissimilarities = push_difference(dissimilarities, where, ...
                with_path(raw_locations(label_hit).detail, raw_locations(label_hit).path));
        elseif isempty(hit)
            dissimilarities = push_difference(dissimilarities, where, ...
                with_path(first_reason, raw_locations(matches(1)).path));
        end
    end
end

for raw_index = 1:numel(raw_locations)
    if raw_ok(raw_index)
        continue
    end
    raw_row = raw_locations(raw_index);
    where = sprintf('NPX_raw_data_folders row %d (%s)', raw_row.row, raw_row.session);
    item_matches = experiment_matches_raw(raw_row, experiments);
    if ~raw_row.parse_ok || ~raw_row.label_ok
        dissimilarities = push_difference(dissimilarities, where, ...
            with_path(raw_row.detail, raw_row.path));
    elseif isempty(item_matches)
        dissimilarities = push_difference(dissimilarities, where, ...
            with_path(sprintf('Session "%s" is not in the experiment info.', raw_row.session), ...
            raw_row.path));
    else
        [~, reason] = path_matches(raw_row.segments, experiments(item_matches(1)));
        dissimilarities = push_difference(dissimilarities, where, ...
            with_path(reason, raw_row.path));
    end
end
end

function dissimilarities = push_difference(dissimilarities, where, detail)
dissimilarities(end + 1).where = where; %#ok<AGROW>
dissimilarities(end).detail = detail;
end

function detail = with_path(detail, path)
if ~isempty(path)
    detail = sprintf('%s\n  %s', detail, path);
end
end

function tf = looks_like_area_limits(headers)
keys = cell(1, numel(headers));
for index = 1:numel(headers)
    keys{index} = header_key(headers{index});
end
tf = any(strcmp(keys, 'depthstart')) && any(strcmp(keys, 'animalname'));
end

function [col, note] = experiment_columns(headers)
roles = cell(1, numel(headers));
for index = 1:numel(headers)
    roles{index} = experiment_role(headers{index});
end
col.animal = find(strcmp(roles, 'animal'), 1);
col.date = find(strcmp(roles, 'date'), 1);
col.recording_folder = find(strcmp(roles, 'recording_folder'), 1);
col.probe = find(strcmp(roles, 'probe'), 1);
col.recording_number = find(strcmp(roles, 'recording_number'), 1);
note = '';
if isempty(col.animal) || isempty(col.date) || isempty(col.recording_folder) ...
        || isempty(col.probe) || isempty(col.recording_number)
    if numel(headers) < 5
        error('Check_NPX_folders_and_load_waveforms:columns', ...
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
        col.session = [];
        for index = 1:numel(headers)
            if strcmp(raw_role(headers{index}), 'session')
                col.session = index;
                break
            end
        end
        if isempty(col.session)
            error('Check_NPX_folders_and_load_waveforms:columns', ...
                'The first raw-folder column looks like a path, and no session column was found.');
        end
    end
end
col.path = [];
for index = 1:numel(headers)
    if strcmp(raw_role(headers{index}), 'path')
        col.path = index;
        break
    end
end
if isempty(col.path)
    col.path = first_path_column(rows, col.session);
end
if isempty(col.path)
    error('Check_NPX_folders_and_load_waveforms:columns', ...
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

function matches = matching_raw(item, raw_locations)
matches = [];
for raw_index = 1:numel(raw_locations)
    if raw_locations(raw_index).parse_ok ...
            && strcmpi(raw_locations(raw_index).animal, item.animal) ...
            && strcmp(raw_locations(raw_index).mmdd, item.mmdd)
        matches(end + 1) = raw_index; %#ok<AGROW>
    end
end
end

function matches = experiment_matches_raw(raw_row, experiments)
matches = [];
if ~raw_row.parse_ok
    return
end
for item_index = 1:numel(experiments)
    if strcmpi(experiments(item_index).animal, raw_row.animal) ...
            && strcmp(experiments(item_index).mmdd, raw_row.mmdd)
        matches(end + 1) = item_index; %#ok<AGROW>
    end
end
end

function [ok, reason] = path_matches(segments, item)
ok = false;
reason = '';
date_at = 0;
for segment_index = 1:numel(segments)
    if segment_matches_date(segments{segment_index}, item)
        date_at = segment_index;
        break
    end
end
if date_at == 0
    reason = sprintf('Path has no date folder for %s.', date_label(item));
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
for segment_index = date_at + 1:max(date_at, last_needed)
    last_part = segment_index + numel(folder_parts) - 1;
    if segment_index < 1 || last_part > numel(segments)
        continue
    end
    if all(strcmpi(segments(segment_index:last_part), folder_parts))
        folder_at = last_part;
        break
    end
end
if folder_at == 0
    reason = sprintf('Path has no NPX recording folder "%s" after the date.', ...
        item.recording_folder);
    return
end
if find_token_after(segments, folder_at, item.probe_tokens) == 0
    reason = sprintf('Path has no probe folder for "%s" (%s) after "%s".', ...
        item.probe, strjoin(item.probe_tokens, ' or '), item.recording_folder);
    return
end
if find_token_after(segments, folder_at, item.recording_tokens) == 0
    reason = sprintf('Path has no recording-number folder for "%s" (%s) after "%s".', ...
        item.recording_number, strjoin(item.recording_tokens, ' or '), item.recording_folder);
    return
end
ok = true;
end

function index = find_token_after(segments, folder_at, tokens)
index = 0;
for segment_index = folder_at + 1:numel(segments)
    for token_index = 1:numel(tokens)
        if token_in_segment(segments{segment_index}, tokens{token_index})
            index = segment_index;
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
if exist('isdatetime', 'builtin') && isdatetime(value) && numel(value) == 1
    info.yyyymmdd = datestr(value, 'yyyymmdd');
    info.mmdd = datestr(value, 'mmdd');
    info.exact = {info.yyyymmdd, datestr(value, 'yyyy-mm-dd')};
    return
end
number = [];
if isnumeric(value) && isscalar(value)
    number = double(value);
end
if ~isempty(number) && number == round(number)
    if number >= 1e7 && number < 3e7
        info = date_from_yyyymmdd(sprintf('%.0f', number));
        return
    end
    if number >= 20000 && number <= 80000
        stamp = number + datenum(1899, 12, 30);
        info.yyyymmdd = datestr(stamp, 'yyyymmdd');
        info.mmdd = info.yyyymmdd(5:8);
        info.exact = {info.yyyymmdd, datestr(stamp, 'yyyy-mm-dd')};
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
for index = 1:numel(info.exact)
    if strcmpi(segment, info.exact{index})
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
segments = segments(:)';
end

function text = date_label(info)
forms = info.exact;
if ~isempty(info.mmdd)
    forms{end + 1} = info.mmdd; %#ok<AGROW>
end
if ~isempty(info.yyyymmdd)
    forms{end + 1} = info.yyyymmdd; %#ok<AGROW>
end
text = strjoin(unique_text(forms), ' / ');
end

function text = session_hint(item)
if isempty(item.mmdd)
    text = item.animal;
else
    text = [item.animal ' ' item.mmdd];
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
if exist('isdatetime', 'builtin') && isdatetime(value) && numel(value) == 1
    text = datestr(value, 'yyyy-mm-dd');
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

function values = unique_text(values)
keep = {};
for index = 1:numel(values)
    if isempty(values{index})
        continue
    end
    if ~any(strcmpi(keep, values{index}))
        keep{end + 1} = values{index}; %#ok<AGROW>
    end
end
values = keep;
end

function paths = paths_for_session(raw_locations, session_id)
paths = {};
for raw_index = 1:numel(raw_locations)
    if ~strcmp(raw_locations(raw_index).session, session_id)
        continue
    end
    candidate = strtrim(raw_locations(raw_index).path);
    if isempty(candidate) || any(strcmp(paths, candidate))
        continue
    end
    paths{end + 1} = candidate; %#ok<AGROW>
end
end

function continuous_file = continuous_file_from_location(location)
continuous_file = location;
if exist(location, 'dir') == 7
    continuous_file = fullfile(location, 'continuous.dat');
end
end

function units = good_units_in_session(session_dir)
units = struct('cluster_id', {}, 'channel', {}, 'spike_times', {});
info_files = dir(fullfile(session_dir, 'cluster_info*.tsv'));
for file_index = 1:numel(info_files)
    info_name = info_files(file_index).name;
    [~, info_base] = fileparts(info_name);
    suffix = '';
    prefix = 'cluster_info';
    if numel(info_base) > numel(prefix)
        suffix = info_base(numel(prefix) + 1:end);
    end
    times_file = fullfile(session_dir, ['spike_times' suffix '.npy']);
    cluster_file = fullfile(session_dir, ['spike_clusters' suffix '.npy']);
    if exist(times_file, 'file') ~= 2 || exist(cluster_file, 'file') ~= 2
        fprintf('  missing spike_times%s.npy or spike_clusters%s.npy\n', suffix, suffix);
        continue
    end
    cluster_info = readtable(fullfile(session_dir, info_name), ...
        'FileType', 'text', 'Delimiter', '\t');
    spike_times = spike_times_to_samples(readNPY(times_file));
    spike_clusters = double(readNPY(cluster_file));
    spike_clusters = spike_clusters(:);
    ids = table_column(cluster_info, {'id', 'cluster_id'});
    channels = table_column(cluster_info, {'ch', 'channel'});
    groups = table_column(cluster_info, {'group', 'KSLabel'});
    for unit_index = 1:numel(ids)
        if ~strcmp(group_label(groups, unit_index), 'good')
            continue
        end
        cluster_id = scalar_number(ids, unit_index);
        channel = scalar_number(channels, unit_index) + 1;
        if ~isfinite(cluster_id) || ~isfinite(channel) || channel ~= round(channel) || channel < 1
            continue
        end
        samples = spike_times(spike_clusters == cluster_id);
        samples = samples(isfinite(samples) & samples == round(samples) & samples >= 1);
        if isempty(samples)
            continue
        end
        units(end + 1).cluster_id = cluster_id; %#ok<AGROW>
        units(end).channel = channel;
        units(end).spike_times = samples(:);
    end
end
end

function spike_samples = spike_times_to_samples(spike_times)
spike_samples = double(spike_times(:));
if isempty(spike_samples)
    return
end
if any(~isfinite(spike_samples)) || any(spike_samples ~= round(spike_samples))
    error('Check_NPX_folders_and_load_waveforms:spikeTime', ...
        'spike_times must be integer sample indices.');
end
if min(spike_samples) == 0
    spike_samples = spike_samples + 1;
elseif any(spike_samples < 0)
    error('Check_NPX_folders_and_load_waveforms:spikeTime', ...
        'spike_times must be non-negative sample indices.');
end
spike_samples = round(spike_samples);
end

function values = table_column(table_data, names)
vars = table_data.Properties.VariableNames;
keys = cell(1, numel(vars));
for index = 1:numel(vars)
    keys{index} = lower(regexprep(vars{index}, '[^A-Za-z0-9]', ''));
end
for name_index = 1:numel(names)
    key = lower(regexprep(names{name_index}, '[^A-Za-z0-9]', ''));
    hit = find(strcmp(keys, key), 1);
    if ~isempty(hit)
        values = table_data.(vars{hit});
        return
    end
end
error('Check_NPX_folders_and_load_waveforms:column', ...
    'cluster_info needs a column named %s.', strjoin(names, ' or '));
end

function label = group_label(groups, index)
if iscell(groups)
    value = groups{index};
elseif ischar(groups)
    value = strtrim(groups(index, :));
elseif exist('isstring', 'builtin') && isstring(groups)
    value = char(groups(index));
elseif exist('iscategorical', 'builtin') && iscategorical(groups)
    value = char(groups(index));
else
    value = groups(index);
end
if exist('isstring', 'builtin') && isstring(value)
    value = char(value);
end
if isnumeric(value)
    label = '';
    return
end
label = lower(strtrim(value));
end

function number = scalar_number(values, index)
if iscell(values)
    number = double(values{index});
else
    number = double(values(index));
end
end

function show_waveform_and_acg(unit, continuous_file, waveforms, channels_kept, ...
    time_samples, fs, acg_window, acg_bin, n_acg_spikes)
[parent_dir, file_name, file_ext] = fileparts(continuous_file);
[~, parent_name] = fileparts(parent_dir);
figure_name = sprintf('%s  cluster %g  ch %g  %s', unit.session_id, ...
    unit.cluster_id, unit.channel, [parent_name filesep file_name file_ext]);
figure('Color', 'w', 'Name', figure_name);

subplot(1, 2, 1);
n_channels = size(waveforms, 2);
n_times = size(waveforms, 3);
mean_wf = nan(n_channels, n_times);
for channel_index = 1:n_channels
    traces = reshape(waveforms(:, channel_index, :), size(waveforms, 1), n_times);
    mean_wf(channel_index, :) = mean_omitnan(traces);
end
peak = max(abs(mean_wf(:)));
if isempty(peak) || ~isfinite(peak) || peak == 0
    peak = 1;
end
gap = 2.2 * peak;
time_ms = time_samples(:)' / fs * 1000;
hold on
y_ticks = (0:n_channels - 1) * gap;
for channel_index = 1:n_channels
    traces = reshape(waveforms(:, channel_index, :), size(waveforms, 1), n_times);
    plot(time_ms, traces' + y_ticks(channel_index), 'Color', [0.75 0.75 0.75]);
end
for channel_index = 1:n_channels
    mean_color = [0 0 0];
    mean_width = 1.5;
    if channels_kept(channel_index) == unit.channel
        mean_color = [0.85 0.1 0.1];
        mean_width = 2.5;
    end
    plot(time_ms, mean_wf(channel_index, :) + y_ticks(channel_index), ...
        'Color', mean_color, 'LineWidth', mean_width);
end
tick_labels = cell(1, n_channels);
for channel_index = 1:n_channels
    tick_labels{channel_index} = sprintf('%g', channels_kept(channel_index));
end
set(gca, 'YTick', y_ticks, 'YTickLabel', tick_labels);
xlabel('Time (ms)')
ylabel('Channel')
title('Waveform')
box off

subplot(1, 2, 2);
spike_seconds = unit.spike_times / fs;
if numel(spike_seconds) > n_acg_spikes
    spike_seconds = spike_seconds(randperm(numel(spike_seconds), n_acg_spikes));
end
lags = autocorrelogram(spike_seconds, acg_window);
edges = acg_window(1):acg_bin:acg_window(2);
if numel(edges) < 2
    edges = [acg_window(1) acg_window(2)];
end
[counts, edges] = histcounts(lags, edges);
centers = (edges(1:end - 1) + edges(2:end)) / 2;
bar(centers, counts, 1, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none');
xlim(acg_window);
xlabel('Lag (s)')
ylabel('Count')
title('Autocorrelogram')
box off
end

function row = mean_omitnan(matrix)
n_times = size(matrix, 2);
row = nan(1, n_times);
for time_index = 1:n_times
    values = matrix(:, time_index);
    values = values(isfinite(values));
    if ~isempty(values)
        row(time_index) = mean(values);
    end
end
end
