function report = check_npx_raw_folder_composition(experiment_info_file, raw_folders_file)
%CHECK_NPX_RAW_FOLDER_COMPOSITION  Compare experiment info with raw-file paths.
%
%   report = check_npx_raw_folder_composition
%   report = check_npx_raw_folder_composition(experiment_info_file, raw_folders_file)
%
% Opens Experiment_Info_GoodLooking.xlsx and NPX_raw_data_folders.xlsx and
% checks that each raw-file location is built from that session's animal
% code, date folder, NPX recording folder, probe, and recording number.
%
% Experiment info columns, by header or in this order:
%   animal code, date folder, NPX recording folder, probe, recording number
%
% The raw-folder sheet's first column is the session name from
% Data\NPX data\NPX raw data (for example 'B2S2 1110 Single2'). Another
% column is the full path.
%
% Prints every mismatch. report.ok is true when the two sheets agree.
% Area_limits_GoodLooking.xlsx is the depth table, not this experiment sheet.

if nargin < 1 || isempty(experiment_info_file)
    experiment_info_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\Experiment_Info_GoodLooking.xlsx'];
end
if nargin < 2 || isempty(raw_folders_file)
    raw_folders_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Data\NPX data\' ...
        'NPX_raw_data_folders.xlsx'];
end

experiment_info = read_excel_as_cell(experiment_info_file);
raw_folders = read_excel_as_cell(raw_folders_file);
report = compare_npx_folder_tables(experiment_info, raw_folders);
print_npx_folder_report(report);
end

function sheet = read_excel_as_cell(path)
if exist(path, 'file') ~= 2
    error('check_npx_raw_folder_composition:file', ...
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
sheet = table_to_sheet(table_data);
end

function sheet = table_to_sheet(table_data)
headers = cellstr(table_data.Properties.VariableNames);
n_rows = height(table_data);
n_cols = width(table_data);
rows = cell(n_rows, n_cols);
for c = 1:n_cols
    rows(:, c) = column_to_cell(table_data.(headers{c}));
end
sheet = [headers; rows];
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

function print_npx_folder_report(report)
for n = 1:numel(report.notes)
    fprintf('%s\n', report.notes{n});
end
if report.ok
    fprintf(['Folder composition agrees for %d experiment rows and ' ...
        '%d raw-folder rows.\n'], report.n_experiment, report.n_raw);
    return
end
fprintf('Folder composition: %d problem(s).\n', report.n_bad);
for k = 1:numel(report.rows)
    row = report.rows(k);
    if strcmp(row.status, 'ok')
        continue
    end
    fprintf('  %s sheet row %d [%s] %s\n    %s\n', ...
        row.source, row.sheet_row, row.status, row.session_name, row.detail);
    if ~isempty(row.actual_path)
        fprintf('    %s\n', row.actual_path);
    end
end
end
