function animal_counts = count_neurons_per_animal(npx_raw_data, csv_file)
%COUNT_NEURONS_PER_ANIMAL  Good units per implanted animal.
%
% animal_counts has Animal, Recording, NNeurons, NSessions, and Sessions.
% A neuron is a cluster labeled good. When cluster_info has both group and
% KSLabel, a non-empty Phy group is the label and KSLabel fills blank
% group rows. A good cluster whose n_spikes column is 0 is not counted.
%
% Session folders are "{animal} {date} {Single|Dual|Single2}". The same
% animal and recording on two dates are added together. Single and Dual
% stay on separate rows: those words are the two animals in the behavior
% file. cluster_info_PAG.tsv and any other suffixed cluster_info are
% counted. The unsuffixed cluster_info.tsv is used only when the session
% has no suffixed file, so the same sorting is not added twice.
%
% Pass csv_file to write that table. The second argument can be omitted.

if nargin < 1 || isempty(npx_raw_data)
    npx_raw_data = default_npx_raw_data();
end
if nargin < 2
    csv_file = '';
end
if exist(npx_raw_data, 'dir') ~= 7
    error('count_neurons_per_animal:folder', ...
        'NPX raw data folder was not found: %s', npx_raw_data);
end

sessions = empty_sessions();
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
    info_files = sorting_info_files(session_dir);
    if isempty(info_files)
        if looks_like_session(session_id)
            warning('count_neurons_per_animal:noSorting', ...
                'No cluster_info file in %s.', session_id);
        end
        continue
    end
    [animal, recording] = animal_from_session(session_id);
    sessions(end + 1).animal = animal; %#ok<AGROW>
    sessions(end).recording = recording;
    sessions(end).session_id = session_id;
    sessions(end).n_neurons = count_good_units(info_files);
end

animal_counts = sessions_to_animal_table(sessions);
if ~isempty(csv_file)
    out_dir = fileparts(csv_file);
    if ~isempty(out_dir) && exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end
    writetable(animal_counts, csv_file);
end
end

function npx_raw_data = default_npx_raw_data()
npx_raw_data = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data\NPX data\NPX raw data'];
end

function sessions = empty_sessions()
sessions = struct('animal', {}, 'recording', {}, 'session_id', {}, ...
    'n_neurons', {});
end

function animal_counts = empty_animal_table()
animal_counts = table(cell(0, 1), cell(0, 1), zeros(0, 1), zeros(0, 1), ...
    cell(0, 1), 'VariableNames', ...
    {'Animal', 'Recording', 'NNeurons', 'NSessions', 'Sessions'});
end

function is_session = looks_like_session(session_id)
tokens = session_tokens(session_id);
is_session = numel(tokens) >= 3;
end

function [animal, recording] = animal_from_session(session_id)
[animal, recording] = session_animal_recording(session_id);
end

function tokens = session_tokens(session_id)
tokens = regexp(strtrim(session_id), '\s+', 'split');
tokens = tokens(~cellfun('isempty', tokens));
end

function files = sorting_info_files(session_dir)
listing = dir(fullfile(session_dir, 'cluster_info*.tsv'));
files = {};
unsuffixed = {};
for index = 1:numel(listing)
    if listing(index).isdir
        continue
    end
    [~, base] = fileparts(listing(index).name);
    path = fullfile(session_dir, listing(index).name);
    if strcmp(base, 'cluster_info')
        unsuffixed{end + 1} = path; %#ok<AGROW>
    elseif strncmp(base, 'cluster_info_', numel('cluster_info_'))
        files{end + 1} = path; %#ok<AGROW>
    end
end
if isempty(files)
    files = unsuffixed;
end
end

function n_neurons = count_good_units(info_files)
n_neurons = 0;
for file_index = 1:numel(info_files)
    n_neurons = n_neurons + count_good_units_in_file(info_files{file_index});
end
end

function n_neurons = count_good_units_in_file(info_file)
cluster_info = readtable(info_file, 'FileType', 'text', 'Delimiter', '\t');
groups = optional_column(cluster_info, {'group'});
labels = optional_column(cluster_info, {'KSLabel'});
spike_counts = optional_column(cluster_info, {'n_spikes', 'nspikes'});
if isempty(groups) && isempty(labels)
    error('count_neurons_per_animal:label', ...
        'cluster_info needs a group or KSLabel column: %s', info_file);
end
n_neurons = 0;
for row_index = 1:height(cluster_info)
    label = '';
    if ~isempty(groups)
        label = group_label(groups, row_index);
    end
    if isempty(label) && ~isempty(labels)
        label = group_label(labels, row_index);
    end
    if ~strcmp(label, 'good')
        continue
    end
    if ~isempty(spike_counts) && ~row_has_spikes(spike_counts, row_index)
        continue
    end
    n_neurons = n_neurons + 1;
end
end

function has_spikes = row_has_spikes(spike_counts, row_index)
count = scalar_number(spike_counts, row_index);
has_spikes = isfinite(count) && count > 0;
end

function animal_counts = sessions_to_animal_table(sessions)
if isempty(sessions)
    animal_counts = empty_animal_table();
    return
end
keys = cell(1, numel(sessions));
for index = 1:numel(sessions)
    keys{index} = [sessions(index).animal, char(0), sessions(index).recording];
end
[~, first] = unique(keys, 'stable');
order = sort(first);
n_animals = numel(order);
animal_col = cell(n_animals, 1);
recording_col = cell(n_animals, 1);
n_neurons = zeros(n_animals, 1);
n_sessions = zeros(n_animals, 1);
session_col = cell(n_animals, 1);
for animal_index = 1:n_animals
    key = keys{order(animal_index)};
    keep = strcmp(keys, key);
    these = sessions(keep);
    [~, session_order] = sort({these.session_id});
    these = these(session_order);
    animal_col{animal_index} = these(1).animal;
    recording_col{animal_index} = these(1).recording;
    n_neurons(animal_index) = sum([these.n_neurons]);
    n_sessions(animal_index) = numel(these);
    session_col{animal_index} = strjoin({these.session_id}, '; ');
end
animal_counts = table(animal_col, recording_col, n_neurons, n_sessions, ...
    session_col, 'VariableNames', ...
    {'Animal', 'Recording', 'NNeurons', 'NSessions', 'Sessions'});
[~, animal_order] = sort(strcat(animal_counts.Animal, ' ', animal_counts.Recording));
animal_counts = animal_counts(animal_order, :);
end

function values = optional_column(table_data, names)
values = [];
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
    value = values{index};
else
    value = values(index);
end
if exist('isstring', 'builtin') && isstring(value)
    value = char(value);
end
if ischar(value)
    number = str2double(value);
else
    number = double(value);
end
end
