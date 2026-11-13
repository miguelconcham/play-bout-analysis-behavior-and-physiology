function plot_counts = count_population_logit_neurons(records)
%COUNT_POPULATION_LOGIT_NEURONS  Neurons kept by Plot_acg_population_logit.
%
% plot_counts has Animal, Recording, NPopulationLogit, NSessions, and
% Sessions. A neuron is counted when align_population_percentiles keeps
% it: the all-runs, high-power, and low-power curves, on both the
% sum-to-1 and rate-ratio scales, share one lag axis. That is the row
% set section 1 of Plot_acg_population_logit stores. Animal and Recording
% come from the session name, the same split as count_neurons_per_animal.
% Dates of the same animal and recording are added together.

if nargin < 1 || isempty(records)
    plot_counts = empty_plot_table();
    return
end
if exist('align_population_percentiles', 'file') ~= 2
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(here), 'Complementary figures', ...
        'ACG power and locking type'));
end
if ~iscell(records)
    records = num2cell(records(:));
end
records = records(:);
[~, ~, ~, record_index] = align_population_percentiles(records);

sessions = struct('animal', {}, 'recording', {}, 'session_id', {}, ...
    'n_neurons', {});
for index = 1:numel(record_index)
    record = records{record_index(index)};
    session_id = record_text(record, 'session');
    [animal, recording] = session_animal_recording(session_id);
    sessions(end + 1).animal = animal; %#ok<AGROW>
    sessions(end).recording = recording;
    sessions(end).session_id = session_id;
    sessions(end).n_neurons = 1;
end
plot_counts = sessions_to_plot_table(sessions);
end


function plot_counts = empty_plot_table()
plot_counts = table(cell(0, 1), cell(0, 1), zeros(0, 1), zeros(0, 1), ...
    cell(0, 1), 'VariableNames', ...
    {'Animal', 'Recording', 'NPopulationLogit', 'NSessions', 'Sessions'});
end


function plot_counts = sessions_to_plot_table(sessions)
if isempty(sessions)
    plot_counts = empty_plot_table();
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
    names = unique({these.session_id}, 'stable');
    n_sessions(animal_index) = numel(names);
    session_col{animal_index} = strjoin(names, '; ');
end
plot_counts = table(animal_col, recording_col, n_neurons, n_sessions, ...
    session_col, 'VariableNames', ...
    {'Animal', 'Recording', 'NPopulationLogit', 'NSessions', 'Sessions'});
[~, animal_order] = sort(strcat(plot_counts.Animal, ' ', plot_counts.Recording));
plot_counts = plot_counts(animal_order, :);
end


function text = record_text(record, field_name)
text = '';
if ~isstruct(record) || ~isfield(record, field_name) || isempty(record.(field_name))
    return
end
value = record.(field_name);
if iscell(value)
    if isempty(value)
        return
    end
    value = value{1};
end
text = strtrim(char(value));
end
