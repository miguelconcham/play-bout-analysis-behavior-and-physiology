function labels = read_oscillatory_labels(label_file)
%READ_OSCILLATORY_LABELS  Session, neuron id, and oscillation class.
%
% The table is oscillatory_neuron_labels.csv. classification is
% 'oscillatory', 'non-oscillatory', or 'not classified'. Rows that were
% not classified stay in the file and are removed later, when the six
% groups are built.

if exist(label_file, 'file') ~= 2
    error('group_ccg:labels', 'Oscillatory label table was not found at %s.', label_file);
end
stored = readtable(label_file);
required = {'session', 'neuron_id', 'oscillatory', 'classification'};
missing = required(~ismember(required, stored.Properties.VariableNames));
if ~isempty(missing)
    error('group_ccg:labels', ...
        'The label table is missing %s.', strjoin(missing, ', '));
end
n = height(stored);
labels.population_row = nan(n, 1);
if ismember('population_row', stored.Properties.VariableNames)
    labels.population_row = double(stored.population_row);
end
labels.session = cell(n, 1);
labels.classification = cell(n, 1);
for i = 1:n
    labels.session{i} = ccg_plain_text(stored.session(i));
    labels.classification{i} = ccg_plain_text(stored.classification(i));
end
labels.neuron_id = double(stored.neuron_id);
labels.oscillatory = double(stored.oscillatory);
if any(~isfinite(labels.neuron_id))
    error('group_ccg:labels', 'A neuron id in %s is not finite.', label_file);
end
end
