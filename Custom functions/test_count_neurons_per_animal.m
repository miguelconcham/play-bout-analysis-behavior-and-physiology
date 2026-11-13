function test_count_neurons_per_animal()
%TEST_COUNT_NEURONS_PER_ANIMAL  Good-unit totals written for each animal.

this_dir = fileparts(mfilename('fullpath'));
if ~isempty(this_dir)
    addpath(this_dir);
end

root = fullfile(tempdir, 'neurons_per_animal_test');
if exist(root, 'dir') == 7
    rmdir(root, 's');
end
mkdir(root);
cleaner = onCleanup(@() rmdir(root, 's')); %#ok<NASGU>

write_session(root, 'B1S1 0101 Single', 'cluster_info.tsv', { ...
    '1', 'good', 'good', '100'; ...
    '2', 'mua', 'mua', '50'; ...
    '3', 'good', 'noise', '10'; ...
    '4', 'good', '', '20'; ...
    '5', 'mua', '', '8'; ...
    '6', 'good', 'good', '0'});
write_session(root, 'B1S1 0202 Single', 'cluster_info_PAG.tsv', { ...
    '1', 'good', 'good', '5'; ...
    '2', 'good', 'good', '6'});
write_tsv(fullfile(root, 'B1S1 0202 Single', 'cluster_info.tsv'), { ...
    '9', 'good', 'good', '100'});
write_session(root, 'B1S1 0303 Dual', 'cluster_info_PAG.tsv', { ...
    '1', 'good', 'good', '4'});
write_tsv(fullfile(root, 'B1S1 0303 Dual', 'cluster_info_CTX.tsv'), { ...
    '2', 'good', 'good', '7'; ...
    '3', 'mua', 'mua', '3'; ...
    '4', 'good', 'good', '9'});
mkdir(fullfile(root, 'B2S1 0101 Single'));
write_kslabel_tsv(fullfile(root, 'B2S1 0101 Single', 'cluster_info.tsv'), { ...
    '1', 'good', '12'; ...
    '2', 'mua', '4'; ...
    '3', 'good', '1'});
mkdir(fullfile(root, 'notes'));

csv_file = fullfile(root, 'neurons_per_animal.csv');
counts = count_neurons_per_animal(root, csv_file);

assert(height(counts) == 3, 'Animals were merged or dropped.');
assert(strcmp(counts.Animal{1}, 'B1S1') && strcmp(counts.Recording{1}, 'Dual'), ...
    'Dual was not its own row.');
assert(counts.NNeurons(1) == 3 && counts.NSessions(1) == 1, ...
    'PAG and CTX good units were not added.');
assert(strcmp(counts.Animal{2}, 'B1S1') && strcmp(counts.Recording{2}, 'Single'), ...
    'The two Single days were not one animal.');
assert(counts.NNeurons(2) == 4 && counts.NSessions(2) == 2, ...
    'Blank Phy groups, noise, mua, empty units, or the unsuffixed copy were mishandled.');
assert(strcmp(counts.Sessions{2}, 'B1S1 0101 Single; B1S1 0202 Single'), ...
    'Session names were not sorted.');
assert(strcmp(counts.Animal{3}, 'B2S1') && counts.NNeurons(3) == 2, ...
    'A file with only KSLabel was counted wrong.');

written = readtable(csv_file, 'FileType', 'text', 'Delimiter', ',');
assert(isequal(written.NNeurons, counts.NNeurons), ...
    'The CSV neuron counts differ from the table.');

empty_root = fullfile(root, 'empty');
mkdir(empty_root);
empty_counts = count_neurons_per_animal(empty_root);
assert(isempty(empty_counts), 'A folder with no sessions returned rows.');

label_root = fullfile(root, 'missing-label');
mkdir(label_root);
mkdir(fullfile(label_root, 'B9S9 0101 Single'));
fid = fopen(fullfile(label_root, 'B9S9 0101 Single', 'cluster_info.tsv'), 'w');
fprintf(fid, 'cluster_id\tch\n1\t0\n');
fclose(fid);
threw = false;
try
    count_neurons_per_animal(label_root);
catch label_error
    threw = strcmp(label_error.identifier, 'count_neurons_per_animal:label');
end
assert(threw, 'A cluster_info file without a label was accepted.');

missing = false;
try
    count_neurons_per_animal(fullfile(root, 'does-not-exist'));
catch folder_error
    missing = strcmp(folder_error.identifier, 'count_neurons_per_animal:folder');
end
assert(missing, 'A missing NPX folder was not reported.');

lag = 0:0.01:1;
curve = ones(size(lag));
kept = population_record('B1S1 0202 Single', 4, 'trough', lag, curve);
also_kept = population_record('B1S1 0101 Single', 8, 'peak', lag, curve);
dropped = population_record('B1S1 0101 Single', 9, 'trough', lag, curve);
dropped.power.low = rmfield(dropped.power.low, 'sum_pctl');
other = population_record('B9S9 0101 Dual', 1, 'nonentrained', lag, curve);
plot_counts = count_population_logit_neurons({kept, also_kept, dropped, other});
assert(height(plot_counts) == 2, 'Population-logit animals were merged or dropped.');
single = plot_counts(strcmp(plot_counts.Recording, 'Single'), :);
dual = plot_counts(strcmp(plot_counts.Recording, 'Dual'), :);
assert(single.NPopulationLogit == 2 && single.NSessions == 2, ...
    'The two Single sessions were not both counted.');
assert(strcmp(single.Sessions{1}, 'B1S1 0101 Single; B1S1 0202 Single'), ...
    'Population-logit session names were not sorted.');
assert(dual.NPopulationLogit == 1, ...
    'A neuron missing a power curve was counted, or Dual was missed.');
combined = add_population_logit_counts(counts, plot_counts);
b1 = combined(strcmp(combined.Animal, 'B1S1') & strcmp(combined.Recording, 'Single'), :);
assert(b1.NNeurons == 4 && b1.NPopulationLogit == 2, ...
    'Good units and population-logit neurons were not kept on one row.');
b9 = combined(strcmp(combined.Animal, 'B9S9'), :);
assert(b9.NNeurons == 0 && b9.NPopulationLogit == 1, ...
    'An animal present only in the population plot was left out.');
unknown = add_population_logit_counts(counts, []);
assert(all(isnan(unknown.NPopulationLogit)), ...
    'A missing checkpoint did not leave the population-logit count empty.');

disp('test_count_neurons_per_animal passed');
end

function write_session(root, session_id, file_name, rows)
session_dir = fullfile(root, session_id);
if exist(session_dir, 'dir') ~= 7
    mkdir(session_dir);
end
write_tsv(fullfile(session_dir, file_name), rows);
end

function write_tsv(path, rows)
fid = fopen(path, 'w');
assert(fid > 0, 'Could not write %s.', path);
fprintf(fid, 'cluster_id\tKSLabel\tgroup\tn_spikes\n');
for row_index = 1:size(rows, 1)
    fprintf(fid, '%s\n', strjoin(rows(row_index, :), '\t'));
end
fclose(fid);
end

function write_kslabel_tsv(path, rows)
fid = fopen(path, 'w');
assert(fid > 0, 'Could not write %s.', path);
fprintf(fid, 'cluster_id\tKSLabel\tn_spikes\n');
for row_index = 1:size(rows, 1)
    fprintf(fid, '%s\n', strjoin(rows(row_index, :), '\t'));
end
fclose(fid);
end

function rec = population_record(session_id, cluster_id, group_name, lag, curve)
rec.session = session_id;
rec.cluster_id = cluster_id;
rec.lock_group = group_name;
rec.all.pctl_lag = lag;
rec.all.sum_pctl = curve;
rec.all.ratio_pctl = curve;
rec.power.high.pctl_lag = lag;
rec.power.high.sum_pctl = curve;
rec.power.high.ratio_pctl = curve;
rec.power.low.pctl_lag = lag;
rec.power.low.sum_pctl = curve;
rec.power.low.ratio_pctl = curve;
end
