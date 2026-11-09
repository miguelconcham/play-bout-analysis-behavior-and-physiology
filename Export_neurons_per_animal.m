%% Export_neurons_per_animal
% Write neurons_per_animal.csv, one row per implanted animal.
%
% Session folders look like "B2S2 1110 Single". Animal is the first word
% and Recording is the third (Single, Dual, or Single2). Dates of the
% same animal and recording are added. A neuron is a cluster labeled good
% in cluster_info. Suffixed files such as cluster_info_PAG.tsv are used
% when they are present.

repo_root = fileparts(mfilename('fullpath'));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];
npx_raw_data = fullfile(data_root, 'NPX data', 'NPX raw data');
csv_file = fullfile(data_root, 'Analysis results', 'neurons_per_animal.csv');

animal_counts = count_neurons_per_animal(npx_raw_data, csv_file);
disp(animal_counts);
fprintf('Wrote %d animals (%d good units) to\n  %s\n', ...
    height(animal_counts), sum(animal_counts.NNeurons), csv_file);
