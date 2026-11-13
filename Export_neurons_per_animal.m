%% Export_neurons_per_animal
% Write neurons_per_animal.csv, one row per implanted animal.
%
% Session folders look like "B2S2 1110 Single". Animal is the first word
% and Recording is the third (Single, Dual, or Single2). Dates of the
% same animal and recording are added. A neuron is a cluster labeled good
% in cluster_info. Suffixed files such as cluster_info_PAG.tsv are used
% when they are present.
%
% NPopulationLogit is the number of those animals' neurons kept by
% Plot_acg_population_logit: all-runs, high-power, and low-power curves
% on both the sum-to-1 and rate-ratio scales. That count is read from
% locking_period_acg_shuffle_n500.mat. When the checkpoint is absent the
% column is left empty.

repo_root = fileparts(mfilename('fullpath'));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

data_root = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
    'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
    'play bout analysis behavior and physiology\Data'];
npx_raw_data = fullfile(data_root, 'NPX data', 'NPX raw data');
csv_file = fullfile(data_root, 'Analysis results', 'neurons_per_animal.csv');
checkpoint = fullfile(data_root, 'Analysis results', 'spike train parsing', ...
    'locking_period_acg_shuffle_n500.mat');

animal_counts = count_neurons_per_animal(npx_raw_data);
if exist(checkpoint, 'file') == 2
    loaded = load(checkpoint, 'shuffle_records');
    plot_counts = count_population_logit_neurons(loaded.shuffle_records);
else
    warning('Export_neurons_per_animal:checkpoint', ...
        ['No population-logit checkpoint at %s. NPopulationLogit is left empty.'], ...
        checkpoint);
    plot_counts = [];
end
animal_counts = add_population_logit_counts(animal_counts, plot_counts);
out_dir = fileparts(csv_file);
if ~isempty(out_dir) && exist(out_dir, 'dir') ~= 7
    mkdir(out_dir);
end
writetable(animal_counts, csv_file);
disp(animal_counts);
n_plot = animal_counts.NPopulationLogit;
n_plot = n_plot(isfinite(n_plot));
fprintf('Wrote %d animals (%d good units, %d population-logit neurons) to\n  %s\n', ...
    height(animal_counts), sum(animal_counts.NNeurons), sum(n_plot), csv_file);
