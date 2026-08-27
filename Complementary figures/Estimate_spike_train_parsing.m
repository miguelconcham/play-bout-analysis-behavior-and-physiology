%% Estimate_spike_train_parsing
% First three sections of Estiamte_phase_copuling, calling GENERATE_SPIKE_TRAIN_PARSING
% (session load + firing-run oscillation phase and power). No plot sections.

%% Paths, animals, and filter
% Local repo Data.
data_root      = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_Raw_Data   = [data_root, '\NPX data\NPX raw data'];
saving_folder  = [data_root, '\Analysis results\spike train parsing'];
run('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');
if ~exist(saving_folder, 'dir')
    mkdir(saving_folder);
end

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];

animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];
animal_names = {};
n_strctut = 1;

sr           = 2500;
filter_order = 2000;
freq_range   = [1 6];   % oscillation band passed to GENERATE (change for delta/gamma)

Hd_freq = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range(1), ...
    'CutoffFrequency2', freq_range(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

%%
tic
for fn = 1:numel(animal_list)

    animal_path = [npx_Raw_Data, '\', animal_list(fn).name];
    if fn == 1
        transt_psth = GENERATE_SPIKE_TRAIN_PARSING(animal_path, Hd_freq);
        parsing_struct = transt_psth;

        n_strctut = n_strctut + numel(parsing_struct);
        animal_names = [animal_names; [repmat({animal_list(fn).name}, numel(parsing_struct), 1) num2cell(1:numel(parsing_struct))']];
    else
        transt_psth = GENERATE_SPIKE_TRAIN_PARSING(animal_path, Hd_freq);

        for sub_j = 1:numel(transt_psth)
            parsing_struct(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut + 1;
        end
        animal_names = [animal_names; [repmat({animal_list(fn).name}, numel(transt_psth), 1) num2cell(1:numel(transt_psth))']];
    end
    toc

end

%% save if needed
disp('saving')
if ~exist(saving_folder, 'dir')
    mkdir(saving_folder);
end
save([saving_folder, '\spike_train_parsing_structure.mat'], 'parsing_struct', '-v7.3');
save([saving_folder, '\spike_train_parsing_animal_names.mat'], 'animal_names');
