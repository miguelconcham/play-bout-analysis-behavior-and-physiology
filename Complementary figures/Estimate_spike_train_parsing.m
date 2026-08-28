%% Estimate_spike_train_parsing
% First three sections of Estiamte_phase_copuling, calling GENERATE_SPIKE_TRAIN_PARSING
% (session load + firing-run oscillation phase and power). No plot sections.
% Each animal is saved as its own .mat (one file per returned struct) so
% later loads do not pull the whole cohort into memory.

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

sr           = 2500;
filter_order = 2000;
freq_range   = [1 6];   % oscillation band passed to GENERATE (change for delta/gamma)
freq_tag     = sprintf('%g-%gHz', freq_range(1), freq_range(2));

Hd_freq = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range(1), ...
    'CutoffFrequency2', freq_range(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

%%
tic
for fn = 5:numel(animal_list)

    animal_path = [npx_Raw_Data, '\', animal_list(fn).name];
    animal_name = animal_list(fn).name;
    disp(['Animal ', num2str(fn), ' of ', num2str(numel(animal_list)), ': ', animal_name])

    transt_psth = GENERATE_SPIKE_TRAIN_PARSING(animal_path, Hd_freq);
    if isempty(transt_psth)
        disp(['  empty result, skipping save: ', animal_name])
        toc
        continue
    end

    for sub_j = 1:numel(transt_psth)
        parsing_struct = transt_psth(sub_j);
        if numel(transt_psth) == 1
            out_file = [animal_name, '_FreqRange_', freq_tag, '_spike_train_parsing.mat'];
        else
            out_file = [animal_name, '_FreqRange_', freq_tag, '_spike_train_parsing_', num2str(sub_j), '.mat'];
        end
        save([saving_folder, '\', out_file], 'parsing_struct', '-v7.3');
        animal_names = [animal_names; {animal_name, sub_j, out_file}]; %#ok<AGROW>
        disp(['  saved ', out_file])
    end
    save([saving_folder, '\spike_train_parsing_animal_names.mat'], 'animal_names');
    clear transt_psth parsing_struct
    toc

end
%%