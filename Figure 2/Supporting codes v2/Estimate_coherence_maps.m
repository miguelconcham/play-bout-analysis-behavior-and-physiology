%% Estimate_coherence_maps
% Driver script: compute PLI/coherence maps between brain areas (PAG recordings).
% Calls GENERATE_COHERNECE_MAPS_STRUCT for the full pairwise maps.
% Plotting is in Analyze_coherence_maps.m. The faster LPAG-center-only
% pipeline is Estimate_coherence_maps_lpag_ref.m.

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];

animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

%% Filter design (delta–theta)
freq_range_1 = [1 5];
freq_range_2 = [6 12];
sr           = 2500;
filter_order = 2000;

Hd_freq1 = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range_1(1), ...
    'CutoffFrequency2', freq_range_1(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

Hd_freq2 = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range_2(1), ...
    'CutoffFrequency2', freq_range_2(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

%% Estimate coherence maps
pli_struct = [];
animal_names = {};

for fn = 1:numel(animal_list)
    disp(fn)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    new_struct = GENERATE_COHERNECE_MAPS_STRUCT(animal_path, Hd_freq1, Hd_freq2, 'PAG');

    if isempty(new_struct)
        continue
    end

    if isempty(pli_struct)
        pli_struct = new_struct;
    else
        start_idx = numel(pli_struct) + 1;
        for sub_j = 1:numel(new_struct)
            pli_struct(start_idx + sub_j - 1) = new_struct(sub_j);
        end
    end

    animal_names = [animal_names; ...
        [repmat({animal_list(fn).name}, numel(new_struct), 1), num2cell(1:numel(new_struct))']];
end

%% Save results
disp('saving')
save(fullfile(saving_folder, 'coherence_pli_struct_2026.mat'), 'pli_struct', '-v7.3');
save(fullfile(saving_folder, 'coherence_pli_structure_animal_names2026.mat'), 'animal_names');
