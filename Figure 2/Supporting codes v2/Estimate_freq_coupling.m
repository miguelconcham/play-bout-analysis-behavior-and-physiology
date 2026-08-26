%% Estimate_freq_coupling
% Driver script: compute cross-frequency coupling between two band-pass signals.
% Calls GENERATE_FREQ_COUPLING_STRUCT for each animal and merges results.
% Outputs: psth_structure, animal_names (save manually when needed).

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

%% Filter design
freq_range_1 = [.1 6];
freq_range_2 = [35 90];
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

%% Estimate frequency coupling
psth_structure = [];
animal_names = {};

for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    new_struct = GENERATE_FREQ_COUPLING_STRUCT(animal_path, Hd_freq1, Hd_freq2);

    if isempty(psth_structure)
        psth_structure = new_struct;
    else
        start_idx = numel(psth_structure) + 1;
        for sub_j = 1:numel(new_struct)
            psth_structure(start_idx + sub_j - 1) = new_struct(sub_j);
        end
    end

    animal_names = [animal_names; ...
        [repmat({animal_list(fn).name}, numel(new_struct), 1), num2cell(1:numel(new_struct))']];
end

%% Save results (uncomment when needed)
% disp('saving')
% save(fullfile(saving_folder, 'couplig_structure_V3.mat'), 'psth_structure', '-v7.3');
% save(fullfile(saving_folder, 'animal_names_coupling_V3.mat'), 'animal_names');
