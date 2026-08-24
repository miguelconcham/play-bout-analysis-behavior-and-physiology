%% Estimate_psth_full_spectrogram
% Driver script: compute full-spectrum PSTH spectrograms for all animals.
% Calls GENERATE_PSTH_SPECTROGRAM for each session and merges results.
% Outputs: psth_structure, animal_names.

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\psth power by frequency and behavior';

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

%% Spectrogram parameters (full spectrum)
wind_length    = 1;
wind_overlap   = .990;
min_separation = .200;
f              = 0.1:.1:30;
freq_pow_range = [35 90];

%% Estimate full-spectrum PSTH
psth_structure = [];
animal_names = {};

for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    new_struct = GENERATE_PSTH_SPECTROGRAM(animal_path, wind_length, wind_overlap, min_separation, f, freq_pow_range);

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

%% Save results
disp('saving')
save(fullfile(saving_folder, 'psth_structure_delta_full_spectrogram.mat'), 'psth_structure', '-v7.3');
save(fullfile(saving_folder, 'animal_names_delta_full_spectrogram.mat'), 'animal_names');
