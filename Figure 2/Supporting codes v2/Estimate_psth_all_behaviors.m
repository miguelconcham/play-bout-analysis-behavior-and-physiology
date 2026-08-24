%% Estimate_psth_all_behaviors
% Driver script: compute delta-band PSTH for all behaviors (self/other/partner).
% Calls GENERATE_PSTH_ALL_BEHAVIORS for each animal and merges results.
% Outputs: psth_structure, animal_names (save manually when needed).

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\psth power by frequency and behavior';

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

%% Spectrogram parameters (delta)
wind_length    = 1;
wind_overlap   = .990;
min_separation = .200;
f              = .1:.5:6;
freq_pow_range = [1 5];

%% Estimate PSTH across all behaviors
psth_structure = [];
animal_names = {};

for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    new_struct = GENERATE_PSTH_ALL_BEHAVIORS(animal_path, wind_length, wind_overlap, min_separation, f, freq_pow_range);

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
% save(fullfile(saving_folder, 'psth_structure_delta_all_behaviors_self_other_partner.mat'), 'psth_structure', '-v7.3');
% save(fullfile(saving_folder, 'animal_names_delta_all_behaviors_self_other_partner.mat'), 'animal_names');
