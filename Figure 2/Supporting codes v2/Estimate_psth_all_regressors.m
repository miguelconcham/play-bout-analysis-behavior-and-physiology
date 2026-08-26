%% Estimate_psth_all_regressors
% Driver script: compute PSTH with call/location regressors across animals and partners.
% Calls GENERATE_PSTH_CALL_LOC_REGRESSORS for each animal-partner pair.
% Outputs: psth_structure, animal_names (save manually when needed).

%% Animal and partner lists
list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single', ...
    'B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};
list_of_paertner = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2 3]};

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

%% Spectrogram parameters (beta/gamma)
f              = {10:.5:35, 35:100};
freq_range     = {[15 30],[40 90]};
wind_length    = {.150, .08};
spect_bin_size = 0.005;

%% Estimate PSTH with all regressors
psth_structure = [];
animal_names = [];

for fn = 1:numel(list_of_animals)
    for pt = list_of_paertner{fn}
        disp([list_of_animals{fn}, ' P', num2str(pt)])

        new_struct = GENERATE_PSTH_CALL_LOC_REGRESSORS(list_of_animals{fn}, pt, f, freq_range, wind_length, spect_bin_size);

        if isempty(psth_structure)
            psth_structure = new_struct;
        else
            start_idx = numel(psth_structure) + 1;
            for sub_j = 1:numel(new_struct)
                psth_structure(start_idx + sub_j - 1) = new_struct(sub_j);
            end
        end

        animal_names = [animal_names; ...
            [repmat({list_of_animals{fn}}, numel(new_struct), 1), ...
             repmat({pt}, numel(new_struct), 1), ...
             num2cell(1:numel(new_struct))']];
    end
end

%% Save results (uncomment when needed)
% save(fullfile(saving_folder, 'psth_structure_all_regressors_beta_gamma.mat'), 'psth_structure', '-v7.3');
% save(fullfile(saving_folder, 'animal_names_all_regressors_beta_gamma.mat'), 'animal_names');
