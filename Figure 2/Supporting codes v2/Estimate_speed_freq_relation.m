%% Estimate_speed_freq_relation
% Driver script: relate locomotion speed to frequency-band power across animals.
% Calls GENERATE_SPEED_FREQ_RELATION for each animal-partner pair.
% Outputs: psth_structure, animal_names.

%% Animal and partner lists
list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single', ...
    'B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};
list_of_paertner = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2 3]};

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\psth power by frequency and behavior';

%% Spectrogram parameters (theta)
wind_length  = .250;
wind_overlap = 0.249;
f            = 5:.1:16;
freq_range   = [6 12];

%% Estimate speed–frequency relation
psth_structure = [];
animal_names = [];

for fn = 1:numel(list_of_animals)
    for pt = list_of_paertner{fn}
        disp([list_of_animals{fn}, ' P', num2str(pt)])

        new_struct = GENERATE_SPEED_FREQ_RELATION(list_of_animals{fn}, pt, freq_range, f, wind_length, wind_overlap);

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

%% Save results
save(fullfile(saving_folder, 'psth_structure_speed_theta_v2.mat'), 'psth_structure')
save(fullfile(saving_folder, 'animal_names_speed_theta_v2.mat'), 'animal_names')
