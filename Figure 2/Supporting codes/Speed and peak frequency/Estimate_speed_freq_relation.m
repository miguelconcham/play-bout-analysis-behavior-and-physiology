%% Estimate_speed_freq_relation
% Rebuild the speed and distance mats, including SpeedCm and DistanceCm.
%
% Run this script again after the pixel-to-cm change. It calls
% GENERATE_SPEED_FREQ_RELATION for each animal-partner pair, once for delta
% and once for theta.
%
% Delta files, read by Analyze_distance_freq_relation.m and
% Analyze_speed_freq_relation.m:
%   psth_structure_speed_delta.mat
%   animal_names_speed_delta.mat
%
% Theta files:
%   psth_structure_speed_theta_v2.mat
%   animal_names_speed_theta_v2.mat
%
% The FOOOF peak-frequency correlations do not use these files. Their speed
% is already pixels per second, and Analyze_speed_peak_frequency.m converts
% that axis to cm/s.

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes\Speed and peak frequency'];
end
addpath(script_dir);

%% Animal and partner lists
list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single', ...
    'B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};
list_of_paertner = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2 3]};

saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

%% One pass per band. Delta is the file the coverage scripts load.
% Defaults are used only when that mat does not already store a grid.
% An existing mat keeps its own f and freq_range, so play minus non-play
% stays on the same power. The centimetre columns do not move those signs.
band_label = {'delta', 'theta'};
fallback_settings = { ...
    struct('freq_range', [1 5], 'f', 1:.1:16, 'wind_length', 0.250, 'wind_overlap', 0.249), ...
    struct('freq_range', [6 12], 'f', 5:.1:16, 'wind_length', 0.250, 'wind_overlap', 0.249)};
structure_files = { ...
    'psth_structure_speed_delta.mat', ...
    'psth_structure_speed_theta_v2.mat'};
names_files = { ...
    'animal_names_speed_delta.mat', ...
    'animal_names_speed_theta_v2.mat'};

fprintf(['Pixel scale %.4f cm/pixel. SpeedCm is pixels per second * %.4f,\n' ...
    'using the time between frames. DistanceCm is pixels * %.4f.\n' ...
    'The saved Speed column stays the within-session z-score.\n'], ...
    pixel_to_cm(), pixel_to_cm(), pixel_to_cm());

for band_n = 1:numel(band_label)
    structure_path = fullfile(saving_folder, structure_files{band_n});
    names_path = fullfile(saving_folder, names_files{band_n});
    band_settings = fallback_settings{band_n};
    if exist(structure_path, 'file') == 2
        previous = load(structure_path, 'psth_structure');
        if isfield(previous, 'psth_structure') && ~isempty(previous.psth_structure)
            band_settings = stored_band_settings(previous.psth_structure(1), band_settings);
        end
    end
    freq_range = band_settings.freq_range;
    f = band_settings.f;
    wind_length = band_settings.wind_length;
    wind_overlap = band_settings.wind_overlap;
    if isfield(band_settings, 'reused') && band_settings.reused
        fprintf(['\n%s reuses the stored grid %.4g-%.4g Hz (%d frequencies), ' ...
            'window %.3f s, overlap %.3f s.\n' ...
            'Play minus non-play is not moved by the centimetre scale.\n'], ...
            band_label{band_n}, freq_range(1), freq_range(2), numel(f), ...
            wind_length, wind_overlap);
    else
        fprintf(['\n%s has no stored grid. Using %.4g-%.4g Hz (%d frequencies).\n'], ...
            band_label{band_n}, freq_range(1), freq_range(2), numel(f));
    end

    psth_structure = [];
    animal_names = [];
    for fn = 1:numel(list_of_animals)
        for pt = list_of_paertner{fn}
            disp([list_of_animals{fn}, ' P', num2str(pt)])
            new_struct = GENERATE_SPEED_FREQ_RELATION( ...
                list_of_animals{fn}, pt, freq_range, f, wind_length, wind_overlap);

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
                 num2cell(1:numel(new_struct))']]; %#ok<AGROW>
        end
    end

    fprintf('Saving %s\n', structure_path);
    save(structure_path, 'psth_structure')
    save(names_path, 'animal_names')
end
