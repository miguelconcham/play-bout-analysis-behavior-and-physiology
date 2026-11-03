%% Estimate_psth_all_behaviors_theta_gamma
% New driver for theta and gamma per-behavior PSTHs. Does not replace
% Estimate_psth_all_behaviors.m or GENERATE_PSTH_ALL_BEHAVIORS.m.
% Calls GENERATE_PSTH_ALL_BEHAVIORS_BANDS once per animal.
%
% Each session stores three alignments, for self and for the partner:
%   1 onset-aligned, 2 offset-aligned, 3 time-wrapped (pre, bout, post).
% Play and non-play behaviors are both included.
% Analyze_psth_all_behaviors_theta_gamma.m splits them and fits the LMEs.
%
% Delta is already saved. This run estimates theta and gamma on the same
% 10 ms spectrogram step (wind_length - wind_overlap) so the onset and
% time-wrapped time bases match the delta per-behavior PSTHs.
% Add 'delta' to bands_to_run to rebuild the delta files (overwrites them).
%
% LFP_PAG is read once per animal. Every band in bands_to_run is then
% estimated from that copy. The spectrogram still runs once per band,
% because delta, theta, and gamma use different window lengths.

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

%% Bands
% Theta is 6-12 Hz (Figure 2, speed, and coherence). A 0.5 s window
% resolves that band at 2 Hz and still steps by 10 ms.
% Gamma is 35-90 Hz (play-bout gamma PSTH, call gamma, and cross-frequency
% coupling). The window is 0.1 s, the same length as the play-bout gamma PSTH.
% Delta parameters are the ones already used for
% psth_structure_delta_all_behaviors_self_other_partner.mat.
bands_to_run = {'theta', 'gamma'};
save_results = true;

bands = struct( ...
    'name', {'delta', 'theta', 'gamma'}, ...
    'wind_length', {1, 0.50, 0.10}, ...
    'wind_overlap', {0.990, 0.490, 0.090}, ...
    'min_separation', {0.200, 0.200, 0.200}, ...
    'f', {.1:.5:6, 4:.5:16, 35:5:90}, ...
    'freq_pow_range', {[1 5], [6 12], [35 90]});

%% Estimate PSTH across all behaviors
% Animals are the outer loop so each LFP_PAG file is opened once.
selected = bands(ismember({bands.name}, bands_to_run));
n_bands = numel(selected);
psth_structure = cell(1, n_bands);
animal_names = cell(1, n_bands);

for band_n = 1:n_bands
    band = selected(band_n);
    step_s = band.wind_length - band.wind_overlap;
    fprintf('Band %s: %g-%g Hz, window %g s, step %g s\n', ...
        band.name, band.freq_pow_range(1), band.freq_pow_range(2), ...
        band.wind_length, step_s);
    psth_structure{band_n} = [];
    animal_names{band_n} = {};
end

for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    fprintf('Loading %s once, then %s\n', animal_list(fn).name, strjoin({selected.name}, ', '));
    band_structs = GENERATE_PSTH_ALL_BEHAVIORS_BANDS(animal_path, selected);

    for band_n = 1:n_bands
        new_struct = band_structs{band_n};
        if ~isempty(new_struct)
            [new_struct.band_name] = deal(selected(band_n).name);
        end

        if isempty(psth_structure{band_n})
            psth_structure{band_n} = new_struct;
        else
            start_idx = numel(psth_structure{band_n}) + 1;
            for sub_j = 1:numel(new_struct)
                psth_structure{band_n}(start_idx + sub_j - 1) = new_struct(sub_j);
            end
        end

        animal_names{band_n} = [animal_names{band_n}; ...
            [repmat({animal_list(fn).name}, numel(new_struct), 1), num2cell(1:numel(new_struct))']];
    end
end

if save_results
    for band_n = 1:n_bands
        band = selected(band_n);
        structure_file = ['psth_structure_' band.name '_all_behaviors_self_other_partner.mat'];
        names_file = ['animal_names_' band.name '_all_behaviors_self_other_partner.mat'];
        disp(['saving ' structure_file])
        save_behavior_psth( ...
            fullfile(saving_folder, structure_file), ...
            fullfile(saving_folder, names_file), ...
            psth_structure{band_n}, animal_names{band_n});
    end
end

function save_behavior_psth(structure_file, names_file, psth_structure, animal_names)
save(structure_file, 'psth_structure', '-v7.3');
save(names_file, 'animal_names');
end
