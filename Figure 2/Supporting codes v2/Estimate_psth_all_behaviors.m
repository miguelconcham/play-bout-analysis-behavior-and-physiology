%% Estimate_psth_all_behaviors
% Driver: band-limited PSTH for every behavior (self and other).
% Calls GENERATE_PSTH_ALL_BEHAVIORS for each animal and merges results.
%
% Each session stores three alignments, for self and for the partner:
%   1 onset-aligned, 2 offset-aligned, 3 time-wrapped (pre, bout, post).
% Play and non-play behaviors are both included. Analyze_psth_all_behaviors.m
% splits them and fits the LMEs.
%
% Delta is already saved. This run estimates theta and gamma on the same
% 10 ms spectrogram step (wind_length - wind_overlap) so the onset and
% time-wrapped time bases match the delta per-behavior PSTHs.
% Add 'delta' to bands_to_run to rebuild the delta files (overwrites them).

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
for band_n = 1:numel(bands)
    band = bands(band_n);
    if ~ismember(band.name, bands_to_run)
        continue
    end

    step_s = band.wind_length - band.wind_overlap;
    fprintf('Band %s: %g-%g Hz, window %g s, step %g s\n', ...
        band.name, band.freq_pow_range(1), band.freq_pow_range(2), ...
        band.wind_length, step_s);

    psth_structure = [];
    animal_names = {};

    for fn = 1:numel(animal_list)
        animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
        fprintf('%s  %s\n', band.name, animal_list(fn).name);
        new_struct = GENERATE_PSTH_ALL_BEHAVIORS(animal_path, ...
            band.wind_length, band.wind_overlap, band.min_separation, ...
            band.f, band.freq_pow_range);

        if ~isempty(new_struct)
            [new_struct.band_name] = deal(band.name);
        end

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

    if ~save_results
        continue
    end

    structure_file = ['psth_structure_' band.name '_all_behaviors_self_other_partner.mat'];
    names_file = ['animal_names_' band.name '_all_behaviors_self_other_partner.mat'];
    disp(['saving ' structure_file])
    save(fullfile(saving_folder, structure_file), 'psth_structure', '-v7.3');
    save(fullfile(saving_folder, names_file), 'animal_names');
end
