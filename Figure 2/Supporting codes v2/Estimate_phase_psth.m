%% Estimate_phase_psth
% Driver: compute filtered LFP phase PSTH around bout onset (band / behavior set in GENERATE).
% Output files are listed under PHASE in:
%   Figure 2/Figure 2 Psth animal names and result combinations.txt
% Calls GENERATE_PHASE_EXPLORATORY_ONSET (change behavior list inside GENERATE for other bout types).

this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Estimate_phase_psth.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
combo_file = fullfile(repo_root, 'Figure 2', 'Figure 2 Psth animal names and result combinations.txt');

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
animal_file_names = cellfun(@(x) ['B', x], strsplit([animal_list.name], 'B'), 'UniformOutput', false)';
animal_file_names(1) = [];
animal2exclude = {''};
animal_list(ismember(animal_file_names, animal2exclude)) = [];

%% Filter design — example: delta (1–5 Hz); see combo file for theta/gamma PSTH pairs
freq_range_1 = [1 5];
sr           = 2500;
filter_order = 2000;

Hd_freq = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range_1(1), ...
    'CutoffFrequency2', freq_range_1(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

%% Estimate phase PSTH across animals
psth_structure = [];
animal_names = {};

for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    new_struct = GENERATE_PHASE_EXPLORATORY_ONSET(animal_path, Hd_freq);

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

%% Save results — filenames from combo file (PHASE section)
disp('saving')
if ~exist(saving_folder, 'dir'), mkdir(saving_folder); end
save(fullfile(saving_folder, 'phase_onset_playbout.mat'), 'psth_structure', '-v7.3');
save(fullfile(saving_folder, 'animal_names_dphase_onset_playbout.mat'), 'animal_names');
