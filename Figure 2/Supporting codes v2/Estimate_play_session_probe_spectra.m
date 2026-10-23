%% Estimate_play_session_probe_spectra
% Channel-wise periodic and aperiodic spectra during annotated partner
% sessions. GENERATE_PLAY_SESSION_PROBE_SPECTRA estimates Welch and direct
% Fourier spectra for every mapped channel and fits FOOOF.
%
% Frequency bands are not averaged here. Each recording is saved on its
% own so a long run keeps completed animals. A folder named
% "B1D1 1007 Dual" is written as
% play_session_probe_periodic_spectra_B1D1_1007_Dual.mat.

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

animal_list = dir(npx_Raw_Data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));

% Add exact Neuropixels folder names here when exclusions are required.
animal2exclude = {''};
animal_list(ismember({animal_list.name}, animal2exclude)) = [];

if ~isfolder(saving_folder)
    mkdir(saving_folder);
end

%% Estimate each recording and save it before moving on
for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    fprintf('Session %d/%d: %s\n', fn, numel(animal_list), ...
        animal_list(fn).name);

    probe_spectra = GENERATE_PLAY_SESSION_PROBE_SPECTRA(animal_path);

    output_file = fullfile(saving_folder, ...
        probe_periodic_spectra_filename(animal_list(fn).name));
    disp(['Saving ', output_file])
    save(output_file, 'probe_spectra', '-v7.3');
end

disp('Channel-wise play-session spectral estimation complete.')
