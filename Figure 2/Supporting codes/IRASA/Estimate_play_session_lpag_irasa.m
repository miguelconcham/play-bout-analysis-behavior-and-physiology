%% Estimate_play_session_lpag_irasa
% Driver for session-wise mid-LPAG IRASA during annotated partner
% sessions. GENERATE_PLAY_SESSION_LPAG_IRASA performs behavior loading,
% LPAG channel selection, and the IRASA separation.
%
% This analysis uses Session_structure / Partners session intervals only.
% It does not identify or analyze play bouts.
%
% Run add_repo_paths once in the MATLAB session, then run this script.
% Analyse_play_session_lpag_irasa loads the file this script saves.

%% Paths and animal selection
npx_Raw_Data = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

animal_list = dir(npx_Raw_Data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));

% Add exact Neuropixels folder names here when exclusions are required.
animal2exclude = {''};
animal_list(ismember({animal_list.name}, animal2exclude)) = [];

%% Estimate all recordings
play_session_lpag_irasa = struct([]);
animal_names = {};

for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    fprintf('Session %d/%d: %s\n', fn, numel(animal_list), ...
        animal_list(fn).name);

    new_struct = GENERATE_PLAY_SESSION_LPAG_IRASA(animal_path);

    if isempty(play_session_lpag_irasa)
        play_session_lpag_irasa = new_struct(:);
    else
        play_session_lpag_irasa(end + 1:end + numel(new_struct), 1) = ...
            new_struct(:);
    end

    animal_names = [animal_names; ...
        [repmat({animal_list(fn).name}, numel(new_struct), 1), ...
         num2cell(1:numel(new_struct))']]; %#ok<AGROW>
end

%% Save the combined session structure
if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
output_file = fullfile(saving_folder, 'play_session_lpag_irasa.mat');
disp(['Saving ', output_file])
save(output_file, 'play_session_lpag_irasa', 'animal_names', '-v7.3');
disp('Play-session LPAG IRASA estimation complete.')
