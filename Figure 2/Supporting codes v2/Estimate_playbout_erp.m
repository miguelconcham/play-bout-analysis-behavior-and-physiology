%% Estimate_playbout_erp
% One raw-LFP ERP per play session, gathered into a single struct.
%
% A play session is one Neuropixels recording. GENERATE_PLAYBOUT_ERP cuts
% the unfiltered mid-LPAG LFP from -2 s to +2 s around every play-bout
% onset and keeps the whole waveform. This script loops the recordings,
% averages the complete epochs in each one, and saves the array playbout_erp.
% That average still contains saturated trials. Analyze_playbout_erp drops
% those trials, refits each session ERP, and estimates the grand average.
%
% Run add_repo_paths.m once per MATLAB session before this script.

%% Paths and epoch
options = playbout_erp_defaults();
window_s = options.window_s;

npx_Raw_Data = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\playbout erp';

if ~isfolder(npx_Raw_Data)
    error('Estimate_playbout_erp:path', ...
        'Neuropixels folder not found: %s', npx_Raw_Data);
end
animal_list = dir(npx_Raw_Data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));
animal2exclude = {''};
animal_list(ismember({animal_list.name}, animal2exclude)) = [];

%% One ERP per recording
playbout_erp = [];
animal_names = cell(0, 2);
n_failed = 0;

for fn = 1:numel(animal_list)
    animal_path = fullfile(npx_Raw_Data, animal_list(fn).name);
    fprintf('Session %d/%d: %s\n', fn, numel(animal_list), animal_list(fn).name);
    try
        new_struct = GENERATE_PLAYBOUT_ERP(animal_path, window_s);
    catch session_error
        n_failed = n_failed + 1;
        warning('Estimate_playbout_erp:session', ...
            'Play-bout ERP failed for %s: %s', animal_list(fn).name, session_error.message);
        continue
    end
    if isempty(new_struct)
        continue
    end

    if isempty(playbout_erp)
        playbout_erp = new_struct(:);
    else
        playbout_erp(end + (1:numel(new_struct)), 1) = new_struct(:); %#ok<SAGROW>
    end
    animal_names = [animal_names; ...
        [repmat({animal_list(fn).name}, numel(new_struct), 1), ...
        num2cell((1:numel(new_struct))')]]; %#ok<AGROW>
end

if isempty(playbout_erp)
    error('Estimate_playbout_erp:empty', ...
        'No play-bout ERP was saved. %d recordings failed.', n_failed);
end

%% Save every recording together
if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
output_file = fullfile(saving_folder, 'playbout_erp.mat');
disp(['Saving ', output_file])
save(output_file, 'playbout_erp', 'animal_names', 'window_s', '-v7.3');
fprintf(['Play-bout ERP estimation complete. %d recording-electrode entries, ' ...
    '%d recordings failed.\n'], numel(playbout_erp), n_failed);
fprintf(['Each erp field is the mean raw waveform of every complete trial. ' ...
    'Saturated-trial rejection is in Analyze_playbout_erp.\n']);
