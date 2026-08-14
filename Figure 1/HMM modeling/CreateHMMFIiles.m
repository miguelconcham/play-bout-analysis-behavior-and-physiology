data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Codes repository\Data';
behavior_folder = [data_root, '\Behavior backups'];
behavior_files = dir(fullfile(behavior_folder, '*.txt'));
dir_output = [data_root, '\HMM data\HMM raw data'];


for j=1:numel(behavior_files)
    animal_code = strrep(behavior_files(j).name, '.txt', '');

    disp('###################################')
    disp(['##### LOADING ', animal_code, '######'])
    disp('###################################')
    Create_HMM_inputs(dir_output,animal_code)
    close all
end

%%


behavior_files = dir(fullfile(behavior_folder, '*.txt'));
dir_output = [data_root, '\Analysis results\HMM 2 and 3 states 2 partners'];

for j=2:numel(behavior_files)
    animal_code = strrep(behavior_files(j).name, '.txt', '');

    disp('###################################')
    disp(['##### LOADING ', animal_code, '######'])
    disp('###################################')
    Load_HMM_outputs_and_analyze(dir_output,animal_code)
    close all
end
