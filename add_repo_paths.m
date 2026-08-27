%% add_repo_paths
% Put shared helpers and GENERATE functions on the MATLAB path.
% Run once per MATLAB session (or add this file to your startup).

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m';
end
repo_root = fileparts(this_file);

addpath(genpath(fullfile(repo_root, 'Custom functions')));
addpath(fullfile(repo_root, 'Figure 1', 'LDA analysis'));
addpath(fullfile(repo_root, 'Figure 1', 'HMM modeling'));
addpath(fullfile(repo_root, 'Figure 2', 'Supporting codes'));
addpath(fullfile(repo_root, 'Figure 2', 'Supporting codes v2'));
addpath(fullfile(repo_root, 'Figure 3', 'Supporting codes'));
addpath(fullfile(repo_root, 'Figure 6', 'Supplementary codes'));
addpath(fullfile(repo_root, 'Complementary figures'));

disp(['Repo paths added from ', repo_root])
