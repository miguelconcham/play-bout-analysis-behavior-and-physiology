%% Estimate_behavior_aligned_csd
% Export raw LFP matrices aligned to behavior/play-bout onset and offset.
%
% Run this script first on a machine with the experimentFS Data share. It
% calls GENERATE_BEHAVIOR_ALIGNED_CSD_LFP once per NPX session and writes
% CSD-ready NumPy files plus MATLAB event/QC metadata.
%
% Output root:
%   Data/Analysis results/behavior aligned CSD/
%
% Main branches:
%   behavior/<type>/onset|offset/<animal>/[probeN]/
%   category/play|aggression|exploration/onset|offset/<animal>/[probeN]/
%   play_role/self_play|other_play/onset|offset/<animal>/[probeN]/
%   play_bout/all/onset|offset/<animal>/[probeN]/
%
% Raw lfp.npy is unfiltered, unreferenced, and not baseline corrected.
% evoked_lfp.npy is the LFP PSTH after applying each event's onset baseline
% to both its onset and offset trials.

data_root     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_raw_data  = fullfile(data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(data_root, 'Analysis results', 'behavior aligned CSD');

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\CSD analysis\Estimate_behavior_aligned_csd.m';
end
repo_root = fileparts(fileparts(this_file));
run(fullfile(repo_root, 'add_repo_paths.m'));
addpath(fullfile(repo_root, 'CSD analysis'));

if ~exist(saving_folder, 'dir')
    mkdir(saving_folder);
end

%% Fixed analysis parameters
% The Python paired loader uses the same defaults. Both snippets are raw;
% baseline subtraction is applied later using only the onset pre-window.
t1 = -0.50;
t2 =  0.50;
recommended_min_events = 5;

%% Optional session exclusions
animal2exclude = {''};

animal_list = dir(npx_raw_data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));
animal_list = animal_list(~ismember({animal_list.name}, animal2exclude));

summary = struct([]);
tic
for fn = 1:numel(animal_list)
    animal_name = animal_list(fn).name;
    npx_dir = fullfile(npx_raw_data, animal_name);
    disp('============================================================')
    disp(sprintf('[%d/%d] Behavior-aligned CSD export: %s', ...
        fn, numel(animal_list), animal_name))

    try
        session_summary = GENERATE_BEHAVIOR_ALIGNED_CSD_LFP( ...
            npx_dir, saving_folder, t1, t2, recommended_min_events);
        if isempty(session_summary)
            warning('Estimate_behavior_aligned_csd:Empty', ...
                '%s produced no datasets', animal_name);
            continue
        end
        if isempty(summary)
            summary = session_summary;
        else
            summary = [summary, session_summary]; %#ok<AGROW>
        end
    catch ME
        warning('Estimate_behavior_aligned_csd:SessionFailed', ...
            '%s failed: %s', animal_name, ME.message);
        error_log_file = fullfile(saving_folder, 'estimate_errors.mat');
        if exist(error_log_file, 'file')
            loaded_errors = load(error_log_file, 'errors');
            errors = loaded_errors.errors;
        else
            errors = struct('animal', {}, 'message', {}, 'identifier', {}, 'stack', {});
        end
        error_entry.animal = animal_name;
        error_entry.message = ME.message;
        error_entry.identifier = ME.identifier;
        error_entry.stack = ME.stack;
        errors(end + 1) = error_entry; %#ok<SAGROW>
        save(error_log_file, 'errors');
    end
end
elapsed_seconds = toc;

save(fullfile(saving_folder, 'estimate_summary.mat'), ...
    'summary', 't1', 't2', 'recommended_min_events', 'elapsed_seconds');

if ~isempty(summary)
    summary_table = struct2table(summary);
    writetable(summary_table, fullfile(saving_folder, 'estimate_summary.csv'));
    disp(sprintf('Saved %d dataset/probe summary rows', height(summary_table)));
    disp(sprintf('Valid exported event pairs: %d', sum(summary_table.n_saved)));
    disp(sprintf('Low-event-count rows (<%d): %d', recommended_min_events, ...
        sum(summary_table.low_event_count)));
end
disp(['Done. Outputs: ', saving_folder])
