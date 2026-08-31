%% Estimate_stim_call_responses
% Loop Acute data sessions (YYYYMMDD_N), call GENERATE_STIM_CALL_RESPONSES,
% collect structs, save next to the acute phase-locking results.

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

acute_root    = fullfile(data_root, 'Acute data');
saving_folder = fullfile(data_root, 'Analysis results', 'acute phase locking');
if ~exist(saving_folder, 'dir')
    mkdir(saving_folder);
end

session_list = dir(acute_root);
session_list = session_list([session_list.isdir]);
session_list = session_list(~ismember({session_list.name}, {'.', '..', 'population analysis'}));
keep = ~cellfun(@isempty, regexp({session_list.name}, '^\d+_\d+$', 'once'));
session_list = session_list(keep);

sr           = 2500;
filter_order = 500;
freq_range   = [1 5];
freq_tag     = sprintf('%g-%gHz', freq_range(1), freq_range(2));
Hd_freq = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range(1), ...
    'CutoffFrequency2', freq_range(2), ...
    'SampleRate', sr);

acute_struct  = struct([]);
session_names = {};
tic
for fn = 1:numel(session_list)
    session_folder = session_list(fn).name;
    npx_data_dir   = fullfile(acute_root, session_folder);
    disp(['Session ', num2str(fn), ' of ', num2str(numel(session_list)), ': ', session_folder])

    if ~exist(fullfile(npx_data_dir, 'continuous.dat'), 'file') ...
            || ~exist(fullfile(npx_data_dir, 'cluster_info.tsv'), 'file')
        disp(['  missing continuous.dat or cluster_info.tsv, skipping ', session_folder])
        continue
    end

    try
        this_struct = GENERATE_STIM_CALL_RESPONSES(npx_data_dir, Hd_freq, '');
    catch ME
        warning('GENERATE_STIM_CALL_RESPONSES failed for %s: %s', session_folder, ME.message);
        toc
        continue
    end
    if isempty(this_struct)
        continue
    end
    if isempty(acute_struct)
        acute_struct = this_struct;
    else
        acute_struct(end + 1) = this_struct; %#ok<AGROW>
    end
    session_names = [session_names; {session_folder}]; %#ok<AGROW>
    toc
end

if isempty(acute_struct)
    error('Estimate_stim_call_responses:empty', 'No sessions were processed.');
end
out_file = ['stim_call_responses_', freq_tag, '.mat'];
save(fullfile(saving_folder, out_file), 'acute_struct', '-v7.3');
save(fullfile(saving_folder, 'stim_call_responses_session_names.mat'), 'session_names');
disp(['saved ', fullfile(saving_folder, out_file)])
toc
