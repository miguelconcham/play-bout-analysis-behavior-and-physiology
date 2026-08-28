%% Estimate_acute_phase_locking
% Loop every Acute data session folder (YYYYMMDD_N), skip "population analysis",
% call GENERATE_ACUTE_PHASE_LOCKING, collect structs, and save.

%% Paths, filter, and output folder
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

acute_root     = fullfile(data_root, 'Acute data');
saving_folder  = fullfile(data_root, 'Analysis results', 'acute phase locking');
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
freq_range   = [1 5];   % LFP band passed to GENERATE (e.g. [1 5] delta, [6 12] theta)
freq_tag     = sprintf('%g-%gHz', freq_range(1), freq_range(2));
Hd_freq = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range(1), ...
    'CutoffFrequency2', freq_range(2), ...
    'SampleRate', sr);

%% Generate all sessions
acute_struct  = struct([]);
session_names = {};
tic
for fn = 1:numel(session_list)
    session_folder = session_list(fn).name;
    npx_data_dir   = fullfile(acute_root, session_folder);
    fig_folder     = fullfile(saving_folder, session_folder);
    disp(['Session ', num2str(fn), ' of ', num2str(numel(session_list)), ': ', session_folder])

    if ~exist(fullfile(npx_data_dir, 'continuous.dat'), 'file') ...
            || ~exist(fullfile(npx_data_dir, 'cluster_info.tsv'), 'file')
        disp(['  missing continuous.dat or cluster_info.tsv, skipping ', session_folder])
        continue
    end

    try
        this_struct = GENERATE_ACUTE_PHASE_LOCKING(npx_data_dir, Hd_freq, fig_folder);
    catch ME
        warning('GENERATE_ACUTE_PHASE_LOCKING failed for %s: %s', session_folder, ME.message);
        toc
        continue
    end
    close all
    if isempty(this_struct)
        disp(['  empty result, skipping ', session_folder])
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

%% Save pooled structure
if isempty(acute_struct)
    error('Estimate_acute_phase_locking:empty', 'No sessions were processed.');
end
out_file = ['acute_phase_locking_', freq_tag, '.mat'];
save(fullfile(saving_folder, out_file), 'acute_struct', '-v7.3');
save(fullfile(saving_folder, 'acute_phase_locking_session_names.mat'), 'session_names');
disp(['saved ', fullfile(saving_folder, out_file)])
toc
