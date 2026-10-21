%% Estimate synchronous-spike waveforms for play-bout coincidences
% Run this first. It loops NPX sessions, calls GENERATE_SYNC_SPIKE_WAVEFORMS,
% and saves one mat file per session. Then run Analyze_sync_spike_waveforms.
%
% Each neuron is read once, and only the spikes inside the play-bout
% windows are read. The run stops at max_pairs (default 80), with at most
% 8 pairs in a session and 2 pairs for a neuron. The first interval read
% is timed. If it exceeds max_first_load_seconds (default 10 minutes),
% that session is not saved and later sessions are not started, so the
% pair cap can be lowered. Re-run with RECOMPUTE_EXISTING = true after
% this pool change; older mats used the session-wide spike pool.
%
% Results:
%   Data/Analysis results/synchronous spike waveforms/
%     {session}_sync_spike_waveforms.mat
%     sync_spike_waveforms_all_sessions.mat

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    run(fullfile(repo_root, 'add_repo_paths.m'));
end

generator_options = sync_analysis_options();
rng('shuffle');
seed_state = rng;
generator_options.random_seed = seed_state.Seed;
RECOMPUTE_EXISTING = false;
sessions_to_skip = {};

data_root = generator_options.data_root;
npx_raw_data = fullfile(data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(data_root, 'Analysis results', 'synchronous spike waveforms');
if exist(npx_raw_data, 'dir') ~= 7
    error('Estimate_sync_spike_waveforms:data', ...
        'NPX raw data folder was not found:\n%s', npx_raw_data);
end
if exist(saving_folder, 'dir') ~= 7
    mkdir(saving_folder);
end

fprintf(['Synchronous-spike waveforms. Seed %d. Up to %d pairs, ' ...
    '%d per session, %d per neuron, %d null draws.\n%s\n'], ...
    generator_options.random_seed, generator_options.max_pairs, ...
    generator_options.max_pairs_per_session, generator_options.max_pairs_per_neuron, ...
    generator_options.n_repeats, saving_folder);

animal_list = dir(npx_raw_data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));
animal_list = animal_list(~ismember({animal_list.name}, sessions_to_skip));

failed_sessions = table('Size', [0 2], ...
    'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'Session', 'Message'});
skipped_sessions = strings(0, 1);
remaining = generator_options.max_pairs;

if ~RECOMPUTE_EXISTING
    for fn = 1:numel(animal_list)
        session_file = fullfile(saving_folder, ...
            [animal_list(fn).name, '_sync_spike_waveforms.mat']);
        if exist(session_file, 'file') ~= 2
            continue
        end
        try
            loaded = load(session_file, 'this_waveforms');
            remaining = remaining - loaded.this_waveforms.n_pairs;
        catch load_error
            warning('Estimate_sync_spike_waveforms:existing', ...
                'Could not read %s: %s', session_file, load_error.message);
        end
    end
end

for fn = 1:numel(animal_list)
    animal_name = animal_list(fn).name;
    session_file = fullfile(saving_folder, [animal_name, '_sync_spike_waveforms.mat']);
    fprintf('\nSession %d/%d: %s\n', fn, numel(animal_list), animal_name);
    if ~RECOMPUTE_EXISTING && exist(session_file, 'file') == 2
        fprintf('Keeping existing file: %s\n', session_file);
        skipped_sessions(end + 1, 1) = string(animal_name); %#ok<SAGROW>
        continue
    end
    if remaining < 1
        fprintf('Pair budget is filled (%d).\n', generator_options.max_pairs);
        break
    end
    call_options = generator_options;
    call_options.max_pairs = remaining;
    try
        this_waveforms = GENERATE_SYNC_SPIKE_WAVEFORMS( ...
            fullfile(npx_raw_data, animal_name), call_options);
        save(session_file, 'this_waveforms', '-v7.3');
        remaining = remaining - this_waveforms.n_pairs;
        fprintf('Saved %s (%d pairs, %d neuron rows). Remaining pair budget %d.\n', ...
            session_file, this_waveforms.n_pairs, numel(this_waveforms.rows), remaining);
    catch exception
        warning('Waveform analysis failed for %s: %s', animal_name, exception.message);
        failed_sessions(end + 1, :) = ...
            {string(animal_name), string(getReport(exception, 'extended'))}; %#ok<SAGROW>
        if strcmp(exception.identifier, 'GENERATE_SYNC_SPIKE_WAVEFORMS:loadTime')
            fprintf(['Stopped because the first neuron read was too slow. ' ...
                'Lower max_pairs or raise max_first_load_seconds, then run again.\n']);
            break
        end
    end
end

session_files = dir(fullfile(saving_folder, '*_sync_spike_waveforms.mat'));
sync_waveform_sessions = struct([]);
for file_index = 1:numel(session_files)
    session_path = fullfile(session_files(file_index).folder, session_files(file_index).name);
    try
        loaded = load(session_path, 'this_waveforms');
    catch load_error
        warning('Estimate_sync_spike_waveforms:combine', ...
            'Skipping %s: %s', session_path, load_error.message);
        continue
    end
    if ~isfield(loaded, 'this_waveforms')
        continue
    end
    if isempty(sync_waveform_sessions)
        sync_waveform_sessions = loaded.this_waveforms;
    else
        sync_waveform_sessions(end + 1) = orderfields( ...
            loaded.this_waveforms, sync_waveform_sessions); %#ok<SAGROW>
    end
end
combined_file = fullfile(saving_folder, 'sync_spike_waveforms_all_sessions.mat');
save(combined_file, 'sync_waveform_sessions', 'failed_sessions', ...
    'skipped_sessions', 'generator_options', '-v7.3');
fprintf('Combined %d sessions in %s. Failed %d. Skipped %d. Remaining budget %d.\n', ...
    numel(sync_waveform_sessions), combined_file, height(failed_sessions), ...
    numel(skipped_sessions), remaining);
