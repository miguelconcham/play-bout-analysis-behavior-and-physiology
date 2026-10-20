%% Estimate raw-LFP power at synchronous play-bout spikes
% Run this first. It loops NPX sessions, calls GENERATE_SYNC_SPIKE_LFP_POWER,
% and saves one mat file per session. Then run Analyze_sync_spike_lfp_power.
%
% Power at a spike is one sample of abs(hilbert(raw LFP)).^2 on that
% neuron's channel. The raw trace is not delta-filtered. The same pair
% caps as the waveform estimate apply (80 pairs, 8 per session, 2 per
% neuron). Each LFP channel is read once per session.
%
% Results:
%   Data/Analysis results/synchronous spike waveforms/
%     {session}_sync_spike_lfp_power.mat
%     sync_spike_lfp_power_all_sessions.mat

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
    error('Estimate_sync_spike_lfp_power:data', ...
        'NPX raw data folder was not found:\n%s', npx_raw_data);
end
if exist(saving_folder, 'dir') ~= 7
    mkdir(saving_folder);
end

fprintf(['Synchronous-spike raw-LFP power. Seed %d. Up to %d pairs, ' ...
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
            [animal_list(fn).name, '_sync_spike_lfp_power.mat']);
        if exist(session_file, 'file') ~= 2
            continue
        end
        try
            loaded = load(session_file, 'this_lfp_power');
            remaining = remaining - loaded.this_lfp_power.n_pairs;
        catch load_error
            warning('Estimate_sync_spike_lfp_power:existing', ...
                'Could not read %s: %s', session_file, load_error.message);
        end
    end
end

for fn = 1:numel(animal_list)
    animal_name = animal_list(fn).name;
    session_file = fullfile(saving_folder, [animal_name, '_sync_spike_lfp_power.mat']);
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
        this_lfp_power = GENERATE_SYNC_SPIKE_LFP_POWER( ...
            fullfile(npx_raw_data, animal_name), call_options);
        save(session_file, 'this_lfp_power', '-v7.3');
        remaining = remaining - this_lfp_power.n_pairs;
        fprintf('Saved %s (%d pairs, %d neuron rows). Remaining pair budget %d.\n', ...
            session_file, this_lfp_power.n_pairs, numel(this_lfp_power.rows), remaining);
    catch exception
        warning('LFP-power analysis failed for %s: %s', animal_name, exception.message);
        failed_sessions(end + 1, :) = ...
            {string(animal_name), string(getReport(exception, 'extended'))}; %#ok<SAGROW>
    end
end

session_files = dir(fullfile(saving_folder, '*_sync_spike_lfp_power.mat'));
sync_lfp_power_sessions = struct([]);
for file_index = 1:numel(session_files)
    session_path = fullfile(session_files(file_index).folder, session_files(file_index).name);
    try
        loaded = load(session_path, 'this_lfp_power');
    catch load_error
        warning('Estimate_sync_spike_lfp_power:combine', ...
            'Skipping %s: %s', session_path, load_error.message);
        continue
    end
    if ~isfield(loaded, 'this_lfp_power')
        continue
    end
    if isempty(sync_lfp_power_sessions)
        sync_lfp_power_sessions = loaded.this_lfp_power;
    else
        sync_lfp_power_sessions(end + 1) = orderfields( ...
            loaded.this_lfp_power, sync_lfp_power_sessions); %#ok<SAGROW>
    end
end
combined_file = fullfile(saving_folder, 'sync_spike_lfp_power_all_sessions.mat');
save(combined_file, 'sync_lfp_power_sessions', 'failed_sessions', ...
    'skipped_sessions', 'generator_options', '-v7.3');
fprintf('Combined %d sessions in %s. Failed %d. Skipped %d. Remaining budget %d.\n', ...
    numel(sync_lfp_power_sessions), combined_file, height(failed_sessions), ...
    numel(skipped_sessions), remaining);
