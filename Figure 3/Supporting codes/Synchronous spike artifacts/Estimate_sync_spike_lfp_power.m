%% Estimate raw-LFP power at synchronous play-bout spikes
% Run this first. It loops NPX sessions, calls GENERATE_SYNC_SPIKE_LFP_POWER,
% and saves one mat file per session. Then run Analyze_sync_spike_lfp_power.
%
% Power at a synchronous spike is the nearest sample of
% abs(hilbert(raw LFP)).^2 on that neuron's channel. The raw trace is not
% delta-filtered. Phase, used only to build the null, is the analytic angle
% of the same 1-5 Hz FIR as the coincidence pipeline. Each of 500 nulls is
% the mean of one random LFP sample per synchronous event, drawn inside the
% play-bout windows from the same 10-degree delta-phase bin. Samples within
% 10 ms of that pair's synchronous spikes are left out. Pairs follow the
% waveform sample: at least 90% of neurons with 1000 spikes in each of SC,
% DLPAG, LPAG, VLPAG, and DR, and no extra pair among neurons already
% observed. Each LFP channel is read once per session. Mats saved before
% this phase-matched null still compare against other spikes. Re-run with
% RECOMPUTE_EXISTING = true to rebuild them.
%
% Results:
%   Data/Analysis results/synchronous spike waveforms/
%     {session}_sync_spike_lfp_power.mat
%     sync_spike_lfp_power_all_sessions.mat

repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
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

fprintf(['Synchronous-spike raw-LFP power. Seed %d. Target %.0f%% of neurons ' ...
    'with >= %d spikes in each area. %d null draws in %d delta-phase bins.\n%s\n'], ...
    generator_options.random_seed, 100 * generator_options.min_neuron_fraction, ...
    generator_options.min_spikes, generator_options.n_repeats, ...
    generator_options.n_phase_bins, saving_folder);

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
        fprintf('Saved %s (%d pairs, %d neuron rows, %d/%d spiking neurons).\n', ...
            session_file, this_lfp_power.n_pairs, numel(this_lfp_power.rows), ...
            this_lfp_power.coverage.n_covered, this_lfp_power.coverage.n_with_spikes);
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
    if isfield(loaded.this_lfp_power, 'rows') && ~isempty(loaded.this_lfp_power.rows) ...
            && ~isfield(loaded.this_lfp_power.rows, 'mean_power_phase_matched')
        warning('Estimate_sync_spike_lfp_power:oldFile', ...
            ['%s still compares synchronous spikes with other spikes. ' ...
             'Set RECOMPUTE_EXISTING = true and run again so the null is ' ...
             'phase-matched LFP inside the play-bout windows.'], ...
            session_path);
    end
    if ~isfield(loaded.this_lfp_power, 'coverage') || isempty(loaded.this_lfp_power.coverage)
        [~, ~, loaded.this_lfp_power.coverage] = select_synchronous_pairs( ...
            struct('cluster_id', {}, 'channel', {}, 'area', {}, ...
            'spike_times', {}, 'spike_samples', {}, 'interval_index', {}), [], struct());
    end
    sync_lfp_power_sessions = append_struct_array( ...
        sync_lfp_power_sessions, loaded.this_lfp_power);
end
combined_file = fullfile(saving_folder, 'sync_spike_lfp_power_all_sessions.mat');
save(combined_file, 'sync_lfp_power_sessions', 'failed_sessions', ...
    'skipped_sessions', 'generator_options', '-v7.3');
fprintf('Combined %d sessions in %s. Failed %d. Skipped %d.\n', ...
    numel(sync_lfp_power_sessions), combined_file, height(failed_sessions), ...
    numel(skipped_sessions));
if ~isempty(sync_lfp_power_sessions) && isfield(sync_lfp_power_sessions, 'coverage')
    summarize_sync_coverage([sync_lfp_power_sessions.coverage], generator_options.min_spikes);
end
