%% Estimate bout-resolved PPC and delta-power controls for all sessions
% Run add_repo_paths.m before this script. The output is intentionally
% separate from the original Figure 3 phase_struct files.

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_raw_data = fullfile(data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(data_root, 'Analysis results', 'phase locking power control');
if ~isfolder(saving_folder)
    mkdir(saving_folder);
end

sample_rate = 2500;
filter_order = 2000;
delta_band = [1 5];
Hd_delta = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', delta_band(1), ...
    'CutoffFrequency2', delta_band(2), ...
    'SampleRate', sample_rate, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');

% Hardcoded analysis parameters. Keep these with the saved output.
generator_options.min_bout_duration = 0.25;
generator_options.phase_method = 'uniformized'; % matches Figure 3 pipeline
generator_options.epoch_duration = 1;
RECOMPUTE_EXISTING = false;

animal_list = dir(npx_raw_data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.','..'}));
animal2exclude = {''};
animal_list = animal_list(~ismember({animal_list.name}, animal2exclude));

% This run stores only sessions processed now. Existing per-session files and
% phase_locking_power_control_all_sessions.mat are left unchanged.
power_control_struct = struct([]);
failed_sessions = table('Size', [0 2], ...
    'VariableTypes', {'string','string'}, ...
    'VariableNames', {'Animal','Message'});
skipped_sessions = strings(0,1);

for fn = 1:numel(animal_list)
    animal_name = animal_list(fn).name;
    session_file = fullfile(saving_folder, ...
        [animal_name, '_phase_locking_power_control.mat']);
    fprintf('\nSession %d/%d: %s\n', fn, numel(animal_list), animal_name);
    if ~RECOMPUTE_EXISTING && isfile(session_file)
        fprintf('Keeping existing file: %s\n', session_file);
        skipped_sessions(end+1,1) = string(animal_name); %#ok<SAGROW>
        continue
    end
    try
        this_control = GENERATE_PHASE_LOCKING_POWER_CONTROL( ...
            fullfile(npx_raw_data, animal_name), Hd_delta, generator_options);
        % Save before aggregation. struct([]) has no fields, so the first
        % subscripted assignment would otherwise discard a completed session.
        save(session_file, 'this_control', '-v7.3');
        if isempty(power_control_struct)
            power_control_struct = this_control;
        else
            power_control_struct(end+1) = orderfields( ...
                this_control, power_control_struct); %#ok<SAGROW>
        end
    catch exception
        warning('Power-control analysis failed for %s: %s', ...
            animal_name, exception.message);
        failed_sessions(end+1,:) = ...
            {string(animal_name), string(getReport(exception, 'extended'))};
    end
end

retry_file = fullfile(saving_folder, ...
    'phase_locking_power_control_retry_sessions.mat');
save(retry_file, 'power_control_struct', 'failed_sessions', ...
    'skipped_sessions', 'generator_options', 'delta_band', '-v7.3');
fprintf(['Kept %d existing sessions; added %d; %d failed.\n' ...
    'Original combined file was not modified.\nNew sessions: %s\n'], ...
    numel(skipped_sessions), numel(power_control_struct), ...
    height(failed_sessions), retry_file);
