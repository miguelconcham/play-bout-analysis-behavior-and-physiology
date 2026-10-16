%% Estimate bout-resolved PPC and delta-power controls for all sessions
% Run add_repo_paths.m before this script. The output is intentionally
% separate from the original Figure 3 phase_struct files.

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
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
generator_options.n_rand = 0; % no PPC significance test for this control
generator_options.min_bout_duration = 0.25;
generator_options.phase_method = 'uniformized'; % matches Figure 3 pipeline
generator_options.epoch_duration = 1;

animal_list = dir(npx_raw_data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.','..'}));
animal2exclude = {''};
animal_list = animal_list(~ismember({animal_list.name}, animal2exclude));

power_control_struct = struct([]);
failed_sessions = table('Size', [0 2], ...
    'VariableTypes', {'string','string'}, ...
    'VariableNames', {'Animal','Message'});

for fn = 1:numel(animal_list)
    animal_name = animal_list(fn).name;
    fprintf('\nSession %d/%d: %s\n', fn, numel(animal_list), animal_name);
    try
        this_control = GENERATE_PHASE_LOCKING_POWER_CONTROL( ...
            fullfile(npx_raw_data, animal_name), Hd_delta, generator_options);
        power_control_struct(end+1) = this_control; %#ok<SAGROW>
        save(fullfile(saving_folder, ...
            [animal_name, '_phase_locking_power_control.mat']), ...
            'this_control', '-v7.3');
    catch exception
        warning('Power-control analysis failed for %s: %s', ...
            animal_name, exception.message);
        failed_sessions(end+1,:) = ...
            {string(animal_name), string(getReport(exception, 'extended'))};
    end
end

save(fullfile(saving_folder, 'phase_locking_power_control_all_sessions.mat'), ...
    'power_control_struct', 'failed_sessions', 'generator_options', ...
    'delta_band', '-v7.3');
fprintf('Saved %d sessions; %d failed.\n', ...
    numel(power_control_struct), height(failed_sessions));
