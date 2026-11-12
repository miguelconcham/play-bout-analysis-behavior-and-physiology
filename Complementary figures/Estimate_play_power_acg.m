%% Estimate play-versus-power autocorrelograms for all sessions
% Epoch autocorrelograms during play, during power-matched time outside
% play, and in the low and high delta-power tails of play. Run
% add_repo_paths.m before this script. Analyze_play_power_acg.m reads the
% per-session files written here.
%
% acg_method 'overlap' is the default: pair counts / (T-|lag|) and the
% original file name. 'eligible' keeps one observation per spike that can
% see the whole lag bin and writes {animal}_play_power_acg_eligible.mat,
% so an eligible estimate does not replace an overlap file. Set the same
% value in Analyze_play_power_acg.m before plotting.

data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_raw_data = fullfile(data_root, 'NPX data', 'NPX raw data');
saving_folder = fullfile(data_root, 'Analysis results', 'play power acg');
if exist(saving_folder, 'dir') ~= 7
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

generator_options.epoch_duration = 1;
generator_options.acg_lag = 1;
generator_options.acg_bin = 0.01;
generator_options.sample_rate = sample_rate;
generator_options.delta_band = delta_band;
generator_options.low_percent = 25;
generator_options.high_percent = 75;
generator_options.caliper_sd = 0.25;
generator_options.areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
generator_options.session_type = 'Partners session';
% 'eligible' or 'overlap'. Overlap is the default. Eligible files are not
% written over overlap files.
acg_method = 'overlap';
generator_options.acg_method = acg_method;
RECOMPUTE_EXISTING = false;

animal_list = dir(npx_raw_data);
animal_list = animal_list([animal_list.isdir]);
animal_list = animal_list(~ismember({animal_list.name}, {'.', '..'}));
animal2exclude = {''};
animal_list = animal_list(~ismember({animal_list.name}, animal2exclude));

failed_sessions = table('Size', [0 2], ...
    'VariableTypes', {'string', 'string'}, ...
    'VariableNames', {'Animal', 'Message'});
skipped_sessions = strings(0, 1);
n_added = 0;

for fn = 1:numel(animal_list)
    animal_name = animal_list(fn).name;
    session_file = fullfile(saving_folder, play_power_acg_filename(animal_name, acg_method));
    fprintf('\nSession %d/%d: %s\n', fn, numel(animal_list), animal_name);
    if ~RECOMPUTE_EXISTING && exist(session_file, 'file') == 2
        fprintf('Keeping existing file: %s\n', session_file);
        skipped_sessions(end + 1, 1) = string(animal_name); %#ok<SAGROW>
        continue
    end
    try
        play_power_acg = GENERATE_PLAY_POWER_ACG( ...
            fullfile(npx_raw_data, animal_name), Hd_delta, generator_options);
        save(session_file, 'play_power_acg', '-v7.3');
        n_added = n_added + 1;
        clear play_power_acg
    catch exception
        warning('Play-power ACG failed for %s: %s', animal_name, exception.message);
        failed_sessions(end + 1, :) = ...
            {string(animal_name), string(getReport(exception, 'extended', 'hyperlinks', 'off'))}; %#ok<SAGROW>
    end
end

retry_file = fullfile(saving_folder, play_power_retry_filename(acg_method));
save(retry_file, 'failed_sessions', 'skipped_sessions', ...
    'generator_options', 'delta_band', 'n_added');
fprintf('Kept %d existing sessions; added %d; %d failed.\nLog: %s\n', ...
    numel(skipped_sessions), n_added, height(failed_sessions), retry_file);
