function spectra_struct = GENERATE_PLAY_SESSION_LPAG_SPECTRA(current_dir)
% GENERATE_PLAY_SESSION_LPAG_SPECTRA
% Estimate mid-LPAG spectra during complete partner-play sessions.
%
% Play sessions are defined only by behavior-table rows with:
%   Animal == 'Session_structure' and Type == 'Partners session'
% This function does not detect or analyze play bouts.
%
% Each play-session interval is split into non-overlapping 10 s windows.
% A Welch PSD and a direct Fourier periodogram are estimated for every
% window, averaged in linear-power space, and then independently
% parameterized with FOOOF.

%% Analysis parameters
sampling_rate       = 2500;
window_seconds      = 10;
window_samples      = window_seconds * sampling_rate;
frequency           = (0.1:0.1:40).';
fooof_peak_widths   = [0.5 12];
fooof_aperiodic_mode = 'fixed';
area2analyze        = 'PAG';

%% Locate repository data and identify the recording
this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
repository_data_root = fullfile(repo_root, 'Data');
recording_data_root = fileparts(fileparts(fileparts(current_dir)));
if isfile(fullfile(recording_data_root, 'Area_limits_GoodLooking.xlsx'))
    data_root = recording_data_root;
else
    data_root = repository_data_root;
end

[~, animal_code] = fileparts(current_dir);
animal_code_params = strsplit(animal_code, ' ');
if numel(animal_code_params) < 3
    error('Unexpected Neuropixels folder name: %s', animal_code);
end
animal_batch    = animal_code_params{1};
repeated_animal = animal_code_params{3};

%% Load synchronization model
local_synch_file = fullfile(current_dir, 'synch_model_video2NPX.mat');
data_synch_file = fullfile(data_root, 'Synch data', animal_code, ...
    'synch_model_video2NPX.mat');
if isfile(local_synch_file)
    synch_file = local_synch_file;
elseif isfile(data_synch_file)
    synch_file = data_synch_file;
else
    error('Synchronization model not found for %s.', animal_code);
end
load(synch_file, 'synch_model_video2NPX');

%% Load behavior and obtain the explicitly annotated play sessions
local_behavior_files = dir(fullfile(current_dir, 'ELAN*'));
data_behavior_file = fullfile(data_root, 'Behavior backups', ...
    [animal_code, '.txt']);
if ~isempty(local_behavior_files)
    behavior_file = fullfile(current_dir, local_behavior_files(1).name);
elseif isfile(data_behavior_file)
    behavior_file = data_behavior_file;
else
    error('Behavior file not found for %s.', animal_code);
end

Behavior = readtable(behavior_file);
Behavior(:, 2) = [];
Behavior.Properties.VariableNames = ...
    {'Animal', 'Start', 'End', 'Length', 'Type'};
Behavior.Start = predict(synch_model_video2NPX, Behavior.Start);
Behavior.End   = predict(synch_model_video2NPX, Behavior.End);

session_mask = ismember(Behavior.Animal, 'Session_structure') & ...
    ismember(Behavior.Type, 'Partners session');
play_session_table = Behavior(session_mask, :);
if isempty(play_session_table)
    error(['No Session_structure / Partners session rows found for %s. ' ...
        'No fallback interval is used.'], animal_code);
end

play_session_ranges = play_session_table{:, {'Start', 'End'}};
[play_session_ranges, session_order] = sortrows(play_session_ranges, 1);
play_session_table = play_session_table(session_order, :);
if any(~isfinite(play_session_ranges), 'all') || ...
        any(play_session_ranges(:, 2) <= play_session_ranges(:, 1))
    error('Invalid play-session start/end times for %s.', animal_code);
end

%% Load LFP
disp(['Loading LFP: ', animal_code])
lfp_mat_file = fullfile(current_dir, 'LFP_PAG.mat');
lfp_dat_file = fullfile(current_dir, 'LFP_PAG.dat');
if isfile(lfp_mat_file)
    NPX_Type = 2;
    loaded_lfp = load(lfp_mat_file, 'LFP');
    if ~isfield(loaded_lfp, 'LFP')
        error('LFP_PAG.mat does not contain a variable named LFP.');
    end
    LFP = loaded_lfp.LFP;
elseif isfile(lfp_dat_file)
    NPX_Type = 1;
    file_pointer = fopen(lfp_dat_file, 'r');
    if file_pointer == -1
        error('Could not open %s.', lfp_dat_file);
    end
    cleaner = onCleanup(@() fclose(file_pointer));
    LFP = fread(file_pointer, 'int16=>int16');
    if mod(numel(LFP), 384) ~= 0
        error('LFP_PAG.dat does not contain a whole number of 384-channel samples.');
    end
    LFP = reshape(LFP, 384, []);
    clear cleaner
else
    error('Neither LFP_PAG.mat nor LFP_PAG.dat exists in %s.', current_dir);
end
disp('LFP loaded')

if size(LFP, 1) ~= 384
    error('Expected 384 LFP channels; found %d.', size(LFP, 1));
end

%% Select the middle LPAG channel(s) from probe geometry and area limits
channel_map_file = fullfile(current_dir, ['chann_map_', area2analyze, '.mat']);
if ~isfile(channel_map_file)
    error('Channel map not found: %s', channel_map_file);
end
load(channel_map_file, 'xcoords', 'ycoords', 'chanMap');
xcoords = xcoords(:);
ycoords = ycoords(:);
chanMap = chanMap(:);

area_limit_file = fullfile(data_root, 'Area_limits_GoodLooking.xlsx');
if ~isfile(area_limit_file)
    error('Area-limit table not found: %s', area_limit_file);
end
area_limit = readtable(area_limit_file);

if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    if numel(animal_batch) < 4
        error('Cannot derive area-table animal name from %s.', animal_code);
    end
    this_animal = ...
        ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);
lpag_rows = ismember(area_limit.area, 'LPAG');
Y_Range = area_limit{lpag_rows, {'ProbeNum', 'depth_start', 'depth_end'}};
if isempty(Y_Range)
    error('No LPAG limits found for %s.', this_animal);
end

% Preserve the missing-channel correction used by the existing PAG code.
if NPX_Type == 1 && ~ismember(192, chanMap)
    pos_191 = find(chanMap == 191, 1);
    pos_193 = find(chanMap == 193, 1);
    if isempty(pos_191) || isempty(pos_193) || pos_193 ~= pos_191 + 1
        error('Inconsistent NPX1 channel map around missing channel 192.');
    end
    xcoords = [xcoords; NaN];
    xcoords(pos_193 + 1:end) = xcoords(pos_193:end - 1);
    xcoords(pos_193) = 43;
    ycoords = [ycoords; NaN];
    ycoords(pos_193 + 1:end) = ycoords(pos_193:end - 1);
    ycoords(pos_193) = 1900;
    chanMap = [chanMap; NaN];
    chanMap(pos_193 + 1:end) = chanMap(pos_193:end - 1);
    chanMap(pos_193) = 192;
end

hard_coded_x_coords = [8 40; 258 290; 508 540; 758 790];
mid_lpag_channel = nan(size(Y_Range, 1), 1);
mid_lpag_lfp_row = nan(size(Y_Range, 1), 1);
mid_lpag_map_index = nan(size(Y_Range, 1), 1);

for probe_n = 1:size(Y_Range, 1)
    if NPX_Type == 1
        in_lpag = ycoords >= Y_Range(probe_n, 2) & ...
            ycoords <= Y_Range(probe_n, 3);
    else
        probe_number = Y_Range(probe_n, 1);
        if probe_number < 1 || probe_number > size(hard_coded_x_coords, 1)
            error('Invalid probe number %g in LPAG area limits.', probe_number);
        end
        in_lpag = ycoords >= Y_Range(probe_n, 2) & ...
            ycoords <= Y_Range(probe_n, 3) & ...
            ismember(xcoords, hard_coded_x_coords(probe_number, :));
    end

    lpag_locations = [xcoords(in_lpag), ycoords(in_lpag)];
    if isempty(lpag_locations)
        error('No mapped channels fall within LPAG limits for probe %d.', probe_n);
    end
    mean_location = mean(lpag_locations, 1);
    distance = sum(abs([xcoords, ycoords] - mean_location), 2);
    distance(~in_lpag) = Inf;
    [~, closest_map_index] = min(distance);
    selected_channel = chanMap(closest_map_index);
    selected_lfp_row = selected_channel + 1;
    if ~isfinite(selected_channel) || selected_lfp_row < 1 || ...
            selected_lfp_row > size(LFP, 1)
        error('Selected LPAG channel %g is not a valid LFP channel.', ...
            selected_channel);
    end
    mid_lpag_channel(probe_n) = selected_channel;
    mid_lpag_lfp_row(probe_n) = selected_lfp_row;
    mid_lpag_map_index(probe_n) = closest_map_index;
end

%% Create full 10-second windows within each play-session row
n_lfp_samples = size(LFP, 2);
session_windows = repmat(struct( ...
    'sample_ranges', zeros(0, 2), ...
    'time_ranges', zeros(0, 2), ...
    'n_windows', 0), size(play_session_ranges, 1), 1);

for session_n = 1:size(play_session_ranges, 1)
    first_sample = max(1, ceil(play_session_ranges(session_n, 1) * sampling_rate));
    last_sample = min(n_lfp_samples, ...
        floor(play_session_ranges(session_n, 2) * sampling_rate));
    starts = first_sample:window_samples:(last_sample - window_samples + 1);
    if isempty(starts)
        error('Play-session row %d contains no complete 10-second window.', ...
            session_n);
    end
    ends = starts + window_samples - 1;
    session_windows(session_n).sample_ranges = [starts(:), ends(:)];
    session_windows(session_n).time_ranges = ...
        [starts(:) / sampling_rate, ends(:) / sampling_rate];
    session_windows(session_n).n_windows = numel(starts);
end

%% Estimate and parameterize spectra for each selected LPAG channel
n_frequency = numel(frequency);
spectra_struct = repmat(struct(), numel(mid_lpag_channel), 1);
welch_window = hann(window_samples, 'periodic');

for channel_n = 1:numel(mid_lpag_channel)
    lfp_row = mid_lpag_lfp_row(channel_n);
    spectra_struct(channel_n).animal_name = animal_code;
    spectra_struct(channel_n).animal_area_table_name = this_animal;
    spectra_struct(channel_n).neuropixels_directory = current_dir;
    spectra_struct(channel_n).npx_type = NPX_Type;
    spectra_struct(channel_n).lpag_channel = ...
        mid_lpag_channel(channel_n);
    spectra_struct(channel_n).lpag_lfp_row = lfp_row;
    spectra_struct(channel_n).lpag_channel_map_index = ...
        mid_lpag_map_index(channel_n);
    spectra_struct(channel_n).lpag_probe_number = Y_Range(channel_n, 1);
    spectra_struct(channel_n).lpag_depth_range = Y_Range(channel_n, 2:3);
    spectra_struct(channel_n).play_session_table = play_session_table;
    spectra_struct(channel_n).play_session_ranges = play_session_ranges;
    spectra_struct(channel_n).n_sessions = size(play_session_ranges, 1);
    spectra_struct(channel_n).n_windows_total = ...
        sum([session_windows.n_windows]);
    spectra_struct(channel_n).frequency = frequency;

    parameters = struct();
    parameters.sampling_rate = sampling_rate;
    parameters.window_seconds = window_seconds;
    parameters.window_samples = window_samples;
    parameters.frequency_range = frequency([1 end]).';
    parameters.frequency_resolution = sampling_rate / window_samples;
    parameters.window_overlap_seconds = 0;
    parameters.welch_window = 'periodic Hann';
    parameters.fourier_window = 'rectangular';
    parameters.detrend = 'constant';
    parameters.averaging_domain = 'linear power';
    parameters.fooof_aperiodic_mode = fooof_aperiodic_mode;
    parameters.fooof_peak_width_limits = fooof_peak_widths;
    spectra_struct(channel_n).parameters = parameters;

    session_results = repmat(struct(), size(play_session_ranges, 1), 1);
    for session_n = 1:size(play_session_ranges, 1)
        window_sample_ranges = session_windows(session_n).sample_ranges;
        n_windows = session_windows(session_n).n_windows;
        welch_window_power = nan(n_windows, n_frequency);
        fourier_window_power = nan(n_windows, n_frequency);

        for window_n = 1:n_windows
            sample_range = window_sample_ranges(window_n, 1): ...
                window_sample_ranges(window_n, 2);
            this_lfp = detrend( ...
                double(LFP(lfp_row, sample_range)), 'constant');

            % One 10-second Hann segment gives a true 0.1 Hz bin spacing.
            this_welch = pwelch(this_lfp, welch_window, 0, ...
                window_samples, sampling_rate, 'onesided');
            welch_window_power(window_n, :) = ...
                this_welch(2:n_frequency + 1).';

            % Direct rectangular-window Fourier periodogram.
            fourier_coefficients = fft(this_lfp, window_samples);
            this_fourier = ...
                abs(fourier_coefficients(1:window_samples / 2 + 1)).^2 ...
                / (sampling_rate * window_samples);
            this_fourier(2:end - 1) = 2 * this_fourier(2:end - 1);
            fourier_window_power(window_n, :) = ...
                this_fourier(2:n_frequency + 1);
        end

        mean_welch_power = mean(welch_window_power, 1, 'omitnan').';
        mean_fourier_power = ...
            mean(fourier_window_power, 1, 'omitnan').';
        mean_welch_power = max(mean_welch_power, realmin('double'));
        mean_fourier_power = max(mean_fourier_power, realmin('double'));

        welch_fooof = fooof_matlab(frequency, mean_welch_power, ...
            'AperiodicMode', fooof_aperiodic_mode, ...
            'PeakWidthLimits', fooof_peak_widths);
        fourier_fooof = fooof_matlab(frequency, mean_fourier_power, ...
            'AperiodicMode', fooof_aperiodic_mode, ...
            'PeakWidthLimits', fooof_peak_widths);

        session_results(session_n).session_index = session_n;
        session_results(session_n).behavior_row = ...
            play_session_table(session_n, :);
        session_results(session_n).time_range = ...
            play_session_ranges(session_n, :);
        session_results(session_n).window_sample_ranges = ...
            window_sample_ranges;
        session_results(session_n).window_time_ranges = ...
            session_windows(session_n).time_ranges;
        session_results(session_n).n_windows = n_windows;

        session_results(session_n).welch.window_power = ...
            welch_window_power;
        session_results(session_n).welch.mean_power = mean_welch_power;
        session_results(session_n).welch.log10_mean_power = ...
            log10(mean_welch_power);
        session_results(session_n).welch.fooof = welch_fooof;
        session_results(session_n).welch.gaussian_periodic_fit = ...
            welch_fooof.peak_fit;
        session_results(session_n).welch.raw_aperiodic_subtracted = ...
            log10(mean_welch_power) - welch_fooof.aperiodic_fit;

        session_results(session_n).fourier.window_power = ...
            fourier_window_power;
        session_results(session_n).fourier.mean_power = ...
            mean_fourier_power;
        session_results(session_n).fourier.log10_mean_power = ...
            log10(mean_fourier_power);
        session_results(session_n).fourier.fooof = fourier_fooof;
        session_results(session_n).fourier.gaussian_periodic_fit = ...
            fourier_fooof.peak_fit;
        session_results(session_n).fourier.raw_aperiodic_subtracted = ...
            log10(mean_fourier_power) - fourier_fooof.aperiodic_fit;
    end
    spectra_struct(channel_n).sessions = session_results;
end
end
