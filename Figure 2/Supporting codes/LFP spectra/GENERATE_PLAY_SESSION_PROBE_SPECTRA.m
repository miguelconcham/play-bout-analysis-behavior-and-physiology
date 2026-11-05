function spectra = GENERATE_PLAY_SESSION_PROBE_SPECTRA(current_dir)
% GENERATE_PLAY_SESSION_PROBE_SPECTRA
% Periodic and aperiodic spectra for every mapped channel.
%
% Play sessions are Session_structure rows of type Partners session.
% Each session is split into non-overlapping 10 s windows. A Welch PSD and
% a direct Fourier periodogram are averaged in linear power, then each
% channel is parameterized with FOOOF (fixed aperiodic mode).
%
% Band means are not computed here. The frequency axis runs through 45 Hz
% so a later plot can show delta (1-5 Hz), theta (6-12 Hz), and gamma
% (30-45 Hz) from the saved spectra.

%% Analysis parameters
sampling_rate = 2500;
window_seconds = 10;
window_samples = window_seconds * sampling_rate;
frequency = (0.1:0.1:45).';
fooof_peak_widths = [0.5 12];
fooof_aperiodic_mode = 'fixed';
area2analyze = 'PAG';
n_frequency = numel(frequency);

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
animal_batch = animal_code_params{1};
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
Behavior.End = predict(synch_model_video2NPX, Behavior.End);

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

%% Channel map, area labels, and the mid-LPAG site used by the original fit
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
if isempty(area_limit)
    error('No area limits found for %s.', this_animal);
end

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
channel_table = channel_table_from_map( ...
    chanMap, xcoords, ycoords, NPX_Type, hard_coded_x_coords, ...
    area_limit, size(LFP, 1));
center_rows = mid_lpag_lfp_rows( ...
    ycoords, xcoords, chanMap, area_limit, NPX_Type, hard_coded_x_coords, ...
    size(LFP, 1));
channel_table.IsLPAGCenter = ismember(channel_table.LFPRow, center_rows);

%% Create full 10-second windows within each play-session row
n_lfp_samples = size(LFP, 2);
n_sessions = size(play_session_ranges, 1);
session_windows = repmat(struct( ...
    'sample_ranges', zeros(0, 2), ...
    'n_windows', 0), n_sessions, 1);

for session_n = 1:n_sessions
    first_sample = max(1, ceil(play_session_ranges(session_n, 1) * sampling_rate));
    last_sample = min(n_lfp_samples, ...
        floor(play_session_ranges(session_n, 2) * sampling_rate));
    starts = first_sample:window_samples:(last_sample - window_samples + 1);
    if isempty(starts)
        error('Play-session row %d contains no complete 10-second window.', ...
            session_n);
    end
    session_windows(session_n).sample_ranges = [starts(:), starts(:) + window_samples - 1];
    session_windows(session_n).n_windows = numel(starts);
end

bin_frequency = (1:n_frequency).' * sampling_rate / window_samples;
if max(abs(bin_frequency - frequency)) > 1e-9
    error('Frequency grid does not match the 10 s Fourier bins.');
end

%% Recording metadata
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
parameters.band_averages = 'not computed; derived in Analyse_play_session_probe_spectra';

spectra = struct();
spectra.animal_name = animal_code;
spectra.animal_area_table_name = this_animal;
spectra.neuropixels_directory = current_dir;
spectra.npx_type = NPX_Type;
spectra.frequency = frequency;
spectra.channel_table = channel_table;
spectra.parameters = parameters;
spectra.play_session_table = play_session_table;
spectra.play_session_ranges = play_session_ranges;
spectra.n_sessions = n_sessions;
spectra.sessions = repmat(struct( ...
    'session_index', [], ...
    'time_range', [], ...
    'n_windows', [], ...
    'welch', [], ...
    'fourier', []), n_sessions, 1);

welch_window = hann(window_samples, 'periodic');
rectangular_window = ones(window_samples, 1);
lfp_rows = channel_table.LFPRow;
n_channels = numel(lfp_rows);

%% Spectra and FOOOF fits for every channel
for session_n = 1:n_sessions
    window_sample_ranges = session_windows(session_n).sample_ranges;
    n_windows = session_windows(session_n).n_windows;
    fprintf(['%s session %d/%d: %d channels, %d windows of %g s\n'], ...
        animal_code, session_n, n_sessions, n_channels, n_windows, ...
        window_seconds);

    welch_sum = zeros(n_channels, n_frequency);
    fourier_sum = zeros(n_channels, n_frequency);
    for window_n = 1:n_windows
        sample_index = window_sample_ranges(window_n, 1): ...
            window_sample_ranges(window_n, 2);
        segment = double(LFP(lfp_rows, sample_index));
        segment = segment - mean(segment, 2);
        welch_sum = welch_sum + onesided_psd( ...
            segment, welch_window, sampling_rate, n_frequency);
        fourier_sum = fourier_sum + onesided_psd( ...
            segment, rectangular_window, sampling_rate, n_frequency);
    end

    spectra.sessions(session_n).session_index = session_n;
    spectra.sessions(session_n).time_range = play_session_ranges(session_n, :);
    spectra.sessions(session_n).n_windows = n_windows;
    spectra.sessions(session_n).welch = fit_power_matrix( ...
        frequency, welch_sum / n_windows, ...
        fooof_aperiodic_mode, fooof_peak_widths, 'Welch');
    spectra.sessions(session_n).fourier = fit_power_matrix( ...
        frequency, fourier_sum / n_windows, ...
        fooof_aperiodic_mode, fooof_peak_widths, 'Fourier');
end
end

function channel_table = channel_table_from_map( ...
        chan_map, xcoords, ycoords, npx_type, shank_x, area_limit, n_lfp_rows)
n_sites = numel(chan_map);
channel = nan(n_sites, 1);
lfp_row = nan(n_sites, 1);
map_index = nan(n_sites, 1);
x_um = nan(n_sites, 1);
y_um = nan(n_sites, 1);
probe_number = nan(n_sites, 1);
area_name = cell(n_sites, 1);
kept = false(n_sites, 1);

for site_n = 1:n_sites
    this_channel = chan_map(site_n);
    this_row = this_channel + 1;
    if ~isfinite(this_channel) || this_row < 1 || this_row > n_lfp_rows
        continue
    end
    kept(site_n) = true;
    channel(site_n) = this_channel;
    lfp_row(site_n) = this_row;
    map_index(site_n) = site_n;
    x_um(site_n) = xcoords(site_n);
    y_um(site_n) = ycoords(site_n);
    probe_number(site_n) = probe_number_from_x( ...
        xcoords(site_n), npx_type, shank_x);
    area_name{site_n} = area_name_at_site( ...
        ycoords(site_n), probe_number(site_n), npx_type, area_limit);
end

if ~any(kept)
    error('Channel map contains no site that indexes the LFP matrix.');
end
if numel(unique(lfp_row(kept))) ~= sum(kept)
    error('Channel map assigns more than one site to the same LFP row.');
end

channel_table = table(channel(kept), lfp_row(kept), map_index(kept), ...
    x_um(kept), y_um(kept), probe_number(kept), area_name(kept), ...
    'VariableNames', {'Channel', 'LFPRow', 'MapIndex', 'Xum', 'Yum', ...
    'ProbeNumber', 'Area'});
end

function probe_number = probe_number_from_x(x_coordinate, npx_type, shank_x)
if npx_type == 1
    probe_number = 1;
    return
end
probe_number = NaN;
if ~isfinite(x_coordinate)
    return
end
for probe_n = 1:size(shank_x, 1)
    if ismember(x_coordinate, shank_x(probe_n, :))
        probe_number = probe_n;
        return
    end
end
end

function label = area_name_at_site(y_coordinate, probe_number, npx_type, area_limit)
label = '';
if ~isfinite(y_coordinate)
    return
end
% Same depth test as the channel-wise area maps: the end limit is inclusive
% for integer micrometre coordinates.
in_depth = y_coordinate >= area_limit.depth_start & ...
    y_coordinate < area_limit.depth_end + 1;
if npx_type == 1
    if ~ismember('Probe_Area', area_limit.Properties.VariableNames)
        error('Area-limit table has no Probe_Area column for NPX1 labels.');
    end
    in_area = in_depth & ismember(area_limit.Probe_Area, 'PAG');
elseif isfinite(probe_number)
    in_area = in_depth & area_limit.ProbeNum == probe_number;
else
    return
end
if ~any(in_area)
    return
end
names = cellstr(string(area_limit.area(in_area)));
names = unique(names(~strcmp(names, '')), 'stable');
if isempty(names)
    return
end
label = strjoin(names, ', ');
end

function center_rows = mid_lpag_lfp_rows( ...
        ycoords, xcoords, chan_map, area_limit, npx_type, shank_x, n_lfp_rows)
lpag_rows = ismember(area_limit.area, 'LPAG');
depth_range = area_limit{lpag_rows, {'ProbeNum', 'depth_start', 'depth_end'}};
center_rows = zeros(0, 1);
if isempty(depth_range)
    disp('No LPAG limits; no central LPAG channel is flagged.')
    return
end

for probe_n = 1:size(depth_range, 1)
    if npx_type == 1
        in_lpag = ycoords >= depth_range(probe_n, 2) & ...
            ycoords <= depth_range(probe_n, 3);
    else
        probe_number = depth_range(probe_n, 1);
        if ~isfinite(probe_number) || probe_number < 1 || ...
                probe_number > size(shank_x, 1)
            error('Invalid probe number %g in LPAG area limits.', probe_number);
        end
        in_lpag = ycoords >= depth_range(probe_n, 2) & ...
            ycoords <= depth_range(probe_n, 3) & ...
            ismember(xcoords, shank_x(probe_number, :));
    end
    if ~any(in_lpag)
        fprintf('No mapped sites inside LPAG probe row %d; center not flagged.\n', ...
            probe_n);
        continue
    end
    mean_location = [mean(xcoords(in_lpag)), mean(ycoords(in_lpag))];
    distance = sum(abs([xcoords, ycoords] - mean_location), 2);
    distance(~in_lpag) = Inf;
    [~, closest_map_index] = min(distance);
    selected_row = chan_map(closest_map_index) + 1;
    if isfinite(selected_row) && selected_row >= 1 && selected_row <= n_lfp_rows
        center_rows(end + 1, 1) = selected_row; %#ok<AGROW>
    end
end
end

function power = onesided_psd(segment, window, sampling_rate, n_frequency)
% One-sided PSD matching pwelch / the rectangular periodogram.
% segment is channels x samples and is already mean-centered.
window = window(:).';
tapered = segment .* window;
n_samples = size(segment, 2);
coefficients = fft(tapered, n_samples, 2);
n_onesided = floor(n_samples / 2) + 1;
power = abs(coefficients(:, 1:n_onesided)).^2;
power = power ./ (sampling_rate * sum(window.^2));
if mod(n_samples, 2) == 0
    power(:, 2:end - 1) = 2 * power(:, 2:end - 1);
else
    power(:, 2:end) = 2 * power(:, 2:end);
end
power = power(:, 2:(n_frequency + 1));
end

function fitted = fit_power_matrix( ...
        frequency, linear_power, aperiodic_mode, peak_widths, method_name)
n_channels = size(linear_power, 1);
n_frequency = numel(frequency);
linear_power = max(linear_power, realmin('double'));

fitted = struct();
fitted.log10_mean_power = log10(linear_power);
fitted.aperiodic_fit = nan(n_channels, n_frequency);
fitted.gaussian_periodic_fit = nan(n_channels, n_frequency);
fitted.raw_aperiodic_subtracted = nan(n_channels, n_frequency);
fitted.aperiodic_params = nan(n_channels, 2);
fitted.r_squared = nan(n_channels, 1);
fitted.error = nan(n_channels, 1);
fitted.fit_ok = false(n_channels, 1);
fitted.peak_params = cell(n_channels, 1);
fitted.gaussian_params = cell(n_channels, 1);

n_failed = 0;
first_error = '';
for channel_n = 1:n_channels
    if mod(channel_n, 50) == 1
        fprintf('  %s FOOOF %d/%d\n', method_name, channel_n, n_channels);
    end
    this_power = linear_power(channel_n, :).';
    if any(~isfinite(this_power))
        fitted.peak_params{channel_n} = zeros(0, 3);
        fitted.gaussian_params{channel_n} = zeros(0, 3);
        n_failed = n_failed + 1;
        continue
    end
    try
        out = fooof_matlab(frequency, this_power, ...
            'AperiodicMode', aperiodic_mode, ...
            'PeakWidthLimits', peak_widths);
    catch fit_error
        fitted.peak_params{channel_n} = zeros(0, 3);
        fitted.gaussian_params{channel_n} = zeros(0, 3);
        n_failed = n_failed + 1;
        if isempty(first_error)
            first_error = fit_error.message;
        end
        continue
    end

    aperiodic_params = out.aperiodic_params(:).';
    if numel(aperiodic_params) ~= 2
        error('Fixed FOOOF mode returned %d aperiodic parameters.', ...
            numel(aperiodic_params));
    end
    fitted.aperiodic_fit(channel_n, :) = out.aperiodic_fit(:).';
    fitted.gaussian_periodic_fit(channel_n, :) = out.peak_fit(:).';
    fitted.raw_aperiodic_subtracted(channel_n, :) = ...
        fitted.log10_mean_power(channel_n, :) - fitted.aperiodic_fit(channel_n, :);
    fitted.aperiodic_params(channel_n, :) = aperiodic_params;
    fitted.r_squared(channel_n) = out.r_squared;
    fitted.error(channel_n) = out.error;
    fitted.fit_ok(channel_n) = true;
    fitted.peak_params{channel_n} = out.peak_params;
    fitted.gaussian_params{channel_n} = out.gaussian_params;
end

fprintf('  %s FOOOF finished: %d/%d channels fitted', ...
    method_name, n_channels - n_failed, n_channels);
if n_failed > 0 && ~isempty(first_error)
    fprintf('. First failure: %s', first_error);
end
fprintf('\n');
end
