function [recordings, parameters] = GENERATE_SPEED_PEAK_FREQUENCY(current_dir, options)
%GENERATE_SPEED_PEAK_FREQUENCY  Delta and theta peak frequency versus speed.
%
% One Neuropixels folder. For each partner tracking file, every overlapping
% Session_structure / Partners session row is analyzed on its own. Speed is
% the focal animal's frame-to-frame displacement in pixels per second.
% Each non-overlapping 10 s window is labeled by its median speed. Those
% window speeds are split into equal-count bins of 5-10% of the windows.
% Spectra in a bin are averaged in linear power and parameterized with
% FOOOF exactly as the partner-session spectra are:
% Welch and direct Fourier spectra, aperiodic fit, Gaussian periodic fit,
% and the raw aperiodic-subtracted spectrum.
%
% The returned recordings keep every window spectrum, every bin spectrum,
% and every FOOOF decomposition. parameters is the struct stored with them.

if nargin < 2 || isempty(options)
    options = struct();
end
parameters = speed_peak_options(options);
options = parameters;

[animal_code, table_name] = recording_names(current_dir);
data_root = resolve_data_root(current_dir, options);
[xcoords, ycoords, chanMap] = load_channel_map(current_dir);
npx_type = neuropixels_type(xcoords);
selected = select_lpag(data_root, table_name, xcoords, ycoords, chanMap, ...
    npx_type, options, animal_code);
if isempty(selected)
    recordings = empty_recording();
    recordings(1) = [];
    return
end

[channel_lfp, n_samples] = read_selected_lfp(current_dir, [selected.lfp_row], ...
    options.n_lfp_channels);
[synch_model, partner_ranges] = load_partner_sessions( ...
    current_dir, data_root, animal_code);
tracking = tracking_files(current_dir, data_root, animal_code);
if isempty(tracking)
    warning('GENERATE_SPEED_PEAK_FREQUENCY:NoTracking', ...
        '%s: no partner tracking files. Skipping.', animal_code);
    recordings = empty_recording();
    recordings(1) = [];
    return
end

n_channels = numel(selected);
recordings = repmat(empty_recording(), n_channels, 1);
session_lists = cell(n_channels, 1);
notes = {};
frequency = options.frequency;
freq_index = (1:options.n_frequency).' + 1;
welch_window = hann(options.window_samples, 'periodic');

for file_n = 1:numel(tracking)
    tracking_path = fullfile(tracking(file_n).folder, tracking(file_n).name);
    try
        [speed, speed_time, track_range] = load_focal_speed( ...
            tracking_path, synch_model, options);
    catch exception
        notes{end + 1} = sprintf('P%d tracking was not used: %s', ...
            tracking(file_n).partner, exception.message); %#ok<AGROW>
        warning('GENERATE_SPEED_PEAK_FREQUENCY:Tracking', ...
            '%s P%d: %s', animal_code, tracking(file_n).partner, exception.message);
        continue
    end

    intervals = session_intervals(partner_ranges, track_range, ...
        options.window_seconds);
    if isempty(intervals)
        notes{end + 1} = sprintf(['P%d has no Partners session overlap ' ...
            'of at least one %g s window. The tracked interval is used.'], ...
            tracking(file_n).partner, options.window_seconds); %#ok<AGROW>
        intervals = struct( ...
            'source', 'tracking', ...
            'behavior_start', NaN, ...
            'behavior_end', NaN, ...
            'start', track_range(1), ...
            'end', track_range(2), ...
            'session_index', 1);
    end

    for interval_n = 1:numel(intervals)
        interval = intervals(interval_n);
        in_session = speed_time >= interval.start & speed_time <= interval.end;
        session_speed = speed(in_session);
        session_speed_time = speed_time(in_session);
        windows = analysis_windows( ...
            [interval.start, interval.end], n_samples, options);
        if windows.n_windows == 0
            notes{end + 1} = sprintf(['P%d session %d has no complete %g s window.'], ...
                tracking(file_n).partner, interval.session_index, ...
                options.window_seconds); %#ok<AGROW>
            continue
        end

        [window_speed, n_speed_samples, frame_bin_time] = window_median_speeds( ...
            session_speed, session_speed_time, windows.time_ranges);
        window_speed(n_speed_samples < options.min_speed_samples_per_window) = NaN;
        bin_definition = speed_quantile_bins(window_speed, ...
            options.speed_bin_fraction_limits);
        if bin_definition.n_bins == 0
            notes{end + 1} = sprintf(['P%d session %d has no window with enough ' ...
                'speed samples to bin.'], tracking(file_n).partner, ...
                interval.session_index); %#ok<AGROW>
            continue
        end

        [welch_power, fourier_power] = channel_window_spectra( ...
            channel_lfp, windows.sample_ranges, welch_window, ...
            options, freq_index);

        for channel_n = 1:n_channels
            session = make_session( ...
                tracking(file_n), tracking_path, interval, track_range, ...
                session_speed, session_speed_time, bin_definition, windows, ...
                window_speed, n_speed_samples, frame_bin_time, ...
                squeeze(welch_power(channel_n, :, :)), ...
                squeeze(fourier_power(channel_n, :, :)), ...
                frequency, options);
            if isempty(session_lists{channel_n})
                session_lists{channel_n} = session;
            else
                session_lists{channel_n}(end + 1, 1) = session; %#ok<AGROW>
            end
        end
    end
end

for channel_n = 1:n_channels
    recording = empty_recording();
    recording.animal_name = animal_code;
    recording.animal_area_table_name = table_name;
    recording.neuropixels_directory = current_dir;
    recording.npx_type = npx_type;
    recording.lpag_channel = selected(channel_n).channel;
    recording.lpag_lfp_row = selected(channel_n).lfp_row;
    recording.lpag_map_index = selected(channel_n).map_index;
    recording.lpag_probe_number = selected(channel_n).probe_number;
    recording.lpag_x_um = selected(channel_n).x_um;
    recording.lpag_y_um = selected(channel_n).y_um;
    recording.lpag_depth_range = selected(channel_n).depth_range;
    recording.parameters = options;
    recording.frequency = frequency;
    if isempty(notes)
        recording.notes = '';
    else
        recording.notes = strjoin(notes(:), newline);
    end
    if ~isempty(session_lists{channel_n})
        recording.sessions = session_lists{channel_n};
    end
    recordings(channel_n) = recording;
end
end

function session = make_session(tracking, tracking_path, interval, track_range, ...
        session_speed, session_speed_time, bin_definition, windows, ...
        window_speed, n_speed_samples, frame_bin_time, ...
        welch_power, fourier_power, frequency, options)

welch_power = reshape(welch_power, windows.n_windows, []);
fourier_power = reshape(fourier_power, windows.n_windows, []);
speed_bins = describe_speed_bins(bin_definition, bin_definition.bin_id, window_speed);
welch = parameterize_binned_spectra(welch_power, bin_definition.bin_id, ...
    bin_definition.n_bins, frequency, options);
fourier = parameterize_binned_spectra(fourier_power, bin_definition.bin_id, ...
    bin_definition.n_bins, frequency, options);
spectrum_speed = [speed_bins.spectrum_speed].';
welch_delta = [welch.delta_peak.center_frequency].';
welch_theta = [welch.theta_peak.center_frequency].';
fourier_delta = [fourier.delta_peak.center_frequency].';
fourier_theta = [fourier.theta_peak.center_frequency].';

session = empty_session();
session.partner_index = tracking.partner;
session.tracking_file = tracking_path;
session.session_index = interval.session_index;
session.session_source = interval.source;
session.behavior_time_range = [interval.behavior_start, interval.behavior_end];
session.tracking_time_range = track_range(:).';
session.analyzed_time_range = [interval.start, interval.end];
session.speed_pixels_per_second = session_speed;
session.speed_time = session_speed_time;
session.speed_bin_id = frame_speed_bins( ...
    frame_bin_time, bin_definition.bin_id);
session.speed_bin_edges = bin_definition.edges;
session.speed_bin_rule = bin_definition.rule;
session.bins_within_upper_fraction = bin_definition.within_upper_fraction;
session.bins_within_lower_fraction = bin_definition.within_lower_fraction;
session.bins_within_requested_fractions = bin_definition.within_requested_fractions;
session.speed_bins = speed_bins;
session.window_sample_ranges = windows.sample_ranges;
session.window_time_ranges = windows.time_ranges;
session.window_speed = window_speed;
session.window_bin_id = bin_definition.bin_id;
session.window_bin_occupancy = n_speed_samples / ...
    (options.window_seconds * options.video_fps);
session.window_n_speed_samples = n_speed_samples;
session.n_windows = windows.n_windows;
session.welch = welch;
session.fourier = fourier;
session.correlation_welch = speed_peak_correlation( ...
    spectrum_speed, welch_delta, welch_theta, options.min_bins_for_correlation);
session.correlation_fourier = speed_peak_correlation( ...
    spectrum_speed, fourier_delta, fourier_theta, options.min_bins_for_correlation);
end

function speed_bins = describe_speed_bins(bin_definition, window_bin, window_speed)
n_bins = bin_definition.n_bins;
speed_bins = repmat(struct( ...
    'n_speed_samples', 0, ...
    'fraction_of_speed_samples', NaN, ...
    'speed_min', NaN, ...
    'speed_max', NaN, ...
    'speed_median', NaN, ...
    'speed_mean', NaN, ...
    'n_windows', 0, ...
    'window_index', zeros(0, 1), ...
    'spectrum_speed', NaN), n_bins, 1);
for bin_n = 1:n_bins
    rows = find(window_bin == bin_n);
    speed_bins(bin_n).n_speed_samples = bin_definition.count(bin_n);
    speed_bins(bin_n).fraction_of_speed_samples = bin_definition.fraction(bin_n);
    speed_bins(bin_n).speed_min = bin_definition.speed_min(bin_n);
    speed_bins(bin_n).speed_max = bin_definition.speed_max(bin_n);
    speed_bins(bin_n).speed_median = bin_definition.speed_median(bin_n);
    speed_bins(bin_n).speed_mean = bin_definition.speed_mean(bin_n);
    speed_bins(bin_n).n_windows = numel(rows);
    speed_bins(bin_n).window_index = rows(:);
    if ~isempty(rows)
        speed_bins(bin_n).spectrum_speed = median(window_speed(rows), 'omitnan');
    end
end
end

function [window_speed, n_speed_samples, frame_in_window] = window_median_speeds( ...
        speed, speed_time, window_time_ranges)

n_windows = size(window_time_ranges, 1);
window_speed = nan(n_windows, 1);
n_speed_samples = zeros(n_windows, 1);
frame_in_window = false(numel(speed_time), n_windows);
for window_n = 1:n_windows
    if window_n == n_windows
        in_window = speed_time >= window_time_ranges(window_n, 1) & ...
            speed_time <= window_time_ranges(window_n, 2);
    else
        in_window = speed_time >= window_time_ranges(window_n, 1) & ...
            speed_time < window_time_ranges(window_n, 2);
    end
    finite = in_window & isfinite(speed);
    n_speed_samples(window_n) = nnz(finite);
    frame_in_window(:, window_n) = in_window;
    if any(finite)
        window_speed(window_n) = median(speed(finite));
    end
end
end

function frame_bin_id = frame_speed_bins(frame_in_window, window_bin)
frame_bin_id = nan(size(frame_in_window, 1), 1);
for window_n = 1:numel(window_bin)
    if ~isfinite(window_bin(window_n))
        continue
    end
    frame_bin_id(frame_in_window(:, window_n)) = window_bin(window_n);
end
end

function [welch_power, fourier_power] = channel_window_spectra( ...
        channel_lfp, sample_ranges, welch_window, options, freq_index)

n_channels = size(channel_lfp, 1);
n_windows = size(sample_ranges, 1);
welch_power = nan(n_channels, n_windows, options.n_frequency);
fourier_power = nan(n_channels, n_windows, options.n_frequency);
for window_n = 1:n_windows
    sample_index = sample_ranges(window_n, 1):sample_ranges(window_n, 2);
    for channel_n = 1:n_channels
        [welch_row, fourier_row] = one_window_spectra( ...
            channel_lfp(channel_n, sample_index), welch_window, ...
            options.window_samples, options.sampling_rate, freq_index);
        welch_power(channel_n, window_n, :) = welch_row;
        fourier_power(channel_n, window_n, :) = fourier_row;
    end
end
end

function [welch_power, fourier_power] = one_window_spectra( ...
        segment, welch_window, window_samples, sampling_rate, freq_index)

segment = detrend(double(segment(:)), 'constant');
welch_full = pwelch(segment, welch_window, 0, window_samples, ...
    sampling_rate, 'onesided');
fourier_coefficients = fft(segment, window_samples);
fourier_full = abs(fourier_coefficients(1:window_samples / 2 + 1)).^2 ...
    / (sampling_rate * window_samples);
fourier_full(2:end - 1) = 2 * fourier_full(2:end - 1);
welch_power = welch_full(freq_index).';
fourier_power = fourier_full(freq_index).';
end

function method = parameterize_binned_spectra( ...
        window_power, window_bin, n_bins, frequency, options)

n_frequency = numel(frequency);
method = struct();
method.window_power = window_power;
method.mean_power = nan(n_frequency, n_bins);
method.log10_mean_power = nan(n_frequency, n_bins);
method.fooof = cell(n_bins, 1);
method.gaussian_periodic_fit = nan(n_frequency, n_bins);
method.raw_aperiodic_subtracted = nan(n_frequency, n_bins);
method.delta_peak = repmat(select_band_peak(zeros(0, 3), options.delta_band), n_bins, 1);
method.theta_peak = repmat(select_band_peak(zeros(0, 3), options.theta_band), n_bins, 1);
method.r_squared = nan(n_bins, 1);
method.fit_error = nan(n_bins, 1);
method.fit_ok = false(n_bins, 1);
method.fit_message = repmat({''}, n_bins, 1);

for bin_n = 1:n_bins
    rows = find(window_bin == bin_n);
    if isempty(rows)
        method.fit_message{bin_n} = 'no window was assigned to this speed bin';
        continue
    end
    mean_power = mean(window_power(rows, :), 1, 'omitnan').';
    if isempty(mean_power) || any(~isfinite(mean_power)) || any(mean_power < 0)
        method.fit_message{bin_n} = 'mean spectrum was not finite and non-negative';
        continue
    end
    mean_power = max(mean_power, realmin('double'));
    method.mean_power(:, bin_n) = mean_power;
    method.log10_mean_power(:, bin_n) = log10(mean_power);
    try
        fit = fooof_matlab(frequency, mean_power, ...
            'AperiodicMode', options.fooof_aperiodic_mode, ...
            'PeakWidthLimits', options.fooof_peak_width_limits, ...
            'MaxNPeaks', options.fooof_max_n_peaks, ...
            'MinPeakHeight', options.fooof_min_peak_height, ...
            'PeakThreshold', options.fooof_peak_threshold);
    catch exception
        method.fit_message{bin_n} = exception.message;
        continue
    end
    method.fooof{bin_n} = fit;
    method.gaussian_periodic_fit(:, bin_n) = fit.peak_fit(:);
    method.raw_aperiodic_subtracted(:, bin_n) = ...
        log10(mean_power) - fit.aperiodic_fit(:);
    method.delta_peak(bin_n) = select_band_peak(fit.peak_params, options.delta_band);
    method.theta_peak(bin_n) = select_band_peak(fit.peak_params, options.theta_band);
    method.r_squared(bin_n) = fit.r_squared;
    method.fit_error(bin_n) = fit.error;
    method.fit_ok(bin_n) = true;
end
end

function selected = select_lpag(data_root, table_name, xcoords, ycoords, ...
        chanMap, npx_type, options, animal_code)

area_limit_file = fullfile(data_root, 'Area_limits_GoodLooking.xlsx');
if ~isfile(area_limit_file)
    error('GENERATE_SPEED_PEAK_FREQUENCY:areas', ...
        'Area-limit table not found: %s', area_limit_file);
end
area_limit = readtable(area_limit_file);
area_limit = area_limit(ismember(area_limit.AnimalName, table_name), :);
lpag_limit = area_limit(ismember(area_limit.area, options.area_name), :);
if isempty(lpag_limit)
    warning('GENERATE_SPEED_PEAK_FREQUENCY:NoLPAG', ...
        '%s: no %s rows for %s. Skipping.', animal_code, options.area_name, table_name);
    selected = [];
    return
end
limit_matrix = lpag_limit{:, {'ProbeNum', 'depth_start', 'depth_end'}};
selected = select_mid_lpag_channels( ...
    xcoords, ycoords, chanMap, limit_matrix, npx_type);
end

function intervals = session_intervals(partner_ranges, track_range, window_seconds)
intervals = struct('source', {}, 'behavior_start', {}, 'behavior_end', {}, ...
    'start', {}, 'end', {}, 'session_index', {});
if isempty(partner_ranges)
    return
end
session_index = 0;
for row_n = 1:size(partner_ranges, 1)
    start_s = max(track_range(1), partner_ranges(row_n, 1));
    end_s = min(track_range(2), partner_ranges(row_n, 2));
    if ~isfinite(start_s) || ~isfinite(end_s) || end_s - start_s < window_seconds
        continue
    end
    session_index = session_index + 1;
    intervals(end + 1).source = 'partners session'; %#ok<AGROW>
    intervals(end).behavior_start = partner_ranges(row_n, 1);
    intervals(end).behavior_end = partner_ranges(row_n, 2);
    intervals(end).start = start_s;
    intervals(end).end = end_s;
    intervals(end).session_index = session_index;
end
end

function windows = analysis_windows(time_range, n_lfp_samples, options)
first_sample = max(1, ceil(time_range(1) * options.sampling_rate));
last_sample = min(n_lfp_samples, floor(time_range(2) * options.sampling_rate));
starts = first_sample:options.window_samples:(last_sample - options.window_samples + 1);
starts = starts(:);
ends = starts + options.window_samples - 1;
windows = struct();
windows.sample_ranges = [starts, ends];
windows.time_ranges = [starts, ends] / options.sampling_rate;
windows.n_windows = numel(starts);
end

function [speed, speed_time, track_range] = load_focal_speed( ...
        tracking_file, synch_model, options)

loaded = load(tracking_file, 'traking_structure');
if ~isfield(loaded, 'traking_structure')
    error('GENERATE_SPEED_PEAK_FREQUENCY:tracking', ...
        'traking_structure is missing from %s.', tracking_file);
end
track = loaded.traking_structure;
if ~isfield(track, 'frames2stract') || ~isfield(track, 'animal_pos')
    error('GENERATE_SPEED_PEAK_FREQUENCY:tracking', ...
        'Tracking file %s has no frames2stract and animal_pos.', tracking_file);
end
frames = round(double(track.frames2stract(:)));
position = double(track.animal_pos);
if size(position, 2) ~= 2 || size(position, 1) ~= numel(frames)
    error('GENERATE_SPEED_PEAK_FREQUENCY:tracking', ...
        'animal_pos must be one [x y] row per frame in %s.', tracking_file);
end
keep = isfinite(frames) & all(isfinite(position), 2);
frames = frames(keep);
position = position(keep, :);
if numel(frames) < 4
    error('GENERATE_SPEED_PEAK_FREQUENCY:tracking', ...
        'Fewer than 4 finite tracked frames in %s.', tracking_file);
end
[frames, order] = sort(frames);
position = position(order, :);
[frames, unique_index] = unique(frames, 'last');
position = position(unique_index, :);
if numel(frames) >= options.tracking_smooth_samples
    position(:, 1) = smoothdata(position(:, 1), options.tracking_smooth, ...
        options.tracking_smooth_samples);
    position(:, 2) = smoothdata(position(:, 2), options.tracking_smooth, ...
        options.tracking_smooth_samples);
end

video_time = frames / options.video_fps;
npx_time = double(predict(synch_model, video_time(:)));
npx_time = npx_time(:);
if any(diff(npx_time) <= 0)
    error('GENERATE_SPEED_PEAK_FREQUENCY:tracking', ...
        'Synchronized tracking times are not increasing in %s.', tracking_file);
end
frame_index = (frames(1):frames(end)).';
full_video_time = frame_index / options.video_fps;
full_npx_time = double(predict(synch_model, full_video_time(:)));
full_npx_time = full_npx_time(:);
if any(diff(full_npx_time) <= 0)
    error('GENERATE_SPEED_PEAK_FREQUENCY:tracking', ...
        'The frame grid is not increasing after synchronization.');
end
full_position = [ ...
    interp1(npx_time, position(:, 1), full_npx_time, options.position_interpolation), ...
    interp1(npx_time, position(:, 2), full_npx_time, options.position_interpolation)];
if any(~isfinite(full_position), 'all')
    error('GENERATE_SPEED_PEAK_FREQUENCY:tracking', ...
        'Interpolated positions are not finite in %s.', tracking_file);
end
[speed, speed_time] = focal_animal_speed(full_position, full_npx_time);
track_range = [full_npx_time(1), full_npx_time(end)];
end

function files = tracking_files(current_dir, data_root, animal_code)
folders = {current_dir, fullfile(data_root, 'Traking backups')};
files = struct('name', {}, 'folder', {}, 'partner', {});
for folder_n = 1:numel(folders)
    listing = dir(fullfile(folders{folder_n}, ...
        [animal_code, ' P* traking_structure.mat']));
    for file_n = 1:numel(listing)
        token = regexp(listing(file_n).name, ...
            ' P(\d+) traking_structure\.mat$', 'tokens', 'once');
        if isempty(token)
            continue
        end
        partner = str2double(token{1});
        if any([files.partner] == partner)
            continue
        end
        files(end + 1).name = listing(file_n).name; %#ok<AGROW>
        files(end).folder = listing(file_n).folder;
        files(end).partner = partner;
    end
end
if ~isempty(files)
    [~, order] = sort([files.partner]);
    files = files(order);
end
end

function [synch_model, ranges] = load_partner_sessions(current_dir, data_root, animal_code)
local_synch = fullfile(current_dir, 'synch_model_video2NPX.mat');
data_synch = fullfile(data_root, 'Synch data', animal_code, 'synch_model_video2NPX.mat');
if isfile(local_synch)
    synch_file = local_synch;
elseif isfile(data_synch)
    synch_file = data_synch;
else
    error('GENERATE_SPEED_PEAK_FREQUENCY:synch', ...
        'Synchronization model not found for %s.', animal_code);
end
loaded_synch = load(synch_file, 'synch_model_video2NPX');
synch_model = loaded_synch.synch_model_video2NPX;

local_behavior = dir(fullfile(current_dir, 'ELAN*'));
data_behavior = fullfile(data_root, 'Behavior backups', [animal_code, '.txt']);
if ~isempty(local_behavior)
    behavior_file = fullfile(current_dir, local_behavior(1).name);
elseif isfile(data_behavior)
    behavior_file = data_behavior;
else
    error('GENERATE_SPEED_PEAK_FREQUENCY:behavior', ...
        'Behavior file not found for %s.', animal_code);
end

behavior_table = readtable(behavior_file);
behavior_table(:, 2) = [];
behavior_table.Properties.VariableNames = ...
    {'Animal', 'Start', 'End', 'Length', 'Type'};
behavior_table.Start = predict(synch_model, behavior_table.Start);
behavior_table.End = predict(synch_model, behavior_table.End);
session_mask = ismember(behavior_table.Animal, 'Session_structure') & ...
    ismember(behavior_table.Type, 'Partners session');
behavior_table = behavior_table(session_mask, :);
if isempty(behavior_table)
    ranges = zeros(0, 2);
    return
end
ranges = behavior_table{:, {'Start', 'End'}};
[ranges, order] = sortrows(ranges, 1);
if any(~isfinite(ranges), 'all') || any(ranges(:, 2) <= ranges(:, 1))
    error('GENERATE_SPEED_PEAK_FREQUENCY:behavior', ...
        'Invalid partner-session times for %s.', animal_code);
end
end

function [channel_lfp, n_samples] = read_selected_lfp(current_dir, rows, n_lfp_channels)
rows = rows(:);
mat_file = fullfile(current_dir, 'LFP_PAG.mat');
dat_file = fullfile(current_dir, 'LFP_PAG.dat');
if isfile(mat_file)
    loaded = load(mat_file, 'LFP');
    if ~isfield(loaded, 'LFP')
        error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', ...
            'LFP_PAG.mat does not contain LFP.');
    end
    if any(rows < 1 | rows > size(loaded.LFP, 1) | rows ~= round(rows))
        error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', ...
            'A selected LPAG row is outside LFP_PAG.mat.');
    end
    channel_lfp = double(loaded.LFP(rows, :));
    n_samples = size(channel_lfp, 2);
    return
end
if ~isfile(dat_file)
    error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', ...
        'Neither LFP_PAG.mat nor LFP_PAG.dat exists in %s.', current_dir);
end
file_info = dir(dat_file);
n_samples = file_info.bytes / (n_lfp_channels * 2);
if n_samples ~= round(n_samples) || n_samples < 1
    error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', ...
        'LFP_PAG.dat is not a multiple of %d int16 channels.', n_lfp_channels);
end
if any(rows < 1 | rows > n_lfp_channels | rows ~= round(rows))
    error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', ...
        'A selected LPAG row is outside the %d-channel LFP file.', n_lfp_channels);
end
channel_lfp = zeros(numel(rows), n_samples);
file_id = fopen(dat_file, 'r', 'ieee-le');
if file_id < 0
    error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', 'Could not open %s.', dat_file);
end
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
for row_n = 1:numel(rows)
    if fseek(file_id, 2 * (rows(row_n) - 1), 'bof') ~= 0
        error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', ...
            'Could not seek to LFP row %d.', rows(row_n));
    end
    [samples, count] = fread(file_id, n_samples, 'int16=>double', ...
        2 * (n_lfp_channels - 1));
    if count ~= n_samples
        error('GENERATE_SPEED_PEAK_FREQUENCY:lfp', ...
            'Read %d of %d samples on row %d.', count, n_samples, rows(row_n));
    end
    channel_lfp(row_n, :) = samples(:).';
end
end

function [xcoords, ycoords, chanMap] = load_channel_map(current_dir)
map_file = fullfile(current_dir, 'chann_map_PAG.mat');
if ~isfile(map_file)
    error('GENERATE_SPEED_PEAK_FREQUENCY:map', ...
        'Channel map not found: %s', map_file);
end
loaded = load(map_file, 'xcoords', 'ycoords', 'chanMap');
xcoords = loaded.xcoords;
ycoords = loaded.ycoords;
chanMap = loaded.chanMap;
end

function npx_type = neuropixels_type(xcoords)
hard_coded_x = [8 40; 258 290; 508 540; 758 790];
if any(ismember(xcoords(:), hard_coded_x(:)))
    npx_type = 2;
else
    npx_type = 1;
end
end

function data_root = resolve_data_root(current_dir, options)
recording_data_root = fileparts(fileparts(fileparts(current_dir)));
if isfile(fullfile(recording_data_root, 'Area_limits_GoodLooking.xlsx'))
    data_root = recording_data_root;
else
    data_root = options.data_root;
end
end

function [animal_code, table_name] = recording_names(current_dir)
[~, animal_code] = fileparts(current_dir);
parts = strsplit(strtrim(animal_code), ' ');
parts = parts(~cellfun(@isempty, parts));
if numel(parts) < 3
    error('GENERATE_SPEED_PEAK_FREQUENCY:name', ...
        'Unexpected Neuropixels folder name: %s', animal_code);
end
animal_batch = parts{1};
repeated_animal = parts{3};
if strcmp(repeated_animal, 'Single2')
    table_name = ['Batch', animal_batch(2), repeated_animal];
else
    if numel(animal_batch) < 4
        error('GENERATE_SPEED_PEAK_FREQUENCY:name', ...
            'Cannot derive the area-table animal name from %s.', animal_code);
    end
    table_name = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
end

function recording = empty_recording()
recording = struct( ...
    'animal_name', '', ...
    'animal_area_table_name', '', ...
    'neuropixels_directory', '', ...
    'npx_type', NaN, ...
    'lpag_channel', NaN, ...
    'lpag_lfp_row', NaN, ...
    'lpag_map_index', NaN, ...
    'lpag_probe_number', NaN, ...
    'lpag_x_um', NaN, ...
    'lpag_y_um', NaN, ...
    'lpag_depth_range', [NaN NaN], ...
    'parameters', struct(), ...
    'frequency', zeros(0, 1), ...
    'notes', '', ...
    'sessions', empty_session());
recording.sessions(1) = [];
end

function session = empty_session()
session = struct( ...
    'partner_index', NaN, ...
    'tracking_file', '', ...
    'session_index', NaN, ...
    'session_source', '', ...
    'behavior_time_range', [NaN NaN], ...
    'tracking_time_range', [NaN NaN], ...
    'analyzed_time_range', [NaN NaN], ...
    'speed_pixels_per_second', zeros(0, 1), ...
    'speed_time', zeros(0, 1), ...
    'speed_bin_id', zeros(0, 1), ...
    'speed_bin_edges', zeros(0, 1), ...
    'speed_bin_rule', '', ...
    'bins_within_upper_fraction', false, ...
    'bins_within_lower_fraction', false, ...
    'bins_within_requested_fractions', false, ...
    'speed_bins', struct([]), ...
    'window_sample_ranges', zeros(0, 2), ...
    'window_time_ranges', zeros(0, 2), ...
    'window_speed', zeros(0, 1), ...
    'window_bin_id', zeros(0, 1), ...
    'window_bin_occupancy', zeros(0, 1), ...
    'window_n_speed_samples', zeros(0, 1), ...
    'n_windows', 0, ...
    'welch', struct(), ...
    'fourier', struct(), ...
    'correlation_welch', struct(), ...
    'correlation_fourier', struct());
end
