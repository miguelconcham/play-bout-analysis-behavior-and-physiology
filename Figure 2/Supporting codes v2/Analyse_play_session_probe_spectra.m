%% Analyse_play_session_probe_spectra
% Single-session view of the channel-wise periodic/aperiodic spectra.
%
% The image is the full spectrum at every channel. The kilosort y coordinate
% is distance to the probe tip, not depth: larger values are further from the
% tip (upper), smaller values are closer to the tip (deeper).
%
% An optional figure pools sessions onto one spectrogram. Two switches
% choose the grouping. Partner sessions are the rows inside one recording
% file, so that figure is one recording day. Recording days pull every
% selected file whose folder starts with the same animal code, for example
% B1S3 in "B1S3 1007 Dual", and that figure also includes the partner
% sessions of those files. Each channel is log10 power with the aperiodic
% fit and the theta Gaussians removed, as in the LPAG analysis. Sites are
% averaged only when the shank and the x and y coordinates match. Channel
% number is not used, because the channel map can change between days and
% the same channel number can sit at a different distance. Frequency is
% the horizontal axis and distance to the tip is the vertical axis. A
% dotted line sits between neighboring sites when their region labels
% differ and the sessions at those coordinates agree.
%
% Each standard region also gets its own spectrum. The aperiodic fit is
% removed, and so are the periodic Gaussians whose centers fall in
% periodic_ranges_to_subtract. At every frequency the estimate is a mixed
% model, Power ~ 1 + (1|Session) + (1|Session:Channel), so channels are
% repeated measures nested in the session. The central LPAG channel is
% drawn as a dotted reference. A second figure tests region at each
% frequency and marks the region with the strongest estimate.
%
% Selected sessions are stacked with mid-LPAG at zero. The figures keep saved
% power, with no rescaling. A table written beside the spectra has the same
% distance, the area, the saved delta/theta/gamma power, that power minus the
% mid-LPAG channel, and that difference divided by the within-shank standard
% deviation. Before that z-score, a channel whose band power is more than
% artifact_trend_sd robust deviations from the local shank trend is removed.
% The scale is a median absolute deviation, and channels at the ends of the
% probe are judged against the trend of the sites inward from those ends.
% SupCol and VLPAG channels stay in the analysis only at distances from
% the LPAG-facing border of that area that more than half of the sessions
% with that region also reach. Those channels, and values that leave the
% shank trend, are removed before the distance correlation and the
% distance mixed model. A proximal/distal comparison can split each
% region into the nearer and farther halves of that available distance.
% It stays off unless plot_proximal_distal is set true.
% Each band and spectral method then has its own scatter of that
% z-scored power against distance to mid-LPAG, with transparent markers.
% The matching distance figure keeps those scatter values. At each
% distance bin, a session contributes one power (the mean of its channels
% in the bin). The model at that distance is Power ~ 1 + (1|Animal):
% a fixed intercept, with the animal's sessions as repeated measures.
% Spearman rho and its p-value use the same channel values.
% A separate figure compares the selected regions with a nested
% repeated-measures model of that same z-scored power. Every channel
% remains an observation nested in its partner session. Region is the
% fixed effect, and session and session-by-region are random. The figure
% shows every channel, a line for each session, and mixed-model brackets
% with Bonferroni correction.
% Figures use Fourier unless spectral_method is changed. The scatter,
% distance-model, region-comparison, and region-distance figures
% draw regions_to_display, which starts as LPAG, VLPAG, and SupCol. PAG sections,
% an upper side that passes through SupCol, and a deeper side that passes
% through DR are colored. Other territory stays gray. Those colored tracks
% are one figure. Its spectrum has each theta Gaussian removed, and any
% Gaussian whose center lies near twice that theta center (the first
% harmonic). The 95% confidence intervals for the three track classes
% are a separate figure, drawn when at least two tracks are available.
%
% Separate figures center each shank on the regions in regions_to_display.
% By default each ventral or dorsal curve is z-scored within the shank and
% set to 0 at its first channel. Set region_distance_scale to 'raw' for the
% saved power. Each panel reports a Spearman correlation of distance and
% power, with its p-value, after the same outlier screens. The center is
% the midpoint of that region's labeled span on the shank. bands_to_display
% selects delta, theta, and gamma. Delta is the default, and the points are
% a scatter. Curves run toward the tip (ventral) and toward the brain
% surface (dorsal). Another figure draws the same values and, at each
% distance, the intercept from Power ~ 1 + (1|Animal). A rectangle marks a
% run of distances where that interval excludes 0.

%% Which recording and session to inspect
saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
area_limit_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Area_limits_GoodLooking.xlsx';

% recording_index selects files from the numbered list printed below.
% session_index selects partner-session rows inside each of those files.
% Leave either one empty to take every recording, or every session.
%
% Every session of every saved recording:
%   recording_index = [];
%   session_index = [];
% Every session of recordings 1 through 12:
%   recording_index = 1:12;
%   session_index = [];
% recording_index = 1:12 and session_index = 1 plots only partner session 1
% of each file, so every figure is labeled session 1.
recording_index = [];
session_index = [];

% The stack is always drawn for the recordings and sessions selected above.
% Single-session figures are optional.
plot_each_session = false;

% Pooled spectrogram. Off until this switch is set true. The color is
% log10 power minus the aperiodic fit minus theta Gaussians whose centers
% fall inside theta_band, the same subtraction as the LPAG analysis.
% The two switches below choose which sessions share a figure. Sites are
% averaged only at a matching shank and matching x and y. Channel number
% is ignored, because the channel map can change between days.
plot_pooled_animal_spectra = false;

% One figure per recording folder. Partner sessions are the rows inside
% that file. Sites that share a shank and the same x and y are averaged.
pool_partner_sessions = true;

% One figure per recorded animal. The animal is the first word of the
% folder name, for example B1S3 in "B1S3 1007 Dual". Every selected
% recording of that animal is included, and so are its partner sessions.
% The same physical site is averaged across days only when shank, x, and
% y match. A channel that moved to new coordinates stays a separate row.
pool_recording_days = true;

% log10 power minus the fitted aperiodic component. Switch to
% 'log10_mean_power' to draw the session-mean spectrum itself.
spectrum_field = 'raw_aperiodic_subtracted';

delta_band = [1 5];
theta_band = [6 12];
gamma_band = [30 45];

% Power and correlation figures use the Fourier spectra. Welch is available:
%   spectral_method = {'Welch'};
%   spectral_method = {'Welch', 'Fourier'};
spectral_method = {'Fourier'};

% Bands on the power and correlation figures. Delta is the default.
%   bands_to_display = {'delta', 'theta'};
%   bands_to_display = {'delta', 'theta', 'gamma'};
bands_to_display = {'delta'};

% Regions drawn on the scatter, distance-model, region-comparison,
% and region-distance figures. SupCol is the superior colliculus. Add 'DR'
% to include it:
%   regions_to_display = {'LPAG', 'VLPAG', 'SupCol', 'DR'};
regions_to_display = {'LPAG', 'VLPAG', 'SupCol'};

% Spectra of log10 power after the aperiodic fit and the periodic ranges
% below are removed. One subplot per region. The list is the initial set.
spectrum_regions = {'SupCol', 'DLPAG', 'LPAG', 'VLPAG', 'DR'};

% Each cell is one region-spectrum analysis. Gaussians whose centers fall
% inside that cell are removed, together with the aperiodic fit. Peaks
% outside the cell stay, so {[1 5]} removes delta and leaves theta.
%   periodic_ranges_to_subtract = {[6 12]};
%   periodic_ranges_to_subtract = {[1 5]};
%   periodic_ranges_to_subtract = {[1 5; 6 12]};
periodic_ranges_to_subtract = {[6 12], [1 5]};

% Ventral and dorsal curves. 'normalized' z-scores each band within the
% shank and subtracts the power of the first channel on that side, so the
% curve starts at 0. 'raw' draws the saved band power.
region_distance_scale = 'normalized';

% Drop a band value when it leaves the robust local trend of its shank by
% more than this many median-absolute-deviation units. 5 is far beyond the
% smooth change between neighboring sites. End-of-probe channels use the
% trend estimated inward of those ends.
artifact_trend_sd = 5;

% SupCol and VLPAG coverage is measured from the border of that area that
% faces LPAG, in bins of this width. A bin stays when more than half of
% the sessions with that region also have a channel there. Near-border
% sites stay in even when LPAG itself is thin. The far tail, reached by
% only a few sessions, is left out. Only those channels are removed.
region_border_bin_um = 100;

% The proximal/distal comparison is a second figure. It stays off until
% this switch is set true.
plot_proximal_distal = false;

% Distance figures keep the scatter. Channels that fall in the same bin
% and the same session become one session power. At that distance the
% model is an intercept with a random effect of animal.
distance_model_bin_um = 20;

%% Distance-figure appearance
% Scatter dots, the confidence shade around each intercept, and the
% rectangles that mark a significant stretch of distance.
scatter_marker_size = 6;
scatter_face_alpha = 0.08;
ci_face_alpha = 0.45;
significance_face_alpha = 0.16;
% A rectangle is drawn only when the intercept interval excludes 0 for a
% contiguous run of at least this many micrometres. The rectangle covers
% the upper 90% of the y-axis.
significance_min_um = 100;

plot_style = struct( ...
    'marker_size', scatter_marker_size, ...
    'marker_alpha', scatter_face_alpha, ...
    'ci_alpha', ci_face_alpha, ...
    'significance_alpha', significance_face_alpha, ...
    'significance_min_um', significance_min_um);

%% Load the requested recording files
file_patterns = probe_periodic_spectra_filename();
file_list = dir(fullfile(saving_folder, file_patterns{1}));
if ~isempty(file_list)
    keep = false(numel(file_list), 1);
    for file_n = 1:numel(file_list)
        [~, ~, extension] = fileparts(file_list(file_n).name);
        keep(file_n) = ~file_list(file_n).isdir && ...
            (isempty(extension) || strcmpi(extension, '.mat'));
    end
    file_list = file_list(keep);
end
if isempty(file_list)
    error(['No channel-wise spectra found in %s. Expected names such as ' ...
        'play_session_probe_periodic_spectra_B1D1_1007_Dual.mat.'], ...
        saving_folder);
end
[~, order] = sort({file_list.name});
file_list = file_list(order);

if isempty(recording_index)
    recording_index = 1:numel(file_list);
end
if any(recording_index < 1) || any(recording_index > numel(file_list))
    error('recording_index must be between 1 and %d.', numel(file_list));
end

fprintf('Saved recordings:\n');
for file_n = 1:numel(file_list)
    fprintf('  %d  %s\n', file_n, file_list(file_n).name);
end

%% Stack every selected shank at its mid-LPAG channel
stacked_frequency = [];
tracks = repmat(empty_track(), 0, 1);
placed = repmat(empty_placed(), 0, 1);
power_pieces = {};
pooled_pieces = {};
range_sets = periodic_range_sets(periodic_ranges_to_subtract);
region_pieces = cell(1, numel(range_sets));
for set_n = 1:numel(range_sets)
    region_pieces{set_n} = {};
end

%% Plot each requested session on its own
for file_n = recording_index(:).'
    input_file = fullfile(saving_folder, file_list(file_n).name);
    fprintf('Loading %s\n', input_file);
    loaded = load(input_file, 'probe_spectra');
    probe_spectra = loaded.probe_spectra;
    frequency = probe_spectra.frequency(:).';
    if isempty(stacked_frequency)
        stacked_frequency = frequency;
    elseif numel(stacked_frequency) ~= numel(frequency) || ...
            max(abs(stacked_frequency - frequency)) > 1e-9
        error('Frequency grids differ between saved recordings.');
    end
    channel_table = probe_spectra.channel_table;
    if frequency(end) + 1e-9 < gamma_band(2)
        error('Spectra in %s end at %g Hz, before the gamma band upper edge.', ...
            probe_spectra.animal_name, frequency(end));
    end

    if ~isfield(probe_spectra.sessions(1).welch, spectrum_field)
        error('Saved spectra have no field "%s".', spectrum_field);
    end

    area_borders = area_borders_for_recording(probe_spectra, area_limit_file);

    if isempty(session_index)
        sessions_in_file = 1:probe_spectra.n_sessions;
    else
        sessions_in_file = session_index;
    end
    if any(sessions_in_file < 1) || any(sessions_in_file > probe_spectra.n_sessions)
        error('%s has %d partner sessions.', probe_spectra.animal_name, ...
            probe_spectra.n_sessions);
    end

    for session_n = sessions_in_file(:).'
        this_session = probe_spectra.sessions(session_n);
        welch_spectrum = this_session.welch.(spectrum_field);
        fourier_spectrum = this_session.fourier.(spectrum_field);
        session_label = sprintf( ...
            'recording %d | %s | session %d of %d | %.0f-%.0f s | %d windows', ...
            file_n, probe_spectra.animal_name, session_n, ...
            probe_spectra.n_sessions, ...
            this_session.time_range(1), this_session.time_range(2), ...
            this_session.n_windows);
        fprintf('Including %s\n', session_label);
        tracks = [tracks; tracks_from_session(channel_table, ...
            welch_spectrum, fourier_spectrum, frequency, ...
            delta_band, theta_band, gamma_band, session_label, ...
            artifact_trend_sd, ...
            gaussian_params_for(this_session.welch, size(welch_spectrum, 1)), ...
            gaussian_params_for(this_session.fourier, size(fourier_spectrum, 1)))]; %#ok<AGROW>
        placed = [placed; placed_from_session(channel_table, ...
            welch_spectrum, fourier_spectrum, frequency, ...
            delta_band, theta_band, gamma_band, session_label, ...
            artifact_trend_sd, probe_spectra.animal_name, session_n)]; %#ok<AGROW>
        power_pieces{end + 1} = lpag_distance_table(channel_table, ...
            welch_spectrum, fourier_spectrum, frequency, ...
            delta_band, theta_band, gamma_band, ...
            probe_spectra.animal_name, session_n, spectrum_field, ...
            artifact_trend_sd); %#ok<SAGROW>
        if plot_pooled_animal_spectra && ...
                (pool_partner_sessions || pool_recording_days)
            pooled_pieces{end + 1} = residual_sites_from_session( ...
                channel_table, this_session, frequency, theta_band, ...
                probe_spectra.animal_name, session_n, spectral_method); %#ok<SAGROW>
        end
        for set_n = 1:numel(range_sets)
            region_pieces{set_n}{end + 1} = region_channel_records( ...
                channel_table, this_session, frequency, range_sets{set_n}, ...
                probe_spectra.animal_name, session_n, spectral_method); %#ok<SAGROW>
        end

        if plot_each_session
            plot_depth_spectrum(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, spectrum_field, ...
                session_label, area_borders, spectral_method);
            plot_band_profiles(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, ...
                delta_band, theta_band, gamma_band, ...
                spectrum_field, session_label, area_borders, spectral_method, ...
                artifact_trend_sd, bands_to_display);
            drawnow
        end
    end
end

if isempty(placed)
    error('No probe channels were found in the selected recordings.');
end
if plot_pooled_animal_spectra
    if ~pool_partner_sessions && ~pool_recording_days
        fprintf(['Pooled residual spectra were skipped because both ' ...
            'pool_partner_sessions and pool_recording_days are false.\n']);
    else
        if isempty(pooled_pieces)
            pooled_sites = repmat(empty_residual_site(), 0, 1);
        else
            pooled_sites = vertcat(pooled_pieces{:});
        end
        plot_pooled_residual_spectra(pooled_sites, stacked_frequency, ...
            spectral_method, theta_band, pool_partner_sessions, ...
            pool_recording_days);
    end
end
power_table = write_lpag_distance_table(power_pieces, saving_folder, ...
    spectrum_field, region_border_bin_um);
placed = apply_screened_power_to_placed(placed, power_table);
plot_normalized_power_scatters(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display, bands_to_display, plot_style);
plot_normalized_power_intervals(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display, distance_model_bin_um, ...
    bands_to_display, plot_style);
plot_region_power_comparison(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display, bands_to_display, plot_style);
if plot_proximal_distal
    plot_region_half_comparison(power_table, delta_band, theta_band, gamma_band, ...
        spectral_method, regions_to_display, bands_to_display, plot_style);
end
if isempty(tracks)
    fprintf(['No shank with a mid-LPAG channel was found; the mid-LPAG ' ...
        'stack is skipped.\n']);
else
    fprintf(['Aligned tracks: %d. Stay in PAG: %d. Upper side through SupCol: ' ...
        '%d. Deeper side through DR: %d.\n'], numel(tracks), ...
        sum([tracks.stays_pag]), sum([tracks.goes_sc]), sum([tracks.goes_dr]));
    plot_aligned_spectra(tracks, stacked_frequency, spectrum_field, ...
        spectral_method);
    plot_aligned_profiles(tracks, spectrum_field, ...
        delta_band, theta_band, gamma_band, spectral_method, bands_to_display);
    plot_aligned_intervals(tracks, spectrum_field, ...
        delta_band, theta_band, gamma_band, spectral_method, bands_to_display);
end
plot_region_distances(placed, spectrum_field, ...
    delta_band, theta_band, gamma_band, spectral_method, ...
    regions_to_display, region_distance_scale, distance_model_bin_um, ...
    bands_to_display, plot_style);
for set_n = 1:numel(range_sets)
    if isempty(region_pieces{set_n})
        region_records = repmat(empty_region_record(), 0, 1);
    else
        region_records = vertcat(region_pieces{set_n}{:});
    end
    plot_region_residual_spectra(region_records, stacked_frequency, ...
        spectral_method, range_sets{set_n}, spectrum_regions, plot_style);
end

function plot_depth_spectrum(frequency, channel_table, ...
        welch_spectrum, fourier_spectrum, spectrum_field, ...
        session_label, area_borders, spectral_method)
probe_groups = probe_groups_by_depth(channel_table);
n_groups = numel(probe_groups);
spectra = {welch_spectrum, fourier_spectrum};
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);

figure('Name', ['Spectrum by depth | ', session_label], 'Color', 'w', ...
    'Position', [80 80 max(360, 280 * n_groups) 720]);
colors = diverging_or_sequential(spectrum_field);

for method_slot = 1:n_methods
    method_n = method_ids(method_slot);
    panel_values = cell(n_groups, 1);
    finite_values = [];
    for group_n = 1:n_groups
        panel_values{group_n} = ...
            spectra{method_n}(probe_groups(group_n).channel_index, :);
        finite_values = [finite_values; ...
            panel_values{group_n}(isfinite(panel_values{group_n}))]; %#ok<AGROW>
    end
    color_limits = spectrum_color_limits(finite_values, spectrum_field);

    for group_n = 1:n_groups
        subplot(n_methods, n_groups, (method_slot - 1) * n_groups + group_n)
        this_group = probe_groups(group_n);
        depth = this_group.depth_um(:);
        use_depth_axis = numel(depth) > 1 && all(isfinite(depth)) && ...
            all(diff(depth) > 0);
        if use_depth_axis
            imagesc(frequency, depth, panel_values{group_n});
            center_y = depth(this_group.is_lpag_center);
        else
            imagesc(frequency, 1:numel(depth), panel_values{group_n});
            depth_ticks(depth);
            center_y = find(this_group.is_lpag_center);
        end
        axis tight
        colormap(colors)
        if all(isfinite(color_limits))
            caxis(color_limits);
        end
        these_borders = borders_for_probe(area_borders, this_group.probe_number);
        if use_depth_axis
            ylim(depth_limits_with_borders(depth, these_borders))
        end
        mark_depth_lines(center_y)
        draw_area_borders(these_borders, depth, use_depth_axis)
        set(gca, 'YDir', 'normal')
        xlabel('Frequency (Hz)')
        if group_n == 1
            ylabel('Distance to tip (µm)')
        end
        title(sprintf('%s | %s', method_names{method_slot}, this_group.label))
        box off
    end
    colorbar
end
sgtitle({session_label, spectrum_title(spectrum_field)}, 'Interpreter', 'none')
end

function item = empty_residual_site()
item = struct( ...
    'recording', '', ...
    'animal_code', '', ...
    'session', NaN, ...
    'shank', NaN, ...
    'x_um', NaN, ...
    'y_um', NaN, ...
    'area_key', '', ...
    'welch', [], ...
    'fourier', []);
end

function sites = residual_sites_from_session(channel_table, session, ...
        frequency, theta_band, animal_name, session_n, spectral_method)
[~, method_ids] = selected_spectral_methods(spectral_method);
method_fields = {'welch', 'fourier'};
n_channels = height(channel_table);
sites = repmat(empty_residual_site(), n_channels, 1);
recording = char(string(animal_name));
animal_code = animal_identity(recording);
residuals = cell(1, 2);
for slot_n = 1:numel(method_ids)
    method_n = method_ids(slot_n);
    residuals{method_n} = residual_matrix( ...
        session.(method_fields{method_n}), frequency, theta_band, n_channels);
end
for channel_n = 1:n_channels
    sites(channel_n).recording = recording;
    sites(channel_n).animal_code = animal_code;
    sites(channel_n).session = session_n;
    sites(channel_n).shank = channel_table.ProbeNumber(channel_n);
    sites(channel_n).x_um = channel_table.Xum(channel_n);
    sites(channel_n).y_um = channel_table.Yum(channel_n);
    sites(channel_n).area_key = region_border_key( ...
        area_label_at(channel_table.Area, channel_n));
    if ~isempty(residuals{1})
        sites(channel_n).welch = residuals{1}(channel_n, :);
    end
    if ~isempty(residuals{2})
        sites(channel_n).fourier = residuals{2}(channel_n, :);
    end
end
end

function code = animal_identity(name)
% The recorded animal is the first chunk of the folder name.
% "B1S3 1007 Dual" and "B1S3_1007_Dual" both give B1S3.
name = strtrim(char(string(name)));
if isempty(name)
    code = '';
    return
end
space_parts = strsplit(name);
if numel(space_parts) > 1
    code = space_parts{1};
    return
end
underscore_parts = strsplit(name, '_');
code = underscore_parts{1};
end

function residual = residual_matrix(fitted, frequency, ranges, n_channels)
required = {'log10_mean_power', 'aperiodic_fit', 'gaussian_params'};
for field_n = 1:numel(required)
    if ~isfield(fitted, required{field_n})
        error(['Saved spectra have no %s, so the periodic ranges and the ' ...
            'aperiodic fit cannot be removed.'], required{field_n});
    end
end
log_power = fitted.log10_mean_power;
aperiodic = fitted.aperiodic_fit;
n_frequency = numel(frequency);
if size(log_power, 1) ~= n_channels || size(aperiodic, 1) ~= n_channels || ...
        size(log_power, 2) ~= n_frequency || size(aperiodic, 2) ~= n_frequency
    error(['Saved log10 power and the aperiodic fit do not match the ' ...
        'channel table or the frequency grid.']);
end
if ~iscell(fitted.gaussian_params) || numel(fitted.gaussian_params) ~= n_channels
    error('gaussian_params must hold one peak matrix per channel.');
end
ranges = periodic_range_matrix(ranges);
residual = nan(n_channels, n_frequency);
for channel_n = 1:n_channels
    residual(channel_n, :) = subtract_periodic_and_aperiodic( ...
        log_power(channel_n, :), aperiodic(channel_n, :), ...
        fitted.gaussian_params{channel_n}, frequency, ranges);
end
end

function residual = subtract_periodic_and_aperiodic( ...
        log_power, aperiodic, gaussians, frequency, ranges)
% FOOOF components are additive in log10 power. A peak is removed when
% its center lies inside one of the requested ranges. Other peaks stay.
residual = log_power(:).' - aperiodic(:).';
if isempty(gaussians) || size(gaussians, 2) < 3 || size(gaussians, 1) < 1
    return
end
center = gaussians(:, 1);
peak_height = gaussians(:, 2);
width = gaussians(:, 3);
keep = isfinite(center) & isfinite(peak_height) & isfinite(width) & ...
    width ~= 0 & center_in_ranges(center, ranges);
center = center(keep);
peak_height = peak_height(keep);
width = width(keep);
residual = residual - sum_gaussians(center, peak_height, width, frequency);
end

function sets = periodic_range_sets(ranges)
if iscell(ranges)
    sets = cell(1, numel(ranges));
    for set_n = 1:numel(ranges)
        sets{set_n} = periodic_range_matrix(ranges{set_n});
    end
    if isempty(sets)
        error('periodic_ranges_to_subtract needs at least one range set.');
    end
    return
end
sets = {periodic_range_matrix(ranges)};
end

function ranges = periodic_range_matrix(ranges)
if isempty(ranges)
    ranges = zeros(0, 2);
    return
end
ranges = double(ranges);
if isvector(ranges)
    if numel(ranges) ~= 2
        error('Each periodic range needs a low and a high frequency.');
    end
    ranges = ranges(:).';
end
if size(ranges, 2) ~= 2
    error(['periodic_ranges_to_subtract must be [low high], a matrix ' ...
        'with those two columns, or a cell of those.']);
end
if any(~isfinite(ranges(:))) || any(ranges(:, 1) > ranges(:, 2))
    error('Each periodic range needs finite edges with low <= high.');
end
end

function inside = center_in_ranges(center, ranges)
inside = false(size(center));
for range_n = 1:size(ranges, 1)
    inside = inside | (center >= ranges(range_n, 1) & ...
        center <= ranges(range_n, 2));
end
end

function label = periodic_range_label(ranges)
if isempty(ranges)
    label = 'no periodic range';
    return
end
label = sprintf('%g-%g Hz', ranges(1, 1), ranges(1, 2));
for range_n = 2:size(ranges, 1)
    label = sprintf('%s, %g-%g Hz', label, ranges(range_n, 1), ranges(range_n, 2));
end
end

function plot_pooled_residual_spectra(sites, frequency, spectral_method, ...
        theta_band, pool_partner_sessions, pool_recording_days)
if isempty(sites)
    fprintf('No channels were available for the pooled residual spectra.\n');
    return
end
[method_names, method_ids] = selected_spectral_methods(spectral_method);
method_fields = {'welch', 'fourier'};
if pool_partner_sessions
    recordings = unique({sites.recording}, 'stable');
    for recording_n = 1:numel(recordings)
        recording_name = recordings{recording_n};
        these_sites = sites(strcmp({sites.recording}, recording_name));
        plot_one_animal_residual(these_sites, recording_name, frequency, ...
            method_names, method_ids, method_fields, theta_band, 'partner');
    end
end
if pool_recording_days
    animal_codes = unique({sites.animal_code}, 'stable');
    for code_n = 1:numel(animal_codes)
        animal_code = animal_codes{code_n};
        these_sites = sites(strcmp({sites.animal_code}, animal_code));
        plot_one_animal_residual(these_sites, animal_code, frequency, ...
            method_names, method_ids, method_fields, theta_band, 'days');
    end
end
end

function plot_one_animal_residual(sites, group_name, frequency, ...
        method_names, method_ids, method_fields, theta_band, pool_mode)
shank_id = reshape([sites.shank], [], 1);
shank_list = unique(shank_id(isfinite(shank_id)));
shank_list = sort(shank_list);
if any(~isfinite(shank_id))
    shank_list(end + 1, 1) = NaN; %#ok<AGROW>
end
n_shanks = numel(shank_list);
n_methods = numel(method_ids);
panels = cell(n_methods, n_shanks);
finite_values = cell(n_methods, 1);
for shank_n = 1:n_shanks
    if isfinite(shank_list(shank_n))
        in_shank = shank_id == shank_list(shank_n);
        shank_label = sprintf('Shank %g', shank_list(shank_n));
    else
        in_shank = ~isfinite(shank_id);
        shank_label = 'Unassigned';
    end
    shank_sites = sites(in_shank);
    for slot_n = 1:n_methods
        method_n = method_ids(slot_n);
        pooled = pool_sites_by_coordinate( ...
            shank_sites, method_fields{method_n});
        pooled.label = shank_label;
        panels{slot_n, shank_n} = pooled;
        if ~isempty(pooled.spectrum)
            finite_values{slot_n} = [finite_values{slot_n}; ...
                pooled.spectrum(isfinite(pooled.spectrum))]; %#ok<AGROW>
        end
    end
end

n_rows = 0;
n_merged = 0;
for shank_n = 1:n_shanks
    n_rows = n_rows + numel(panels{1, shank_n}.depth_um);
    n_merged = n_merged + nnz(panels{1, shank_n}.n_averaged > 1);
end
[n_recordings, n_partner] = pooled_group_counts(sites);
if strcmp(pool_mode, 'days')
    figure_prefix = 'Pooled recording days';
    match_line = ['Sites averaged across days only where shank, x, and y ' ...
        'match. Channel number is not used.'];
else
    figure_prefix = 'Pooled partner sessions';
    match_line = ['Partner sessions of one recording day, averaged only ' ...
        'where shank, x, and y match. Channel number is not used.'];
end
grouping_line = sprintf('%d recording folders | %d partner sessions', ...
    n_recordings, n_partner);
fprintf(['%s | %s | %s | %d sites from %d session channels | ' ...
    '%d sites combined matching coordinates\n'], ...
    figure_prefix, group_name, grouping_line, n_rows, numel(sites), n_merged);

figure('Name', [figure_prefix, ' | ', group_name], 'Color', 'w', ...
    'Position', [80 80 max(360, 280 * max(n_shanks, 1)) 720]);
colors = diverging_or_sequential('raw_aperiodic_subtracted');
for slot_n = 1:n_methods
    color_limits = spectrum_color_limits( ...
        finite_values{slot_n}, 'raw_aperiodic_subtracted');
    for shank_n = 1:n_shanks
        subplot(n_methods, n_shanks, (slot_n - 1) * n_shanks + shank_n)
        pooled = panels{slot_n, shank_n};
        depth = pooled.depth_um(:);
        if isempty(depth)
            axis off
            title(sprintf('%s | %s | no sites', ...
                method_names{slot_n}, pooled.label), 'Interpreter', 'none')
            continue
        end
        use_depth_axis = numel(depth) > 1 && all(isfinite(depth)) && ...
            all(diff(depth) > 0);
        if use_depth_axis
            imagesc(frequency, depth, pooled.spectrum);
        else
            imagesc(frequency, 1:numel(depth), pooled.spectrum);
            coordinate_ticks(depth, pooled.x_um);
        end
        axis tight
        colormap(colors)
        if all(isfinite(color_limits))
            caxis(color_limits);
        end
        set(gca, 'YDir', 'normal')
        mark_region_divisions(depth, pooled.area_key, use_depth_axis)
        xlabel('Frequency (Hz)')
        if shank_n == 1
            ylabel('Distance to tip (µm)')
        end
        title(sprintf('%s | %s | %d sites', method_names{slot_n}, ...
            pooled.label, numel(depth)), 'Interpreter', 'none')
        box off
    end
    color_bar = colorbar;
    ylabel(color_bar, {'log10 power', 'minus aperiodic and theta'})
end
sgtitle({[figure_prefix, ' | ', group_name], grouping_line, ...
    sprintf('log10 power minus aperiodic fit and theta (%.0f-%.0f Hz)', ...
    theta_band(1), theta_band(2)), match_line}, 'Interpreter', 'none')
drawnow
end

function [n_recordings, n_partner] = pooled_group_counts(sites)
n_recordings = numel(unique({sites.recording}, 'stable'));
session_keys = cell(numel(sites), 1);
for site_n = 1:numel(sites)
    session_keys{site_n} = sprintf('%s|%g', sites(site_n).recording, ...
        sites(site_n).session);
end
n_partner = numel(unique(session_keys, 'stable'));
end

function pooled = pool_sites_by_coordinate(sites, method_field)
pooled = struct('label', '', 'depth_um', [], 'x_um', [], ...
    'spectrum', [], 'n_averaged', [], 'area_key', {{}});
if isempty(sites)
    return
end
n_sites = numel(sites);
x_um = reshape([sites.x_um], [], 1);
y_um = reshape([sites.y_um], [], 1);
spectrum = vertcat(sites.(method_field));
if isempty(spectrum)
    return
end
% Neuropixels coordinates are integer micrometres. Rounding keeps the
% same physical site together when a saved value carries a fraction.
x_key = round(x_um);
y_key = round(y_um);
row_id = zeros(n_sites, 1);
missing = ~isfinite(x_key) | ~isfinite(y_key);
row_id(missing) = find(missing);
x_key(~isfinite(x_key)) = 0;
y_key(~isfinite(y_key)) = 0;
[~, ~, group] = unique([x_key, y_key, row_id], 'rows', 'stable');
n_groups = max(group);
depth = nan(n_groups, 1);
lateral = nan(n_groups, 1);
n_averaged = zeros(n_groups, 1);
area_key = cell(n_groups, 1);
pooled_spectrum = nan(n_groups, size(spectrum, 2));
site_areas = {sites.area_key};
site_areas = site_areas(:);
for group_n = 1:n_groups
    in_group = group == group_n;
    depth(group_n) = mean(y_um(in_group), 'omitnan');
    lateral(group_n) = mean(x_um(in_group), 'omitnan');
    n_averaged(group_n) = nnz(in_group);
    area_key{group_n} = agreed_area_key(site_areas(in_group));
    pooled_spectrum(group_n, :) = mean(spectrum(in_group, :), 1, 'omitnan');
end
order_y = depth;
order_y(~isfinite(order_y)) = inf;
order_x = lateral;
order_x(~isfinite(order_x)) = inf;
[~, order] = sortrows([order_y, order_x], [1 2]);
pooled.depth_um = depth(order);
pooled.x_um = lateral(order);
pooled.spectrum = pooled_spectrum(order, :);
pooled.n_averaged = n_averaged(order);
pooled.area_key = area_key(order);
end

function key = region_border_key(label)
names = area_to_region_keys(label);
if isempty(names)
    key = '';
    return
end
names = sort(names);
key = names{1};
for name_n = 2:numel(names)
    key = [key, ',', names{name_n}]; %#ok<AGROW>
end
end

function key = agreed_area_key(keys)
keys = keys(~cellfun('isempty', keys));
if isempty(keys)
    key = '';
    return
end
unique_keys = unique(keys, 'stable');
if numel(unique_keys) == 1
    key = unique_keys{1};
else
    key = '';
end
end

function mark_region_divisions(depth_um, area_keys, use_depth_axis)
if numel(depth_um) < 2 || numel(area_keys) ~= numel(depth_um)
    return
end
hold on
x_limits = xlim;
for row_n = 1:numel(depth_um) - 1
    below = area_keys{row_n};
    above = area_keys{row_n + 1};
    if isempty(below) || isempty(above) || strcmp(below, above)
        continue
    end
    if use_depth_axis
        y_value = (depth_um(row_n) + depth_um(row_n + 1)) / 2;
    else
        y_value = row_n + 0.5;
    end
    if ~isfinite(y_value)
        continue
    end
    yline(y_value, ':', 'Color', [0.15 0.15 0.15], ...
        'LineWidth', 1.1, 'HandleVisibility', 'off');
    text(x_limits(1), y_value, ['  ', ...
        region_border_label(below), ' | ', region_border_label(above)], ...
        'Color', [0.15 0.15 0.15], 'FontSize', 7, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'Interpreter', 'none', 'Clipping', 'on');
end
end

function label = region_border_label(key)
parts = strsplit(key, ',');
label = region_display_name(parts{1});
for part_n = 2:numel(parts)
    label = [label, ', ', region_display_name(parts{part_n})]; %#ok<AGROW>
end
end

function coordinate_ticks(depth_um, x_um)
n_rows = numel(depth_um);
n_ticks = min(8, n_rows);
if n_ticks < 1
    return
end
tick_index = unique(round(linspace(1, n_rows, n_ticks)));
y_key = round(depth_um);
tick_labels = cell(numel(tick_index), 1);
for tick_n = 1:numel(tick_index)
    row_n = tick_index(tick_n);
    if ~isfinite(depth_um(row_n))
        tick_labels{tick_n} = '';
        continue
    end
    shared = isfinite(y_key(row_n)) && nnz(y_key == y_key(row_n)) > 1;
    if shared && isfinite(x_um(row_n))
        tick_labels{tick_n} = sprintf('%.0f (x %.0f)', ...
            depth_um(row_n), x_um(row_n));
    else
        tick_labels{tick_n} = sprintf('%.0f', depth_um(row_n));
    end
end
set(gca, 'YTick', tick_index, 'YTickLabel', tick_labels)
end

function plot_band_profiles(frequency, channel_table, ...
        welch_spectrum, fourier_spectrum, ...
        delta_band, theta_band, gamma_band, spectrum_field, ...
        session_label, area_borders, spectral_method, artifact_trend_sd, ...
        bands_to_display)
probe_groups = probe_groups_by_depth(channel_table);
n_groups = numel(probe_groups);
spectra = {welch_spectrum, fourier_spectrum};
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);
selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band);

figure('Name', ['Band profiles | ', session_label], 'Color', 'w', ...
    'Position', [120 120 max(360, 280 * n_groups) 720]);

for method_slot = 1:n_methods
    method_n = method_ids(method_slot);
    for group_n = 1:n_groups
        subplot(n_methods, n_groups, (method_slot - 1) * n_groups + group_n)
        this_group = probe_groups(group_n);
        this_spectrum = spectra{method_n}(this_group.channel_index, :);
        band_power = [ ...
            mean_in_band(this_spectrum, frequency, delta_band), ...
            mean_in_band(this_spectrum, frequency, theta_band), ...
            mean_in_band(this_spectrum, frequency, gamma_band)];
        band_power = reject_trend_outliers(this_group.depth_um, ...
            band_power, artifact_trend_sd);
        hold on
        line_handles = gobjects(numel(selection.index), 1);
        for slot_n = 1:numel(selection.index)
            band_n = selection.index(slot_n);
            line_handles(slot_n) = plot(band_power(:, band_n), ...
                this_group.depth_um, ...
                'Color', selection.colors(slot_n, :), 'LineWidth', 1.4);
        end
        if strcmp(spectrum_field, 'raw_aperiodic_subtracted')
            xline(0, ':k');
        end
        finite_depth = this_group.depth_um(isfinite(this_group.depth_um));
        center_depth = this_group.depth_um(this_group.is_lpag_center);
        center_depth = center_depth(isfinite(center_depth));
        for center_n = 1:numel(center_depth)
            yline(center_depth(center_n), '--k', 'LPAG center', ...
                'LabelHorizontalAlignment', 'left', 'FontSize', 8, ...
                'HandleVisibility', 'off');
        end
        set(gca, 'YDir', 'normal')
        these_borders = borders_for_probe(area_borders, this_group.probe_number);
        if ~isempty(finite_depth)
            ylim(depth_limits_with_borders(finite_depth, these_borders))
        end
        draw_area_borders(these_borders, this_group.depth_um, true)
        xlabel(profile_xlabel(spectrum_field))
        if group_n == 1
            ylabel('Distance to tip (µm)')
        end
        title(sprintf('%s | %s', method_names{method_slot}, this_group.label))
        if method_slot == 1 && group_n == 1
            legend(line_handles, selection.labels, 'Location', 'best')
        end
        box off
    end
end
sgtitle({session_label, 'Band means along the probe'}, 'Interpreter', 'none')
end

function probe_groups = probe_groups_by_depth(channel_table)
probe_numbers = channel_table.ProbeNumber;
assigned = isfinite(probe_numbers);
probe_ids = unique(probe_numbers(assigned));
labels = arrayfun(@(probe_id) sprintf('Shank %g', probe_id), ...
    probe_ids, 'UniformOutput', false);
group_probe_ids = probe_ids;
if any(~assigned)
    group_probe_ids(end + 1, 1) = NaN; %#ok<AGROW>
    labels{end + 1} = 'Unassigned'; %#ok<AGROW>
end

probe_groups = repmat(struct( ...
    'label', '', ...
    'probe_number', [], ...
    'channel_index', [], ...
    'depth_um', [], ...
    'areas', {{}}, ...
    'is_lpag_center', []), numel(group_probe_ids), 1);

for group_n = 1:numel(group_probe_ids)
    if isfinite(group_probe_ids(group_n))
        in_group = probe_numbers == group_probe_ids(group_n);
    else
        in_group = ~assigned;
    end
    channel_index = find(in_group);
    depth_um = channel_table.Yum(channel_index);
    depth_key = depth_um;
    depth_key(~isfinite(depth_key)) = Inf;
    x_key = channel_table.Xum(channel_index);
    x_key(~isfinite(x_key)) = Inf;
    [~, order] = sortrows([depth_key, x_key], [1 2]);
    channel_index = channel_index(order);
    probe_groups(group_n).label = labels{group_n};
    probe_groups(group_n).probe_number = group_probe_ids(group_n);
    probe_groups(group_n).channel_index = channel_index;
    probe_groups(group_n).depth_um = channel_table.Yum(channel_index);
    probe_groups(group_n).areas = cellstr(string(channel_table.Area(channel_index)));
    probe_groups(group_n).is_lpag_center = ...
        logical(channel_table.IsLPAGCenter(channel_index));
end
end

function band_power = mean_in_band(spectrum, frequency, band)
in_band = frequency >= band(1) & frequency <= band(2);
if ~any(in_band)
    band_power = nan(size(spectrum, 1), 1);
    return
end
band_power = mean(spectrum(:, in_band), 2, 'omitnan');
end

function [bands, n_removed, n_extreme] = reject_trend_outliers( ...
        distance_um, bands, n_fence)
n_removed = 0;
n_extreme = 0;
if isempty(bands) || ~isfinite(n_fence) || n_fence <= 0
    return
end
distance_um = distance_um(:);
for band_n = 1:size(bands, 2)
    [bands(:, band_n), removed_n, extreme_n] = reject_one_band( ...
        distance_um, bands(:, band_n), n_fence);
    n_removed = n_removed + removed_n;
    n_extreme = n_extreme + extreme_n;
end
end

function [power, n_removed, n_extreme] = reject_one_band( ...
        distance_um, power, n_fence)
n_removed = 0;
n_extreme = 0;
power = power(:);
finite = find(isfinite(distance_um) & isfinite(power));
if numel(finite) < 8
    return
end
[~, order] = sort(distance_um(finite));
idx = finite(order);
y = power(idx);
n_edge = edge_count(numel(y));
bad = false(numel(y), 1);
trend = y;
for pass_n = 1:4
    trend = trend_holding_out(y, bad);
    hold_edge = bad;
    if n_edge >= 1
        hold_edge(1:n_edge) = true;
        hold_edge(end - n_edge + 1:end) = true;
        trend_inside = trend_holding_out(y, hold_edge);
        trend(1:n_edge) = trend_inside(1:n_edge);
        trend(end - n_edge + 1:end) = trend_inside(end - n_edge + 1:end);
    end
    sigma_row = ~hold_edge & isfinite(y) & isfinite(trend);
    if nnz(sigma_row) < 4
        sigma_row = isfinite(y) & isfinite(trend);
    end
    sigma = robust_sigma(y(sigma_row) - trend(sigma_row));
    data_sigma = robust_sigma(y(sigma_row));
    if isfinite(data_sigma)
        sigma = max(sigma, 0.05 * data_sigma);
    end
    if ~isfinite(sigma) || sigma == 0
        break
    end
    newly = isfinite(y) & isfinite(trend) & abs(y - trend) > n_fence * sigma;
    if ~any(newly & ~bad)
        bad = newly;
        break
    end
    bad = bad | newly;
end
n_removed = nnz(bad);
if n_removed == 0
    return
end
edge = false(numel(y), 1);
if n_edge >= 1
    edge(1:n_edge) = true;
    edge(end - n_edge + 1:end) = true;
end
power(idx(bad)) = NaN;
n_extreme = nnz(bad & edge);
end

function n_edge = edge_count(n_sites)
n_edge = max(3, ceil(0.05 * n_sites));
n_edge = min(n_edge, floor((n_sites - 4) / 2));
n_edge = max(n_edge, 0);
end

function trend = trend_holding_out(y, hold_out)
y_fit = y;
y_fit(hold_out) = NaN;
trend = local_robust_trend(y_fit);
if all(isfinite(trend))
    return
end
known = find(isfinite(trend));
if numel(known) < 2
    return
end
trend = interp1(known, trend(known), (1:numel(y)).', 'linear', 'extrap');
end

function trend = local_robust_trend(y)
y = y(:);
finite_n = nnz(isfinite(y));
trend = y;
if finite_n < 3
    return
end
window = max(5, 2 * floor(finite_n / 8) + 1);
window = min(window, finite_n);
if mod(window, 2) == 0
    window = max(window - 1, 3);
end
try
    trend = smoothdata(y, 1, 'rlowess', window, 'omitnan');
catch
    trend = y;
end
end

function sigma = robust_sigma(residual)
residual = residual(isfinite(residual));
if numel(residual) < 4
    sigma = NaN;
    return
end
center = median(residual);
sigma = 1.4826 * median(abs(residual - center));
end

function color_limits = spectrum_color_limits(values, spectrum_field)
if isempty(values)
    color_limits = [NaN NaN];
    return
end
if strcmp(spectrum_field, 'raw_aperiodic_subtracted')
    limit = prctile(abs(values), 98);
    if ~isfinite(limit) || limit == 0
        limit = max(abs(values));
    end
    if ~isfinite(limit) || limit == 0
        limit = 1;
    end
    color_limits = [-limit limit];
else
    color_limits = prctile(values, [2 98]);
    if color_limits(1) == color_limits(2)
        color_limits = color_limits + [-1 1];
    end
end
end

function colors = diverging_or_sequential(spectrum_field)
if ~strcmp(spectrum_field, 'raw_aperiodic_subtracted')
    colors = parula(256);
    return
end
n_colors = 256;
half = floor(n_colors / 2);
blue = [0.12 0.27 0.63];
white = [1 1 1];
red = [0.70 0.13 0.13];
low = [linspace(blue(1), white(1), half)', ...
    linspace(blue(2), white(2), half)', ...
    linspace(blue(3), white(3), half)'];
n_high = n_colors - half;
high = [linspace(white(1), red(1), n_high)', ...
    linspace(white(2), red(2), n_high)', ...
    linspace(white(3), red(3), n_high)'];
colors = [low; high];
end

function title_text = spectrum_title(spectrum_field)
if strcmp(spectrum_field, 'raw_aperiodic_subtracted')
    title_text = 'log10 power minus aperiodic fit';
else
    title_text = 'Session-mean log10 power';
end
end

function label = profile_xlabel(spectrum_field)
if strcmp(spectrum_field, 'raw_aperiodic_subtracted')
    label = 'Mean log10 power above aperiodic fit';
else
    label = 'Mean log10 power';
end
end

function depth_ticks(depth_um)
n_channels = numel(depth_um);
n_ticks = min(8, n_channels);
if n_ticks < 1
    return
end
tick_index = unique(round(linspace(1, n_channels, n_ticks)));
tick_labels = strings(numel(tick_index), 1);
for tick_n = 1:numel(tick_index)
    this_depth = depth_um(tick_index(tick_n));
    if isfinite(this_depth)
        tick_labels(tick_n) = sprintf('%.0f', this_depth);
    else
        tick_labels(tick_n) = '';
    end
end
    set(gca, 'YTick', tick_index, 'YTickLabel', cellstr(tick_labels))
end

function mark_depth_lines(y_values)
y_values = y_values(isfinite(y_values));
for line_n = 1:numel(y_values)
    yline(y_values(line_n), '--k', 'LPAG center', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 8, ...
        'HandleVisibility', 'off');
end
end

function borders = area_borders_for_recording(probe_spectra, area_limit_file)
borders = repmat(struct( ...
    'probe_number', [], 'name', '', 'depth_start', [], ...
    'depth_end', [], 'color', []), 0, 1);
if ~isfield(probe_spectra, 'animal_area_table_name') || ...
        ~isfield(probe_spectra, 'npx_type')
    fprintf('Saved spectra have no area-table name. Area borders are omitted.\n');
    return
end
if ~isfile(area_limit_file)
    fprintf('Area-limit table not found. Area borders are omitted.\n%s\n', ...
        area_limit_file);
    return
end

area_limit = readtable(area_limit_file);
area_limit = area_limit(ismember(area_limit.AnimalName, ...
    probe_spectra.animal_area_table_name), :);
if isempty(area_limit)
    fprintf('No area limits for %s.\n', probe_spectra.animal_area_table_name);
    return
end

names = cellstr(string(area_limit.area));
depth_start = double(area_limit.depth_start);
depth_end = double(area_limit.depth_end);
probe_number = double(area_limit.ProbeNum);
keep = isfinite(depth_start) & isfinite(depth_end) & ~strcmp(names, '');
if probe_spectra.npx_type == 1
    if ~ismember('Probe_Area', area_limit.Properties.VariableNames)
        fprintf('Area table has no Probe_Area column. NPX1 borders are omitted.\n');
        return
    end
    keep = keep & ismember(area_limit.Probe_Area, 'PAG');
    probe_number(:) = 1;
else
    keep = keep & isfinite(probe_number);
end
names = names(keep);
depth_start = depth_start(keep);
depth_end = depth_end(keep);
probe_number = probe_number(keep);
if isempty(names)
    fprintf('No area borders for %s.\n', probe_spectra.animal_area_table_name);
    return
end

palette = lines(numel(unique(names, 'stable')));
name_list = unique(names, 'stable');
borders = repmat(struct( ...
    'probe_number', [], 'name', '', 'depth_start', [], ...
    'depth_end', [], 'color', []), numel(names), 1);
for border_n = 1:numel(names)
    color_index = find(strcmp(name_list, names{border_n}), 1);
    borders(border_n).probe_number = probe_number(border_n);
    borders(border_n).name = names{border_n};
    borders(border_n).depth_start = depth_start(border_n);
    borders(border_n).depth_end = depth_end(border_n);
    borders(border_n).color = palette(color_index, :);
end
fprintf('Area borders for %s: %d\n', ...
    probe_spectra.animal_area_table_name, numel(borders));
end

function these = borders_for_probe(area_borders, probe_number)
if isempty(area_borders) || ~isfinite(probe_number)
    these = area_borders([]);
    return
end
these = area_borders([area_borders.probe_number] == probe_number);
end

function limits = depth_limits_with_borders(depth_um, borders)
finite_depth = depth_um(isfinite(depth_um));
limits = [min(finite_depth), max(finite_depth)];
if isempty(borders)
    return
end
edges = [[borders.depth_start].'; [borders.depth_end].'];
edges = edges(isfinite(edges));
near = edges >= limits(1) - 50 & edges <= limits(2) + 50;
if any(near)
    limits(1) = min([limits(1); edges(near)]);
    limits(2) = max([limits(2); edges(near)]);
end
end

function draw_area_borders(borders, depth_um, use_depth_axis)
if isempty(borders)
    return
end
hold on
x_limits = xlim;
for border_n = 1:numel(borders)
    edges = [borders(border_n).depth_start, borders(border_n).depth_end];
    [~, shallow] = min(edges);
    for edge_n = 1:2
        if use_depth_axis
            y_value = edges(edge_n);
            if ~isfinite(y_value)
                continue
            end
            y_limits = ylim;
            if y_value < y_limits(1) || y_value > y_limits(2)
                continue
            end
        else
            y_value = depth_to_row(depth_um, edges(edge_n));
            if ~isfinite(y_value)
                continue
            end
        end
        yline(y_value, ':', 'Color', borders(border_n).color, ...
            'LineWidth', 1, 'HandleVisibility', 'off');
        if edge_n == shallow
            text(x_limits(1), y_value, ['  ' borders(border_n).name], ...
                'Color', borders(border_n).color, 'FontSize', 7, ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
                'Interpreter', 'none', 'Clipping', 'on');
        end
    end
end
end

function row = depth_to_row(depth_um, border_depth)
row = NaN;
if ~isfinite(border_depth)
    return
end
finite = isfinite(depth_um(:));
if nnz(finite) < 2
    return
end
depth_values = depth_um(finite);
row_index = find(finite);
if border_depth < depth_values(1) || border_depth > depth_values(end)
    return
end
upper = find(depth_values >= border_depth, 1, 'first');
if upper == 1
    row = row_index(1);
    return
end
span = depth_values(upper) - depth_values(upper - 1);
if span == 0
    row = row_index(upper);
    return
end
fraction = (border_depth - depth_values(upper - 1)) / span;
row = row_index(upper - 1) + fraction * (row_index(upper) - row_index(upper - 1));
end

function item = empty_placed()
item = struct( ...
    'label', '', ...
    'animal', '', ...
    'session', NaN, ...
    'shank', NaN, ...
    'channel', [], ...
    'tip_um', [], ...
    'areas', {{}}, ...
    'welch_bands', [], ...
    'fourier_bands', []);
end

function placed = placed_from_session(channel_table, welch_spectrum, ...
        fourier_spectrum, frequency, delta_band, theta_band, gamma_band, ...
        session_label, artifact_trend_sd, animal_name, session_n)
groups = probe_groups_by_depth(channel_table);
placed = repmat(empty_placed(), 0, 1);
for group_n = 1:numel(groups)
    this_group = groups(group_n);
    tip_um = this_group.depth_um(:);
    if ~any(isfinite(tip_um))
        continue
    end
    welch_rows = welch_spectrum(this_group.channel_index, :);
    fourier_rows = fourier_spectrum(this_group.channel_index, :);
    item = empty_placed();
    item.label = sprintf('%s | %s', session_label, this_group.label);
    item.animal = char(string(animal_name));
    item.session = session_n;
    item.shank = this_group.probe_number;
    item.channel = channel_table.Channel(this_group.channel_index);
    item.channel = item.channel(:);
    item.tip_um = tip_um;
    item.areas = this_group.areas;
    item.welch_bands = reject_trend_outliers(tip_um, [ ...
        mean_in_band(welch_rows, frequency, delta_band), ...
        mean_in_band(welch_rows, frequency, theta_band), ...
        mean_in_band(welch_rows, frequency, gamma_band)], artifact_trend_sd);
    item.fourier_bands = reject_trend_outliers(tip_um, [ ...
        mean_in_band(fourier_rows, frequency, delta_band), ...
        mean_in_band(fourier_rows, frequency, theta_band), ...
        mean_in_band(fourier_rows, frequency, gamma_band)], artifact_trend_sd);
    placed(end + 1, 1) = item; %#ok<AGROW>
end
end

function placed = apply_screened_power_to_placed(placed, power_table)
if isempty(placed) || isempty(power_table) || ...
        ~ismember('DeltaPower', power_table.Properties.VariableNames)
    return
end
band_fields = {'DeltaPower', 'ThetaPower', 'GammaPower'};
nan_bands = containers.Map('KeyType', 'char', 'ValueType', 'any');
row_key = site_key(power_table.Animal, power_table.Session, ...
    power_table.Shank, power_table.Channel, power_table.Method);
for row_n = 1:height(power_table)
    mask = false(1, 3);
    for band_n = 1:3
        mask(band_n) = ~isfinite(power_table.(band_fields{band_n})(row_n));
    end
    key = row_key{row_n};
    if isKey(nan_bands, key)
        nan_bands(key) = nan_bands(key) | mask;
    else
        nan_bands(key) = mask;
    end
end
n_blanked = 0;
method_names = {'Welch', 'Fourier'};
for placed_n = 1:numel(placed)
    item = placed(placed_n);
    for site_n = 1:numel(item.channel)
        for method_n = 1:2
            key = site_key(item.animal, item.session, item.shank, ...
                item.channel(site_n), method_names{method_n});
            key = key{1};
            if ~isKey(nan_bands, key)
                continue
            end
            mask = nan_bands(key);
            if method_n == 1
                bands = item.welch_bands;
            else
                bands = item.fourier_bands;
            end
            for band_n = 1:3
                if mask(band_n) && isfinite(bands(site_n, band_n))
                    bands(site_n, band_n) = NaN;
                    n_blanked = n_blanked + 1;
                end
            end
            if method_n == 1
                item.welch_bands = bands;
            else
                item.fourier_bands = bands;
            end
        end
    end
    placed(placed_n) = item;
end
if n_blanked > 0
    fprintf(['Region-distance plots leave out %d further band values that ' ...
        'the distance table left blank. Those channels failed the shank-trend ' ...
        'screen or sit beyond the SupCol or VLPAG border distance reached by ' ...
        'most sessions. The session itself stays.\n'], n_blanked);
end
end

function key = site_key(animal, session, shank, channel, method)
key = strcat(as_key_text(animal), '|', ...
    as_key_text(session), '|', ...
    as_key_text(shank), '|', ...
    as_key_text(channel), '|', ...
    as_key_text(method));
end

function text = as_key_text(value)
if ischar(value)
    text = {value};
    return
end
if isnumeric(value)
    value = value(:);
    text = cell(numel(value), 1);
    for value_n = 1:numel(value)
        if isfinite(value(value_n))
            text{value_n} = sprintf('%.12g', value(value_n));
        else
            text{value_n} = 'NaN';
        end
    end
    return
end
text = cellstr(string(value(:)));
end

function power_table = lpag_distance_table(channel_table, welch_spectrum, ...
        fourier_spectrum, frequency, delta_band, theta_band, gamma_band, ...
        animal_name, session_n, spectrum_field, artifact_trend_sd)
groups = probe_groups_by_depth(channel_table);
method_names = {'Welch', 'Fourier'};
spectra = {welch_spectrum, fourier_spectrum};
variable_names = { ...
    'DistanceToMidLPAG_um', ...
    'DeltaPower', 'ThetaPower', 'GammaPower', 'Area', ...
    'DeltaPower_centeredAtLPAG', 'ThetaPower_centeredAtLPAG', ...
    'GammaPower_centeredAtLPAG', ...
    'DeltaPower_zWithinShank', 'ThetaPower_zWithinShank', ...
    'GammaPower_zWithinShank', ...
    'Animal', 'Session', 'Shank', 'Channel', 'Method', 'Spectrum'};
pieces = cell(1, 0);
n_removed = 0;
n_extreme = 0;
for group_n = 1:numel(groups)
    this_group = groups(group_n);
    center_mask = this_group.is_lpag_center(:);
    center_depth = this_group.depth_um(center_mask);
    center_depth = center_depth(isfinite(center_depth));
    if isempty(center_depth)
        continue
    end
    distance_um = this_group.depth_um(:) - mean(center_depth);
    channel = channel_table.Channel(this_group.channel_index);
    n_sites = numel(distance_um);
    for method_n = 1:2
        rows = spectra{method_n}(this_group.channel_index, :);
        bands = [ ...
            mean_in_band(rows, frequency, delta_band), ...
            mean_in_band(rows, frequency, theta_band), ...
            mean_in_band(rows, frequency, gamma_band)];
        [bands, removed_n, extreme_n] = reject_trend_outliers( ...
            distance_um, bands, artifact_trend_sd);
        n_removed = n_removed + removed_n;
        n_extreme = n_extreme + extreme_n;
        reference = mean(bands(center_mask, :), 1, 'omitnan');
        centered = bands - reference;
        band_std = std(bands, 0, 1, 'omitnan');
        band_std(band_std == 0) = NaN;
        z_power = centered ./ band_std;
        pieces{end + 1} = table( ...
            distance_um, ...
            bands(:, 1), bands(:, 2), bands(:, 3), this_group.areas(:), ...
            centered(:, 1), centered(:, 2), centered(:, 3), ...
            z_power(:, 1), z_power(:, 2), z_power(:, 3), ...
            repmat({animal_name}, n_sites, 1), ...
            repmat(session_n, n_sites, 1), ...
            repmat(this_group.probe_number, n_sites, 1), ...
            channel, ...
            repmat(method_names(method_n), n_sites, 1), ...
            repmat({spectrum_field}, n_sites, 1), ...
            'VariableNames', variable_names); %#ok<AGROW>
    end
end
if n_removed > 0
    fprintf(['  Removed %d band values more than %g robust SD from the ' ...
        'shank trend, %d of them at the probe ends.\n'], ...
        n_removed, artifact_trend_sd, n_extreme);
end
if isempty(pieces)
    power_table = table();
    return
end
power_table = vertcat(pieces{:});
end

function power_table = write_lpag_distance_table(power_pieces, ...
        saving_folder, spectrum_field, region_border_bin_um)
kept = power_pieces(~cellfun(@isempty, power_pieces));
if isempty(kept)
    fprintf(['No shank with a mid-LPAG channel was found; ' ...
        'the distance table was not written.\n']);
    power_table = table();
    return
end
power_table = vertcat(kept{:});
power_table = sortrows(power_table, ...
    {'Animal', 'Session', 'Shank', 'Method', 'DistanceToMidLPAG_um'});
power_table = restrict_region_border_distances(power_table, region_border_bin_um);
table_file = fullfile(saving_folder, sprintf( ...
    'play_session_probe_power_by_lpag_distance_%s.csv', spectrum_field));
writetable(power_table, table_file);
fprintf(['Wrote %d channel rows to %s\n' ...
    'Distance is positive dorsal (farther from the tip). ' ...
    'centeredAtLPAG subtracts the mid-LPAG channel. ' ...
    'zWithinShank divides that difference by the standard deviation ' ...
    'of the same band on that shank. Band values left blank failed the ' ...
    'shank-trend screen or lie beyond the SupCol or VLPAG border distance ' ...
    'reached by most sessions.\n'], height(power_table), table_file);
end

function power_table = restrict_region_border_distances(power_table, bin_um)
if isempty(power_table) || ~isfinite(bin_um) || bin_um <= 0
    return
end
region_keys = region_keys_for_rows(power_table.Area);
distance_um = power_table.DistanceToMidLPAG_um(:);
border_um = distance_to_area_border(power_table, region_keys, distance_um);
power_table.DistanceToAreaBorder_um = border_um;
screen_regions = {'VLPAG', 'SC'};
drop = false(height(power_table), 1);
for region_n = 1:numel(screen_regions)
    drop = drop | uncommon_border_rows(power_table, region_keys, ...
        screen_regions{region_n}, bin_um);
end
if ~any(drop)
    fprintf(['Every VLPAG and SupCol channel falls within a border distance ' ...
        'shared by most sessions.\n']);
    return
end
power_fields = { ...
    'DeltaPower', 'ThetaPower', 'GammaPower', ...
    'DeltaPower_centeredAtLPAG', 'ThetaPower_centeredAtLPAG', ...
    'GammaPower_centeredAtLPAG', ...
    'DeltaPower_zWithinShank', 'ThetaPower_zWithinShank', ...
    'GammaPower_zWithinShank'};
for field_n = 1:numel(power_fields)
    power_table.(power_fields{field_n})(drop) = NaN;
end
end

function border_um = distance_to_area_border(power_table, region_keys, distance_um)
border_um = nan(height(power_table), 1);
regions = {'VLPAG', 'SC', 'LPAG'};
for region_n = 1:numel(regions)
    this_border = border_for_region(power_table, region_keys, ...
        regions{region_n}, distance_um);
    unset = isnan(border_um) & isfinite(this_border);
    border_um(unset) = this_border(unset);
end
end

function border_um = border_for_region(power_table, region_keys, ...
        region_name, distance_um)
border_um = nan(height(power_table), 1);
group_key = strcat(cellstr(string(power_table.Animal)), '|', ...
    cellstr(string(power_table.Session)), '|', ...
    cellstr(string(power_table.Shank)), '|', ...
    cellstr(string(power_table.Method)));
[~, ~, group_id] = unique(group_key, 'stable');
for group_n = 1:max(group_id)
    rows = find(group_id == group_n);
    in_region = region_row_has(region_keys(rows), region_name);
    these = rows(in_region);
    these = these(isfinite(distance_um(these)));
    if isempty(these)
        continue
    end
    if strcmp(region_name, 'LPAG')
        border_um(these) = abs(distance_um(these));
    else
        border_um(these) = depth_from_lpag_border(distance_um(these));
    end
end
end

function depth = depth_from_lpag_border(distance_um)
distance_um = distance_um(:);
if median(distance_um, 'omitnan') >= 0
    border = min(distance_um);
    depth = distance_um - border;
else
    border = max(distance_um);
    depth = border - distance_um;
end
end

function drop = uncommon_border_rows(power_table, region_keys, ...
        region_name, bin_um)
drop = false(height(power_table), 1);
distance_um = power_table.DistanceToMidLPAG_um(:);
border_um = border_for_region(power_table, region_keys, region_name, distance_um);
methods = unique(cellstr(string(power_table.Method)), 'stable');
if isempty(methods)
    return
end
in_region = region_row_has(region_keys, region_name);
in_count = in_region & strcmp(power_table.Method, methods{1}) & ...
    isfinite(border_um);
if ~any(in_count)
    fprintf('No %s channels with a border distance.\n', ...
        region_display_name(region_name));
    return
end
session_key = strcat(cellstr(string(power_table.Animal)), '|', ...
    cellstr(string(power_table.Session)));
n_sessions = numel(unique(session_key(in_count)));
bin_id = floor(border_um(:) / bin_um);
bins = sort(unique(bin_id(in_count)));
n_at = zeros(numel(bins), 1);
for bin_n = 1:numel(bins)
    in_bin = in_count & bin_id == bins(bin_n);
    n_at(bin_n) = numel(unique(session_key(in_bin)));
end
need = floor(n_sessions / 2) + 1;
kept_bin = n_at >= need;
if ~any(kept_bin)
    kept_bin = n_at == max(n_at);
    fprintf(['No %s border-distance bin reached %d sessions. ' ...
        'Keeping bins with the largest coverage (%d sessions).\n'], ...
        region_display_name(region_name), need, max(n_at));
end
fprintf(['%s distance to its LPAG-facing border. ' ...
    '%d sessions have %s. Bins are %g µm. ' ...
    'A bin stays when at least %d sessions have a channel there.\n'], ...
    region_display_name(region_name), n_sessions, ...
    region_display_name(region_name), bin_um, need);
for bin_n = 1:numel(bins)
    bin_start = bins(bin_n) * bin_um;
    if kept_bin(bin_n)
        status = 'kept';
    else
        status = 'set aside';
    end
    fprintf('  %8g to %8g µm: %3d sessions  %s\n', ...
        bin_start, bin_start + bin_um, n_at(bin_n), status);
end
drop = in_region & isfinite(border_um) & ~ismember(bin_id, bins(kept_bin));
if ~any(drop)
    fprintf('Every %s channel is within a border distance shared by most sessions.\n', ...
        region_display_name(region_name));
    return
end
kept_ids = bins(kept_bin);
fprintf(['Set aside %d %s rows beyond %g µm from the border, ' ...
    'a distance reached by only a minority of sessions.\n'], ...
    nnz(drop), region_display_name(region_name), ...
    (max(kept_ids) + 1) * bin_um);
end

function plot_normalized_power_scatters(power_table, ...
        delta_band, theta_band, gamma_band, spectral_method, ...
        regions_to_display, bands_to_display, plot_style)
if isempty(power_table)
    return
end
[method_names, ~] = selected_spectral_methods(spectral_method);
selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band);
region_keys = region_keys_for_rows(power_table.Area);
region_names = selected_regions(regions_to_display);
figure_n = 0;
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    distance_um = power_table.DistanceToMidLPAG_um;
    for slot_n = 1:numel(selection.index)
        figure_n = figure_n + 1;
        power = power_table.(selection.fields{slot_n});
        figure('Name', sprintf('%s | %s vs distance to mid-LPAG', ...
            method_names{method_n}, selection.labels{slot_n}), ...
            'Color', 'w', 'Position', [60 + 24 * figure_n, 60, 760, 560]);
        hold on
        yline(0, ':k', 'HandleVisibility', 'off');
        xline(0, '--k', 'mid-LPAG', 'HandleVisibility', 'off', ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 8);
        legend_handles = gobjects(0, 1);
        legend_labels = {};
        for region_n = 1:numel(region_names)
            this_region = region_names{region_n};
            in_region = in_method & region_row_has(region_keys, this_region);
            keep = in_region & isfinite(distance_um) & isfinite(power);
            if ~any(keep)
                continue
            end
            this_color = region_color(this_region);
            scatter(distance_um(keep), power(keep), plot_style.marker_size, ...
                this_color, 'filled', ...
                'MarkerFaceAlpha', plot_style.marker_alpha, ...
                'MarkerEdgeAlpha', plot_style.marker_alpha, ...
                'HandleVisibility', 'off');
            legend_handles(end + 1, 1) = plot(nan, nan, 'o', ...
                'MarkerFaceColor', this_color, 'MarkerEdgeColor', this_color, ...
                'MarkerSize', 7); %#ok<AGROW>
            legend_labels{end + 1} = sprintf('%s (n=%d)', ...
                region_display_name(this_region), nnz(keep)); %#ok<AGROW>
        end
        xlabel({'Distance to mid-LPAG (µm)', ...
            'positive = dorsal, farther from the tip'})
        ylabel({'Power z-scored within shank', '0 = mid-LPAG channel'})
        title(sprintf('%s | %s', method_names{method_n}, selection.labels{slot_n}))
        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, ...
                'Location', 'eastoutside', 'FontSize', 8)
        end
        box off
    end
end
end

function keys = region_keys_for_rows(area_labels)
area_labels = cellstr(string(area_labels));
keys = cell(numel(area_labels), 1);
for row_n = 1:numel(area_labels)
    tokens = split_area_label(area_labels{row_n});
    if isempty(tokens)
        keys{row_n} = {'Other'};
        continue
    end
    names = cell(1, numel(tokens));
    for token_n = 1:numel(tokens)
        names{token_n} = canonical_region_name(tokens{token_n});
    end
    names = unique(names, 'stable');
    colored = names(~strcmp(names, 'Other'));
    if isempty(colored)
        keys{row_n} = {'Other'};
    else
        keys{row_n} = colored(:).';
    end
end
end

function [method_names, method_ids] = selected_spectral_methods(spectral_method)
if ischar(spectral_method) || isstring(spectral_method)
    requested = cellstr(spectral_method);
else
    requested = cellstr(string(spectral_method));
end
requested = requested(:).';
catalog = {'Welch', 'Fourier'};
method_names = cell(1, numel(requested));
method_ids = zeros(1, numel(requested));
for item_n = 1:numel(requested)
    match = strcmpi(catalog, strtrim(requested{item_n}));
    if ~any(match)
        error('spectral_method must be Welch, Fourier, or both.');
    end
    method_ids(item_n) = find(match, 1);
    method_names{item_n} = catalog{method_ids(item_n)};
end
end

function selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band)
catalog = {'delta', 'theta', 'gamma'};
labels = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
fields = { ...
    'DeltaPower_zWithinShank', ...
    'ThetaPower_zWithinShank', ...
    'GammaPower_zWithinShank'};
colors = [0.12 0.35 0.75; 0.85 0.40 0.10; 0.45 0.62 0.20];
if ischar(bands_to_display) || isstring(bands_to_display)
    requested = cellstr(bands_to_display);
else
    requested = cellstr(string(bands_to_display));
end
requested = requested(:).';
if isempty(requested)
    error('bands_to_display must include delta, theta, gamma, or a combination.');
end
index = zeros(1, numel(requested));
for item_n = 1:numel(requested)
    match = strcmpi(catalog, strtrim(requested{item_n}));
    if ~any(match)
        error('bands_to_display must include delta, theta, gamma, or a combination.');
    end
    index(item_n) = find(match, 1);
end
selection = struct( ...
    'index', index, ...
    'labels', {labels(index)}, ...
    'fields', {fields(index)}, ...
    'colors', colors(index, :));
end

function region_names = selected_regions(regions_to_display)
if ischar(regions_to_display) || isstring(regions_to_display)
    requested = cellstr(regions_to_display);
else
    requested = cellstr(string(regions_to_display));
end
requested = requested(:).';
region_names = cell(1, numel(requested));
for item_n = 1:numel(requested)
    token = strtrim(requested{item_n});
    if strcmpi(token, 'Other')
        region_names{item_n} = 'Other';
        continue
    end
    region_names{item_n} = canonical_region_name(token);
    if strcmp(region_names{item_n}, 'Other')
        error(['Unknown region "%s". Use LPAG, VLPAG, SupCol, DR, ' ...
            'DLPAG, DMPAG, or PAG.'], token);
    end
end
end

function label = region_display_name(name)
[parent, half] = split_region_half(name);
if strcmp(parent, 'SC')
    parent = 'SupCol';
end
if isempty(half)
    label = parent;
else
    label = [parent, ' ', half];
end
end

function [parent, half] = split_region_half(name)
half = '';
parent = name;
suffixes = {' proximal', ' distal'};
for suffix_n = 1:numel(suffixes)
    token = suffixes{suffix_n};
    if numel(name) > numel(token) && ...
            strcmp(name(end - numel(token) + 1:end), token)
        parent = name(1:end - numel(token));
        half = token(2:end);
        return
    end
end
end

function name = canonical_region_name(token)
if name_is_supcol(token)
    name = 'SC';
elseif name_is_dr(token) || name_is_pag(token)
    name = token;
else
    name = 'Other';
end
end

function ordered = region_draw_order(names)
foreground = {'PAG', 'DMPAG', 'DLPAG', 'LPAG', 'VLPAG', 'SC', 'DR'};
background = names(~ismember(names, foreground));
background = background(:).';
present_foreground = foreground(ismember(foreground, names));
ordered = [background, present_foreground];
end

function mask = region_row_has(region_keys, region_name)
mask = false(numel(region_keys), 1);
for row_n = 1:numel(region_keys)
    mask(row_n) = any(strcmp(region_keys{row_n}, region_name));
end
end

function color = region_color(name)
[parent, half] = split_region_half(name);
color = base_region_color(parent);
if strcmp(half, 'proximal')
    color = color * 0.72;
elseif strcmp(half, 'distal')
    color = min(color + 0.28, 1);
end
end

function color = base_region_color(name)
switch name
    case 'SC'
        color = [0.90 0.45 0.05];
    case 'DLPAG'
        color = [0.16 0.42 0.78];
    case 'LPAG'
        color = [0.13 0.55 0.36];
    case 'VLPAG'
        color = [0.55 0.72 0.16];
    case 'DR'
        color = [0.45 0.28 0.72];
    case 'DMPAG'
        color = [0.05 0.55 0.58];
    case 'PAG'
        color = [0.20 0.38 0.30];
    case {'Other', 'Unlabeled'}
        color = [0 0 0];
    otherwise
        color = [0 0 0];
end
end

function plot_normalized_power_intervals(power_table, ...
        delta_band, theta_band, gamma_band, spectral_method, ...
        regions_to_display, distance_model_bin_um, bands_to_display, plot_style)
if isempty(power_table)
    return
end
[method_names, ~] = selected_spectral_methods(spectral_method);
selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band);
region_keys = region_keys_for_rows(power_table.Area);
region_names = selected_regions(regions_to_display);
region_names = region_names(~strcmp(region_names, 'Other'));
animal = cellstr(string(power_table.Animal));
session = power_table.Session;
figure_n = 0;
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    distance_um = power_table.DistanceToMidLPAG_um;
    for slot_n = 1:numel(selection.index)
        figure_n = figure_n + 1;
        power = power_table.(selection.fields{slot_n});
        figure('Name', sprintf( ...
            '%s | %s intercept by distance to mid-LPAG', ...
            method_names{method_n}, selection.labels{slot_n}), ...
            'Color', 'w', 'Position', [80 + 24 * figure_n, 100, 980, 560]);
        hold on
        yline(0, ':k', 'HandleVisibility', 'off');
        xline(0, '--k', 'mid-LPAG', 'HandleVisibility', 'off', ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 8);
        legend_handles = gobjects(0, 1);
        legend_labels = {};
        x_span = [];
        drawn_profiles = cell(0, 1);
        drawn_colors = zeros(0, 3);
        for region_n = 1:numel(region_names)
            this_region = region_names{region_n};
            keep = in_method & region_row_has(region_keys, this_region) & ...
                isfinite(distance_um) & isfinite(power);
            if ~any(keep)
                continue
            end
            this_distance = distance_um(keep);
            this_power = power(keep);
            this_color = region_color(this_region);
            scatter(this_distance, this_power, plot_style.marker_size, ...
                this_color, 'filled', ...
                'MarkerFaceAlpha', plot_style.marker_alpha, ...
                'MarkerEdgeAlpha', plot_style.marker_alpha, ...
                'HandleVisibility', 'off');
            profile = fit_intercepts_by_distance(this_distance, this_power, ...
                animal(keep), session(keep), distance_model_bin_um);
            [rho, p_rho] = rank_correlation(this_distance, this_power);
            report_intercept_profile(sprintf('%s | %s | %s', ...
                method_names{method_n}, selection.labels{slot_n}, ...
                region_display_name(this_region)), profile, rho, p_rho);
            draw_distance_intercepts(profile, this_color, plot_style.ci_alpha);
            x_span = [x_span; this_distance(:)]; %#ok<AGROW>
            drawn_profiles{end + 1, 1} = profile; %#ok<AGROW>
            drawn_colors(end + 1, :) = this_color; %#ok<AGROW>
            legend_handles(end + 1, 1) = plot(nan, nan, '-o', ...
                'Color', this_color, 'MarkerFaceColor', this_color, ...
                'LineWidth', 2, 'MarkerSize', 4); %#ok<AGROW>
            legend_labels{end + 1} = intercept_legend( ...
                this_region, profile, rho, p_rho); %#ok<AGROW>
        end
        if numel(x_span) > 1
            pad = 0.02 * max(max(x_span) - min(x_span), 1);
            xlim([min(x_span) - pad, max(x_span) + pad]);
        end
        y_limits = ylim;
        for drawn_n = 1:numel(drawn_profiles)
            draw_distance_significance(drawn_profiles{drawn_n}, ...
                drawn_colors(drawn_n, :), plot_style.significance_alpha, ...
                plot_style.significance_min_um, y_limits);
        end
        ylim(y_limits);
        xlabel({'Distance to mid-LPAG (µm)', ...
            'positive = dorsal, farther from the tip'})
        ylabel({'Power z-scored within shank', '0 = mid-LPAG channel'})
        title(sprintf('%s | %s | intercept at each %g µm', ...
            method_names{method_n}, selection.labels{slot_n}, distance_model_bin_um))
        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, ...
                'Location', 'eastoutside', 'FontSize', 8, 'Interpreter', 'tex')
        end
        box off
    end
end
fprintf(['At each distance bin, one session power is the mean of that ' ...
    'session''s channels. Power ~ 1 + (1|Animal). The marker is the ' ...
    'intercept and the shade is its 95%% interval. A rectangle covers the ' ...
    'upper 90%% of the axis where that interval excludes 0 for at least ' ...
    '%g µm. Spearman rho and p use the channel scatter, after extreme ' ...
    'values and uncommon border distances are removed.\n'], ...
    plot_style.significance_min_um);
end

function profile = empty_distance_profile(bin_um)
if nargin < 1
    bin_um = NaN;
end
profile = struct( ...
    'bin_um', bin_um, ...
    'distance_um', [], ...
    'intercept', [], ...
    'lower', [], ...
    'upper', [], ...
    'n_animals', [], ...
    'n_sessions', [], ...
    'n_fit', 0, ...
    'n_fallback', 0, ...
    'n_skipped', 0);
end

function profile = fit_intercepts_by_distance( ...
        distance, power, animal, session, bin_um)
profile = empty_distance_profile(bin_um);
distance = distance(:);
power = power(:);
animal = cellstr(string(animal(:)));
session = cellstr(string(session(:)));
keep = isfinite(distance) & isfinite(power);
distance = distance(keep);
power = power(keep);
animal = animal(keep);
session = session(keep);
if isempty(power) || ~isfinite(bin_um) || bin_um <= 0
    return
end
bin_id = round(distance / bin_um);
bins = unique(bin_id);
bins = sort(bins);
for bin_n = 1:numel(bins)
    in_bin = bin_id == bins(bin_n);
    [session_power, session_animal] = session_mean_power( ...
        power(in_bin), animal(in_bin), session(in_bin));
    n_animals = numel(unique(session_animal));
    n_sessions = numel(session_power);
    if n_animals < 2 || n_sessions < 2
        profile.n_skipped = profile.n_skipped + 1;
        continue
    end
    [intercept, lower, upper, used_fallback] = intercept_mixed_model( ...
        session_power, session_animal);
    if ~isfinite(intercept)
        profile.n_skipped = profile.n_skipped + 1;
        continue
    end
    profile.distance_um(end + 1, 1) = bins(bin_n) * bin_um; %#ok<AGROW>
    profile.intercept(end + 1, 1) = intercept; %#ok<AGROW>
    profile.lower(end + 1, 1) = lower; %#ok<AGROW>
    profile.upper(end + 1, 1) = upper; %#ok<AGROW>
    profile.n_animals(end + 1, 1) = n_animals; %#ok<AGROW>
    profile.n_sessions(end + 1, 1) = n_sessions; %#ok<AGROW>
    profile.n_fit = profile.n_fit + 1;
    if used_fallback
        profile.n_fallback = profile.n_fallback + 1;
    end
end
end

function [session_power, session_animal] = session_mean_power( ...
        power, animal, session)
group_key = strcat(animal(:), '|', session(:));
[~, ~, group_id] = unique(group_key, 'stable');
n_groups = max(group_id);
session_power = accumarray(group_id, power(:), [n_groups, 1], @mean);
session_animal = cell(n_groups, 1);
for group_n = 1:n_groups
    session_animal{group_n} = animal{find(group_id == group_n, 1, 'first')};
end
end

function [intercept, lower, upper, used_fallback] = intercept_mixed_model( ...
        power, animal)
intercept = NaN;
lower = NaN;
upper = NaN;
used_fallback = false;
power = power(:);
animal = categorical(animal(:));
try
    model_table = table(power, animal, 'VariableNames', {'Power', 'Animal'});
    lme = fitlme(model_table, 'Power ~ 1 + (1|Animal)', 'FitMethod', 'REML');
    [beta, ~, stats] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    intercept = beta(1);
    try
        lower = stats.Lower(1);
        upper = stats.Upper(1);
    catch
        [lower, upper] = interval_from_se(intercept, stats.SE(1), stats.DF(1));
    end
catch
    try
        model_table = table(power, animal, 'VariableNames', {'Power', 'Animal'});
        lme = fitlme(model_table, 'Power ~ 1 + (1|Animal)', 'FitMethod', 'REML');
        beta = fixedEffects(lme);
        intercept = beta(1);
        se = sqrt(max(lme.CoefficientCovariance(1, 1), 0));
        [lower, upper] = interval_from_se(intercept, se, max(numel(unique(animal)) - 1, 1));
    catch
        intercept = NaN;
    end
end
if isfinite(intercept) && isfinite(lower) && isfinite(upper)
    return
end
used_fallback = true;
[animal_id, ~, group_id] = unique(animal, 'stable');
animal_mean = accumarray(group_id, power, [numel(animal_id), 1], @mean);
if numel(animal_mean) < 2
    intercept = NaN;
    lower = NaN;
    upper = NaN;
    return
end
intercept = mean(animal_mean);
se = std(animal_mean, 0) / sqrt(numel(animal_mean));
[lower, upper] = interval_from_se(intercept, se, numel(animal_mean) - 1);
end

function [lower, upper] = interval_from_se(intercept, se, df)
if ~isfinite(intercept) || ~isfinite(se) || ~isfinite(df) || df <= 0
    lower = NaN;
    upper = NaN;
    return
end
half = tinv(0.975, df) * se;
lower = intercept - half;
upper = intercept + half;
end

function draw_distance_intercepts(profile, color, ci_alpha)
if isempty(profile.distance_um)
    return
end
[x_um, order] = sort(profile.distance_um(:));
y = profile.intercept(order);
y_low = profile.lower(order);
y_high = profile.upper(order);
adjacent = [true; diff(x_um) <= profile.bin_um * 1.5];
run_start = find([true; ~adjacent(2:end)]);
run_stop = [run_start(2:end) - 1; numel(x_um)];
for run_n = 1:numel(run_start)
    index = run_start(run_n):run_stop(run_n);
    if numel(index) >= 2 && all(isfinite(y_low(index))) && all(isfinite(y_high(index)))
        segment_x = x_um(index);
        fill([segment_x; flipud(segment_x)], ...
            [y_low(index); flipud(y_high(index))], color, ...
            'FaceAlpha', ci_alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(segment_x, y(index), '-', 'Color', color, 'LineWidth', 1.6, ...
            'HandleVisibility', 'off');
    end
    plot(x_um(index), y(index), 'o', 'Color', color, ...
        'MarkerFaceColor', color, 'MarkerSize', 4, 'HandleVisibility', 'off');
end
end

function draw_distance_significance(profile, color, face_alpha, min_um, y_limits)
if isempty(profile.distance_um) || ~isfinite(min_um) || numel(y_limits) < 2
    return
end
span = y_limits(2) - y_limits(1);
if ~isfinite(span) || span <= 0
    return
end
[x_um, order] = sort(profile.distance_um(:));
lower = profile.lower(order);
upper = profile.upper(order);
significant = isfinite(lower) & isfinite(upper) & (upper < 0 | lower > 0);
y0 = y_limits(1) + 0.10 * span;
y1 = y_limits(2);
bin_um = profile.bin_um;
run_start = 0;
for point_n = 1:numel(x_um)
    continues = significant(point_n) && run_start > 0 && ...
        (x_um(point_n) - x_um(point_n - 1) <= bin_um * 1.5);
    if continues
        continue
    end
    if run_start > 0
        mark_significant_span(x_um(run_start), x_um(point_n - 1), ...
            bin_um, min_um, y0, y1, color, face_alpha);
        run_start = 0;
    end
    if significant(point_n)
        run_start = point_n;
    end
end
if run_start > 0
    mark_significant_span(x_um(run_start), x_um(end), ...
        bin_um, min_um, y0, y1, color, face_alpha);
end
end

function mark_significant_span(x_first, x_last, bin_um, min_um, y0, y1, color, face_alpha)
width = x_last - x_first + bin_um;
if width < min_um
    return
end
x0 = x_first - bin_um / 2;
x1 = x_last + bin_um / 2;
rectangle_handle = fill([x0, x1, x1, x0], [y0, y0, y1, y1], color, ...
    'FaceAlpha', face_alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack(rectangle_handle, 'bottom');
end

function label = intercept_legend(region_name, profile, rho, p_rho)
if profile.n_fit > 0
    fit_text = sprintf('intercept, %d distances', profile.n_fit);
else
    fit_text = 'no distance with 2 animals';
end
if isfinite(rho)
    rho_text = sprintf('\\rho = %.2f, p = %s', rho, format_p_value(p_rho));
else
    rho_text = '\rho = n/a';
end
label = sprintf('%s | %s | %s', region_display_name(region_name), ...
    fit_text, rho_text);
end

function report_intercept_profile(context, profile, rho, p_rho)
fprintf('%s\n', context);
fprintf(['  Power ~ 1 + (1|Animal) at each %g µm bin. ' ...
    'One power per session, sessions repeated within the animal.\n' ...
    '  %d distances fit, %d skipped.\n'], ...
    profile.bin_um, profile.n_fit, profile.n_skipped);
if profile.n_fallback > 0
    fprintf(['  %d distances used the mean of animal means because the ' ...
        'mixed model did not return an intercept.\n'], profile.n_fallback);
end
if isfinite(rho)
    fprintf('  Spearman rho = %.2f, p = %s on the channel values.\n', ...
        rho, format_p_value(p_rho));
else
    fprintf('  Spearman rho is undefined.\n');
end
end

function [rho, p_value] = rank_correlation(distance, power)
rho = NaN;
p_value = NaN;
distance = distance(:);
power = power(:);
keep = isfinite(distance) & isfinite(power);
if nnz(keep) < 3
    return
end
distance = distance(keep);
power = power(keep);
if numel(unique(distance)) < 2 || numel(unique(power)) < 2
    return
end
[rho, p_value] = corr(distance, power, 'Type', 'Spearman');
end

function plot_region_power_comparison(power_table, ...
        delta_band, theta_band, gamma_band, spectral_method, ...
        regions_to_display, bands_to_display, plot_style)
if isempty(power_table)
    return
end
region_keys = region_keys_for_rows(power_table.Area);
region_names = selected_regions(regions_to_display);
region_names = region_names(~strcmp(region_names, 'Other'));
region_names = unique(region_names, 'stable');
plot_nested_region_bars(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, region_keys, region_names, ...
    'nested region comparison', 'nested repeated measures', ...
    ['Channels stay in the test, nested in the session. ' ...
    'Lines join one session. Brackets are mixed-model contrasts, Bonferroni.'], ...
    bands_to_display, plot_style);
end

function plot_region_half_comparison(power_table, ...
        delta_band, theta_band, gamma_band, spectral_method, ...
        regions_to_display, bands_to_display, plot_style)
if isempty(power_table) || ~ismember('DistanceToAreaBorder_um', power_table.Properties.VariableNames)
    return
end
region_names = selected_regions(regions_to_display);
region_names = region_names(~strcmp(region_names, 'Other'));
region_names = unique(region_names, 'stable');
[half_keys, half_names] = proximal_distal_keys(power_table, region_names);
if numel(half_names) < 2
    fprintf(['Proximal and distal comparison needs two halves with channels ' ...
        'in the selected regions.\n']);
    return
end
plot_nested_region_bars(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, half_keys, half_names, ...
    'proximal and distal comparison', 'proximal and distal', ...
    ['Proximal is the nearer half of the available distance. ' ...
    'Distal is the farther half. Lines join one session. ' ...
    'Brackets are mixed-model contrasts, Bonferroni.'], ...
    bands_to_display, plot_style);
end

function [half_keys, half_names] = proximal_distal_keys(power_table, region_names)
region_keys = region_keys_for_rows(power_table.Area);
distance_um = power_table.DistanceToMidLPAG_um(:);
alive = isfinite(power_table.DeltaPower_zWithinShank) | ...
    isfinite(power_table.ThetaPower_zWithinShank) | ...
    isfinite(power_table.GammaPower_zWithinShank);
half_keys = cell(height(power_table), 1);
for row_n = 1:height(power_table)
    half_keys{row_n} = {};
end
half_names = cell(1, 0);
methods = unique(cellstr(string(power_table.Method)), 'stable');
if isempty(methods)
    return
end
in_reference = strcmp(power_table.Method, methods{1});
for region_n = 1:numel(region_names)
    this_region = region_names{region_n};
    border_um = border_for_region(power_table, region_keys, ...
        this_region, distance_um);
    in_region = region_row_has(region_keys, this_region) & alive & ...
        isfinite(border_um);
    span_rows = in_region & in_reference;
    span = border_um(span_rows);
    if isempty(span)
        span = border_um(in_region);
    end
    if isempty(span)
        fprintf('%s has no channels left for a proximal and distal split.\n', ...
            region_display_name(this_region));
        continue
    end
    midpoint = (min(span) + max(span)) / 2;
    proximal_name = [this_region, ' proximal'];
    distal_name = [this_region, ' distal'];
    if strcmp(this_region, 'LPAG')
        reference = 'mid-LPAG';
    else
        reference = 'the LPAG-facing border';
    end
    fprintf(['%s available distance %g to %g µm from %s. ' ...
        'Proximal is through %g µm. Distal is beyond that.\n'], ...
        region_display_name(this_region), min(span), max(span), ...
        reference, midpoint);
    half_names{end + 1} = proximal_name; %#ok<AGROW>
    half_names{end + 1} = distal_name; %#ok<AGROW>
    for row_n = find(in_region).'
        if border_um(row_n) <= midpoint
            half_keys{row_n}{end + 1} = proximal_name; %#ok<AGROW>
        else
            half_keys{row_n}{end + 1} = distal_name; %#ok<AGROW>
        end
    end
end
end

function plot_nested_region_bars(power_table, ...
        delta_band, theta_band, gamma_band, spectral_method, ...
        region_keys, region_names, figure_name, title_name, note, ...
        bands_to_display, plot_style)
if numel(region_names) < 2
    fprintf(['Region comparison needs at least two regions in ' ...
        'regions_to_display.\n']);
    return
end
[method_names, ~] = selected_spectral_methods(spectral_method);
selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band);
n_bands = numel(selection.index);
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    band_obs = cell(1, n_bands);
    enough = false(1, n_bands);
    for slot_n = 1:n_bands
        band_obs{slot_n} = nested_channel_observations(power_table, ...
            in_method, region_keys, region_names, ...
            power_table.(selection.fields{slot_n}));
        enough(slot_n) = band_obs{slot_n}.n_sessions >= 2 && ...
            numel(band_obs{slot_n}.region_names) >= 2;
    end
    if ~any(enough)
        fprintf(['%s %s needs channels from at least two sessions in ' ...
            'at least two regions.\n'], method_names{method_n}, figure_name);
        continue
    end
    figure('Name', sprintf('%s | %s', method_names{method_n}, figure_name), ...
        'Color', 'w', ...
        'Position', [110, 80, max(720, 360 * n_bands), 560]);
    for slot_n = 1:n_bands
        subplot(1, n_bands, slot_n)
        obs = band_obs{slot_n};
        if ~enough(slot_n)
            axis off
            text(0.5, 0.5, sprintf('%s\nfewer than 2 sessions', ...
                selection.labels{slot_n}), 'HorizontalAlignment', 'center');
            fprintf(['%s | %s | %s: fewer than 2 sessions with channels in ' ...
                'at least two regions.\n'], method_names{method_n}, ...
                figure_name, selection.labels{slot_n});
            continue
        end
        result = fit_nested_region_model(obs);
        report_nested_region_tests(method_names{method_n}, selection.labels{slot_n}, ...
            region_names, obs, result);
        draw_nested_region_comparison(obs, result, selection.labels{slot_n}, plot_style);
    end
    sgtitle({sprintf('%s | %s', method_names{method_n}, title_name), note}, ...
        'FontSize', 10)
end
end

function obs = nested_channel_observations(power_table, row_mask, ...
        region_keys, region_names, power)
obs = empty_nested_observations();
session_index = session_groups(power_table, row_mask);
n_regions = numel(region_names);
power_cells = cell(1, 0);
region_cells = cell(1, 0);
session_cells = cell(1, 0);
for session_n = 1:numel(session_index)
    rows = session_index{session_n};
    these_keys = region_keys(rows);
    these_power = power(rows);
    for region_n = 1:n_regions
        in_region = region_row_has(these_keys, region_names{region_n});
        samples = these_power(in_region);
        samples = samples(isfinite(samples));
        if isempty(samples)
            continue
        end
        power_cells{end + 1} = samples(:); %#ok<AGROW>
        region_cells{end + 1} = repmat(region_n, numel(samples), 1); %#ok<AGROW>
        session_cells{end + 1} = repmat(session_n, numel(samples), 1); %#ok<AGROW>
    end
end
if isempty(power_cells)
    return
end
region_id = vertcat(region_cells{:});
present = false(1, n_regions);
present(unique(region_id)) = true;
kept = find(present);
remap = zeros(1, n_regions);
remap(kept) = 1:numel(kept);
obs.power = vertcat(power_cells{:});
obs.region_id = remap(region_id);
obs.region_id = obs.region_id(:);
obs.session_id = vertcat(session_cells{:});
[~, ~, obs.session_id] = unique(obs.session_id, 'stable');
obs.region_names = region_names(kept);
obs.n_sessions = numel(unique(obs.session_id));
obs.session_mean = session_region_means(obs);
end

function obs = empty_nested_observations()
obs = struct( ...
    'power', [], ...
    'region_id', [], ...
    'session_id', [], ...
    'region_names', {{}}, ...
    'n_sessions', 0, ...
    'session_mean', []);
end

function session_mean = session_region_means(obs)
% Anchors for the session lines. The mixed model keeps every channel.
n_sessions = obs.n_sessions;
n_regions = numel(obs.region_names);
session_mean = nan(n_sessions, n_regions);
for session_n = 1:n_sessions
    for region_n = 1:n_regions
        samples = obs.power(obs.session_id == session_n & ...
            obs.region_id == region_n);
        if isempty(samples)
            continue
        end
        session_mean(session_n, region_n) = mean(samples);
    end
end
end

function group_index = session_groups(power_table, row_mask)
rows = find(row_mask);
animal = cellstr(string(power_table.Animal(rows)));
session_key = cellstr(string(power_table.Session(rows)));
group_key = strcat(animal, '|', session_key);
[~, ~, group_id] = unique(group_key, 'stable');
n_groups = max([group_id; 0]);
group_index = cell(n_groups, 1);
for group_n = 1:n_groups
    group_index{group_n} = rows(group_id == group_n);
end
end

function result = fit_nested_region_model(obs)
n_regions = numel(obs.region_names);
n_pairs = n_regions * (n_regions - 1) / 2;
result = struct( ...
    'ok', false, ...
    'p_omnibus', NaN, ...
    'f_stat', NaN, ...
    'df1', NaN, ...
    'df2', NaN, ...
    'p_raw', nan(n_pairs, 1), ...
    'p_adjusted', nan(n_pairs, 1), ...
    'pairs', nchoosek(1:n_regions, 2), ...
    'emm', nan(1, n_regions), ...
    'emm_se', nan(1, n_regions), ...
    'session_sd', NaN, ...
    'cell_sd', NaN, ...
    'channel_sd', NaN, ...
    'message', '');
region_label = obs.region_names(obs.region_id);
region_label = region_label(:);
session_label = cellstr(string(obs.session_id));
try
    model_table = table(obs.power(:), ...
        categorical(region_label, obs.region_names), ...
        categorical(session_label), ...
        'VariableNames', {'Power', 'Region', 'Session'});
    lme = fitlme(model_table, ...
        'Power ~ Region + (1|Session) + (1|Session:Region)', ...
        'FitMethod', 'REML', 'DummyVarCoding', 'reference');
catch ME
    result.message = ME.message;
    return
end
try
    stats = anova(lme, 'DFMethod', 'satterthwaite');
    [result.f_stat, result.df1, result.df2, result.p_omnibus] = ...
        region_anova_row(stats);
catch ME
    result.message = ME.message;
end
coef_names = coefficient_names(lme);
beta = fixedEffects(lme);
covB = lme.CoefficientCovariance;
for region_n = 1:n_regions
    L = level_row(coef_names, obs.region_names{region_n}, obs.region_names{1});
    if any(~isfinite(L))
        continue
    end
    result.emm(region_n) = L * beta(:);
    result.emm_se(region_n) = sqrt(max(L * covB * L', 0));
end
for pair_n = 1:size(result.pairs, 1)
    L = contrast_vector(coef_names, obs.region_names, ...
        result.pairs(pair_n, 1), result.pairs(pair_n, 2));
    if any(~isfinite(L)) || ~any(L)
        continue
    end
    try
        result.p_raw(pair_n) = coefTest(lme, L, 0, 'DFMethod', 'satterthwaite');
    catch
        try
            result.p_raw(pair_n) = coefTest(lme, L, 0);
        catch
            result.p_raw(pair_n) = NaN;
        end
    end
end
result.p_adjusted = min(result.p_raw * size(result.pairs, 1), 1);
try
    [psi, channel_var] = covarianceParameters(lme);
    result.session_sd = sqrt(max(psi{1}(1, 1), 0));
    result.cell_sd = sqrt(max(psi{2}(1, 1), 0));
    result.channel_sd = sqrt(max(channel_var, 0));
catch
end
result.ok = isfinite(result.p_omnibus) || any(isfinite(result.p_adjusted));
if ~result.ok && isempty(result.message)
    result.message = 'The region effect could not be tested.';
end
end

function names = coefficient_names(lme)
names = lme.CoefficientNames;
if isstring(names)
    names = cellstr(names);
end
names = names(:).';
end

function [f_stat, df1, df2, p_value] = region_anova_row(stats)
f_stat = NaN;
df1 = NaN;
df2 = NaN;
p_value = NaN;
terms = column_as_cells(stats.Term);
row = find(strcmp(strtrim(terms), 'Region'), 1);
if isempty(row)
    return
end
f_stat = stats.FStat(row);
df1 = stats.DF1(row);
df2 = stats.DF2(row);
p_value = stats.pValue(row);
end

function terms = column_as_cells(column)
if iscell(column)
    terms = column(:);
    return
end
terms = cellstr(string(column));
terms = terms(:);
end

function L = level_row(coef_names, region_name, reference_name)
L = zeros(1, numel(coef_names));
intercept = find(strcmp(coef_names, '(Intercept)'), 1);
if isempty(intercept)
    L(:) = NaN;
    return
end
L(intercept) = 1;
if strcmp(region_name, reference_name)
    return
end
idx = coefficient_index(coef_names, region_name);
if isempty(idx)
    L(:) = NaN;
    return
end
L(idx) = 1;
end

function L = contrast_vector(coef_names, region_names, index_a, index_b)
reference = region_names{1};
weight_a = level_row(coef_names, region_names{index_a}, reference);
weight_b = level_row(coef_names, region_names{index_b}, reference);
if any(~isfinite(weight_a)) || any(~isfinite(weight_b))
    L = nan(1, numel(coef_names));
    return
end
L = weight_a - weight_b;
end

function idx = coefficient_index(coef_names, region_name)
candidates = {['Region_' region_name], ...
    ['Region_' strrep(region_name, ' ', '_')]};
for candidate_n = 1:numel(candidates)
    idx = find(strcmp(coef_names, candidates{candidate_n}), 1);
    if ~isempty(idx)
        return
    end
end
idx = [];
flat_name = strrep(region_name, ' ', '_');
for name_n = 1:numel(coef_names)
    name = coef_names{name_n};
    prefix = 'Region_';
    if numel(name) <= numel(prefix) || ~strcmp(name(1:numel(prefix)), prefix)
        continue
    end
    level = name(numel(prefix) + 1:end);
    if strcmp(level, region_name) || strcmp(level, flat_name)
        idx = name_n;
        return
    end
end
end

function report_nested_region_tests(method_name, band_name, ...
        requested_names, obs, result)
region_list = strjoin(cellfun(@region_display_name, obs.region_names, ...
    'UniformOutput', false), ', ');
fprintf(['%s | %s: %d channels nested in %d sessions (%s).\n' ...
    '  Power ~ Region + (1|Session) + (1|Session:Region)\n'], ...
    method_name, band_name, numel(obs.power), obs.n_sessions, region_list);
missing = requested_names(~ismember(requested_names, obs.region_names));
if ~isempty(missing)
    missing_list = strjoin(cellfun(@region_display_name, missing, ...
        'UniformOutput', false), ', ');
    fprintf('  No channels in %s.\n', missing_list);
end
if ~result.ok
    if isempty(result.message)
        detail = 'The region effect could not be tested.';
    else
        detail = result.message;
    end
    fprintf('  Nested mixed model was not fit: %s\n', detail);
    return
end
fprintf('  Region effect F(%s, %s) = %s, p = %s\n', ...
    format_p_value(result.df1), format_p_value(result.df2), ...
    format_p_value(result.f_stat), format_p_value(result.p_omnibus));
if isfinite(result.session_sd)
    fprintf(['  SD: session %s, session-by-region %s, channel %s\n'], ...
        format_p_value(result.session_sd), format_p_value(result.cell_sd), ...
        format_p_value(result.channel_sd));
    if result.cell_sd < 1e-4
        fprintf(['  The session-by-region variance is at the boundary, ' ...
            'so this p-value leans on the channel residual.\n']);
    end
end
for pair_n = 1:size(result.pairs, 1)
    fprintf('  %s vs %s: p = %s, Bonferroni p = %s %s\n', ...
        region_display_name(obs.region_names{result.pairs(pair_n, 1)}), ...
        region_display_name(obs.region_names{result.pairs(pair_n, 2)}), ...
        format_p_value(result.p_raw(pair_n)), ...
        format_p_value(result.p_adjusted(pair_n)), ...
        significance_label(result.p_adjusted(pair_n)));
end
end

function draw_nested_region_comparison(obs, result, band_name, plot_style)
n_sessions = obs.n_sessions;
n_regions = numel(obs.region_names);
if result.ok && all(isfinite(result.emm))
    means = result.emm;
    standard_error = result.emm_se;
else
    means = mean(obs.session_mean, 1, 'omitnan');
    n_per_region = sum(isfinite(obs.session_mean), 1);
    standard_error = std(obs.session_mean, 0, 1, 'omitnan') ./ sqrt(n_per_region);
end
hold on
yline(0, ':k', 'HandleVisibility', 'off');
jitter = ((1:n_sessions).' - (n_sessions + 1) / 2) / n_sessions * 0.22;
for session_n = 1:n_sessions
    y = obs.session_mean(session_n, :);
    keep = isfinite(y);
    if nnz(keep) < 2
        continue
    end
    plot(find(keep) + jitter(session_n), y(keep), '-', ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 0.6, ...
        'HandleVisibility', 'off');
end
bar_half = 0.275;
for region_n = 1:n_regions
    if ~isfinite(means(region_n))
        continue
    end
    fill(region_n + [-bar_half, bar_half, bar_half, -bar_half], ...
        [0, 0, means(region_n), means(region_n)], ...
        region_color(obs.region_names{region_n}), ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
errorbar(1:n_regions, means, standard_error, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
colors = zeros(numel(obs.power), 3);
for region_n = 1:n_regions
    in_region = obs.region_id == region_n;
    colors(in_region, :) = repmat(region_color(obs.region_names{region_n}), ...
        nnz(in_region), 1);
end
scatter(obs.region_id + jitter(obs.session_id), obs.power, ...
    plot_style.marker_size, colors, 'filled', ...
    'MarkerFaceAlpha', plot_style.marker_alpha, 'MarkerEdgeColor', 'none', ...
    'HandleVisibility', 'off');
y_extent = [0; obs.power(:); (means + standard_error).'; (means - standard_error).'];
y_extent = y_extent(isfinite(y_extent));
y_low = min(y_extent);
y_high = max(y_extent);
span = y_high - y_low;
if ~isfinite(span) || span < 0.5
    span = 0.5;
end
if result.ok
    draw_significance_brackets(result.pairs, result.p_adjusted, y_high, span);
    n_pairs = size(result.pairs, 1);
else
    n_pairs = 0;
end
margin = 0.08 * span;
ylim([y_low - margin, y_high + margin + n_pairs * (0.12 * span) + 0.12 * span]);
xlim([0.4, n_regions + 0.6]);
set(gca, 'XTick', 1:n_regions, 'XTickLabel', ...
    cellfun(@region_display_name, obs.region_names, 'UniformOutput', false));
if n_regions > 4
    xtickangle(30);
end
ylabel({'Power z-scored within shank', '0 = mid-LPAG channel'})
if result.ok && isfinite(result.p_omnibus)
    omnibus_label = sprintf('nested p = %.3g', result.p_omnibus);
else
    omnibus_label = 'nested model unavailable';
end
title(sprintf('%s | %s | %d sessions', band_name, omnibus_label, n_sessions))
box off
end

function label = format_p_value(p_value)
if ~isfinite(p_value)
    label = 'n/a';
else
    label = sprintf('%.4g', p_value);
end
end

function draw_significance_brackets(pairs, p_adjusted, y_high, span)
if isempty(pairs)
    return
end
span_width = pairs(:, 2) - pairs(:, 1);
[~, draw_order] = sort(span_width, 'ascend');
step = 0.12 * span;
tick = 0.025 * span;
base = y_high + 0.08 * span;
for draw_n = 1:numel(draw_order)
    pair_n = draw_order(draw_n);
    x1 = pairs(pair_n, 1);
    x2 = pairs(pair_n, 2);
    y = base + (draw_n - 1) * step;
    plot([x1, x1, x2, x2], [y - tick, y, y, y - tick], 'k-', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    text((x1 + x2) / 2, y, significance_label(p_adjusted(pair_n)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9, 'Interpreter', 'none');
end
end

function label = significance_label(p_value)
if ~isfinite(p_value)
    label = 'n/a';
elseif p_value < 0.001
    label = '***';
elseif p_value < 0.01
    label = '**';
elseif p_value < 0.05
    label = '*';
else
    label = 'n.s.';
end
end

function track = empty_track()
track = struct( ...
    'label', '', ...
    'rel_um', [], ...
    'welch_bands', [], ...
    'fourier_bands', [], ...
    'welch_spectrum', [], ...
    'fourier_spectrum', [], ...
    'segment_color', [], ...
    'goes_sc', false, ...
    'goes_dr', false, ...
    'stays_pag', false);
end

function gaussians = gaussian_params_for(fitted, n_channels)
if ~isfield(fitted, 'gaussian_params') || ~iscell(fitted.gaussian_params) || ...
        numel(fitted.gaussian_params) ~= n_channels
    error(['Saved spectra have no per-channel gaussian_params, so theta ' ...
        'and its harmonic cannot be removed from the mid-LPAG spectra.']);
end
gaussians = fitted.gaussian_params;
end

function tracks = tracks_from_session(channel_table, welch_spectrum, ...
        fourier_spectrum, frequency, delta_band, theta_band, gamma_band, ...
        session_label, artifact_trend_sd, welch_gaussians, fourier_gaussians)
groups = probe_groups_by_depth(channel_table);
tracks = repmat(empty_track(), 0, 1);
for group_n = 1:numel(groups)
    this_group = groups(group_n);
    center = this_group.depth_um(this_group.is_lpag_center);
    center = center(isfinite(center));
    if isempty(center)
        fprintf(['  %s, %s has no mid-LPAG channel and is left out of ' ...
            'the stack.\n'], session_label, this_group.label);
        continue
    end
    rel = this_group.depth_um(:) - mean(center);
    [segment_color, goes_sc, goes_dr, stays_pag] = ...
        classify_track(rel, this_group.areas);
    welch_rows = welch_spectrum(this_group.channel_index, :);
    fourier_rows = fourier_spectrum(this_group.channel_index, :);
    track = empty_track();
    track.label = sprintf('%s | %s', session_label, this_group.label);
    track.rel_um = rel;
    track.welch_bands = reject_trend_outliers(rel, [ ...
        mean_in_band(welch_rows, frequency, delta_band), ...
        mean_in_band(welch_rows, frequency, theta_band), ...
        mean_in_band(welch_rows, frequency, gamma_band)], artifact_trend_sd);
    track.fourier_bands = reject_trend_outliers(rel, [ ...
        mean_in_band(fourier_rows, frequency, delta_band), ...
        mean_in_band(fourier_rows, frequency, theta_band), ...
        mean_in_band(fourier_rows, frequency, gamma_band)], artifact_trend_sd);
    track.welch_spectrum = welch_rows - theta_harmonic_image( ...
        welch_gaussians, this_group.channel_index, frequency, theta_band);
    track.fourier_spectrum = fourier_rows - theta_harmonic_image( ...
        fourier_gaussians, this_group.channel_index, frequency, theta_band);
    track.segment_color = segment_color;
    track.goes_sc = goes_sc;
    track.goes_dr = goes_dr;
    track.stays_pag = stays_pag;
    tracks(end + 1, 1) = track; %#ok<AGROW>
end
end

function [segment_color, goes_sc, goes_dr, stays_pag] = ...
        classify_track(rel_um, area_labels)
palette = track_colors();
n_sites = numel(rel_um);
segment_color = repmat(palette.other, n_sites, 1);
is_pag = false(n_sites, 1);
is_sc = false(n_sites, 1);
is_dr = false(n_sites, 1);
has_label = false(n_sites, 1);
has_other = false(n_sites, 1);
for site_n = 1:n_sites
    tokens = split_area_label(area_labels{site_n});
    if isempty(tokens)
        continue
    end
    has_label(site_n) = true;
    pag_token = cellfun(@name_is_pag, tokens);
    sc_token = cellfun(@name_is_supcol, tokens);
    dr_token = cellfun(@name_is_dr, tokens);
    is_pag(site_n) = any(pag_token);
    is_sc(site_n) = any(sc_token);
    is_dr(site_n) = any(dr_token);
    has_other(site_n) = any(~pag_token & ~sc_token & ~dr_token);
end
goes_sc = any(is_sc & rel_um(:) > 0);
goes_dr = any(is_dr & rel_um(:) < 0);
stays_pag = any(has_label) && ~goes_sc && ~goes_dr && ...
    ~any(has_other) && ~any(is_sc) && ~any(is_dr);
for site_n = 1:n_sites
    if is_sc(site_n) && rel_um(site_n) > 0 && goes_sc
        segment_color(site_n, :) = palette.upper;
    elseif is_dr(site_n) && rel_um(site_n) < 0 && goes_dr
        segment_color(site_n, :) = palette.deeper;
    elseif is_pag(site_n) && ~has_other(site_n)
        segment_color(site_n, :) = palette.pag;
    end
end
end

function tokens = split_area_label(label)
if isempty(label) || (isstring(label) && ismissing(label))
    tokens = {};
    return
end
tokens = strtrim(strsplit(char(label), ','));
tokens = tokens(~cellfun('isempty', tokens));
end

function tf = name_is_pag(name)
tf = numel(name) >= 3 && strcmp(name(end - 2:end), 'PAG');
end

function tf = name_is_supcol(name)
tf = strcmpi(name, 'SupCol') || strcmpi(name, 'SC');
end

function tf = name_is_dr(name)
tf = strcmp(name, 'DR');
end

function palette = track_colors()
palette.pag = [0.13 0.55 0.36];
palette.upper = [0.90 0.45 0.05];
palette.deeper = [0.36 0.28 0.72];
palette.other = [0.72 0.72 0.72];
end

function plot_aligned_profiles(tracks, spectrum_field, ...
        delta_band, theta_band, gamma_band, spectral_method, bands_to_display)
palette = track_colors();
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);
selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band);
n_bands = numel(selection.index);
figure('Name', 'Profiles aligned on mid-LPAG', 'Color', 'w', ...
    'Position', [60 60 max(420, 360 * n_bands) 720]);
for method_slot = 1:n_methods
    method_n = method_ids(method_slot);
    for slot_n = 1:n_bands
        band_n = selection.index(slot_n);
        subplot(n_methods, n_bands, (method_slot - 1) * n_bands + slot_n)
        hold on
        yline(0, '--k', 'mid-LPAG', 'HandleVisibility', 'off', ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 8);
        if strcmp(spectrum_field, 'raw_aperiodic_subtracted')
            xline(0, ':k', 'HandleVisibility', 'off');
        end
        for track_n = 1:numel(tracks)
            power = band_column(tracks(track_n), method_n, band_n);
            draw_colored_profile(tracks(track_n).rel_um, power, ...
                tracks(track_n).segment_color);
        end
        guide = [ ...
            plot(nan, nan, 'Color', palette.pag, 'LineWidth', 1.5); ...
            plot(nan, nan, 'Color', palette.upper, 'LineWidth', 1.5); ...
            plot(nan, nan, 'Color', palette.deeper, 'LineWidth', 1.5); ...
            plot(nan, nan, 'Color', palette.other, 'LineWidth', 1.5)];
        guide_labels = {'Within PAG', 'Upper side through SupCol', ...
            'Deeper side through DR', 'Other'};
        xlabel(profile_xlabel(spectrum_field))
        if slot_n == 1
            ylabel({'Distance to tip relative to mid-LPAG (µm)', ...
                'positive = upper, further from tip'})
        end
        title(sprintf('%s | %s', method_names{method_slot}, selection.labels{slot_n}))
        if method_slot == 1 && slot_n == n_bands
            legend(guide, guide_labels, 'Location', 'best', 'FontSize', 7)
        end
        box off
    end
end
sgtitle('Sessions stacked on the mid-LPAG channel', 'Interpreter', 'none')
end

function plot_aligned_intervals(tracks, spectrum_field, ...
        delta_band, theta_band, gamma_band, spectral_method, bands_to_display)
palette = track_colors();
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);
selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band);
n_bands = numel(selection.index);
figure('Name', 'Confidence intervals aligned on mid-LPAG', 'Color', 'w', ...
    'Position', [100 100 max(420, 360 * n_bands) 720]);
for method_slot = 1:n_methods
    method_n = method_ids(method_slot);
    for slot_n = 1:n_bands
        band_n = selection.index(slot_n);
        subplot(n_methods, n_bands, (method_slot - 1) * n_bands + slot_n)
        hold on
        yline(0, '--k', 'mid-LPAG', 'HandleVisibility', 'off', ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 8);
        if strcmp(spectrum_field, 'raw_aperiodic_subtracted')
            xline(0, ':k', 'HandleVisibility', 'off');
        end
        ci_handles = gobjects(0, 1);
        ci_labels = {};
        [ci_handles, ci_labels] = add_class_interval(ci_handles, ci_labels, ...
            tracks, method_n, band_n, 'stays_pag', palette.pag, 'Stays in PAG');
        [ci_handles, ci_labels] = add_class_interval(ci_handles, ci_labels, ...
            tracks, method_n, band_n, 'goes_sc', palette.upper, ...
            'Through SupCol');
        [ci_handles, ci_labels] = add_class_interval(ci_handles, ci_labels, ...
            tracks, method_n, band_n, 'goes_dr', palette.deeper, ...
            'Through DR');
        xlabel(profile_xlabel(spectrum_field))
        if slot_n == 1
            ylabel({'Distance to tip relative to mid-LPAG (µm)', ...
                'positive = upper, further from tip'})
        end
        title(sprintf('%s | %s', method_names{method_slot}, selection.labels{slot_n}))
        if method_slot == 1 && slot_n == n_bands && ~isempty(ci_handles)
            legend(ci_handles, ci_labels, 'Location', 'best', 'FontSize', 7)
        end
        box off
    end
end
sgtitle('95% confidence intervals, centered on mid-LPAG', 'Interpreter', 'none')
end

function draw_colored_profile(rel_um, power, segment_color)
[rel_um, order] = sort(rel_um(:));
power = power(order);
segment_color = segment_color(order, :);
for site_n = 1:numel(rel_um) - 1
    if ~all(isfinite(rel_um(site_n:site_n + 1))) || ...
            ~all(isfinite(power(site_n:site_n + 1)))
        continue
    end
    plot(power(site_n:site_n + 1), rel_um(site_n:site_n + 1), ...
        'Color', segment_color(site_n, :), 'LineWidth', 1, ...
        'HandleVisibility', 'off');
end
end

function [handles, labels] = add_class_interval(handles, labels, ...
        tracks, method_n, band_n, field_name, color, class_label)
selected = tracks([tracks.(field_name)]);
n_tracks = numel(selected);
if n_tracks == 0
    return
end
[grid, mean_line, lower_line, upper_line] = class_interval( ...
    selected, method_n, band_n);
if strcmp(field_name, 'goes_sc')
    off_side = grid < 0;
    mean_line(off_side) = NaN;
    lower_line(off_side) = NaN;
    upper_line(off_side) = NaN;
elseif strcmp(field_name, 'goes_dr')
    off_side = grid > 0;
    mean_line(off_side) = NaN;
    lower_line(off_side) = NaN;
    upper_line(off_side) = NaN;
end
if ~any(isfinite(mean_line))
    fprintf('%s: %d tracks, no confidence interval (need at least 2).\n', ...
        class_label, n_tracks);
    return
end
handles(end + 1, 1) = plot(mean_line, grid, 'Color', color, ...
    'LineWidth', 2.4); %#ok<AGROW>
plot(lower_line, grid, 'Color', color, 'LineWidth', 1, ...
    'HandleVisibility', 'off');
plot(upper_line, grid, 'Color', color, 'LineWidth', 1, ...
    'HandleVisibility', 'off');
labels{end + 1} = sprintf('%s 95%% CI (n=%d)', class_label, n_tracks); %#ok<AGROW>
end

function power = band_column(track, method_n, band_n)
if method_n == 1
    power = track.welch_bands(:, band_n);
else
    power = track.fourier_bands(:, band_n);
end
end

function [grid, mean_line, lower_line, upper_line] = class_interval( ...
        tracks, method_n, band_n)
[grid, resampled] = resample_tracks(tracks, method_n, band_n);
n_at_grid = sum(isfinite(resampled), 1);
mean_line = mean(resampled, 1, 'omitnan');
deviation = std(resampled, 0, 1, 'omitnan');
standard_error = deviation ./ sqrt(n_at_grid);
critical = nan(size(n_at_grid));
enough = n_at_grid >= 2;
critical(enough) = tinv(0.975, n_at_grid(enough) - 1);
half_width = critical .* standard_error;
lower_line = mean_line - half_width;
upper_line = mean_line + half_width;
mean_line(~enough) = NaN;
lower_line(~enough) = NaN;
upper_line(~enough) = NaN;
end

function [grid, resampled] = resample_tracks(tracks, method_n, band_n)
grid = alignment_grid(tracks);
resampled = nan(numel(tracks), numel(grid));
for track_n = 1:numel(tracks)
    resampled(track_n, :) = resample_one(tracks(track_n).rel_um, ...
        band_column(tracks(track_n), method_n, band_n), grid);
end
end

function values = resample_one(rel_um, samples, grid)
[rel_um, order] = sort(rel_um(:));
samples = samples(order);
keep = isfinite(rel_um) & isfinite(samples);
rel_um = rel_um(keep);
samples = samples(keep);
[rel_um, unique_index] = unique(rel_um, 'stable');
samples = samples(unique_index);
values = nan(size(grid));
if numel(rel_um) < 2
    return
end
values = interp1(rel_um, samples, grid, 'linear', NaN);
end

function grid = alignment_grid(tracks)
samples = [];
for track_n = 1:numel(tracks)
    samples = [samples; tracks(track_n).rel_um(:)]; %#ok<AGROW>
end
samples = samples(isfinite(samples));
if isempty(samples)
    grid = 0;
    return
end
step = 20;
grid = (floor(min(samples) / step) * step):step:(ceil(max(samples) / step) * step);
end

function plot_aligned_spectra(tracks, frequency, spectrum_field, spectral_method)
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);
class_fields = {'stays_pag', 'goes_sc', 'goes_dr'};
class_titles = {'Stays in PAG', 'Through SupCol', 'Through DR'};
figure('Name', 'Spectra aligned on mid-LPAG', 'Color', 'w', ...
    'Position', [80 80 1100 720]);
colors = diverging_or_sequential(spectrum_field);
for method_slot = 1:n_methods
    method_n = method_ids(method_slot);
    class_means = cell(1, 3);
    finite_values = [];
    for class_n = 1:3
        selected = tracks([tracks.(class_fields{class_n})]);
        [grid, mean_spectrum] = mean_aligned_spectrum(selected, method_n, frequency);
        class_means{class_n} = struct('grid', grid, 'spectrum', mean_spectrum, ...
            'n', numel(selected));
        if ~isempty(mean_spectrum)
            finite_values = [finite_values; ...
                mean_spectrum(isfinite(mean_spectrum))]; %#ok<AGROW>
        end
    end
    color_limits = spectrum_color_limits(finite_values, spectrum_field);
    for class_n = 1:3
        subplot(n_methods, 3, (method_slot - 1) * 3 + class_n)
        this_class = class_means{class_n};
        if this_class.n == 0 || isempty(this_class.spectrum)
            axis off
            text(0.5, 0.5, sprintf('%s\n%s\nn = 0', ...
                method_names{method_slot}, class_titles{class_n}), ...
                'HorizontalAlignment', 'center', 'Interpreter', 'none');
            continue
        end
        imagesc(frequency, this_class.grid, this_class.spectrum);
        set(gca, 'YDir', 'normal')
        axis tight
        colormap(colors)
        if all(isfinite(color_limits))
            caxis(color_limits);
        end
        yline(0, '--k', 'mid-LPAG', 'LineWidth', 1, ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 8);
        xlabel('Frequency (Hz)')
        if class_n == 1
            ylabel({'Distance to tip relative to mid-LPAG (µm)', ...
                'positive = upper, further from tip'})
        end
        title(sprintf('%s | %s (n=%d)', method_names{method_slot}, ...
            class_titles{class_n}, this_class.n))
        box off
    end
    colorbar
end
sgtitle({spectrum_title(spectrum_field), ...
    'Theta Gaussians and their first harmonic removed', ...
    'Mean of tracks in each class, centered on mid-LPAG'}, ...
    'Interpreter', 'none')
end

function image = theta_harmonic_image(gaussians, channel_index, frequency, theta_band)
n_rows = numel(channel_index);
image = zeros(n_rows, numel(frequency));
for row_n = 1:n_rows
    image(row_n, :) = reconstruct_theta_and_harmonic( ...
        gaussians{channel_index(row_n)}, frequency, theta_band);
end
end

function fit = reconstruct_theta_and_harmonic(gaussians, frequency, theta_band)
% Theta peaks are Gaussians with a center inside theta_band. The first
% harmonic is another Gaussian whose center lies near twice that theta
% center. The window is the larger of 1 Hz and 10% of the doubled
% frequency. Each peak is included once.
fit = zeros(1, numel(frequency));
if isempty(gaussians) || size(gaussians, 2) < 3 || size(gaussians, 1) < 1
    return
end
center = gaussians(:, 1);
peak_height = gaussians(:, 2);
width = gaussians(:, 3);
keep = isfinite(center) & isfinite(peak_height) & isfinite(width) & width ~= 0;
center = center(keep);
peak_height = peak_height(keep);
width = width(keep);
if isempty(center)
    return
end
is_theta = center >= theta_band(1) & center <= theta_band(2);
is_harmonic = false(size(center));
theta_index = find(is_theta);
for theta_n = theta_index(:).'
    expected = 2 * center(theta_n);
    tolerance = max(1, 0.1 * expected);
    is_harmonic = is_harmonic | (abs(center - expected) <= tolerance);
end
is_harmonic(is_theta) = false;
selected = is_theta | is_harmonic;
fit = sum_gaussians(center(selected), peak_height(selected), ...
    width(selected), frequency);
end

function fit = sum_gaussians(center, peak_height, width, frequency)
fit = zeros(1, numel(frequency));
for gaussian_n = 1:numel(center)
    fit = fit + peak_height(gaussian_n) .* exp( ...
        -(frequency - center(gaussian_n)).^2 ./ (2 * width(gaussian_n).^2));
end
end

function [grid, mean_spectrum] = mean_aligned_spectrum(tracks, method_n, frequency)
if isempty(tracks)
    grid = [];
    mean_spectrum = [];
    return
end
grid = alignment_grid(tracks);
n_frequency = numel(frequency);
resampled = nan(numel(tracks), numel(grid), n_frequency);
for track_n = 1:numel(tracks)
    if method_n == 1
        spectrum = tracks(track_n).welch_spectrum;
    else
        spectrum = tracks(track_n).fourier_spectrum;
    end
    for frequency_n = 1:n_frequency
        resampled(track_n, :, frequency_n) = resample_one( ...
            tracks(track_n).rel_um, spectrum(:, frequency_n), grid);
    end
end
mean_spectrum = squeeze(mean(resampled, 1, 'omitnan'));
if isvector(mean_spectrum)
    mean_spectrum = reshape(mean_spectrum, numel(grid), n_frequency);
end
end

function plot_region_distances(placed, spectrum_field, ...
        delta_band, theta_band, gamma_band, spectral_method, ...
        regions_to_display, region_distance_scale, distance_model_bin_um, ...
        bands_to_display, plot_style)
if ~any(strcmp(region_distance_scale, {'normalized', 'raw'}))
    error('region_distance_scale must be ''normalized'' or ''raw''.');
end
region_names = selected_regions(regions_to_display);
directions = {'ventral', 'dorsal'};
direction_titles = { ...
    'Ventral, toward the tip', ...
    'Dorsal, toward the brain surface'};
direction_xlabels = { ...
    'Distance from region center toward the tip (µm)', ...
    'Distance from region center toward the brain surface (µm)'};
selection = selected_bands(bands_to_display, delta_band, theta_band, gamma_band);
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_regions = numel(region_names);
n_placed = numel(placed);
centers = nan(n_placed, n_regions);
for placed_n = 1:n_placed
    for region_n = 1:n_regions
        centers(placed_n, region_n) = region_center_um( ...
            placed(placed_n).tip_um, placed(placed_n).areas, ...
            region_names{region_n});
    end
end
for region_n = 1:n_regions
    fprintf('%s center available on %d of %d shanks.\n', ...
        region_display_name(region_names{region_n}), ...
        nnz(isfinite(centers(:, region_n))), n_placed);
end

panel_kinds = {'tracks', 'model'};
for kind_n = 1:numel(panel_kinds)
    panel_kind = panel_kinds{kind_n};
    for method_slot = 1:numel(method_ids)
        method_n = method_ids(method_slot);
        if strcmp(panel_kind, 'tracks')
            figure_name = sprintf('%s power vs distance to region center', ...
                method_names{method_slot});
            figure_position = [40 40 980 max(480, 240 * n_regions)];
            figure_note = region_distance_note(region_distance_scale, false);
        else
            figure_name = sprintf('%s intercept by distance to region center', ...
                method_names{method_slot});
            figure_position = [90 40 980 max(480, 240 * n_regions)];
            figure_note = region_distance_note(region_distance_scale, true);
            figure_note = sprintf(['%s Rectangles cover the upper 90%% of the ' ...
                'axis where the intercept interval excludes 0 for at least %g µm.'], ...
                figure_note, plot_style.significance_min_um);
        end
        figure('Name', figure_name, 'Color', 'w', 'Position', figure_position);
        legend_drawn = false;
        for region_n = 1:n_regions
            for direction_n = 1:2
                subplot(n_regions, 2, (region_n - 1) * 2 + direction_n)
                available = find(isfinite(centers(:, region_n)));
                if isempty(available)
                    axis off
                    text(0.5, 0.5, sprintf('%s\nnot available', ...
                        region_display_name(region_names{region_n})), ...
                        'HorizontalAlignment', 'center', 'Interpreter', 'none');
                    continue
                end
                hold on
                side_tracks = repmat(empty_side_track(), 0, 1);
                farthest = 0;
                for placed_n = available(:).'
                    [distance, welch_power, fourier_power] = side_samples( ...
                        placed(placed_n), centers(placed_n, region_n), ...
                        directions{direction_n});
                    if isempty(distance)
                        continue
                    end
                    welch_std = std(placed(placed_n).welch_bands, 0, 1, 'omitnan');
                    fourier_std = std(placed(placed_n).fourier_bands, 0, 1, 'omitnan');
                    [welch_power, fourier_power] = scale_side_power( ...
                        distance, welch_power, fourier_power, ...
                        welch_std, fourier_std, region_distance_scale);
                    if strcmp(panel_kind, 'tracks')
                        draw_distance_scatter(distance, ...
                            power_for_method(welch_power, fourier_power, method_n), ...
                            selection, plot_style);
                    end
                    farthest = max(farthest, max(distance(isfinite(distance))));
                    side_track = empty_side_track();
                    side_track.rel_um = distance;
                    side_track.welch_bands = welch_power;
                    side_track.fourier_bands = fourier_power;
                    side_track.animal = placed(placed_n).animal;
                    side_track.session = placed(placed_n).session;
                    side_tracks(end + 1, 1) = side_track; %#ok<AGROW>
                end
                if isempty(side_tracks)
                    axis off
                    text(0.5, 0.5, sprintf('%s\n%s\nno channels on this side', ...
                        region_display_name(region_names{region_n}), ...
                        direction_titles{direction_n}), ...
                        'HorizontalAlignment', 'center', 'Interpreter', 'none');
                    continue
                end
                x_max = farthest;
                model_ok = false;
                drawn_profiles = cell(0, 1);
                drawn_colors = zeros(0, 3);
                if strcmp(panel_kind, 'model')
                    for slot_n = 1:numel(selection.index)
                        band_n = selection.index(slot_n);
                        [point_distance, point_power, point_animal, point_session] = ...
                            side_track_points(side_tracks, method_n, band_n);
                        finite_point = isfinite(point_distance) & isfinite(point_power);
                        if any(finite_point)
                            scatter(point_distance(finite_point), ...
                                point_power(finite_point), plot_style.marker_size, ...
                                selection.colors(slot_n, :), 'filled', ...
                                'MarkerFaceAlpha', plot_style.marker_alpha, ...
                                'MarkerEdgeAlpha', plot_style.marker_alpha, ...
                                'HandleVisibility', 'off');
                        end
                        profile = fit_intercepts_by_distance( ...
                            point_distance, point_power, point_animal, ...
                            point_session, distance_model_bin_um);
                        [rho, p_rho] = rank_correlation(point_distance, point_power);
                        report_intercept_profile(sprintf('%s | %s | %s | %s', ...
                            method_names{method_slot}, selection.labels{slot_n}, ...
                            region_display_name(region_names{region_n}), ...
                            direction_titles{direction_n}), ...
                            profile, rho, p_rho);
                        if profile.n_fit == 0
                            continue
                        end
                        model_ok = true;
                        draw_distance_intercepts(profile, selection.colors(slot_n, :), ...
                            plot_style.ci_alpha);
                        drawn_profiles{end + 1, 1} = profile; %#ok<AGROW>
                        drawn_colors(end + 1, :) = selection.colors(slot_n, :); %#ok<AGROW>
                        x_max = max(x_max, max(profile.distance_um));
                    end
                end
                if strcmp(region_distance_scale, 'normalized') || ...
                        strcmp(spectrum_field, 'raw_aperiodic_subtracted')
                    yline(0, ':k', 'HandleVisibility', 'off');
                end
                xline(0, '--k', 'HandleVisibility', 'off');
                if ~isfinite(x_max) || x_max <= 0
                    x_max = 20;
                end
                xlim([0, x_max]);
                if strcmp(panel_kind, 'model')
                    y_limits = ylim;
                    for drawn_n = 1:numel(drawn_profiles)
                        draw_distance_significance(drawn_profiles{drawn_n}, ...
                            drawn_colors(drawn_n, :), plot_style.significance_alpha, ...
                            plot_style.significance_min_um, y_limits);
                    end
                    ylim(y_limits);
                end
                xlabel(direction_xlabels{direction_n})
                if direction_n == 1
                    ylabel(region_distance_ylabel(spectrum_field, region_distance_scale))
                end
                title(sprintf('%s | %s (n=%d)', ...
                    region_display_name(region_names{region_n}), ...
                    direction_titles{direction_n}, numel(side_tracks)))
                draw_rank_correlations(side_tracks, method_n, selection);
                if strcmp(panel_kind, 'model') && ~model_ok
                    text(0.5, 0.5, 'No distance had 2 animals', ...
                        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                        'FontSize', 8);
                end
                if ~legend_drawn
                    guide = gobjects(numel(selection.index), 1);
                    for slot_n = 1:numel(selection.index)
                        guide(slot_n) = plot(nan, nan, 'o', ...
                            'Color', selection.colors(slot_n, :), ...
                            'MarkerFaceColor', selection.colors(slot_n, :), ...
                            'MarkerSize', 5);
                    end
                    legend(guide, selection.labels, 'Location', 'best', 'FontSize', 7)
                    legend_drawn = true;
                end
                box off
            end
        end
        sgtitle({sprintf('%s | %s', method_names{method_slot}, ...
            spectrum_title(spectrum_field)), figure_note}, ...
            'Interpreter', 'none')
    end
end
end

function note = region_distance_note(scale_mode, with_model)
if strcmp(scale_mode, 'normalized')
    note = ['Each shank is z-scored and set to 0 at its first channel. ' ...
        'Spearman rho and p compare distance with that power ' ...
        'after extreme values and uncommon border distances are removed.'];
else
    note = ['Saved band power. Spearman rho and p compare distance with power ' ...
        'after extreme values and uncommon border distances are removed.'];
end
if with_model
    note = [note, [' Dots are the channel values. At each distance, one ' ...
        'session power enters Power ~ 1 + (1|Animal). The marker is the ' ...
        'intercept and the shade is its 95% interval.']];
end
end

function label = region_distance_ylabel(spectrum_field, scale_mode)
if strcmp(scale_mode, 'normalized')
    label = {'Z-scored power', '0 = first channel on this side'};
    return
end
label = profile_xlabel(spectrum_field);
end

function [welch_power, fourier_power] = scale_side_power( ...
        distance, welch_power, fourier_power, welch_std, fourier_std, scale_mode)
if ~strcmp(scale_mode, 'normalized')
    return
end
welch_power = center_and_scale_side(distance, welch_power, welch_std);
fourier_power = center_and_scale_side(distance, fourier_power, fourier_std);
end

function power = center_and_scale_side(distance, power, band_std)
distance = distance(:);
finite_distance = distance(isfinite(distance));
if isempty(finite_distance)
    return
end
start = isfinite(distance) & distance == min(finite_distance);
reference = mean(power(start, :), 1, 'omitnan');
band_std = band_std(:).';
band_std(band_std == 0 | ~isfinite(band_std)) = NaN;
power = (power - reference) ./ band_std;
end

function draw_rank_correlations(side_tracks, method_n, selection)
for slot_n = 1:numel(selection.index)
    band_n = selection.index(slot_n);
    [rho, p_value] = panel_rank_correlation(side_tracks, method_n, band_n);
    if isfinite(rho)
        label = sprintf('\\rho = %.2f, p = %s', rho, format_p_value(p_value));
    else
        label = '\rho = n/a';
    end
    text(0.98, 0.98 - (slot_n - 1) * 0.09, label, ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'Color', selection.colors(slot_n, :), ...
        'FontSize', 8, 'Interpreter', 'tex');
end
end

function [rho, p_value] = panel_rank_correlation(side_tracks, method_n, band_n)
[distance, power] = side_track_points(side_tracks, method_n, band_n);
[rho, p_value] = rank_correlation(distance, power);
end

function [distance, power, animal, session] = side_track_points( ...
        side_tracks, method_n, band_n)
distance = [];
power = [];
animal = {};
session = [];
for track_n = 1:numel(side_tracks)
    this_distance = side_tracks(track_n).rel_um(:);
    this_power = band_column(side_tracks(track_n), method_n, band_n);
    n_points = numel(this_distance);
    distance = [distance; this_distance]; %#ok<AGROW>
    power = [power; this_power(:)]; %#ok<AGROW>
    animal = [animal; repmat({side_tracks(track_n).animal}, n_points, 1)]; %#ok<AGROW>
    session = [session; repmat(side_tracks(track_n).session, n_points, 1)]; %#ok<AGROW>
end
end

function center = region_center_um(tip_um, areas, region_name)
center = NaN;
tip_um = tip_um(:);
in_region = false(numel(tip_um), 1);
for site_n = 1:numel(tip_um)
    if ~isfinite(tip_um(site_n))
        continue
    end
    tokens = split_area_label(areas{site_n});
    for token_n = 1:numel(tokens)
        if token_is_region(tokens{token_n}, region_name)
            in_region(site_n) = true;
            break
        end
    end
end
span = tip_um(in_region);
if isempty(span)
    return
end
center = (min(span) + max(span)) / 2;
end

function tf = token_is_region(name, region_name)
switch region_name
    case 'SC'
        tf = name_is_supcol(name);
    case 'DR'
        tf = name_is_dr(name);
    case {'DLPAG', 'LPAG', 'VLPAG'}
        tf = strcmp(name, region_name);
    otherwise
        tf = false;
end
end

function [distance, welch_power, fourier_power] = side_samples(item, center, direction)
tip_um = item.tip_um(:);
if strcmp(direction, 'ventral')
    keep = isfinite(tip_um) & tip_um <= center;
    distance = center - tip_um(keep);
else
    keep = isfinite(tip_um) & tip_um >= center;
    distance = tip_um(keep) - center;
end
welch_power = item.welch_bands(keep, :);
fourier_power = item.fourier_bands(keep, :);
finite_row = isfinite(distance) & ...
    (any(isfinite(welch_power), 2) | any(isfinite(fourier_power), 2));
distance = distance(finite_row);
welch_power = welch_power(finite_row, :);
fourier_power = fourier_power(finite_row, :);
end

function power = power_for_method(welch_power, fourier_power, method_n)
if method_n == 1
    power = welch_power;
else
    power = fourier_power;
end
end

function item = empty_side_track()
item = struct('rel_um', [], 'welch_bands', [], 'fourier_bands', [], ...
    'animal', '', 'session', NaN);
end

function draw_distance_scatter(distance, power, selection, plot_style)
distance = distance(:);
for slot_n = 1:numel(selection.index)
    band_n = selection.index(slot_n);
    this_power = power(:, band_n);
    finite = isfinite(distance) & isfinite(this_power);
    if ~any(finite)
        continue
    end
    scatter(distance(finite), this_power(finite), plot_style.marker_size, ...
        selection.colors(slot_n, :), 'filled', ...
        'MarkerFaceAlpha', plot_style.marker_alpha, ...
        'MarkerEdgeAlpha', plot_style.marker_alpha, ...
        'HandleVisibility', 'off');
end
end

function item = empty_region_record()
item = struct( ...
    'animal', '', ...
    'session', NaN, ...
    'shank', NaN, ...
    'channel', NaN, ...
    'method_id', NaN, ...
    'regions', {{}}, ...
    'is_center', false, ...
    'spectrum', []);
end

function records = region_channel_records(channel_table, session, ...
        frequency, ranges, animal_name, session_n, spectral_method)
[~, method_ids] = selected_spectral_methods(spectral_method);
method_fields = {'welch', 'fourier'};
if ~ismember('IsLPAGCenter', channel_table.Properties.VariableNames)
    error(['Channel table has no IsLPAGCenter column, so the central ' ...
        'LPAG reference cannot be drawn.']);
end
n_channels = height(channel_table);
animal = char(string(animal_name));
records = repmat(empty_region_record(), 0, 1);
for slot_n = 1:numel(method_ids)
    method_n = method_ids(slot_n);
    residual = residual_matrix(session.(method_fields{method_n}), ...
        frequency, ranges, n_channels);
    these = repmat(empty_region_record(), n_channels, 1);
    for channel_n = 1:n_channels
        these(channel_n).animal = animal;
        these(channel_n).session = session_n;
        these(channel_n).shank = channel_table.ProbeNumber(channel_n);
        these(channel_n).channel = channel_table.Channel(channel_n);
        these(channel_n).method_id = method_n;
        these(channel_n).regions = area_to_region_keys( ...
            area_label_at(channel_table.Area, channel_n));
        these(channel_n).is_center = ...
            channel_table.IsLPAGCenter(channel_n) ~= 0;
        these(channel_n).spectrum = residual(channel_n, :);
    end
    records = [records; these]; %#ok<AGROW>
end
end

function label = area_label_at(area_column, index)
if iscell(area_column)
    label = area_column{index};
elseif isstring(area_column)
    label = area_column(index);
else
    label = area_column(index, :);
end
end

function keys = area_to_region_keys(label)
tokens = split_area_label(label);
keys = cell(1, 0);
for token_n = 1:numel(tokens)
    name = canonical_region_name(tokens{token_n});
    if strcmp(name, 'Other') || any(strcmp(keys, name))
        continue
    end
    keys{end + 1} = name; %#ok<AGROW>
end
end

function plot_region_residual_spectra(records, frequency, spectral_method, ...
        ranges, spectrum_regions, plot_style)
if isempty(records) || isempty(frequency)
    fprintf('No channels were available for the region spectra.\n');
    return
end
range_label = periodic_range_label(ranges);
region_names = selected_regions(spectrum_regions);
region_names = region_names(~strcmp(region_names, 'Other'));
region_names = unique(region_names, 'stable');
region_names = region_names(:).';
if isempty(region_names)
    fprintf('spectrum_regions is empty, so the region spectra are skipped.\n');
    return
end
[method_names, method_ids] = selected_spectral_methods(spectral_method);
for slot_n = 1:numel(method_ids)
    method_n = method_ids(slot_n);
    method_records = records(reshape([records.method_id], [], 1) == method_n);
    if isempty(method_records)
        continue
    end
    fprintf(['Region spectra | %s | aperiodic fit and %s removed | ' ...
        '%d channels\n'], method_names{slot_n}, range_label, ...
        numel(method_records));
    [fits, center_fit] = fit_region_spectrum_means( ...
        method_records, frequency, region_names);
    effect = fit_region_effect_by_frequency( ...
        method_records, frequency, region_names);
    draw_region_spectrum_panels(frequency, region_names, fits, center_fit, ...
        method_names{slot_n}, range_label, plot_style);
    draw_region_effect_figure(frequency, region_names, effect, ...
        method_names{slot_n}, range_label);
    report_region_effect(method_names{slot_n}, region_names, effect);
end
end

function [fits, center_fit] = fit_region_spectrum_means( ...
        records, frequency, region_names)
n_regions = numel(region_names);
fits = repmat(empty_spectrum_fit(numel(frequency)), 1, n_regions);
for region_n = 1:n_regions
    subset = records(record_in_region(records, region_names{region_n}));
    fprintf('  %s\n', region_display_name(region_names{region_n}));
    fits(region_n) = fit_spectrum_intercepts(subset, frequency);
end
center_records = records(reshape([records.is_center], [], 1));
fprintf('  Central LPAG\n');
center_fit = fit_spectrum_intercepts(center_records, frequency);
end

function fit = empty_spectrum_fit(n_frequency)
fit = struct( ...
    'intercept', nan(1, n_frequency), ...
    'lower', nan(1, n_frequency), ...
    'upper', nan(1, n_frequency), ...
    'n_channels', 0, ...
    'n_sessions', 0, ...
    'formula', '', ...
    'n_fit', 0);
end

function fit = fit_spectrum_intercepts(records, frequency)
n_frequency = numel(frequency);
fit = empty_spectrum_fit(n_frequency);
if isempty(records)
    fprintf('    no channels\n');
    return
end
[power, session_key, channel_key] = record_columns(records);
[formula, ~] = nested_formula(session_key, channel_key, false);
fit.formula = formula;
fit.n_channels = size(power, 1);
fit.n_sessions = numel(unique(session_key));
fprintf('    %d channels in %d sessions | %s\n', ...
    fit.n_channels, fit.n_sessions, formula);
if isempty(strfind(formula, 'Session:Channel')) %#ok<STREMP>
    fprintf(['    Each channel occurs once in a session, so the channel ' ...
        'term is the residual and the random effect is the session.\n']);
end
if fit.n_sessions < 2
    fprintf('    fewer than 2 sessions, so the mixed model is skipped\n');
    return
end
for frequency_n = 1:n_frequency
    [intercept, lower, upper] = fit_one_intercept( ...
        power(:, frequency_n), session_key, channel_key, formula);
    fit.intercept(frequency_n) = intercept;
    fit.lower(frequency_n) = lower;
    fit.upper(frequency_n) = upper;
    if isfinite(intercept)
        fit.n_fit = fit.n_fit + 1;
    end
end
fprintf('    fitted %d of %d frequencies\n', fit.n_fit, n_frequency);
end

function [formula, nested] = nested_formula(session_key, channel_key, with_region)
% One row per channel leaves Session:Channel identical to the residual.
% The fitted model then keeps every channel nested in its session through
% the session random intercept. A repeated channel inside one session
% uses the explicit nested term.
pair_key = strcat(session_key(:), '|', channel_key(:));
nested = numel(unique(pair_key)) < numel(pair_key);
if with_region
    if nested
        formula = 'Power ~ Region + (1|Session) + (1|Session:Channel)';
    else
        formula = 'Power ~ Region + (1|Session)';
    end
else
    if nested
        formula = 'Power ~ 1 + (1|Session) + (1|Session:Channel)';
    else
        formula = 'Power ~ 1 + (1|Session)';
    end
end
end

function [intercept, lower, upper] = fit_one_intercept( ...
        power, session_key, channel_key, formula)
intercept = NaN;
lower = NaN;
upper = NaN;
[intercept, lower, upper] = fit_intercept_table( ...
    power, session_key, channel_key, formula);
if isfinite(intercept)
    return
end
if ~isempty(strfind(formula, 'Session:Channel')) %#ok<STREMP>
    [intercept, lower, upper] = fit_intercept_table( ...
        power, session_key, channel_key, 'Power ~ 1 + (1|Session)');
end
end

function [intercept, lower, upper] = fit_intercept_table( ...
        power, session_key, channel_key, formula)
intercept = NaN;
lower = NaN;
upper = NaN;
power = power(:);
keep = isfinite(power);
if nnz(keep) < 2 || numel(unique(session_key(keep))) < 2
    return
end
model_table = table(power(keep), ...
    categorical(session_key(keep)), categorical(channel_key(keep)), ...
    'VariableNames', {'Power', 'Session', 'Channel'});
warning_state = warning;
warning('off', 'all');
try
    lme = fitlme(model_table, formula, 'FitMethod', 'REML');
    try
        [beta, ~, stats] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
        intercept = beta(1);
        try
            lower = stats.Lower(1);
            upper = stats.Upper(1);
        catch
            [lower, upper] = interval_from_se( ...
                intercept, stats.SE(1), stats.DF(1));
        end
    catch
        beta = fixedEffects(lme);
        intercept = beta(1);
        se = sqrt(max(lme.CoefficientCovariance(1, 1), 0));
        [lower, upper] = interval_from_se( ...
            intercept, se, max(lme.DFE, 1));
    end
catch
    intercept = NaN;
    lower = NaN;
    upper = NaN;
end
warning(warning_state);
end

function mask = record_in_region(records, region_name)
mask = false(numel(records), 1);
for record_n = 1:numel(records)
    mask(record_n) = any(strcmp(records(record_n).regions, region_name));
end
end

function [power, session_key, channel_key] = record_columns(records)
power = vertcat(records.spectrum);
n_records = numel(records);
session_key = cell(n_records, 1);
channel_key = cell(n_records, 1);
for record_n = 1:n_records
    session_key{record_n} = [records(record_n).animal, '|', ...
        sprintf('%.12g', records(record_n).session)];
    channel_key{record_n} = nested_channel_key( ...
        records(record_n).shank, records(record_n).channel);
end
end

function key = nested_channel_key(shank, channel)
key = [number_key(shank), '|', number_key(channel)];
end

function key = number_key(value)
if ~isfinite(value)
    key = 'NaN';
else
    key = sprintf('%.12g', value);
end
end

function effect = fit_region_effect_by_frequency(records, frequency, region_names)
n_frequency = numel(frequency);
n_regions = numel(region_names);
effect = struct( ...
    'emm', nan(n_regions, n_frequency), ...
    'p_omnibus', nan(1, n_frequency), ...
    'winner', zeros(1, n_frequency), ...
    'formula', '', ...
    'n_fit', 0);
[power, region_id, session_key, channel_key] = stack_region_rows( ...
    records, region_names);
if isempty(power)
    fprintf('  Region test: no channels in the selected regions\n');
    return
end
[formula, ~] = nested_formula(session_key, channel_key, true);
effect.formula = formula;
fprintf('  Region test | %d channel rows in %d sessions | %s\n', ...
    size(power, 1), numel(unique(session_key)), formula);
if isempty(strfind(formula, 'Session:Channel')) %#ok<STREMP>
    fprintf(['  Each channel occurs once in a session, so the channel ' ...
        'term is the residual and the random effect is the session.\n']);
end
for frequency_n = 1:n_frequency
    [emm, p_value, winner] = fit_one_region_test( ...
        power(:, frequency_n), region_id, region_names, ...
        session_key, channel_key, formula);
    effect.emm(:, frequency_n) = emm;
    effect.p_omnibus(frequency_n) = p_value;
    effect.winner(frequency_n) = winner;
    if isfinite(p_value)
        effect.n_fit = effect.n_fit + 1;
    end
end
fprintf('  Region test fitted %d of %d frequencies\n', ...
    effect.n_fit, n_frequency);
end

function [power, region_id, session_key, channel_key] = stack_region_rows( ...
        records, region_names)
power_cells = cell(1, 0);
region_cells = cell(1, 0);
session_cells = cell(1, 0);
channel_cells = cell(1, 0);
for region_n = 1:numel(region_names)
    subset = records(record_in_region(records, region_names{region_n}));
    if isempty(subset)
        continue
    end
    [this_power, this_session, this_channel] = record_columns(subset);
    for row_n = 1:numel(this_channel)
        this_channel{row_n} = [this_channel{row_n}, '|', region_names{region_n}];
    end
    power_cells{end + 1} = this_power; %#ok<AGROW>
    region_cells{end + 1} = repmat(region_n, size(this_power, 1), 1); %#ok<AGROW>
    session_cells{end + 1} = this_session; %#ok<AGROW>
    channel_cells{end + 1} = this_channel; %#ok<AGROW>
end
if isempty(power_cells)
    power = [];
    region_id = [];
    session_key = {};
    channel_key = {};
    return
end
power = vertcat(power_cells{:});
region_id = vertcat(region_cells{:});
session_key = vertcat(session_cells{:});
channel_key = vertcat(channel_cells{:});
end

function [emm, p_value, winner] = fit_one_region_test( ...
        power, region_id, region_names, session_key, channel_key, formula)
n_regions = numel(region_names);
emm = nan(n_regions, 1);
p_value = NaN;
winner = 0;
power = power(:);
keep = isfinite(power);
if nnz(keep) < 4
    return
end
power = power(keep);
region_id = region_id(keep);
session_key = session_key(keep);
channel_key = channel_key(keep);
present = unique(region_id, 'stable');
present = present(present >= 1 & present <= n_regions);
if numel(present) < 2 || numel(unique(session_key)) < 2
    return
end
present_names = region_names(present);
present_names = present_names(:).';
if any(strcmp(present_names, 'LPAG'))
    present_names = [{'LPAG'}, present_names(~strcmp(present_names, 'LPAG'))];
end
row_names = region_names(region_id);
row_names = row_names(:);
model_table = table(power, ...
    categorical(row_names, present_names), ...
    categorical(session_key), categorical(channel_key), ...
    'VariableNames', {'Power', 'Region', 'Session', 'Channel'});
warning_state = warning;
warning('off', 'all');
lme = [];
try
    lme = fitlme(model_table, formula, ...
        'FitMethod', 'REML', 'DummyVarCoding', 'reference');
catch
    if ~isempty(strfind(formula, 'Session:Channel')) %#ok<STREMP>
        try
            lme = fitlme(model_table, 'Power ~ Region + (1|Session)', ...
                'FitMethod', 'REML', 'DummyVarCoding', 'reference');
        catch
            lme = [];
        end
    end
end
warning(warning_state);
if isempty(lme)
    return
end
try
    stats = anova(lme, 'DFMethod', 'satterthwaite');
    [~, ~, ~, p_value] = region_anova_row(stats);
catch
    try
        stats = anova(lme);
        [~, ~, ~, p_value] = region_anova_row(stats);
    catch
        p_value = NaN;
    end
end
coef_names = coefficient_names(lme);
try
    beta = fixedEffects(lme);
catch
    return
end
reference_name = present_names{1};
for region_n = 1:n_regions
    if ~any(present == region_n)
        continue
    end
    L = level_row(coef_names, region_names{region_n}, reference_name);
    if any(~isfinite(L)) || numel(L) ~= numel(beta)
        continue
    end
    emm(region_n) = L * beta(:);
end
finite_region = find(isfinite(emm));
if isempty(finite_region) || ~isfinite(p_value) || p_value >= 0.05
    return
end
[~, best] = max(emm(finite_region));
winner = finite_region(best);
end

function draw_region_spectrum_panels(frequency, region_names, fits, ...
        center_fit, method_name, range_label, plot_style)
n_regions = numel(region_names);
figure('Name', sprintf('Region spectra | %s | %s removed', ...
    method_name, range_label), 'Color', 'w', ...
    'Position', [70 30 860 180 * n_regions]);
y_limits = spectrum_panel_limits(fits, center_fit);
for region_n = 1:n_regions
    subplot(n_regions, 1, region_n)
    hold on
    yline(0, ':', 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
    this_fit = fits(region_n);
    draw_spectrum_interval(frequency, this_fit.intercept, ...
        this_fit.lower, this_fit.upper, ...
        region_color(region_names{region_n}), plot_style.ci_alpha);
    center_line = plot(frequency, center_line_values(center_fit.intercept), ...
        'k:', 'LineWidth', 1.3);
    region_line = plot(nan, nan, '-', ...
        'Color', region_color(region_names{region_n}), 'LineWidth', 1.6);
    ci_patch = plot(nan, nan, 's', ...
        'MarkerFaceColor', region_color(region_names{region_n}), ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 8);
    xlim([frequency(1), frequency(end)])
    if all(isfinite(y_limits))
        ylim(y_limits)
    end
    ylabel('log10 power')
    title(sprintf('%s | %d sessions', ...
        region_display_name(region_names{region_n}), this_fit.n_sessions))
    if region_n == 1
        legend([region_line, ci_patch, center_line], ...
            {'Region estimate', '95% CI', 'Central LPAG'}, ...
            'Location', 'northeast', 'Box', 'off');
    end
    if region_n == n_regions
        xlabel('Frequency (Hz)')
    else
        set(gca, 'XTickLabel', [])
    end
    box off
end
sgtitle({sprintf('%s | aperiodic fit and %s removed', method_name, range_label), ...
    'Channels nested in session. Dotted line is central LPAG.'}, ...
    'FontSize', 10)
drawnow
end

function values = center_line_values(intercept)
values = intercept(:).';
values(~isfinite(values)) = NaN;
end

function limits = spectrum_panel_limits(fits, center_fit)
values = center_fit.intercept(:);
for region_n = 1:numel(fits)
    values = [values; fits(region_n).intercept(:); ...
        fits(region_n).lower(:); fits(region_n).upper(:)]; %#ok<AGROW>
end
values = values(isfinite(values));
if numel(values) < 2
    limits = [NaN NaN];
    return
end
low = prctile(values, 1);
high = prctile(values, 99);
if low == high
    low = low - 0.1;
    high = high + 0.1;
end
pad = 0.08 * (high - low);
limits = [low - pad, high + pad];
end

function draw_spectrum_interval(frequency, intercept, lower, upper, color, ci_alpha)
frequency = frequency(:);
intercept = intercept(:);
lower = lower(:);
upper = upper(:);
finite = isfinite(frequency) & isfinite(intercept) & ...
    isfinite(lower) & isfinite(upper);
if ~any(finite)
    return
end
    index = find(finite);
    run_break = [1; find(diff(index) > 1) + 1; numel(index) + 1];
    for run_n = 1:numel(run_break) - 1
        selected = index(run_break(run_n):run_break(run_n + 1) - 1);
        x = frequency(selected);
        if numel(selected) >= 2
            fill([x; flipud(x)], ...
                [lower(selected); flipud(upper(selected))], color, ...
                'FaceAlpha', ci_alpha, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        plot(x, intercept(selected), '-', 'Color', color, 'LineWidth', 1.6, ...
            'HandleVisibility', 'off');
    end
end

function draw_region_effect_figure(frequency, region_names, effect, ...
        method_name, range_label)
n_regions = numel(region_names);
labels = cellfun(@region_display_name, region_names, 'UniformOutput', false);
figure('Name', sprintf('Region effect by frequency | %s | %s removed', ...
    method_name, range_label), ...
    'Color', 'w', 'Position', [100 40 960 760]);
heatmap_axis = subplot(3, 1, 1);
finite_values = effect.emm(isfinite(effect.emm));
imagesc(heatmap_axis, frequency, 1:n_regions, effect.emm);
set(heatmap_axis, 'YTick', 1:n_regions, 'YTickLabel', labels, 'YDir', 'reverse');
axis(heatmap_axis, 'tight')
colormap(heatmap_axis, diverging_or_sequential('raw_aperiodic_subtracted'))
color_limits = spectrum_color_limits(finite_values, 'raw_aperiodic_subtracted');
if all(isfinite(color_limits))
    set(heatmap_axis, 'CLim', color_limits)
end
colorbar(heatmap_axis)
ylabel(heatmap_axis, 'Region')
title(heatmap_axis, 'Estimated log10 power')
box(heatmap_axis, 'off')

subplot(3, 1, 2)
p_values = effect.p_omnibus;
p_values(p_values <= 0) = realmin;
plot(frequency, -log10(p_values), 'k', 'LineWidth', 1.2)
hold on
yline(-log10(0.05), '--', 'p = 0.05', 'LabelHorizontalAlignment', 'left');
xlim([frequency(1), frequency(end)])
ylabel('-log10 p')
title('Test of region at each frequency')
box off

subplot(3, 1, 3)
hold on
step = median(diff(frequency));
if ~isfinite(step) || step <= 0
    step = 1;
end
for frequency_n = 1:numel(frequency)
    x0 = frequency(frequency_n) - step / 2;
    winner = effect.winner(frequency_n);
    if winner >= 1 && winner <= n_regions
        color = region_color(region_names{winner});
    else
        color = [0.86 0.86 0.86];
    end
    patch([x0, x0 + step, x0 + step, x0], [0, 0, 1, 1], color, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end
guides = gobjects(n_regions + 1, 1);
guide_labels = [labels, {'Region effect p >= 0.05'}];
for region_n = 1:n_regions
    guides(region_n) = plot(nan, nan, 's', ...
        'MarkerFaceColor', region_color(region_names{region_n}), ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 8);
end
guides(end) = plot(nan, nan, 's', ...
    'MarkerFaceColor', [0.86 0.86 0.86], ...
    'MarkerEdgeColor', 'none', 'MarkerSize', 8);
legend(guides, guide_labels, 'Location', 'eastoutside', 'Box', 'off')
xlim([frequency(1) - step / 2, frequency(end) + step / 2])
ylim([0, 1])
set(gca, 'YTick', [])
xlabel('Frequency (Hz)')
title('Strongest region where regions differ')
box off
sgtitle({sprintf('%s | regional effect by frequency', method_name), ...
    sprintf(['%s removed with the aperiodic fit. ' ...
    'Gray frequencies have no region effect at p < 0.05.'], range_label)}, ...
    'FontSize', 10)
drawnow
end

function report_region_effect(method_name, region_names, effect)
p_values = effect.p_omnibus;
significant = isfinite(p_values) & p_values < 0.05;
q_values = benjamini_hochberg(p_values);
fdr = isfinite(q_values) & q_values < 0.05;
fprintf(['%s | region effect at %d frequencies, %d with p < 0.05, ' ...
    '%d with FDR q < 0.05\n'], method_name, effect.n_fit, ...
    nnz(significant), nnz(fdr));
if ~any(significant)
    return
end
for region_n = 1:numel(region_names)
    fprintf('  %s is strongest at %d of those frequencies\n', ...
        region_display_name(region_names{region_n}), ...
        nnz(effect.winner == region_n));
end
end

function q_values = benjamini_hochberg(p_values)
q_values = nan(size(p_values));
p_column = p_values(:);
valid = isfinite(p_column);
if ~any(valid)
    return
end
[sorted, order] = sort(p_column(valid));
n_tests = numel(sorted);
q_sorted = sorted(:) .* n_tests ./ (1:n_tests).';
q_sorted = min(q_sorted, 1);
for rank_n = n_tests - 1:-1:1
    q_sorted(rank_n) = min(q_sorted(rank_n), q_sorted(rank_n + 1));
end
q_valid = nan(n_tests, 1);
q_valid(order) = q_sorted;
q_column = nan(numel(p_column), 1);
q_column(valid) = q_valid;
q_values(:) = q_column;
end
