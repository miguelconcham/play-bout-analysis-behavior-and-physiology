%% Analyse_play_session_probe_spectra
% Single-session view of the channel-wise periodic/aperiodic spectra.
%
% The image is the full spectrum at every channel. The kilosort y coordinate
% is distance to the probe tip, not depth: larger values are further from the
% tip (upper), smaller values are closer to the tip (deeper).
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
% The matching distance figure fits a mixed model of that power on
% distance. Every channel stays an observation, and sessions are nested
% in the animal. The line is the fixed effect. Spearman rho and its
% p-value use the same channels.
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
% are one figure. The 95% confidence intervals for the three track classes
% are a separate figure, drawn when at least two tracks are available.
%
% Separate figures center each shank on the regions in regions_to_display.
% By default each ventral or dorsal curve is z-scored within the shank and
% set to 0 at its first channel. Set region_distance_scale to 'raw' for the
% saved power. Each panel reports a Spearman correlation of distance and
% power, with its p-value, after the same outlier screens. The center is
% the midpoint of that region's labeled span on the shank. Delta, theta,
% and gamma are drawn against distance from that center, once toward the
% tip (ventral) and once toward the brain surface (dorsal). Individual
% shanks are one figure. The mixed-model fit of power on distance is another.

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

% log10 power minus the fitted aperiodic component. Switch to
% 'log10_mean_power' to draw the session-mean spectrum itself.
spectrum_field = 'raw_aperiodic_subtracted';

delta_band = [1 5];
theta_band = [6 12];
gamma_band = [30 45];

% Figures use the Fourier spectra. Welch is available when requested:
%   spectral_method = {'Welch'};
%   spectral_method = {'Welch', 'Fourier'};
spectral_method = {'Fourier'};

% Regions drawn on the scatter, distance-model, region-comparison,
% and region-distance figures. SupCol is the superior colliculus. Add 'DR'
% to include it:
%   regions_to_display = {'LPAG', 'VLPAG', 'SupCol', 'DR'};
regions_to_display = {'LPAG', 'VLPAG', 'SupCol'};

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
            artifact_trend_sd)]; %#ok<AGROW>
        placed = [placed; placed_from_session(channel_table, ...
            welch_spectrum, fourier_spectrum, frequency, ...
            delta_band, theta_band, gamma_band, session_label, ...
            artifact_trend_sd, probe_spectra.animal_name, session_n)]; %#ok<AGROW>
        power_pieces{end + 1} = lpag_distance_table(channel_table, ...
            welch_spectrum, fourier_spectrum, frequency, ...
            delta_band, theta_band, gamma_band, ...
            probe_spectra.animal_name, session_n, spectrum_field, ...
            artifact_trend_sd); %#ok<SAGROW>

        if plot_each_session
            plot_depth_spectrum(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, spectrum_field, ...
                session_label, area_borders, spectral_method);
            plot_band_profiles(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, ...
                delta_band, theta_band, gamma_band, ...
                spectrum_field, session_label, area_borders, spectral_method, ...
                artifact_trend_sd);
            drawnow
        end
    end
end

if isempty(placed)
    error('No probe channels were found in the selected recordings.');
end
power_table = write_lpag_distance_table(power_pieces, saving_folder, ...
    spectrum_field, region_border_bin_um);
placed = apply_screened_power_to_placed(placed, power_table);
plot_normalized_power_scatters(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display);
plot_normalized_power_intervals(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display);
plot_region_power_comparison(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display);
if plot_proximal_distal
    plot_region_half_comparison(power_table, delta_band, theta_band, gamma_band, ...
        spectral_method, regions_to_display);
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
        delta_band, theta_band, gamma_band, spectral_method);
    plot_aligned_intervals(tracks, spectrum_field, ...
        delta_band, theta_band, gamma_band, spectral_method);
end
plot_region_distances(placed, spectrum_field, ...
    delta_band, theta_band, gamma_band, spectral_method, ...
    regions_to_display, region_distance_scale);

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

function plot_band_profiles(frequency, channel_table, ...
        welch_spectrum, fourier_spectrum, ...
        delta_band, theta_band, gamma_band, spectrum_field, ...
        session_label, area_borders, spectral_method, artifact_trend_sd)
probe_groups = probe_groups_by_depth(channel_table);
n_groups = numel(probe_groups);
spectra = {welch_spectrum, fourier_spectrum};
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
band_colors = [0.12 0.35 0.75; 0.85 0.40 0.10; 0.45 0.62 0.20];

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
        line_handles = gobjects(3, 1);
        for band_n = 1:3
            line_handles(band_n) = plot(band_power(:, band_n), ...
                this_group.depth_um, ...
                'Color', band_colors(band_n, :), 'LineWidth', 1.4);
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
            legend(line_handles, band_names, 'Location', 'best')
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
        delta_band, theta_band, gamma_band, spectral_method, regions_to_display)
if isempty(power_table)
    return
end
[method_names, ~] = selected_spectral_methods(spectral_method);
band_fields = { ...
    'DeltaPower_zWithinShank', ...
    'ThetaPower_zWithinShank', ...
    'GammaPower_zWithinShank'};
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
region_keys = region_keys_for_rows(power_table.Area);
region_names = selected_regions(regions_to_display);
figure_n = 0;
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    distance_um = power_table.DistanceToMidLPAG_um;
    for band_n = 1:numel(band_fields)
        figure_n = figure_n + 1;
        power = power_table.(band_fields{band_n});
        figure('Name', sprintf('%s | %s vs distance to mid-LPAG', ...
            method_names{method_n}, band_names{band_n}), ...
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
            scatter(distance_um(keep), power(keep), 22, this_color, 'filled', ...
                'MarkerFaceAlpha', 0.22, 'MarkerEdgeAlpha', 0.22, ...
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
        title(sprintf('%s | %s', method_names{method_n}, band_names{band_n}))
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
        delta_band, theta_band, gamma_band, spectral_method, regions_to_display)
if isempty(power_table)
    return
end
[method_names, ~] = selected_spectral_methods(spectral_method);
band_fields = { ...
    'DeltaPower_zWithinShank', ...
    'ThetaPower_zWithinShank', ...
    'GammaPower_zWithinShank'};
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
region_keys = region_keys_for_rows(power_table.Area);
region_names = selected_regions(regions_to_display);
region_names = region_names(~strcmp(region_names, 'Other'));
animal = cellstr(string(power_table.Animal));
session = power_table.Session;
figure_n = 0;
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    distance_um = power_table.DistanceToMidLPAG_um;
    for band_n = 1:numel(band_fields)
        figure_n = figure_n + 1;
        power = power_table.(band_fields{band_n});
        figure('Name', sprintf('%s | %s mixed model vs distance to mid-LPAG', ...
            method_names{method_n}, band_names{band_n}), ...
            'Color', 'w', 'Position', [80 + 24 * figure_n, 100, 980, 560]);
        hold on
        yline(0, ':k', 'HandleVisibility', 'off');
        xline(0, '--k', 'mid-LPAG', 'HandleVisibility', 'off', ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 8);
        legend_handles = gobjects(0, 1);
        legend_labels = {};
        x_span = [];
        for region_n = 1:numel(region_names)
            this_region = region_names{region_n};
            keep = in_method & region_row_has(region_keys, this_region) & ...
                isfinite(distance_um) & isfinite(power);
            if ~any(keep)
                continue
            end
            this_distance = distance_um(keep);
            this_power = power(keep);
            result = fit_distance_mixed_model(this_distance, this_power, ...
                animal(keep), session(keep));
            [rho, p_rho] = rank_correlation(this_distance, this_power);
            report_distance_model(sprintf('%s | %s | %s', ...
                method_names{method_n}, band_names{band_n}, ...
                region_display_name(this_region)), result, rho, p_rho);
            this_color = region_color(this_region);
            if result.ok
                draw_distance_model_line(this_distance, result, this_color);
                x_span = [x_span; this_distance(:)]; %#ok<AGROW>
            end
            legend_handles(end + 1, 1) = plot(nan, nan, '-', ...
                'Color', this_color, 'LineWidth', 2); %#ok<AGROW>
            legend_labels{end + 1} = distance_model_legend( ...
                this_region, result, rho, p_rho); %#ok<AGROW>
        end
        if numel(x_span) > 1
            pad = 0.02 * max(max(x_span) - min(x_span), 1);
            xlim([min(x_span) - pad, max(x_span) + pad]);
        end
        xlabel({'Distance to mid-LPAG (µm)', ...
            'positive = dorsal, farther from the tip'})
        ylabel({'Power z-scored within shank', '0 = mid-LPAG channel'})
        title(sprintf('%s | %s | mixed model across animals', ...
            method_names{method_n}, band_names{band_n}))
        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, ...
                'Location', 'eastoutside', 'FontSize', 8, 'Interpreter', 'tex')
        end
        box off
    end
end
fprintf(['Distance fits are mixed models of channel power on distance, ' ...
    'with sessions nested in the animal. Spearman rho and p use the same ' ...
    'channels, after extreme values and uncommon border distances are removed.\n']);
end

function result = fit_distance_mixed_model(distance, power, animal, session)
result = struct( ...
    'ok', false, ...
    'intercept', NaN, ...
    'slope', NaN, ...
    'p_slope', NaN, ...
    'n_channels', 0, ...
    'n_animals', 0, ...
    'n_sessions', 0, ...
    'formula', '', ...
    'message', '');
distance = distance(:);
power = power(:);
animal = cellstr(string(animal(:)));
session = cellstr(string(session(:)));
keep = isfinite(distance) & isfinite(power);
distance = distance(keep);
power = power(keep);
animal = animal(keep);
session = session(keep);
result.n_channels = numel(power);
if result.n_channels == 0
    result.message = 'No finite channels.';
    return
end
session_id = strcat(animal, '|', session);
result.n_animals = numel(unique(animal));
result.n_sessions = numel(unique(session_id));
if result.n_channels < 4 || numel(unique(distance)) < 2 || ...
        numel(unique(power)) < 2
    result.message = 'Not enough variation to fit distance.';
    return
end
formulas = {};
if result.n_animals >= 2 && result.n_sessions > result.n_animals
    formulas{end + 1} = 'Power ~ Distance + (1|Animal) + (1|Session)'; %#ok<AGROW>
end
if result.n_animals >= 2
    formulas{end + 1} = 'Power ~ Distance + (1|Animal)'; %#ok<AGROW>
end
if result.n_sessions >= 2
    formulas{end + 1} = 'Power ~ Distance + (1|Session)'; %#ok<AGROW>
end
if isempty(formulas)
    result.message = ['A mixed model needs repeated measures in more than ' ...
        'one animal or session.'];
    return
end
model_table = table(power, distance, categorical(animal), ...
    categorical(session_id), ...
    'VariableNames', {'Power', 'Distance', 'Animal', 'Session'});
last_error = '';
for formula_n = 1:numel(formulas)
    try
        lme = fitlme(model_table, formulas{formula_n}, 'FitMethod', 'REML');
    catch ME
        last_error = ME.message;
        continue
    end
    coef_names = coefficient_names(lme);
    slope_index = find(strcmp(coef_names, 'Distance'), 1);
    intercept_index = find(strcmp(coef_names, '(Intercept)'), 1);
    if isempty(slope_index) || isempty(intercept_index)
        last_error = 'The distance slope was not in the model.';
        continue
    end
    beta = fixedEffects(lme);
    contrast = zeros(1, numel(coef_names));
    contrast(slope_index) = 1;
    try
        p_slope = coefTest(lme, contrast, 0, 'DFMethod', 'satterthwaite');
    catch
        try
            p_slope = coefTest(lme, contrast, 0);
        catch ME
            last_error = ME.message;
            continue
        end
    end
    result.ok = isfinite(p_slope) && isfinite(beta(slope_index)) && ...
        isfinite(beta(intercept_index));
    if ~result.ok
        last_error = 'The distance slope could not be tested.';
        continue
    end
    result.intercept = beta(intercept_index);
    result.slope = beta(slope_index);
    result.p_slope = p_slope;
    result.formula = formulas{formula_n};
    result.message = '';
    return
end
if isempty(last_error)
    last_error = 'The distance model could not be fit.';
end
result.message = last_error;
end

function draw_distance_model_line(distance, result, color)
distance = distance(isfinite(distance));
if ~result.ok || numel(distance) < 2
    return
end
x_line = [min(distance), max(distance)];
if x_line(1) == x_line(2)
    return
end
plot(x_line, result.intercept + result.slope * x_line, '-', ...
    'Color', color, 'LineWidth', 2.0, 'HandleVisibility', 'off');
end

function label = distance_model_legend(region_name, result, rho, p_rho)
if result.ok
    slope_text = sprintf('slope p = %s, %d animals', ...
        format_p_value(result.p_slope), result.n_animals);
else
    slope_text = 'model unavailable';
end
if isfinite(rho)
    rho_text = sprintf('\\rho = %.2f, p = %s', rho, format_p_value(p_rho));
else
    rho_text = '\rho = n/a';
end
label = sprintf('%s | %s | %s', region_display_name(region_name), ...
    slope_text, rho_text);
end

function report_distance_model(context, result, rho, p_rho)
fprintf('%s\n', context);
if result.ok
    fprintf('  %s\n', result.formula);
    if contains(result.formula, 'Session')
        fprintf('  Session is the partner session within the animal.\n');
    end
    fprintf(['  slope = %s per µm, p = %s ' ...
        '(%d channels, %d animals, %d sessions)\n'], ...
        format_p_value(result.slope), format_p_value(result.p_slope), ...
        result.n_channels, result.n_animals, result.n_sessions);
else
    if isempty(result.message)
        detail = 'The distance model could not be fit.';
    else
        detail = result.message;
    end
    fprintf('  Mixed model was not fit: %s\n', detail);
end
if isfinite(rho)
    fprintf('  Spearman rho = %.2f, p = %s\n', rho, format_p_value(p_rho));
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
        delta_band, theta_band, gamma_band, spectral_method, regions_to_display)
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
    'Lines join one session. Brackets are mixed-model contrasts, Bonferroni.']);
end

function plot_region_half_comparison(power_table, ...
        delta_band, theta_band, gamma_band, spectral_method, regions_to_display)
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
    'Brackets are mixed-model contrasts, Bonferroni.']);
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
        region_keys, region_names, figure_name, title_name, note)
if numel(region_names) < 2
    fprintf(['Region comparison needs at least two regions in ' ...
        'regions_to_display.\n']);
    return
end
[method_names, ~] = selected_spectral_methods(spectral_method);
band_fields = { ...
    'DeltaPower_zWithinShank', ...
    'ThetaPower_zWithinShank', ...
    'GammaPower_zWithinShank'};
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    band_obs = cell(1, numel(band_fields));
    enough = false(1, numel(band_fields));
    for band_n = 1:numel(band_fields)
        band_obs{band_n} = nested_channel_observations(power_table, ...
            in_method, region_keys, region_names, ...
            power_table.(band_fields{band_n}));
        enough(band_n) = band_obs{band_n}.n_sessions >= 2 && ...
            numel(band_obs{band_n}.region_names) >= 2;
    end
    if ~any(enough)
        fprintf(['%s %s needs channels from at least two sessions in ' ...
            'at least two regions.\n'], method_names{method_n}, figure_name);
        continue
    end
    figure('Name', sprintf('%s | %s', method_names{method_n}, figure_name), ...
        'Color', 'w', ...
        'Position', [110, 80, max(1080, 160 * numel(region_names)), 560]);
    for band_n = 1:numel(band_fields)
        subplot(1, numel(band_fields), band_n)
        obs = band_obs{band_n};
        if ~enough(band_n)
            axis off
            text(0.5, 0.5, sprintf('%s\nfewer than 2 sessions', ...
                band_names{band_n}), 'HorizontalAlignment', 'center');
            fprintf(['%s | %s | %s: fewer than 2 sessions with channels in ' ...
                'at least two regions.\n'], method_names{method_n}, ...
                figure_name, band_names{band_n});
            continue
        end
        result = fit_nested_region_model(obs);
        report_nested_region_tests(method_names{method_n}, band_names{band_n}, ...
            region_names, obs, result);
        draw_nested_region_comparison(obs, result, band_names{band_n});
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

function draw_nested_region_comparison(obs, result, band_name)
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
scatter(obs.region_id + jitter(obs.session_id), obs.power, 16, colors, ...
    'filled', 'MarkerFaceAlpha', 0.28, 'MarkerEdgeColor', 'none', ...
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

function tracks = tracks_from_session(channel_table, welch_spectrum, ...
        fourier_spectrum, frequency, delta_band, theta_band, gamma_band, ...
        session_label, artifact_trend_sd)
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
    track.welch_spectrum = welch_rows;
    track.fourier_spectrum = fourier_rows;
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
        delta_band, theta_band, gamma_band, spectral_method)
palette = track_colors();
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
figure('Name', 'Profiles aligned on mid-LPAG', 'Color', 'w', ...
    'Position', [60 60 1100 720]);
for method_slot = 1:n_methods
    method_n = method_ids(method_slot);
    for band_n = 1:3
        subplot(n_methods, 3, (method_slot - 1) * 3 + band_n)
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
        if band_n == 1
            ylabel({'Distance to tip relative to mid-LPAG (µm)', ...
                'positive = upper, further from tip'})
        end
        title(sprintf('%s | %s', method_names{method_slot}, band_names{band_n}))
        if method_slot == 1 && band_n == 3
            legend(guide, guide_labels, 'Location', 'best', 'FontSize', 7)
        end
        box off
    end
end
sgtitle('Sessions stacked on the mid-LPAG channel', 'Interpreter', 'none')
end

function plot_aligned_intervals(tracks, spectrum_field, ...
        delta_band, theta_band, gamma_band, spectral_method)
palette = track_colors();
[method_names, method_ids] = selected_spectral_methods(spectral_method);
n_methods = numel(method_ids);
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
figure('Name', 'Confidence intervals aligned on mid-LPAG', 'Color', 'w', ...
    'Position', [100 100 1100 720]);
for method_slot = 1:n_methods
    method_n = method_ids(method_slot);
    for band_n = 1:3
        subplot(n_methods, 3, (method_slot - 1) * 3 + band_n)
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
        if band_n == 1
            ylabel({'Distance to tip relative to mid-LPAG (µm)', ...
                'positive = upper, further from tip'})
        end
        title(sprintf('%s | %s', method_names{method_slot}, band_names{band_n}))
        if method_slot == 1 && band_n == 3 && ~isempty(ci_handles)
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
    'Mean of tracks in each class, centered on mid-LPAG'}, ...
    'Interpreter', 'none')
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
        regions_to_display, region_distance_scale)
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
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
band_colors = [0.12 0.35 0.75; 0.85 0.40 0.10; 0.45 0.62 0.20];
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
            figure_name = sprintf('%s mixed model vs distance to region center', ...
                method_names{method_slot});
            figure_position = [90 40 980 max(480, 240 * n_regions)];
            figure_note = region_distance_note(region_distance_scale, true);
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
                        draw_distance_track(distance, ...
                            power_for_method(welch_power, fourier_power, method_n), ...
                            band_colors);
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
                slope_p = nan(1, 3);
                model_ok = false;
                if strcmp(panel_kind, 'model')
                    for band_n = 1:3
                        [point_distance, point_power, point_animal, point_session] = ...
                            side_track_points(side_tracks, method_n, band_n);
                        result = fit_distance_mixed_model(point_distance, ...
                            point_power, point_animal, point_session);
                        [rho, p_rho] = rank_correlation(point_distance, point_power);
                        report_distance_model(sprintf('%s | %s | %s | %s', ...
                            method_names{method_slot}, band_names{band_n}, ...
                            region_display_name(region_names{region_n}), ...
                            direction_titles{direction_n}), ...
                            result, rho, p_rho);
                        slope_p(band_n) = result.p_slope;
                        if ~result.ok
                            continue
                        end
                        model_ok = true;
                        finite_distance = point_distance( ...
                            isfinite(point_distance) & isfinite(point_power));
                        draw_distance_model_line(finite_distance, result, ...
                            band_colors(band_n, :));
                        if ~isempty(finite_distance)
                            x_max = max(x_max, max(finite_distance));
                        end
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
                xlabel(direction_xlabels{direction_n})
                if direction_n == 1
                    ylabel(region_distance_ylabel(spectrum_field, region_distance_scale))
                end
                title(sprintf('%s | %s (n=%d)', ...
                    region_display_name(region_names{region_n}), ...
                    direction_titles{direction_n}, numel(side_tracks)))
                draw_rank_correlations(side_tracks, method_n, band_colors, ...
                    strcmp(panel_kind, 'model'), slope_p);
                if strcmp(panel_kind, 'model') && ~model_ok
                    text(0.5, 0.5, 'Mixed model was not fit', ...
                        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                        'FontSize', 8);
                end
                if ~legend_drawn
                    guide = gobjects(3, 1);
                    for band_n = 1:3
                        guide(band_n) = plot(nan, nan, ...
                            'Color', band_colors(band_n, :), 'LineWidth', 2);
                    end
                    legend(guide, band_names, 'Location', 'best', 'FontSize', 7)
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
    note = [note, [' The line is a mixed model of power on distance, ' ...
        'with sessions nested in the animal.']];
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

function draw_rank_correlations(side_tracks, method_n, band_colors, ...
        with_slope, slope_p)
if nargin < 4
    with_slope = false;
    slope_p = nan(1, 3);
end
if with_slope
    gap = 0.16;
else
    gap = 0.09;
end
for band_n = 1:3
    [rho, p_value] = panel_rank_correlation(side_tracks, method_n, band_n);
    if isfinite(rho)
        label = sprintf('\\rho = %.2f, p = %s', rho, format_p_value(p_value));
    else
        label = '\rho = n/a';
    end
    if with_slope
        label = sprintf('%s\nslope p = %s', label, format_p_value(slope_p(band_n)));
    end
    text(0.98, 0.98 - (band_n - 1) * gap, label, ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'Color', band_colors(band_n, :), ...
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

function draw_distance_track(distance, power, band_colors)
[distance, order] = sort(distance(:));
power = power(order, :);
for band_n = 1:size(power, 2)
    this_color = band_colors(band_n, :);
    this_power = power(:, band_n);
    finite = isfinite(distance) & isfinite(this_power);
    if nnz(finite) == 1
        plot(distance(finite), this_power(finite), 'o', 'Color', this_color, ...
            'MarkerSize', 4, 'HandleVisibility', 'off');
        continue
    end
    plot(distance, this_power, 'Color', this_color, 'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end
end
