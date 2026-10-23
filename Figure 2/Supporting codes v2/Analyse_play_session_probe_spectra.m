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
% deviation. Each band and spectral method then has its own scatter of that
% z-scored power against distance to mid-LPAG, with transparent markers.
% A separate figure compares the selected regions with repeated measures
% on that same z-scored power: one bar per region, the shank means as
% connected points, and sign-rank brackets with Bonferroni correction.
% Figures use Fourier unless spectral_method is changed. The scatter,
% confidence-interval, region-comparison, and region-distance figures
% draw regions_to_display, which starts as LPAG, VLPAG, and SupCol. PAG sections,
% an upper side that passes through SupCol, and a deeper side that passes
% through DR are colored. Other territory stays gray. Those colored tracks
% are one figure. The 95% confidence intervals for the three track classes
% are a separate figure, drawn when at least two tracks are available.
%
% Separate figures center each shank on the regions in regions_to_display.
% By default each ventral or dorsal curve is z-scored within the shank and
% set to 0 at its first channel. Set region_distance_scale to 'raw' for the
% saved power. Each panel reports a Spearman correlation of distance and power.
% The center is the midpoint of that region's
% labeled span on the shank. Delta, theta, and gamma are drawn against
% distance from that center, once toward the tip (ventral) and once toward
% the brain surface (dorsal). Individual shanks and the confidence intervals
% are separate figures.

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

% Regions drawn on the scatter, confidence-interval, region-comparison,
% and region-distance figures. SupCol is the superior colliculus. Add 'DR'
% to include it:
%   regions_to_display = {'LPAG', 'VLPAG', 'SupCol', 'DR'};
regions_to_display = {'LPAG', 'VLPAG', 'SupCol'};

% Ventral and dorsal curves. 'normalized' z-scores each band within the
% shank and subtracts the power of the first channel on that side, so the
% curve starts at 0. 'raw' draws the saved band power.
region_distance_scale = 'normalized';

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
            delta_band, theta_band, gamma_band, session_label)]; %#ok<AGROW>
        placed = [placed; placed_from_session(channel_table, ...
            welch_spectrum, fourier_spectrum, frequency, ...
            delta_band, theta_band, gamma_band, session_label)]; %#ok<AGROW>
        power_pieces{end + 1} = lpag_distance_table(channel_table, ...
            welch_spectrum, fourier_spectrum, frequency, ...
            delta_band, theta_band, gamma_band, ...
            probe_spectra.animal_name, session_n, spectrum_field); %#ok<SAGROW>

        if plot_each_session
            plot_depth_spectrum(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, spectrum_field, ...
                session_label, area_borders, spectral_method);
            plot_band_profiles(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, ...
                delta_band, theta_band, gamma_band, ...
                spectrum_field, session_label, area_borders, spectral_method);
            drawnow
        end
    end
end

if isempty(placed)
    error('No probe channels were found in the selected recordings.');
end
power_table = write_lpag_distance_table(power_pieces, saving_folder, spectrum_field);
plot_normalized_power_scatters(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display);
plot_normalized_power_intervals(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display);
plot_region_power_comparison(power_table, delta_band, theta_band, gamma_band, ...
    spectral_method, regions_to_display);
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
        session_label, area_borders, spectral_method)
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
    'tip_um', [], ...
    'areas', {{}}, ...
    'welch_bands', [], ...
    'fourier_bands', []);
end

function placed = placed_from_session(channel_table, welch_spectrum, ...
        fourier_spectrum, frequency, delta_band, theta_band, gamma_band, ...
        session_label)
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
    item.tip_um = tip_um;
    item.areas = this_group.areas;
    item.welch_bands = [ ...
        mean_in_band(welch_rows, frequency, delta_band), ...
        mean_in_band(welch_rows, frequency, theta_band), ...
        mean_in_band(welch_rows, frequency, gamma_band)];
    item.fourier_bands = [ ...
        mean_in_band(fourier_rows, frequency, delta_band), ...
        mean_in_band(fourier_rows, frequency, theta_band), ...
        mean_in_band(fourier_rows, frequency, gamma_band)];
    placed(end + 1, 1) = item; %#ok<AGROW>
end
end

function power_table = lpag_distance_table(channel_table, welch_spectrum, ...
        fourier_spectrum, frequency, delta_band, theta_band, gamma_band, ...
        animal_name, session_n, spectrum_field)
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
if isempty(pieces)
    power_table = table();
    return
end
power_table = vertcat(pieces{:});
end

function power_table = write_lpag_distance_table(power_pieces, saving_folder, spectrum_field)
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
table_file = fullfile(saving_folder, sprintf( ...
    'play_session_probe_power_by_lpag_distance_%s.csv', spectrum_field));
writetable(power_table, table_file);
fprintf(['Wrote %d channel rows to %s\n' ...
    'Distance is positive dorsal (farther from the tip). ' ...
    'centeredAtLPAG subtracts the mid-LPAG channel. ' ...
    'zWithinShank divides that difference by the standard deviation ' ...
    'of the same band on that shank.\n'], height(power_table), table_file);
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
if strcmp(name, 'SC')
    label = 'SupCol';
else
    label = name;
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
distance_step_um = 5;
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
figure_n = 0;
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    distance_um = power_table.DistanceToMidLPAG_um;
    colored_row = false(height(power_table), 1);
    for region_n = 1:numel(region_names)
        colored_row = colored_row | region_row_has(region_keys, region_names{region_n});
    end
    colored_row = colored_row & in_method;
    grid = fine_distance_grid(distance_um(colored_row), distance_step_um);
    for band_n = 1:numel(band_fields)
        figure_n = figure_n + 1;
        power = power_table.(band_fields{band_n});
        figure('Name', sprintf('%s | %s 95%% CI vs distance to mid-LPAG', ...
            method_names{method_n}, band_names{band_n}), ...
            'Color', 'w', 'Position', [80 + 24 * figure_n, 100, 760, 560]);
        hold on
        yline(0, ':k', 'HandleVisibility', 'off');
        xline(0, '--k', 'mid-LPAG', 'HandleVisibility', 'off', ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 8);
        legend_handles = gobjects(0, 1);
        legend_labels = {};
        for region_n = 1:numel(region_names)
            this_region = region_names{region_n};
            in_region = in_method & region_row_has(region_keys, this_region);
            [mean_line, lower_line, upper_line, n_profiles] = ...
                interpolated_session_interval( ...
                power_table, in_region, distance_um, power, grid);
            if n_profiles < 2 || ~any(isfinite(mean_line))
                continue
            end
            this_color = region_color(this_region);
            draw_band_interval(grid, mean_line, lower_line, upper_line, this_color);
            legend_handles(end + 1, 1) = plot(nan, nan, '-', ...
                'Color', this_color, 'LineWidth', 2); %#ok<AGROW>
            legend_labels{end + 1} = sprintf('%s (n=%d shanks)', ...
                region_display_name(this_region), n_profiles); %#ok<AGROW>
        end
        if ~isempty(grid) && numel(grid) > 1
            xlim([grid(1), grid(end)]);
        end
        xlabel({'Distance to mid-LPAG (µm)', ...
            'positive = dorsal, farther from the tip'})
        ylabel({'Power z-scored within shank', '0 = mid-LPAG channel'})
        title(sprintf('%s | %s | 95%% CI on a %g µm grid', ...
            method_names{method_n}, band_names{band_n}, distance_step_um))
        if ~isempty(legend_handles)
            legend(legend_handles, legend_labels, ...
                'Location', 'eastoutside', 'FontSize', 8)
        end
        box off
    end
end
fprintf(['Region confidence intervals use a %g µm distance axis. ' ...
    'Each shank is interpolated on its own, and the interval is across shanks.\n'], ...
    distance_step_um);
end

function plot_region_power_comparison(power_table, ...
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
region_names = unique(region_names, 'stable');
if numel(region_names) < 2
    fprintf(['Region comparison needs at least two regions in ' ...
        'regions_to_display.\n']);
    return
end
region_list = strjoin(cellfun(@region_display_name, region_names, ...
    'UniformOutput', false), ', ');
for method_n = 1:numel(method_names)
    in_method = strcmp(power_table.Method, method_names{method_n});
    band_values = cell(1, numel(band_fields));
    for band_n = 1:numel(band_fields)
        band_values{band_n} = region_mean_matrix(power_table, in_method, ...
            region_keys, region_names, power_table.(band_fields{band_n}));
    end
    n_complete = cellfun(@(values) size(values, 1), band_values);
    if ~any(n_complete >= 2)
        fprintf(['%s region comparison needs at least two shanks that ' ...
            'contain %s.\n'], method_names{method_n}, region_list);
        continue
    end
    figure('Name', sprintf('%s | region comparison', method_names{method_n}), ...
        'Color', 'w', 'Position', [110, 80, 1080, 480]);
    for band_n = 1:numel(band_fields)
        subplot(1, numel(band_fields), band_n)
        if n_complete(band_n) < 2
            axis off
            text(0.5, 0.5, sprintf('%s\nfewer than 2 complete shanks', ...
                band_names{band_n}), 'HorizontalAlignment', 'center');
            fprintf(['%s | %s: fewer than 2 shanks contain %s.\n'], ...
                method_names{method_n}, band_names{band_n}, region_list);
            continue
        end
        draw_region_comparison(band_values{band_n}, region_names, ...
            band_names{band_n}, method_names{method_n}, region_list);
    end
    sgtitle({sprintf('%s | repeated-measures region comparison', ...
        method_names{method_n}), ...
        ['Mean z-scored power in each region (0 = mid-LPAG). ' ...
        'Lines join one shank. Brackets: sign-rank, Bonferroni.']}, ...
        'FontSize', 10)
end
end

function values = region_mean_matrix(power_table, row_mask, ...
        region_keys, region_names, power)
group_index = session_shank_groups(power_table, row_mask);
n_groups = numel(group_index);
n_regions = numel(region_names);
values = nan(n_groups, n_regions);
for group_n = 1:n_groups
    rows = group_index{group_n};
    these_keys = region_keys(rows);
    for region_n = 1:n_regions
        in_region = region_row_has(these_keys, region_names{region_n});
        samples = power(rows);
        samples = samples(in_region);
        samples = samples(isfinite(samples));
        if isempty(samples)
            continue
        end
        values(group_n, region_n) = mean(samples);
    end
end
values = values(all(isfinite(values), 2), :);
end

function draw_region_comparison(values, region_names, band_name, ...
        method_name, region_list)
n_shanks = size(values, 1);
n_regions = size(values, 2);
means = mean(values, 1);
standard_error = std(values, 0, 1) ./ sqrt(n_shanks);
[p_raw, p_adjusted, pairs] = paired_region_tests(values);
p_omnibus = friedman_p(values);
report_region_tests(method_name, band_name, region_list, region_names, ...
    n_shanks, p_omnibus, pairs, p_raw, p_adjusted);

hold on
yline(0, ':k', 'HandleVisibility', 'off');
jitter_span = 0.22;
jitter = ((1:n_shanks).' - (n_shanks + 1) / 2) / n_shanks * jitter_span;
for shank_n = 1:n_shanks
    plot((1:n_regions) + jitter(shank_n), values(shank_n, :), '-', ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 0.7, ...
        'HandleVisibility', 'off');
end
bar_half = 0.275;
for region_n = 1:n_regions
    fill(region_n + [-bar_half, bar_half, bar_half, -bar_half], ...
        [0, 0, means(region_n), means(region_n)], ...
        region_color(region_names{region_n}), ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
errorbar(1:n_regions, means, standard_error, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1, 'CapSize', 8, 'HandleVisibility', 'off');
for region_n = 1:n_regions
    scatter(region_n + jitter, values(:, region_n), 28, ...
        region_color(region_names{region_n}), 'filled', ...
        'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'MarkerEdgeAlpha', 0.45, 'HandleVisibility', 'off');
end
y_extent = [0; values(:); (means + standard_error).'; (means - standard_error).'];
y_extent = y_extent(isfinite(y_extent));
y_low = min(y_extent);
y_high = max(y_extent);
span = y_high - y_low;
if ~isfinite(span) || span < 0.5
    span = 0.5;
end
draw_significance_brackets(pairs, p_adjusted, y_high, span);
margin = 0.08 * span;
n_pairs = size(pairs, 1);
ylim([y_low - margin, y_high + margin + n_pairs * (0.12 * span) + 0.12 * span]);
xlim([0.4, n_regions + 0.6]);
set(gca, 'XTick', 1:n_regions, 'XTickLabel', ...
    cellfun(@region_display_name, region_names, 'UniformOutput', false));
ylabel({'Mean power z-scored within shank', '0 = mid-LPAG channel'})
if n_regions >= 3
    if isfinite(p_omnibus)
        omnibus_label = sprintf('Friedman p = %.3g', p_omnibus);
    else
        omnibus_label = 'Friedman p = n/a';
    end
    title(sprintf('%s | %s | n = %d shanks', band_name, omnibus_label, n_shanks))
else
    title(sprintf('%s | n = %d shanks', band_name, n_shanks))
end
box off
end

function [p_raw, p_adjusted, pairs] = paired_region_tests(values)
n_regions = size(values, 2);
pairs = nchoosek(1:n_regions, 2);
n_pairs = size(pairs, 1);
p_raw = nan(n_pairs, 1);
for pair_n = 1:n_pairs
    sample_a = values(:, pairs(pair_n, 1));
    sample_b = values(:, pairs(pair_n, 2));
    if all(sample_a - sample_b == 0)
        p_raw(pair_n) = 1;
        continue
    end
    try
        p_raw(pair_n) = signrank(sample_a, sample_b);
    catch
        p_raw(pair_n) = NaN;
    end
end
p_adjusted = min(p_raw * n_pairs, 1);
end

function p_value = friedman_p(values)
p_value = NaN;
if size(values, 2) < 3 || size(values, 1) < 2
    return
end
try
    p_value = friedman(values, 1, 'off');
catch
    p_value = NaN;
end
end

function report_region_tests(method_name, band_name, region_list, ...
        region_names, n_shanks, p_omnibus, pairs, p_raw, p_adjusted)
fprintf('%s | %s: %d shanks contain %s.\n', ...
    method_name, band_name, n_shanks, region_list);
if numel(region_names) >= 3
    if isfinite(p_omnibus)
        fprintf('  Friedman p = %.4g\n', p_omnibus);
    else
        fprintf('  Friedman p = n/a\n');
    end
end
for pair_n = 1:size(pairs, 1)
    fprintf('  %s vs %s: sign-rank p = %s, Bonferroni p = %s %s\n', ...
        region_display_name(region_names{pairs(pair_n, 1)}), ...
        region_display_name(region_names{pairs(pair_n, 2)}), ...
        format_p_value(p_raw(pair_n)), format_p_value(p_adjusted(pair_n)), ...
        significance_label(p_adjusted(pair_n)));
end
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

function grid = fine_distance_grid(distance_um, step_um)
distance_um = distance_um(isfinite(distance_um));
if isempty(distance_um)
    grid = 0;
    return
end
grid = (floor(min(distance_um) / step_um) * step_um):step_um: ...
    (ceil(max(distance_um) / step_um) * step_um);
end

function [mean_line, lower_line, upper_line, n_profiles] = ...
        interpolated_session_interval(power_table, row_mask, ...
        distance_um, power, grid)
group_index = session_shank_groups(power_table, row_mask);
n_groups = numel(group_index);
resampled = nan(n_groups, numel(grid));
n_profiles = 0;
for group_n = 1:n_groups
    rows = group_index{group_n};
    resampled(group_n, :) = resample_one( ...
        distance_um(rows), power(rows), grid);
    finite = isfinite(distance_um(rows)) & isfinite(power(rows));
    if nnz(finite) >= 2
        n_profiles = n_profiles + 1;
    end
end
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

function group_index = session_shank_groups(power_table, row_mask)
rows = find(row_mask);
animal = cellstr(string(power_table.Animal(rows)));
session_key = cellstr(string(power_table.Session(rows)));
shank_key = cellstr(string(power_table.Shank(rows)));
group_key = strcat(animal, '|', session_key, '|', shank_key);
[~, ~, group_id] = unique(group_key, 'stable');
n_groups = max([group_id; 0]);
group_index = cell(n_groups, 1);
for group_n = 1:n_groups
    group_index{group_n} = rows(group_id == group_n);
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
        session_label)
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
    track.welch_bands = [ ...
        mean_in_band(welch_rows, frequency, delta_band), ...
        mean_in_band(welch_rows, frequency, theta_band), ...
        mean_in_band(welch_rows, frequency, gamma_band)];
    track.fourier_bands = [ ...
        mean_in_band(fourier_rows, frequency, delta_band), ...
        mean_in_band(fourier_rows, frequency, theta_band), ...
        mean_in_band(fourier_rows, frequency, gamma_band)];
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

panel_kinds = {'tracks', 'interval'};
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
            figure_name = sprintf('%s 95%% CI vs distance to region center', ...
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
                if strcmp(panel_kind, 'interval') && numel(side_tracks) >= 2
                    for band_n = 1:3
                        [grid, mean_line, lower_line, upper_line] = ...
                            class_interval(side_tracks, method_n, band_n);
                        draw_band_interval(grid, mean_line, lower_line, ...
                            upper_line, band_colors(band_n, :));
                        finite_grid = grid(isfinite(mean_line));
                        if ~isempty(finite_grid)
                            x_max = max(x_max, max(finite_grid));
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
                draw_rank_correlations(side_tracks, method_n, band_colors);
                if strcmp(panel_kind, 'interval') && numel(side_tracks) < 2
                    text(0.5, 0.5, 'CI needs at least 2 shanks', ...
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

function note = region_distance_note(scale_mode, with_interval)
if strcmp(scale_mode, 'normalized')
    note = ['Each shank is z-scored and set to 0 at its first channel. ' ...
        'Spearman rho compares distance with that power.'];
else
    note = ['Saved band power. Spearman rho compares distance with power.'];
end
if with_interval
    note = [note, ' Thick line and shade are the mean and 95% CI where n >= 2.'];
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

function draw_rank_correlations(side_tracks, method_n, band_colors)
for band_n = 1:3
    rho = panel_rank_correlation(side_tracks, method_n, band_n);
    if isfinite(rho)
        label = sprintf('\\rho = %.2f', rho);
    else
        label = '\rho = n/a';
    end
    text(0.98, 0.98 - (band_n - 1) * 0.09, label, ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'Color', band_colors(band_n, :), ...
        'FontSize', 8, 'Interpreter', 'tex');
end
end

function rho = panel_rank_correlation(side_tracks, method_n, band_n)
distance = [];
power = [];
for track_n = 1:numel(side_tracks)
    distance = [distance; side_tracks(track_n).rel_um(:)]; %#ok<AGROW>
    power = [power; band_column(side_tracks(track_n), method_n, band_n)]; %#ok<AGROW>
end
rho = NaN;
keep = isfinite(distance) & isfinite(power);
if nnz(keep) < 3
    return
end
distance = distance(keep);
power = power(keep);
if numel(unique(distance)) < 2 || numel(unique(power)) < 2
    return
end
rho = corr(distance, power, 'Type', 'Spearman');
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
item = struct('rel_um', [], 'welch_bands', [], 'fourier_bands', []);
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

function draw_band_interval(grid, mean_line, lower_line, upper_line, color)
x = grid(:);
y_low = lower_line(:);
y_high = upper_line(:);
finite = isfinite(x) & isfinite(y_low) & isfinite(y_high);
if ~any(finite)
    return
end
breaks = diff([false; finite; false]);
starts = find(breaks == 1);
stops = find(breaks == -1) - 1;
for segment_n = 1:numel(starts)
    index = starts(segment_n):stops(segment_n);
    if numel(index) < 2
        continue
    end
    segment_x = x(index).';
    fill([segment_x, fliplr(segment_x)], ...
        [y_low(index).', fliplr(y_high(index).')], color, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
plot(grid, mean_line, 'Color', color, 'LineWidth', 2.0, ...
    'HandleVisibility', 'off');
end
