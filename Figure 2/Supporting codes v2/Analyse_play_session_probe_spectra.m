%% Analyse_play_session_probe_spectra
% Single-session view of the channel-wise periodic/aperiodic spectra.
%
% The image is the full spectrum at every channel. The kilosort y coordinate
% is distance to the probe tip, not depth: larger values are further from the
% tip (upper), smaller values are closer to the tip (deeper).
%
% Selected sessions are stacked on one pair of figures, with mid-LPAG at zero.
% PAG sections, an upper side that passes through SupCol, and a deeper side
% that passes through DR are colored. Other territory stays gray. A 95%
% confidence interval is drawn for those three track classes when at least
% two tracks are available.
%
% Separate figures center each shank, when the region is present, on SC,
% DLPAG, LPAG, VLPAG, and DR. The center is the midpoint of that region's
% labeled span on the shank. Delta, theta, and gamma are drawn against
% distance from that center, once toward the tip (ventral) and once toward
% the brain surface (dorsal).

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

%% Load the requested recording files
file_list = dir(fullfile(saving_folder, ...
    'play_session_probe_periodic_aperiodic_spectra_*.mat'));
if isempty(file_list)
    error('No channel-wise spectra found in %s.', saving_folder);
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

        if plot_each_session
            plot_depth_spectrum(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, spectrum_field, ...
                session_label, area_borders);
            plot_band_profiles(frequency, channel_table, ...
                welch_spectrum, fourier_spectrum, ...
                delta_band, theta_band, gamma_band, ...
                spectrum_field, session_label, area_borders);
            drawnow
        end
    end
end

if isempty(placed)
    error('No probe channels were found in the selected recordings.');
end
if isempty(tracks)
    fprintf(['No shank with a mid-LPAG channel was found; the mid-LPAG ' ...
        'stack is skipped.\n']);
else
    fprintf(['Aligned tracks: %d. Stay in PAG: %d. Upper side through SupCol: ' ...
        '%d. Deeper side through DR: %d.\n'], numel(tracks), ...
        sum([tracks.stays_pag]), sum([tracks.goes_sc]), sum([tracks.goes_dr]));
    plot_aligned_spectra(tracks, stacked_frequency, spectrum_field);
    plot_aligned_profiles(tracks, spectrum_field, ...
        delta_band, theta_band, gamma_band);
end
plot_region_distances(placed, spectrum_field, ...
    delta_band, theta_band, gamma_band);

function plot_depth_spectrum(frequency, channel_table, ...
        welch_spectrum, fourier_spectrum, spectrum_field, ...
        session_label, area_borders)
probe_groups = probe_groups_by_depth(channel_table);
n_groups = numel(probe_groups);
spectra = {welch_spectrum, fourier_spectrum};
method_names = {'Welch', 'Fourier'};

figure('Name', ['Spectrum by depth | ', session_label], 'Color', 'w', ...
    'Position', [80 80 max(360, 280 * n_groups) 720]);
colors = diverging_or_sequential(spectrum_field);

for method_n = 1:2
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
        subplot(2, n_groups, (method_n - 1) * n_groups + group_n)
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
        title(sprintf('%s | %s', method_names{method_n}, this_group.label))
        box off
    end
    colorbar
end
sgtitle({session_label, spectrum_title(spectrum_field)}, 'Interpreter', 'none')
end

function plot_band_profiles(frequency, channel_table, ...
        welch_spectrum, fourier_spectrum, ...
        delta_band, theta_band, gamma_band, spectrum_field, ...
        session_label, area_borders)
probe_groups = probe_groups_by_depth(channel_table);
n_groups = numel(probe_groups);
spectra = {welch_spectrum, fourier_spectrum};
method_names = {'Welch', 'Fourier'};
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
band_colors = [0.12 0.35 0.75; 0.85 0.40 0.10; 0.45 0.62 0.20];

figure('Name', ['Band profiles | ', session_label], 'Color', 'w', ...
    'Position', [120 120 max(360, 280 * n_groups) 720]);

for method_n = 1:2
    for group_n = 1:n_groups
        subplot(2, n_groups, (method_n - 1) * n_groups + group_n)
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
        title(sprintf('%s | %s', method_names{method_n}, this_group.label))
        if method_n == 1 && group_n == 1
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
        delta_band, theta_band, gamma_band)
palette = track_colors();
method_names = {'Welch', 'Fourier'};
band_names = { ...
    sprintf('Delta %g-%g Hz', delta_band(1), delta_band(2)), ...
    sprintf('Theta %g-%g Hz', theta_band(1), theta_band(2)), ...
    sprintf('Gamma %g-%g Hz', gamma_band(1), gamma_band(2))};
figure('Name', 'Profiles aligned on mid-LPAG', 'Color', 'w', ...
    'Position', [60 60 1100 720]);
for method_n = 1:2
    for band_n = 1:3
        subplot(2, 3, (method_n - 1) * 3 + band_n)
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
        title(sprintf('%s | %s', method_names{method_n}, band_names{band_n}))
        if method_n == 1 && band_n == 3
            legend([guide; ci_handles], [guide_labels, ci_labels], ...
                'Location', 'best', 'FontSize', 7)
        end
        box off
    end
end
sgtitle('Sessions stacked on the mid-LPAG channel', 'Interpreter', 'none')
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

function plot_aligned_spectra(tracks, frequency, spectrum_field)
method_names = {'Welch', 'Fourier'};
class_fields = {'stays_pag', 'goes_sc', 'goes_dr'};
class_titles = {'Stays in PAG', 'Through SupCol', 'Through DR'};
figure('Name', 'Spectra aligned on mid-LPAG', 'Color', 'w', ...
    'Position', [80 80 1100 720]);
colors = diverging_or_sequential(spectrum_field);
for method_n = 1:2
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
        subplot(2, 3, (method_n - 1) * 3 + class_n)
        this_class = class_means{class_n};
        if this_class.n == 0 || isempty(this_class.spectrum)
            axis off
            text(0.5, 0.5, sprintf('%s\n%s\nn = 0', ...
                method_names{method_n}, class_titles{class_n}), ...
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
        title(sprintf('%s | %s (n=%d)', method_names{method_n}, ...
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
        delta_band, theta_band, gamma_band)
region_names = {'SC', 'DLPAG', 'LPAG', 'VLPAG', 'DR'};
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
method_names = {'Welch', 'Fourier'};
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
        region_names{region_n}, nnz(isfinite(centers(:, region_n))), n_placed);
end

for method_n = 1:2
    figure('Name', sprintf('%s power vs distance to region center', ...
        method_names{method_n}), 'Color', 'w', ...
        'Position', [40 40 980 1180]);
    legend_drawn = false;
    for region_n = 1:n_regions
        for direction_n = 1:2
            subplot(n_regions, 2, (region_n - 1) * 2 + direction_n)
            available = find(isfinite(centers(:, region_n)));
            if isempty(available)
                axis off
                text(0.5, 0.5, sprintf('%s\nnot available', ...
                    region_names{region_n}), ...
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
                draw_distance_track(distance, ...
                    power_for_method(welch_power, fourier_power, method_n), ...
                    band_colors);
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
                    region_names{region_n}, direction_titles{direction_n}), ...
                    'HorizontalAlignment', 'center', 'Interpreter', 'none');
                continue
            end
            x_max = farthest;
            if numel(side_tracks) >= 2
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
            if strcmp(spectrum_field, 'raw_aperiodic_subtracted')
                yline(0, ':k', 'HandleVisibility', 'off');
            end
            xline(0, '--k', 'HandleVisibility', 'off');
            if ~isfinite(x_max) || x_max <= 0
                x_max = 20;
            end
            xlim([0, x_max]);
            xlabel(direction_xlabels{direction_n})
            if direction_n == 1
                ylabel(profile_xlabel(spectrum_field))
            end
            title(sprintf('%s | %s (n=%d)', region_names{region_n}, ...
                direction_titles{direction_n}, numel(side_tracks)))
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
    sgtitle({sprintf('%s | %s', method_names{method_n}, ...
        spectrum_title(spectrum_field)), ...
        ['Distance 0 is the region center. Pale lines are shanks; ' ...
        'the thick line and shade are the mean and 95% CI where n >= 2.']}, ...
        'Interpreter', 'none')
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
item = struct('rel_um', [], 'welch_bands', [], 'fourier_bands', []);
end

function draw_distance_track(distance, power, band_colors)
[distance, order] = sort(distance(:));
power = power(order, :);
for band_n = 1:size(power, 2)
    pale = 0.35 * band_colors(band_n, :) + 0.65 * [1 1 1];
    this_power = power(:, band_n);
    finite = isfinite(distance) & isfinite(this_power);
    if nnz(finite) == 1
        plot(distance(finite), this_power(finite), 'o', 'Color', pale, ...
            'MarkerSize', 4, 'HandleVisibility', 'off');
        continue
    end
    plot(distance, this_power, 'Color', pale, 'LineWidth', 0.6, ...
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
