%% Analyse_play_session_probe_spectra
% Single-session view of the channel-wise periodic/aperiodic spectra.
%
% The image is the full spectrum at every channel, sorted by depth within
% each shank. The profile figure is a line plot, not an image: mean power
% in delta (1-5 Hz), theta (6-12 Hz), and gamma (30-45 Hz) against depth.
% Both use the same saved spectrum. Area depth limits from the area table
% are drawn on both figures, together with the mid-LPAG channel.
% No cross-session summary is made here.

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
recording_index = 1;
session_index = 1;

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

%% Plot each requested session on its own
for file_n = recording_index(:).'
    input_file = fullfile(saving_folder, file_list(file_n).name);
    fprintf('Loading %s\n', input_file);
    loaded = load(input_file, 'probe_spectra');
    probe_spectra = loaded.probe_spectra;
    frequency = probe_spectra.frequency(:).';
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
        fprintf('Plotting %s\n', session_label);

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
        xlabel('Frequency (Hz)')
        if group_n == 1
            ylabel('Depth (µm)')
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
        set(gca, 'YDir', 'reverse')
        these_borders = borders_for_probe(area_borders, this_group.probe_number);
        if ~isempty(finite_depth)
            ylim(depth_limits_with_borders(finite_depth, these_borders))
        end
        draw_area_borders(these_borders, this_group.depth_um, true)
        xlabel(profile_xlabel(spectrum_field))
        if group_n == 1
            ylabel('Depth (µm)')
        end
        title(sprintf('%s | %s', method_names{method_n}, this_group.label))
        if method_n == 1 && group_n == 1
            legend(line_handles, band_names, 'Location', 'best')
        end
        box off
    end
end
sgtitle({session_label, 'Band means along depth'}, 'Interpreter', 'none')
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
