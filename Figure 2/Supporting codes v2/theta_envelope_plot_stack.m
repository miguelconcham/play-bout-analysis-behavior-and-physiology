function theta_envelope_plot_stack(grouped, actor, measure)
%THETA_ENVELOPE_PLOT_STACK  Behaviors stacked above their mean and interval.
%
% The image has one row per behavior. The color is the mean waveform of
% animals with enough onsets of that behavior. Play rows sit above a white
% divider, then the other behaviors. Below, play and other each have a mean
% across behaviors and a 95% confidence interval.

if isempty(grouped)
    return
end
selected = grouped(strcmp({grouped.actor}, actor) & strcmp({grouped.measure}, measure));
if isempty(selected)
    return
end
if strcmp(actor, 'self')
    actor_label = 'Implanted animal';
else
    actor_label = 'Playing partner';
end
if strcmp(measure, 'envelope')
    measure_label = 'Hilbert envelope';
    y_label = 'Theta envelope (z)';
else
    measure_label = 'Spectrogram power';
    y_label = 'Theta power (z)';
end

time = selected(1).epoch_time(:).';
n_behaviors = numel(selected);
matrix = nan(n_behaviors, numel(time));
labels = cell(n_behaviors, 1);
is_play = false(n_behaviors, 1);
for behavior_n = 1:n_behaviors
    matrix(behavior_n, :) = selected(behavior_n).grand_mean_epoch(:).';
    labels{behavior_n} = selected(behavior_n).behavior;
    is_play(behavior_n) = selected(behavior_n).is_play;
end

figure('Name', sprintf('Stack, %s, %s', measure_label, actor_label), 'Color', 'w');
image_axis = subplot(3, 1, 1);
imagesc(time, 1:n_behaviors, matrix);
set(image_axis, 'YDir', 'reverse', 'TickDir', 'out');
yticks(1:n_behaviors);
yticklabels(labels);
ylim([0.5, n_behaviors + 0.5]);
hold on
n_play = sum(is_play);
if n_play > 0 && n_play < n_behaviors
    plot(time([1 end]), [n_play n_play] + 0.5, 'w', 'LineWidth', 1.5);
end
plot([0 0], [0.5, n_behaviors + 0.5], 'w', 'LineWidth', 1.5);
hold off
finite = matrix(isfinite(matrix));
if isempty(finite) || max(abs(finite)) == 0
    color_limit = 1;
else
    color_limit = max(abs(finite));
end
axes(image_axis);
colormap(blue_white_red(256));
caxis([-color_limit, color_limit]);
colorbar
xlabel('Time from behavior onset (s)');
title(sprintf('%s, %s', measure_label, actor_label), 'Interpreter', 'none');

mean_axis = subplot(3, 1, 2:3);
hold on
span = waveform_span(matrix);
patch([-2 0 0 -2], [span(1) span(1) span(2) span(2)], [0.85 0.9 1], ...
    'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
play_line = plot_group(time, matrix(is_play, :), [0.75 0.15 0.15], 'Play');
other_line = plot_group(time, matrix(~is_play, :), [0.15 0.15 0.15], 'Other');
plot([0 0], span, 'r', 'HandleVisibility', 'off');
plot(time([1 end]), [0 0], 'k:', 'HandleVisibility', 'off');
hold off
xlim([time(1), time(end)]);
xlabel('Time from behavior onset (s)');
ylabel(y_label);
title('Mean across behaviors, 95% confidence interval');
if ~isempty(play_line) && ~isempty(other_line)
    legend([play_line, other_line], 'Location', 'northeast');
elseif ~isempty(play_line)
    legend(play_line, 'Location', 'northeast');
elseif ~isempty(other_line)
    legend(other_line, 'Location', 'northeast');
end
set(mean_axis, 'TickDir', 'out');
drawnow;
image_position = get(image_axis, 'Position');
mean_position = get(mean_axis, 'Position');
mean_position(1) = image_position(1);
mean_position(3) = image_position(3);
set(mean_axis, 'Position', mean_position);
end

function line_handle = plot_group(time, matrix, color, label)
line_handle = [];
if isempty(matrix)
    return
end
[mu, lo, hi] = row_interval(matrix);
finite = isfinite(mu) & isfinite(lo) & isfinite(hi);
if ~any(finite)
    return
end
band_time = [time(finite), fliplr(time(finite))];
band_value = [lo(finite), fliplr(hi(finite))];
fill(band_time, band_value, color, 'FaceAlpha', 0.2, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
line_handle = plot(time, mu, 'Color', color, 'LineWidth', 1.5, ...
    'DisplayName', label);
end

function [mu, lo, hi] = row_interval(matrix)
n_time = size(matrix, 2);
mu = nan(1, n_time);
lo = nan(1, n_time);
hi = nan(1, n_time);
for column = 1:n_time
    values = matrix(:, column);
    values = values(isfinite(values));
    n = numel(values);
    if n == 0
        continue
    end
    mu(column) = mean(values);
    if n < 2
        lo(column) = mu(column);
        hi(column) = mu(column);
        continue
    end
    standard_error = std(values) / sqrt(n);
    half_width = t_critical(n - 1) * standard_error;
    lo(column) = mu(column) - half_width;
    hi(column) = mu(column) + half_width;
end
end

function critical = t_critical(degrees)
if exist('tinv', 'file') == 2
    critical = tinv(0.975, degrees);
    return
end
table_df = [1 2 3 4 5 6 7 8 9 10 12 15 20 30 40 60 120];
table_t = [12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 ...
    2.228 2.179 2.131 2.086 2.042 2.021 2.000 1.980];
if degrees >= table_df(end)
    critical = 1.96;
else
    critical = interp1(table_df, table_t, degrees, 'linear', 1.96);
end
end

function span = waveform_span(matrix)
values = matrix(isfinite(matrix));
if isempty(values)
    span = [-1 1];
    return
end
span = [min(values), max(values)];
pad = 0.05 * (span(2) - span(1));
if pad == 0
    pad = 1;
end
span = span + [-pad, pad];
end

function colors = blue_white_red(n_colors)
half = floor(n_colors / 2);
rising = linspace(0, 1, half).';
falling = linspace(1, 0, n_colors - half).';
red = [rising; ones(n_colors - half, 1)];
green = [rising; falling];
blue = [ones(half, 1); falling];
colors = [red, green, blue];
end
