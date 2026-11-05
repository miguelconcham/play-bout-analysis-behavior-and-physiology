function theta_envelope_plot_aligned(grouped, actor, measure)
%THETA_ENVELOPE_PLOT_ALIGNED  One panel per behavior for one actor and measure.
%
% Gray lines are animal means. The black line is the mean of animals with
% enough onsets of that behavior. The blue band is -2 to 0 s and the red
% line is behavior onset.

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

n_panels = numel(selected);
n_columns = 4;
n_rows = ceil(n_panels / n_columns);
figure('Name', sprintf('%s, %s', measure_label, actor_label), 'Color', 'w');
for panel_n = 1:n_panels
    item = selected(panel_n);
    subplot(n_rows, n_columns, panel_n)
    curves = item.mean_epoch;
    y_values = curves(isfinite(curves));
    if isempty(y_values)
        span = [-1 1];
    else
        span = [min(y_values), max(y_values)];
        if span(1) == span(2)
            span = span + [-1 1];
        end
    end
    hold on
    patch([-2 0 0 -2], [span(1) span(1) span(2) span(2)], [0.85 0.9 1], ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    animal_line = plot(nan, nan, 'Color', [0.75 0.75 0.75], 'DisplayName', 'Animal');
    for animal_n = 1:size(curves, 1)
        plot(item.epoch_time, curves(animal_n, :), 'Color', [0.75 0.75 0.75], ...
            'HandleVisibility', 'off');
    end
    mean_line = plot(item.epoch_time, item.grand_mean_epoch, 'k', ...
        'LineWidth', 1.5, 'DisplayName', 'Mean');
    if panel_n == 1
        onset_line = plot([0 0], span, 'r', 'DisplayName', 'Onset');
    else
        plot([0 0], span, 'r', 'HandleVisibility', 'off');
    end
    hold off
    if item.is_play
        kind = 'play';
    else
        kind = 'other';
    end
    title(sprintf('%s (%s), p = %.3f', item.behavior, kind, item.locked_p), ...
        'Interpreter', 'none')
    if panel_n == 1
        xlabel('Time from behavior onset (s)')
        ylabel(y_label)
        legend([animal_line, mean_line, onset_line], 'Location', 'northeast')
    end
    if ~isempty(item.epoch_time)
        xlim([item.epoch_time(1), item.epoch_time(end)])
    end
end
annotation_line = sprintf('%s, %s. p is the 2-3 Hz onset-shift test.', ...
    measure_label, actor_label);
try
    sgtitle(annotation_line, 'Interpreter', 'none');
catch
    set(gcf, 'Name', annotation_line);
end
end
