function plot_physical_scatter(x, y, at_cut, x_label, y_label, plot_title, below_label, above_label)
%PLOT_PHYSICAL_SCATTER  One session per point, x already in cm or cm/s.

hold on
x = x(:);
y = y(:);
at_cut = at_cut(:);
has_point = isfinite(x) & isfinite(y);
below = has_point & ~at_cut;
above = has_point & at_cut;
scatter(x(below), y(below), 28, [0.6 0.6 0.6], 'filled')
scatter(x(above), y(above), 36, 'k', 'filled')
yline(0, 'k:')
xlabel(x_label)
ylabel(y_label)
title(plot_title)
if any(below) && any(above)
    legend({below_label, above_label}, 'Location', 'best')
end
end
