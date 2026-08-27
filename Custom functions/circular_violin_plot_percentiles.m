function ranges = circular_violin_plot_percentiles(angles, groups, percentiles, sort_by_mean)

if nargin < 4
    sort_by_mean = false;
end

if nargin < 3 || isempty(percentiles)
    percentiles = [25 75]; % default IQR
end

angles = angles(:);
groups = groups(:);

group_ids = unique(groups);
n_groups = length(group_ids);

ranges = nan(n_groups, 3);

% --- Compute circular means ---
group_means = zeros(n_groups,1);
for i = 1:n_groups
    ang = angles(groups == group_ids(i));
    group_means(i) = angle(mean(exp(1i * ang), 'omitmissing'));
end
group_means = mod(group_means, 2*pi);

% --- Sort if requested ---
if sort_by_mean
    [~, order] = sort(group_means, 'ascend');
    group_ids = group_ids(order);
    group_means = group_means(order);
end

% --- Plot ---
figure; hold on;
colors = lines(n_groups);
max_width = 0.4;

for i = 1:n_groups
    g = group_ids(i);
    ang = angles(groups == g);

    % --- Wrap data for KDE ---
    ang_wrap = [ang; ang + 2*pi; ang - 2*pi];

    theta = linspace(0, 2*pi, 200);
    bw = 0.3;

    density = zeros(size(theta));
    for j = 1:length(ang_wrap)
        density = density + exp(cos(theta - ang_wrap(j)) / bw^2);
    end
    density = density / sum(density);
    density = density / max(density) * max_width;

    x0 = i;

    % --- Violin ---
    x_left  = x0 - density;
    x_right = x0 + density;

    fill([x_left, fliplr(x_right)], ...
         [theta, fliplr(theta)], ...
         colors(i,:), ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % --- Density-aware swarm scatter (FIXED interpolation) ---
    % Wrap angles into [0, 2pi]
    ang_wrapped = mod(ang, 2*pi);

    density_interp = interp1(theta, density, ang_wrapped, 'linear');

    % Handle any NaNs at boundaries
    density_interp(isnan(density_interp)) = 0;
    density_interp = density_interp / max(density_interp);

    spread = (rand(size(ang)) - 0.5) * 2;
    x_swarm = x0 + spread .* density_interp * max_width * 0.9;

    scatter([x_swarm;x_swarm], [ang;ang+2*pi], 6, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

    % --- Mean line ---
    mean_ang = group_means(i);
    plot([x0 - max_width, x0 + max_width], ...
         [mean_ang, mean_ang], ...
         'k', 'LineWidth', 2);

    % --- Percentile lines (NEW) ---
    % unwrap around circular mean for stable percentile calc
    % --- Percentile lines (and output) ---
    ang_centered = angle(exp(1i*(ang - mean_ang)));

    p_vals = prctile(ang_centered, percentiles);

    % map back to original space
    p_vals = mod(p_vals + mean_ang, 2*pi);

    % store in output
    ranges(i,1:2) = p_vals;
    ranges(i,3)   = g;

    for p = 1:length(p_vals)
        plot([x0 - max_width*0.8, x0 + max_width*0.8], ...
             [p_vals(p), p_vals(p)], ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 2);
        plot([x0 - max_width*0.8, x0 + max_width*0.8], ...
             [p_vals(p), p_vals(p)], ...
              'Color', [1 0 0], 'LineWidth', 1);
    end
end

% --- Formatting ---
ylim([0 2*pi])
yticks([0 pi/2 pi 3*pi/2 2*pi])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})

xlim([0 n_groups+1])
xticks(1:n_groups)
xticklabels(string(group_ids))

ylabel('Angle (rad)')
xlabel('Group')

title(sprintf('Circular Violin + Swarm + Percentiles [%d %d]', percentiles(1), percentiles(2)))

box on
end