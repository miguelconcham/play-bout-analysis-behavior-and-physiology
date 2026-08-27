function circular_violin_plot(angles, groups, sort_by_mean)

if nargin < 3
    sort_by_mean = false;
end

angles = angles(:);
groups = groups(:);

group_ids = unique(groups);
n_groups = length(group_ids);

% --- Compute circular means ---
group_means = zeros(n_groups,1);
for i = 1:n_groups
    ang = angles(groups == group_ids(i));
    group_means(i) = angle(mean(exp(1i * ang), 'omitmissing'));
end
group_means = mod(group_means, 2*pi);
group_means_mod = mod(group_means, 2*pi);

% --- Sort if requested ---
if sort_by_mean
    [~, order] = sort(group_means_mod, 'ascend');
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

    % --- Wrap data ---
    ang_wrap = [ang; ang + 2*pi; ang - 2*pi];

    % --- KDE ---
    theta = linspace(0, 2*pi, 200);
    bw = 0.3;

    density = zeros(size(theta));
    for j = 1:length(ang_wrap)
        density = density + exp(cos(theta - ang_wrap(j)) / bw^2);
    end
    density = density / sum(density);

    % Normalize width
    density = density / max(density) * max_width;

    x0 = i;

    % --- Violin ---
    x_left  = x0 - density;
    x_right = x0 + density;

    x = [x_left, fliplr(x_right)];
    y = [theta, fliplr(theta)];

    fill(x, y, colors(i,:), ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % --- Scatter (raw data) ---
    % --- Density-aware swarm scatter ---
    jitter_base = max_width * 0.9;

    % Interpolate density at each data point
    double_ang = [ang;ang+2*pi];
    density_interp = interp1(theta, density, double_ang, 'linear', 'extrap');

    % Normalize density (0 → 1)
    density_interp = density_interp / max(density_interp);

    % Generate symmetric spread (swarm-like)
    spread = (rand(size(double_ang)) - 0.5) * 2;

    x_swarm = x0 + spread .* density_interp * jitter_base;

    scatter(x_swarm, double_ang, ...
        2, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.5);
   

    % --- Mean line ---
    mean_ang = group_means(i);
    plot([x0 - max_width, x0 + max_width], ...
         [mean_ang, mean_ang], ...
         'k', 'LineWidth', 2);
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

if sort_by_mean
    title('Circular Violin Plot + Scatter (Sorted)')
else
    title('Circular Violin Plot + Scatter')
end

box on
end