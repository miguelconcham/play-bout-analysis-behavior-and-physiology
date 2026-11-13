%% Analyze pooled cross-correlations between the six cell groups
% Reads the per-session files from Estimate_group_cross_correlation.
% Each session contributes one excess-rate curve per directed pair. The
% session is the unit: sessions are averaged with equal weight.
%
% Excess rate is the target group's conditional rate around a reference
% spike, minus the target's own Partners-session rate. A negative value
% means the target fired less than its own rate while the reference group
% was active. The printed test uses the mean excess inside +/- central_lag
% on the unsmoothed curve. The drawn curves are Gaussian-smoothed first,
% and the ribbon is the 95% t interval of those smoothed session curves.
% Holm-Sidak is applied across the 30 directed pairs.
%
% figure_folder empty keeps the figures on screen and still writes
% group_ccg_summary.txt in the result folder. A non-empty figure_folder
% also writes SVG and JPEG copies there.

%% Paths and parameters
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
run(fullfile(repo_root, 'add_repo_paths.m'));

data_root = group_ccg_data_root();
result_folder = fullfile(data_root, 'Analysis results', 'group cross correlation');
central_lag = 0.05;
smooth_half = 0.03;
alpha = 0.05;
figure_folder = result_folder;

sessions = load_group_ccg_sessions(result_folder);
summary = pool_group_ccg_sessions(sessions, central_lag);
summary_path = fullfile(result_folder, 'group_ccg_summary.txt');
write_group_ccg_summary(summary_path, summary, alpha);
fprintf('Summary: %s\n', summary_path);

draw_central_matrix(summary, alpha, figure_folder);
draw_excess_curves(summary, smooth_half, figure_folder);


function sessions = load_group_ccg_sessions(result_folder)
files = dir(fullfile(result_folder, '*_group_ccg.mat'));
if isempty(files)
    error('group_ccg:analyze', ...
        'No per-session correlograms in %s. Run Estimate_group_cross_correlation first.', ...
        result_folder);
end
sessions = [];
for i = 1:numel(files)
    stored = load(fullfile(files(i).folder, files(i).name), 'group_ccg');
    if ~isfield(stored, 'group_ccg')
        warning('group_ccg:analyze', 'No group_ccg in %s.', files(i).name);
        continue
    end
    session = stored.group_ccg;
    if ~isfield(session, 'version') || session.version ~= 1
        error('group_ccg:analyze', ...
            '%s was saved by a different estimate. Re-run Estimate_group_cross_correlation.', ...
            files(i).name);
    end
    if isempty(sessions)
        sessions = session;
    else
        sessions(end + 1) = session; %#ok<AGROW>
    end
end
if isempty(sessions)
    error('group_ccg:analyze', 'None of the files in %s contained group_ccg.', result_folder);
end
end


function write_group_ccg_summary(path, summary, alpha)
names = group_ccg_names();
lines = cell(numel(summary.pair_reference) + 8, 1);
lines{1} = 'Pooled cross-correlation between lock type and oscillation class.';
lines{2} = ['Excess rate is target spikes/s around a reference spike, minus the target rate. ', ...
    'Negative means the target was below its own rate.'];
lines{3} = sprintf('Central window: +/- %.0f ms, unsmoothed. A star in the matrix is Holm-Sidak p < %.2f.', ...
    1000 * summary.central_lag, alpha);
lines{4} = sprintf('%-28s %-28s %6s %10s %10s %10s %10s %10s', ...
    'Reference', 'Target', 'N', 'Excess', 'CI low', 'CI high', 'Ratio', 'Holm p');
n_lines = 4;
for k = 1:numel(summary.pair_reference)
    ratio_values = summary.central_ratio(:, k);
    ratio_values = ratio_values(isfinite(ratio_values));
    if isempty(ratio_values)
        ratio_mean = NaN;
    else
        ratio_mean = mean(ratio_values);
    end
    n_lines = n_lines + 1;
    lines{n_lines} = sprintf('%-28s %-28s %6d %10.3f %10.3f %10.3f %10.3f %10.4g', ...
        names.title{summary.pair_reference(k)}, names.title{summary.pair_target(k)}, ...
        summary.n_sessions(k), summary.central_mean(k), summary.central_low(k), ...
        summary.central_high(k), ratio_mean, summary.holm_p(k));
end
fid = fopen(path, 'w');
if fid < 0
    error('group_ccg:analyze', 'Could not write %s.', path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', lines{1:n_lines});
clear cleanup
end


function draw_central_matrix(summary, alpha, figure_folder)
names = group_ccg_names();
n_groups = numel(names.group_id);
matrix = nan(n_groups, n_groups);
stars = false(n_groups, n_groups);
for k = 1:numel(summary.pair_reference)
    matrix(summary.pair_reference(k), summary.pair_target(k)) = summary.central_mean(k);
    stars(summary.pair_reference(k), summary.pair_target(k)) = ...
        isfinite(summary.holm_p(k)) && summary.holm_p(k) < alpha;
end
fig = figure('Name', 'Group CCG central excess', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.15 0.15 0.7 0.7]);
imagesc(matrix, 'AlphaData', isfinite(matrix));
set(gca, 'Color', [0.93 0.93 0.93]);
axis square
colormap(gca, blue_red_colormap());
finite_values = matrix(isfinite(matrix));
if isempty(finite_values)
    limit = 1;
else
    limit = max(abs(finite_values));
    if ~(limit > 0)
        limit = 1;
    end
end
caxis([-limit, limit]);
colorbar
set(gca, 'XTick', 1:n_groups, 'YTick', 1:n_groups, ...
    'XTickLabel', names.title, 'YTickLabel', names.title, 'TickLength', [0 0]);
xtickangle(35);
xlabel('Target group');
ylabel('Reference group');
title(sprintf('Mean excess rate inside +/- %.0f ms', 1000 * summary.central_lag));
for i = 1:n_groups
    for j = 1:n_groups
        if ~isfinite(matrix(i, j))
            continue
        end
        label = sprintf('%.2f', matrix(i, j));
        if stars(i, j)
            label = [label, '*']; %#ok<AGROW>
        end
        text(j, i, label, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, 'Color', 'k');
    end
end
save_group_figure(fig, figure_folder, 'group_ccg_central_excess');
end


function draw_excess_curves(summary, smooth_half, figure_folder)
names = group_ccg_names();
base_title = names.title;
bin_w = median(diff(summary.centers));
sigma_bins = smooth_half / bin_w;
fig = figure('Name', 'Group CCG excess curves', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
for reference = 1:numel(names.group_id)
    subplot(2, 3, reference);
    hold on
    plotted = false(numel(names.group_id), 1);
    legend_text = base_title;
    for k = 1:numel(summary.pair_reference)
        if summary.pair_reference(k) ~= reference
            continue
        end
        target = summary.pair_target(k);
        sample = squeeze(summary.excess(:, :, k));
        if size(summary.excess, 1) == 1
            sample = sample(:)';
        end
        smoothed = nan(size(sample));
        for s = 1:size(sample, 1)
            if any(isfinite(sample(s, :)))
                smoothed(s, :) = gaussian_smooth_acg(sample(s, :), sigma_bins);
            end
        end
        [mu, lo, hi, n] = row_mean_ci(smoothed);
        if ~any(isfinite(mu))
            continue
        end
        color = names.color(target, :);
        draw_ci(summary.centers, lo, hi, color);
        plot(summary.centers, mu, 'Color', color, 'LineWidth', 1.5);
        plotted(target) = true;
        legend_text{target} = sprintf('%s (n=%d)', base_title{target}, max(n));
    end
    plot(summary.centers([1 end]), [0 0], 'Color', [0.75 0.75 0.75], ...
        'HandleVisibility', 'off');
    xlim([summary.centers(1), summary.centers(end)]);
    xlabel('Lag (s)');
    ylabel('Excess rate (spikes/s)');
    title(base_title{reference}, 'Color', names.color(reference, :));
    if any(plotted)
        legend(legend_text(plotted), 'Location', 'best', 'Box', 'off');
    end
end
if exist('sgtitle', 'file') == 2 || exist('sgtitle', 'builtin') == 5
    sgtitle(sprintf('Smoothed excess, %.0f ms Gaussian. Negative: target below its own rate.', ...
        1000 * smooth_half));
end
save_group_figure(fig, figure_folder, 'group_ccg_excess_curves');
end


function [mu, lo, hi, n] = row_mean_ci(sample)
n_bins = size(sample, 2);
mu = nan(1, n_bins);
lo = nan(1, n_bins);
hi = nan(1, n_bins);
n = zeros(1, n_bins);
for b = 1:n_bins
    values = sample(:, b);
    values = values(isfinite(values));
    n(b) = numel(values);
    if isempty(values)
        continue
    end
    mu(b) = mean(values);
    if numel(values) < 2
        continue
    end
    half = tinv(0.975, numel(values) - 1) * std(values) / sqrt(numel(values));
    lo(b) = mu(b) - half;
    hi(b) = mu(b) + half;
end
end


function draw_ci(centers, lo, hi, color)
covered = isfinite(lo) & isfinite(hi);
if sum(covered) < 2
    return
end
fill([centers(covered), fliplr(centers(covered))], ...
    [lo(covered), fliplr(hi(covered))], color, ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end


function colors = blue_red_colormap()
n = 256;
half = n / 2;
blue = [linspace(0.15, 1, half)', linspace(0.25, 1, half)', linspace(0.70, 1, half)'];
red = [linspace(1, 0.70, half)', linspace(1, 0.15, half)', linspace(1, 0.15, half)'];
colors = [blue; red];
end


function save_group_figure(fig, figure_folder, name)
if isempty(figure_folder)
    return
end
if exist(figure_folder, 'dir') ~= 7
    mkdir(figure_folder);
end
print(fig, fullfile(figure_folder, [name, '.svg']), '-dsvg');
print(fig, fullfile(figure_folder, [name, '.jpg']), '-djpeg', '-r200');
end
