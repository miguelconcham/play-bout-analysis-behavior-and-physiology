%% Analyze pooled cross-correlations between the six cell groups
% Reads the per-session files from Estimate_group_cross_correlation.
% Each session curve is normalized before the sessions are combined:
% conditional rate divided by that session's target rate, minus 1. Zero
% is independence. A negative value is the fraction by which the target
% group fell below its own rate. The mean and the 95% t interval are the
% spread of those normalized session values, with equal session weight.
%
% The printed test uses the mean normalized excess inside +/- central_lag
% on the unsmoothed curve. Each reference group has its own figure.
% Columns are the other groups. The top image is the sessions sorted by
% the central value, the middle image keeps recording order, and the
% bottom row is the mean of the smoothed correlograms and its 95% interval.
% smooth_ms is the Gaussian standard deviation in milliseconds. A value
% no larger than the bin width leaves every correlogram unsmoothed.
% x_lim is the lag window in seconds. Holm-Sidak is
% applied across the 30 directed pairs.
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
smooth_ms = 30;
x_lim = [-0.5 0.5];
alpha = 0.05;
figure_folder = result_folder;

sessions = load_group_ccg_sessions(result_folder);
summary = pool_group_ccg_sessions(sessions, central_lag);
summary_path = fullfile(result_folder, 'group_ccg_summary.txt');
write_group_ccg_summary(summary_path, summary, alpha);
fprintf('Summary: %s\n', summary_path);

draw_central_matrix(summary, alpha, figure_folder);
draw_excess_curves(summary, smooth_ms, x_lim, figure_folder);


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
    if ~isfield(session, 'version') || session.version ~= 2
        error('group_ccg:analyze', ...
            ['%s was saved at a bin width other than 1 ms. ' ...
            'Re-run Estimate_group_cross_correlation.'], ...
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
lines{2} = ['Each session is normalized by its own target rate before averaging. ', ...
    'Normalized excess is rate/expected - 1. Zero is independence. ', ...
    'Negative means the target was below its own rate. The confidence interval is across sessions.'];
lines{3} = sprintf('Central window: +/- %.0f ms, unsmoothed. A star in the matrix is Holm-Sidak p < %.2f.', ...
    1000 * summary.central_lag, alpha);
lines{4} = sprintf('%-28s %-28s %6s %10s %10s %10s %10s %10s', ...
    'Reference', 'Target', 'N', 'Norm', 'CI low', 'CI high', 'Ratio', 'Holm p');
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
title(sprintf('Mean normalized excess inside +/- %.0f ms', 1000 * summary.central_lag));
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


function draw_excess_curves(summary, smooth_ms, x_lim, figure_folder)
names = group_ccg_names();
for reference = 1:numel(names.group_id)
    columns = repmat(struct( ...
        'title', '', 'top', [], 'middle', [], ...
        'top_label', 'Sessions, sorted', 'middle_label', 'Sessions', ...
        'mean', [], 'lo', [], 'hi', [], 'note', ''), numel(names.group_id) - 1, 1);
    slot = 0;
    for k = 1:numel(summary.pair_reference)
        if summary.pair_reference(k) ~= reference
            continue
        end
        slot = slot + 1;
        sample = squeeze(summary.excess(:, :, k));
        if size(summary.excess, 1) == 1
            sample = sample(:)';
        end
        smoothed = smooth_ccg_matrix(sample, summary.centers, smooth_ms);
        [mu, lo, hi, n] = row_mean_ci(smoothed);
        shown = smoothed(any(isfinite(smoothed), 2), :);
        columns(slot).title = names.title{summary.pair_target(k)};
        columns(slot).top = sort_session_rows(shown, summary.centers, summary.central_lag);
        columns(slot).middle = shown;
        columns(slot).mean = mu;
        columns(slot).lo = lo;
        columns(slot).hi = hi;
        columns(slot).note = sprintf('n = %d', max([n, 0]));
    end
    draw_ccg_columns(figure_folder, ...
        ['group_ccg_', names.group_id{reference}], ...
        names.title{reference}, columns, summary.centers, x_lim);
end
end


function image = sort_session_rows(image, centers, central_lag)
if isempty(image)
    return
end
in_center = abs(centers) <= central_lag + 1e-12;
score = mean(image(:, in_center), 2, 'omitnan');
[~, order] = sort(score, 'descend');
image = image(order, :);
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
