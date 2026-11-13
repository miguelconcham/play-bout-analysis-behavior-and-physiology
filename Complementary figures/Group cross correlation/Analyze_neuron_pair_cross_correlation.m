%% Analyze neuron-by-neuron cross-correlations
% Reads the pair matrices from Estimate_neuron_pair_cross_correlation.
% Each pair is normalized by its own target rate. Pairs in one recording
% are averaged, and the mean and 95% t interval are taken across
% recordings.
%
% One figure per reference group. Each column is one of the other groups.
% The top row is every neuron pair. The middle row is one row per
% session. The bottom row is the mean of the smoothed session
% correlograms and its confidence interval. smooth_ms is the Gaussian
% standard deviation in milliseconds. A value no larger than the bin
% width leaves every correlogram unsmoothed. x_lim is the lag window in
% seconds. Overlaid target groups are not drawn in the same axes.

%% Paths and parameters
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
run(fullfile(repo_root, 'add_repo_paths.m'));

data_root = group_ccg_data_root();
result_folder = fullfile(data_root, 'Analysis results', 'neuron pair cross correlation');
central_lag = 0.05;
smooth_ms = 30;
x_lim = [-0.5 0.5];
alpha = 0.05;
figure_folder = result_folder;

sessions = load_neuron_pair_sessions(result_folder);
summary = pool_neuron_pair_sessions(sessions, central_lag);
summary_path = fullfile(result_folder, 'neuron_pair_ccg_summary.txt');
write_neuron_pair_summary(summary_path, summary, alpha);
fprintf('Summary: %s\n', summary_path);
draw_neuron_pair_figures(summary, alpha, smooth_ms, x_lim, figure_folder);


function sessions = load_neuron_pair_sessions(result_folder)
files = dir(fullfile(result_folder, '*_neuron_pair_ccg.mat'));
if isempty(files)
    error('group_ccg:analyze', ...
        ['No neuron-pair correlograms in %s. ' ...
        'Run Estimate_neuron_pair_cross_correlation first.'], result_folder);
end
sessions = [];
for i = 1:numel(files)
    stored = load(fullfile(files(i).folder, files(i).name), 'neuron_pair_ccg');
    if ~isfield(stored, 'neuron_pair_ccg')
        warning('group_ccg:analyze', 'No neuron_pair_ccg in %s.', files(i).name);
        continue
    end
    session = stored.neuron_pair_ccg;
    if ~isfield(session, 'version') || session.version ~= 2
        error('group_ccg:analyze', ...
            '%s was not saved with the 1 ms neuron-pair estimate.', files(i).name);
    end
    if isempty(sessions)
        sessions = session;
    else
        sessions(end + 1) = session; %#ok<AGROW>
    end
end
if isempty(sessions)
    error('group_ccg:analyze', 'None of the files in %s contained neuron_pair_ccg.', result_folder);
end
end


function write_neuron_pair_summary(path, summary, alpha)
names = group_ccg_names();
lines = cell(numel(summary.pair_reference) + 8, 1);
lines{1} = 'Neuron-by-neuron cross-correlation between lock type and oscillation class.';
lines{2} = ['Each pair is normalized by its own target rate. The session value is the mean of its pairs. ', ...
    'Normalized excess is rate/expected - 1. The confidence interval is across sessions.'];
lines{3} = sprintf('Central window: +/- %.0f ms. A star would be Holm-Sidak p < %.2f.', ...
    1000 * summary.central_lag, alpha);
lines{4} = sprintf('%-28s %-28s %6s %8s %10s %10s %10s %10s', ...
    'Reference', 'Target', 'N ses', 'N pairs', 'Norm', 'CI low', 'CI high', 'Holm p');
n_lines = 4;
for k = 1:numel(summary.pair_reference)
    n_lines = n_lines + 1;
    lines{n_lines} = sprintf('%-28s %-28s %6d %8d %10.3f %10.3f %10.3f %10.4g', ...
        names.title{summary.pair_reference(k)}, names.title{summary.pair_target(k)}, ...
        summary.n_sessions(k), summary.n_pairs(k), summary.central_mean(k), ...
        summary.central_low(k), summary.central_high(k), summary.holm_p(k));
end
fid = fopen(path, 'w');
if fid < 0
    error('group_ccg:analyze', 'Could not write %s.', path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', lines{1:n_lines});
clear cleanup
end


function draw_neuron_pair_figures(summary, alpha, smooth_ms, x_lim, figure_folder)
names = group_ccg_names();
for reference = 1:numel(names.group_id)
    columns = empty_columns(numel(names.group_id) - 1);
    slot = 0;
    for k = 1:numel(summary.pair_reference)
        if summary.pair_reference(k) ~= reference
            continue
        end
        slot = slot + 1;
        pairs = smooth_ccg_matrix(finite_image_rows(summary.pair_excess{k}), ...
            summary.centers, smooth_ms);
        sessions = squeeze(summary.session_excess(:, :, k));
        if size(summary.session_excess, 1) == 1
            sessions = sessions(:)';
        end
        sessions = smooth_ccg_matrix(finite_image_rows(sessions), ...
            summary.centers, smooth_ms);
        [mu, lo, hi] = smoothed_session_ci(sessions, numel(summary.centers));
        columns(slot).title = names.title{summary.pair_target(k)};
        if isfinite(summary.holm_p(k)) && summary.holm_p(k) < alpha
            columns(slot).title = [columns(slot).title, ' *'];
        end
        columns(slot).top = sort_image_rows(pairs, summary.centers, summary.central_lag);
        columns(slot).middle = sort_image_rows(sessions, summary.centers, summary.central_lag);
        columns(slot).top_label = 'Neuron pairs';
        columns(slot).middle_label = 'Sessions';
        columns(slot).mean = mu;
        columns(slot).lo = lo;
        columns(slot).hi = hi;
        columns(slot).note = sprintf('sessions %d, pairs %d', ...
            summary.n_sessions(k), summary.n_pairs(k));
    end
    draw_ccg_columns(figure_folder, ...
        ['neuron_pair_ccg_', names.group_id{reference}], ...
        names.title{reference}, columns, summary.centers, x_lim);
end
end


function columns = empty_columns(n_targets)
columns = repmat(struct( ...
    'title', '', 'top', [], 'middle', [], ...
    'top_label', 'Neuron pairs', 'middle_label', 'Sessions', ...
    'mean', [], 'lo', [], 'hi', [], 'note', ''), n_targets, 1);
end


function image = finite_image_rows(image)
if isempty(image)
    image = zeros(0, 0);
    return
end
image = image(any(isfinite(image), 2), :);
end


function [mu, lo, hi] = smoothed_session_ci(sample, n_bins)
mu = nan(1, n_bins);
lo = mu;
hi = mu;
if isempty(sample)
    return
end
for b = 1:size(sample, 2)
    values = sample(:, b);
    values = values(isfinite(values));
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


function image = sort_image_rows(image, centers, central_lag)
if isempty(image)
    return
end
in_center = abs(centers) <= central_lag + 1e-12;
score = mean(image(:, in_center), 2, 'omitnan');
[~, order] = sort(score, 'descend');
image = image(order, :);
end
