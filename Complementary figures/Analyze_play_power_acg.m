%% Analyze_play_power_acg
% Are autocorrelogram changes during play a play effect or a delta-power
% effect? Reads the per-session files from Estimate_play_power_acg.
%
% Four paired comparisons, each inside a cell type. Cell type is the
% delta-lock label from EntireSession in delta_all_neurons_v2.mat:
% trough-locked, peak-locked, non-entrained, and delta-locked (peak and
% trough together).
%
%   Play, different delta power
%     low 25% versus high 25% of play-epoch delta power
%   Matched delta power, play versus no play
%     every matched play epoch versus its outside partner
%   Low delta power, play versus no play
%   High delta power, play versus no play
%
% The outside partner is the nearest unused epoch outside play whose
% decibel power falls inside the channel caliper. Time-wise, a sign-flip
% cluster test is run on the positive lags of the smoothed autocorrelogram.
% The same neurons are then compared on the locking-period properties:
% burst (0-50 ms), first-peak prominence, second-cycle prominence, folded
% autocorrelogram delta power after 100 ms, and the two peak latencies.
% Sum-to-1 and rate-ratio curves are both tested. A rate-ratio curve is the
% overlap-weighted autocorrelogram divided by rate squared times bin width,
% so a Poisson train sits at 1.
%
% figure_folder empty keeps figures on screen only. Statistics are printed
% either way. Each show_* pair is [save, keep open]. [false false] skips
% that figure and still prints the test.

%% Paths and parameters
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
run('\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\add_repo_paths.m');

result_folder = fullfile(data_root, 'Analysis results', 'play power acg');
phase_folder = fullfile(data_root, 'Analysis results', 'phase locking data');

alpha = 0.01;
non_entrained_p = 0.1;
min_spikes = 20;
delta_lag = [0.2 0.5];
smooth_half = 0.06;
skip_abs_lag = 0.02;
acg_refractory_s = 0.005;
acg_cluster_alpha = 0.05;
acg_cluster_perm = 1000;
plot_rate_scaled = true;
show_mean_acg = [true true];
show_rate_ratio_acg = [true true];
show_paired_tests = [true true];
areas_to_keep = {'SC', 'LPAG', 'VLPAG', 'DLPAG', 'DR'};
figure_folder = '';

group_order = {'trough', 'peak', 'delta', 'nonentrained', 'all'};
group_title = {'Trough-locked', 'Peak-locked', 'Delta-locked', 'Non-entrained', 'All neurons'};
contrasts = play_power_contrasts();
condition_names = play_power_condition_names();
compare_options = struct( ...
    'refractory_s', acg_refractory_s, ...
    'smooth_half', smooth_half, ...
    'skip_abs_lag', skip_abs_lag, ...
    'delta_lag', delta_lag, ...
    'cluster_alpha', acg_cluster_alpha, ...
    'cluster_perm', acg_cluster_perm);

if ~isempty(figure_folder)
    if exist(figure_folder, 'dir') ~= 7
        mkdir(figure_folder);
    end
    diary(fullfile(figure_folder, ['play_power_acg_statistics_', datestr(now, 'yyyy-mm-dd'), '.txt']));
    diary_cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
end

%% Cell types
load(fullfile(phase_folder, 'delta_all_neurons_v2.mat'), 'all_neurons');
id_all = double(cluster_id_column(all_neurons));
n_cells = height(all_neurons);
lock_lookup = cell(n_cells, 1);
for i = 1:n_cells
    lock_lookup{i} = delta_lock_group( ...
        all_neurons.EntireSession.PPCPval(i), ...
        all_neurons.EntireSession.PreferedAngle(i), ...
        alpha, non_entrained_p);
end

%% Stack session autocorrelograms
session_files = dir(fullfile(result_folder, '*_play_power_acg.mat'));
if isempty(session_files)
    error('Analyze_play_power_acg:empty', ...
        'No play-power ACG files in %s. Run Estimate_play_power_acg first.', result_folder);
end
stack = [];
n_missing = 0;
n_unlabeled = 0;
for f = 1:numel(session_files)
    loaded = load(fullfile(result_folder, session_files(f).name), 'play_power_acg');
    if ~isfield(loaded, 'play_power_acg')
        warning('Analyze_play_power_acg:contents', ...
            'No play_power_acg in %s', session_files(f).name);
        continue
    end
    session = loaded.play_power_acg;
    if isempty(stack)
        stack = empty_neuron_stack(numel(session.centers), condition_names);
        stack.centers = session.centers;
        stack.bin_edges = session.bin_edges;
    elseif numel(session.centers) ~= numel(stack.centers) ...
            || max(abs(session.centers - stack.centers)) > 1e-8
        error('Analyze_play_power_acg:bins', ...
            'Lag bins in %s do not match the first session.', session.animal_code);
    end
    for nn = 1:numel(session.neurons)
        neuron = session.neurons(nn);
        hit = find(session_matches(all_neurons.session, session.animal_code) ...
            & id_all == neuron.cluster_id, 1);
        if isempty(hit)
            n_missing = n_missing + 1;
            continue
        end
        this_group = lock_lookup{hit};
        if isempty(this_group)
            n_unlabeled = n_unlabeled + 1;
            continue
        end
        stack.session{end + 1, 1} = session.animal_code; %#ok<SAGROW>
        stack.cluster_id(end + 1, 1) = neuron.cluster_id; %#ok<SAGROW>
        stack.area{end + 1, 1} = neuron.area; %#ok<SAGROW>
        stack.lock_group{end + 1, 1} = this_group; %#ok<SAGROW>
        for c = 1:numel(condition_names)
            name = condition_names{c};
            block = neuron.conditions.(name);
            stack.acg.(name)(end + 1, :) = block.acg; %#ok<SAGROW>
            stack.n_spikes.(name)(end + 1, 1) = block.n_spikes; %#ok<SAGROW>
            stack.rate.(name)(end + 1, 1) = block.rate; %#ok<SAGROW>
            stack.mean_power_db.(name)(end + 1, 1) = block.mean_power_db; %#ok<SAGROW>
            stack.median_pair_gap_db.(name)(end + 1, 1) = block.median_pair_gap_db; %#ok<SAGROW>
        end
    end
    clear loaded session
end
if isempty(stack) || isempty(stack.lock_group)
    error('Analyze_play_power_acg:none', ...
        'No labeled neurons were found in the play-power ACG files.');
end

fprintf('\nPlay versus delta power, autocorrelogram\n')
fprintf('Neurons: %d labeled, %d missing from the delta table, %d between the lock thresholds\n', ...
    numel(stack.lock_group), n_missing, n_unlabeled)
fprintf('Minimum spikes on each side: %d\n', min_spikes)
for g = 1:numel(group_order)
    fprintf('  %s: %d\n', group_title{g}, sum(in_lock_group(stack.lock_group, group_order{g})))
end

%% Paired tests
scales = {'sum'};
scale_labels = {'Sum to 1'};
scale_flags = {show_mean_acg};
if plot_rate_scaled
    scales{end + 1} = 'ratio'; %#ok<SAGROW>
    scale_labels{end + 1} = 'Rate ratio (Poisson = 1)'; %#ok<SAGROW>
    scale_flags{end + 1} = show_rate_ratio_acg; %#ok<SAGROW>
end
n_contrasts = numel(contrasts);
n_scales = numel(scales);
n_groups = numel(group_order);
results = cell(n_contrasts, n_scales);
for s = 1:n_scales
    for c = 1:n_contrasts
        contrast = contrasts(c);
        fprintf('\n%s — %s\n', contrast.question, scale_labels{s})
        groups = struct();
        for g = 1:n_groups
            rows = in_lock_group(stack.lock_group, group_order{g}) ...
                & stack.n_spikes.(contrast.side_a) >= min_spikes ...
                & stack.n_spikes.(contrast.side_b) >= min_spikes;
            if ~any(rows)
                comparison = empty_comparison();
            else
                comparison = compare_play_power_acgs( ...
                    stack.acg.(contrast.side_a)(rows, :), ...
                    stack.acg.(contrast.side_b)(rows, :), ...
                    stack.centers, ...
                    stack.rate.(contrast.side_a)(rows), ...
                    stack.rate.(contrast.side_b)(rows), ...
                    scales{s}, compare_options);
            end
            kept = find(rows);
            kept = kept(comparison.keep);
            groups(g).name = group_order{g};
            groups(g).title = group_title{g};
            groups(g).comparison = comparison;
            groups(g).power_a = stack.mean_power_db.(contrast.side_a)(kept);
            groups(g).power_b = stack.mean_power_db.(contrast.side_b)(kept);
            groups(g).gap_db = stack.median_pair_gap_db.(contrast.side_a)(kept);
            groups(g).n_spike_pairs = sum(rows);
            report_group(groups(g), contrast, scale_labels{s})
        end
        results{c, s} = groups;
    end
end

%% Figures
curve_colors = [0.20 0.45 0.75; 0.85 0.25 0.20];
group_colors = [0.20 0.45 0.80; 0.85 0.25 0.20; 0.10 0.50 0.35; 0.45 0.45 0.45; 0.15 0.15 0.15];
for s = 1:n_scales
    mean_flag = scale_flags{s};
    y_lim = shared_curve_limits(results(:, s), contrasts);
    for c = 1:n_contrasts
        if any(mean_flag)
            plot_mean_acg_figure(results{c, s}, contrasts(c), scale_labels{s}, ...
                stack.centers, delta_lag, y_lim, curve_colors, mean_flag, figure_folder)
        end
        if any(show_paired_tests)
            plot_property_figure(results{c, s}, contrasts(c), scale_labels{s}, ...
                group_colors, show_paired_tests, figure_folder)
        end
    end
end


function id_col = cluster_id_column(ci)
if ismember('cluster_id', ci.Properties.VariableNames)
    id_col = ci.cluster_id;
elseif ismember('id', ci.Properties.VariableNames)
    id_col = ci.id;
else
    error('Analyze_play_power_acg:id', 'Neuron table has neither cluster_id nor id.');
end
end


function tf = session_matches(session_col, animal_code)
if iscell(session_col)
    tf = strcmp(session_col, animal_code);
    empty = cellfun(@isempty, session_col);
    tf(empty) = false;
    nested = cellfun(@iscell, session_col);
    if any(nested)
        tf(nested) = cellfun(@(x) strcmp(x{1}, animal_code), session_col(nested));
    end
else
    tf = strcmp(string(session_col), animal_code);
end
end


function stack = empty_neuron_stack(n_bins, names)
stack.session = {};
stack.cluster_id = zeros(0, 1);
stack.area = {};
stack.lock_group = {};
stack.centers = [];
stack.bin_edges = [];
for i = 1:numel(names)
    stack.acg.(names{i}) = zeros(0, n_bins);
    stack.n_spikes.(names{i}) = zeros(0, 1);
    stack.rate.(names{i}) = zeros(0, 1);
    stack.mean_power_db.(names{i}) = zeros(0, 1);
    stack.median_pair_gap_db.(names{i}) = zeros(0, 1);
end
end


function tf = in_lock_group(lock_group, name)
if strcmp(name, 'all')
    tf = ismember(lock_group, {'trough', 'peak', 'nonentrained'});
elseif strcmp(name, 'delta')
    tf = ismember(lock_group, {'trough', 'peak'});
else
    tf = strcmp(lock_group, name);
end
tf = tf(:);
end


function comparison = empty_comparison()
comparison.keep = false(0, 1);
comparison.n = 0;
comparison.smooth_a = [];
comparison.smooth_b = [];
comparison.scores_a = struct();
comparison.scores_b = struct();
comparison.clusters = struct('clusters', repmat(struct( ...
    'x0', [], 'x1', [], 'mass', [], 'p', [], 'survives', []), 0, 1));
comparison.properties = struct();
end


function report_group(group, contrast, scale_label)
fprintf('  %s, %s\n', group.title, scale_label)
fprintf('    neurons with enough spikes: %d; finite on both sides: %d\n', ...
    group.n_spike_pairs, group.comparison.n)
if group.comparison.n < 1
    return
end
fprintf('    median power, %s %.3g dB, %s %.3g dB, B minus A %.3g dB\n', ...
    contrast.side_a_label, median(group.power_a, 'omitnan'), ...
    contrast.side_b_label, median(group.power_b, 'omitnan'), ...
    median(group.power_b - group.power_a, 'omitnan'))
gap = group.gap_db(isfinite(group.gap_db));
if ~isempty(gap)
    fprintf('    median absolute matched-power gap %.3g dB\n', median(gap))
end
report_clusters(group.comparison.clusters, contrast)
specs = acg_property_specs();
for m = 1:numel(specs)
    if ~isfield(group.comparison.properties, specs(m).name)
        continue
    end
    stats = group.comparison.properties.(specs(m).name);
    if strcmp(specs(m).kind, 'level')
        fprintf('    %s: median B-A %.3g, n %d, two-sided p %.3g, increase p %.3g\n', ...
            stats.label, stats.median_delta, stats.n, stats.p_two, stats.p_increase)
    else
        fprintf('    %s: median B-A %.3g s, n %d, two-sided p %.3g\n', ...
            stats.label, stats.median_delta, stats.n, stats.p_two)
    end
end
end


function report_clusters(result, contrast)
fprintf('    cluster test, %s minus %s, positive lags, uncorrected p < %.3g, %d sign flips\n', ...
    contrast.side_b_label, contrast.side_a_label, result.alpha, result.n_perm)
if ~isfield(result, 'clusters') || isempty(result.clusters)
    fprintf('    No contiguous uncorrected difference.\n')
    return
end
for k = 1:numel(result.clusters)
    cluster = result.clusters(k);
    if cluster.survives
        fate = 'survives';
    else
        fate = 'does not survive';
    end
    fprintf('    Cluster %.3g-%.3g s, cluster p %.3g, %s\n', ...
        cluster.x0, cluster.x1, cluster.p, fate)
end
end


function y_lim = shared_curve_limits(scale_results, ~)
values = [];
for c = 1:numel(scale_results)
    groups = scale_results{c};
    for g = 1:numel(groups)
        comparison = groups(g).comparison;
        if comparison.n < 1
            continue
        end
        values = [values, curve_band(comparison.smooth_a), curve_band(comparison.smooth_b)]; %#ok<AGROW>
    end
end
values = values(isfinite(values));
if isempty(values)
    y_lim = [0 1];
    return
end
pad = 0.08 * max(max(values) - min(values), eps);
y_lim = [min(values) - pad, max(values) + pad];
end


function band = curve_band(curves)
if isempty(curves)
    band = [];
    return
end
mu = mean(curves, 1, 'omitnan');
n = sum(isfinite(curves), 1);
sem = std(curves, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
band = [mu - sem, mu + sem];
end


function plot_mean_acg_figure(groups, contrast, scale_label, centers, delta_lag, y_lim, curve_colors, flag, folder)
fig = figure('Name', sprintf('%s, %s', contrast.question, scale_label), ...
    'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.08 0.08 0.8 0.75]);
for g = 1:numel(groups)
    subplot(2, 3, g)
    hold on
    shade_lags(delta_lag, y_lim, [0.75 0.82 0.90])
    shade_lags(-delta_lag([2 1]), y_lim, [0.75 0.82 0.90])
    comparison = groups(g).comparison;
    if comparison.n >= 1
        plot_mean_sem(centers, comparison.smooth_a, curve_colors(1, :))
        plot_mean_sem(centers, comparison.smooth_b, curve_colors(2, :))
        shade_clusters(comparison.clusters, y_lim)
    end
    xlim([centers(1), centers(end)])
    ylim(y_lim)
    xlabel('Lag (s)')
    ylabel(scale_label)
    title(sprintf('%s (n=%d)', groups(g).title, comparison.n))
    if g == 1
        legend({contrast.side_a_label, contrast.side_b_label}, 'Location', 'northeast')
    end
    box off
end
finish_figure(fig, flag, folder)
end


function plot_property_figure(groups, contrast, scale_label, group_colors, flag, folder)
specs = acg_property_specs();
fig = figure('Name', sprintf('%s properties, %s', contrast.question, scale_label), ...
    'Color', 'w', 'Units', 'normalized', 'OuterPosition', [0.08 0.08 0.8 0.75]);
for m = 1:numel(specs)
    subplot(2, 3, m)
    hold on
    any_point = false;
    for g = 1:numel(groups)
        comparison = groups(g).comparison;
        if comparison.n < 1 || ~isfield(comparison.properties, specs(m).name)
            continue
        end
        ok = acg_measure_mask(comparison.scores_a, comparison.scores_b, specs(m).name);
        if ~any(ok)
            continue
        end
        xa = comparison.scores_a.(specs(m).name);
        xb = comparison.scores_b.(specs(m).name);
        plot(xa(ok), xb(ok), '.', 'Color', group_colors(g, :), ...
            'MarkerSize', 12, 'HandleVisibility', 'off')
        any_point = true;
    end
    if any_point
        pts = [xlim, ylim];
        lo = min(pts);
        hi = max(pts);
        plot([lo hi], [lo hi], 'k-', 'HandleVisibility', 'off')
        xlim([lo hi])
        ylim([lo hi])
    end
    xlabel(contrast.side_a_label)
    ylabel(contrast.side_b_label)
    title(specs(m).label)
    axis square
    box off
end
subplot(2, 3, 1)
hold on
legend_handles = gobjects(numel(groups), 1);
legend_labels = cell(numel(groups), 1);
for g = 1:numel(groups)
    legend_handles(g) = plot(nan, nan, '.', 'Color', group_colors(g, :), 'MarkerSize', 12);
    legend_labels{g} = groups(g).title;
end
legend(legend_handles, legend_labels, 'Location', 'best')
finish_figure(fig, flag, folder)
end


function shade_lags(lag_win, y_lim, color)
patch([lag_win(1) lag_win(2) lag_win(2) lag_win(1)], ...
    [y_lim(1) y_lim(1) y_lim(2) y_lim(2)], color, ...
    'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off')
end


function shade_clusters(result, y_lim)
if ~isfield(result, 'clusters') || isempty(result.clusters)
    return
end
span = y_lim(2) - y_lim(1);
y0 = y_lim(2) - 0.1 * span;
y1 = y_lim(2);
for k = 1:numel(result.clusters)
    if ~result.clusters(k).survives
        continue
    end
    x0 = result.clusters(k).x0;
    x1 = result.clusters(k).x1;
    patch([x0 x1 x1 x0], [y0 y0 y1 y1], [0.25 0.25 0.25], ...
        'FaceAlpha', 0.45, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    patch([-x1 -x0 -x0 -x1], [y0 y0 y1 y1], [0.25 0.25 0.25], ...
        'FaceAlpha', 0.45, 'EdgeColor', 'none', 'HandleVisibility', 'off')
end
end


function h = plot_mean_sem(t, curves, color)
mu = mean(curves, 1, 'omitnan');
n = sum(isfinite(curves), 1);
sem = std(curves, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
lo = mu - sem;
hi = mu + sem;
ok = isfinite(mu) & isfinite(lo) & isfinite(hi);
breaks = diff([false, ok, false]);
i1 = find(breaks == 1);
i2 = find(breaks == -1) - 1;
for k = 1:numel(i1)
    idx = i1(k):i2(k);
    patch([t(idx), fliplr(t(idx))], [lo(idx), fliplr(hi(idx))], color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
end
h = plot(t, mu, 'Color', color, 'LineWidth', 1.5);
end


function finish_figure(fig, flag, folder)
do_save = flag(1) && ~isempty(folder);
keep_open = flag(2);
if do_save
    if exist(folder, 'dir') ~= 7
        mkdir(folder);
    end
    day_folder = fullfile(folder, ['PlayPower_', datestr(now, 'yyyy-mm-dd')]);
    if exist(day_folder, 'dir') ~= 7
        mkdir(day_folder);
    end
    name = regexprep(strtrim(fig.Name), '[^\w]+', '_');
    name = regexprep(name, '^_|_$', '');
    if isempty(name)
        name = sprintf('figure_%d', fig.Number);
    end
    set(fig, 'PaperPositionMode', 'auto')
    print(fig, fullfile(day_folder, [name, '.svg']), '-dsvg')
    print(fig, fullfile(day_folder, [name, '.jpg']), '-djpeg', '-r200')
    fprintf('Saved %s\n', fullfile(day_folder, [name, '.svg']))
end
if ~keep_open
    close(fig)
end
end
