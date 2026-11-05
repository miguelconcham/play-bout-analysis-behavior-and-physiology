%% Analyze_speed_freq_relation
% Band power vs speed during play vs non-play.
% Earlier sections z-score power and bin by RelativeSpeed.
% The last section keeps the saved log10 power and finds sessions that
% still have non-play samples at high locomotion speed.

this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_speed_freq_relation.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

%% Load precomputed speed-power model data
load([saving_folder,'\psth_structure_speed_delta.mat'],'psth_structure')
load([saving_folder,'\animal_names_speed_delta.mat'],'animal_names')

%% Session-level speed-power regression statistics
stasts_table = cell(numel(psth_structure), 6);

for j = 1:numel(psth_structure)
    stasts_table(j, :) = {psth_structure(j).lm.Coefficients.Estimate(4), ...
        psth_structure(j).lm.Coefficients.pValue(4), ...
        psth_structure(j).lm.Coefficients.tStat(4), ...
        animal_names{j, :}};
end

stasts_table = cell2table(stasts_table);
stasts_table.Properties.VariableNames = {'Estimate', 'pValue', 'tStat', 'Animal', 'Partner', 'Electrode'};
stasts_table = stasts_table(stasts_table.Electrode == 1, :);

lme = fitlme(stasts_table, 'Estimate ~ 1 + (1|Animal)');
coef = fixedEffects(lme);
ci = coefCI(lme);
pval = lme.Coefficients.pValue;

fprintf('Mean = %.3f, 95%% CI [%.3f, %.3f], p = %.4f\n', coef, ci(1), ci(2), pval);

%% Bin power by relative speed (play vs non-play)
n_grid = 101;
variable_grid = linspace(-2, 2, n_grid);
data_together = [];
mean_powers_play = nan(numel(psth_structure), n_grid);
mean_powers_noplay = nan(numel(psth_structure), n_grid);

for j = 1:numel(psth_structure)
    tbl = psth_structure(j).model_data;
    tbl = tbl(~isnan(tbl.Speed), :);

    Y = tbl.Power;
    Y = (Y - mean(Y, 'omitmissing')) / std(Y, 'omitmissing');
    X = tbl.RelativeSpeed;
    Xa = X(tbl.Play == 'true');
    Ya = Y(tbl.Play == 'true');
    Xb = X(tbl.Play == 'false');
    Yb = Y(tbl.Play == 'false');

    data_together = [data_together; [X Y]];

    edges = [variable_grid Inf];
    [~, ~, binA] = histcounts(Xa, edges);
    [~, ~, binB] = histcounts(Xb, edges);

    meanA = nan(size(variable_grid));
    meanB = nan(size(variable_grid));
    for i = 1:numel(variable_grid)
        meanA(i) = mean(Ya(binA == i), 'omitnan');
        meanB(i) = mean(Yb(binB == i), 'omitnan');
    end

    mean_powers_play(j, :) = meanA;
    mean_powers_noplay(j, :) = meanB;
end

%% Scatter and binned play-minus-non-play difference
power_diff = mean_powers_play - mean_powers_noplay;
[~, p_diff, ci_diff] = ttest(power_diff);
h_diff = p_diff < 0.05;
mean_diff = mean(power_diff, 1, 'omitnan');

[~, ~, samples2plot] = histcounts(data_together(:, 1), -5:.1:10);

figure
bins_with_values = unique(samples2plot)';
max_count = 500;
min_count = 50;
indexes2plot = [];

for j = bins_with_values
    indexes = find(samples2plot == j);
    if numel(indexes) >= max_count
        indexes = datasample(indexes, max_count);
    end
    if numel(indexes) >= min_count
        indexes2plot = [indexes2plot; indexes];
    end
end

selection_index = ~any(isnan(data_together(indexes2plot, :)), 2) & data_together(indexes2plot, 1) < 1 & abs(data_together(indexes2plot, 2)) < 3;

subplot(1, 2, 1)
plot(data_together(indexes2plot(selection_index), 1), data_together(indexes2plot(selection_index), 2), '.k', 'MarkerSize', .1)
[c, p] = corr(data_together(indexes2plot(selection_index), 1), data_together(indexes2plot(selection_index), 2));
title(num2str([c, p]))
axis xy

subplot(1, 2, 2)
hold on
no_nan = ~any(isnan(ci_diff), 1);
fill([variable_grid(no_nan) fliplr(variable_grid(no_nan))], ...
    [ci_diff(1, no_nan) fliplr(ci_diff(2, no_nan))], ...
    'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
plot(variable_grid, power_diff, ':k');
plot(variable_grid, mean_diff, 'k', 'LineWidth', 1.5);
yline(0, 'k:');
axis tight
y_lim = ylim;

% Significance bars in the top 10% of the axis (above the 90th percentile of ylim)
fill_lo = y_lim(1) + 0.9 * (y_lim(2) - y_lim(1));
fill_hi = y_lim(2);
start_end = [find(diff([0 h_diff 0]) == 1)' (find(diff([0 h_diff 0]) == -1) - 1)'];
for j = 1:size(start_end, 1)
    fill(variable_grid(start_end(j, [1 2 2 1])), [fill_lo fill_lo fill_hi fill_hi], ...
        'k', 'EdgeColor', 'none', 'FaceAlpha', .35)
end
xlabel('Relative speed (z)')
ylabel('Play − non-play power (z)')
title('Binned difference ± 95% CI (ttest); black bar: p < 0.05')

%% Mixed linear model of the play−non-play difference (repeated measures by animal)
% For each relative-speed bin: Difference ~ 1 + (1|Animal)
animal_labels = categorical(animal_names(:, 1));
n_bins = numel(variable_grid);

beta_est = nan(n_bins, 1);
SE = nan(n_bins, 1);
CI_low = nan(n_bins, 1);
CI_high = nan(n_bins, 1);
pValue = nan(n_bins, 1);

for bin_n = 1:n_bins
    y = power_diff(:, bin_n);
    keep = ~isnan(y);
    if nnz(keep) < 3 || numel(unique(animal_labels(keep))) < 2
        continue
    end

    tbl = table(y(keep), animal_labels(keep), 'VariableNames', {'Difference', 'Animal'});
    lme_diff = fitlme(tbl, 'Difference ~ 1 + (1|Animal)');

    coef = lme_diff.Coefficients;
    ci_lme = coefCI(lme_diff, 'Alpha', 0.05);

    beta_est(bin_n) = coef.Estimate(1);
    SE(bin_n) = coef.SE(1);
    CI_low(bin_n) = ci_lme(1, 1);
    CI_high(bin_n) = ci_lme(1, 2);
    pValue(bin_n) = coef.pValue(1);
end

sig_bins = find(pValue < 0.05);
fprintf('Speed bins with LME intercept p < 0.05: %d / %d\n', numel(sig_bins), n_bins);

figure
hold on
no_nan_lme = ~isnan(CI_low) & ~isnan(CI_high);
fill([variable_grid(no_nan_lme) fliplr(variable_grid(no_nan_lme))], ...
    [CI_low(no_nan_lme)' fliplr(CI_high(no_nan_lme)')], ...
    'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
plot(variable_grid, beta_est, 'k', 'LineWidth', 1.5)
yline(0, 'k:');
axis tight
y_lim = ylim;
fill_lo = y_lim(1) + 0.9 * (y_lim(2) - y_lim(1));
fill_hi = y_lim(2);
h_lme = pValue(:)' < 0.05;
h_lme(isnan(pValue(:)')) = false;
start_end = [find(diff([0 h_lme 0]) == 1)' (find(diff([0 h_lme 0]) == -1) - 1)'];
for j = 1:size(start_end, 1)
    fill(variable_grid(start_end(j, [1 2 2 1])), [fill_lo fill_lo fill_hi fill_hi], ...
        'k', 'EdgeColor', 'none', 'FaceAlpha', .35)
end
xlabel('Relative speed (z)')
ylabel('LME intercept (play − non-play)')
title('Mixed model Difference ~ 1 + (1|Animal); black bar: p < 0.05')

%% Full-dataset correlation
selection_index = ~any(isnan(data_together(:, :)), 2) & data_together(:, 1) < 1 & abs(data_together(:, 2)) < 3;
[c, p] = corr(data_together(selection_index, 1), data_together(selection_index, 2));

%% Sessions with enough high-speed non-play (raw log power, electrode 1)
% model_data.Power is log10 spectrogram amplitude. It is not z-scored in the
% mat. The sections above z-score it only for those plots. This section does not.
% model_data.Speed is already a within-session z-score, so "highest speed" is
% the top of that session's own speed distribution. Tracking units are not
% stored. Mats written before the RelativeSpeed z-score fix stored a rescaled
% Distance in that column; the check below detects that.

high_z = 1;
ceiling_bin_width = 0.5;
min_nonplay_s = 2;
min_play_s = 1;
min_ceiling_nonplay_s = 1;
min_ceiling_play_s = 0.5;
curve_bin_width = 0.25;
min_bin_nonplay_s = 0.5;
min_bin_play_s = 0.5;

if size(animal_names, 1) ~= numel(psth_structure)
    error('Analyze_speed_freq_relation:NameMismatch', ...
        'animal_names rows (%d) do not match psth_structure (%d).', ...
        size(animal_names, 1), numel(psth_structure));
end
electrode_id = cellfun(@double, animal_names(:, 3));
session_rows = find(electrode_id == 1);
if isempty(session_rows)
    error('Analyze_speed_freq_relation:NoSessions', ...
        'animal_names has no electrode-1 rows.');
end

if isfield(psth_structure(session_rows(1)), 'freq_range')
    band = psth_structure(session_rows(1)).freq_range;
    fprintf('Band in this file: %.1f–%.1f Hz\n', band(1), band(2));
end
fprintf(['Power is the saved log10 spectrogram amplitude (not z-scored here).\n' ...
    'Speed is the saved within-session z-score.\n']);

[r_distance, r_speed] = local_relative_speed_correlations(psth_structure, session_rows);
fprintf('Median |corr(RelativeSpeed, Distance)| = %.3f\n', median(abs(r_distance), 'omitnan'));
fprintf('Median |corr(RelativeSpeed, Speed)| = %.3f\n', median(abs(r_speed), 'omitnan'));
relative_is_distance = median(abs(r_distance), 'omitnan') > 0.95;
if relative_is_distance
    fprintf(['RelativeSpeed in this file tracks Distance.\n' ...
        'The generator used to write rescaled distance into that column.\n' ...
        'High-speed coverage below uses Speed. Re-run Estimate_speed_freq_relation\n' ...
        'before treating RelativeSpeed as the change in inter-animal distance.\n']);
end

local_report_speed_coverage(psth_structure, animal_names, session_rows, ...
    'Speed', high_z, ceiling_bin_width, min_nonplay_s, min_play_s, ...
    min_ceiling_nonplay_s, min_ceiling_play_s, ...
    curve_bin_width, min_bin_nonplay_s, min_bin_play_s);

if ~relative_is_distance
    fprintf(['RelativeSpeed is the within-session z-score of the change in\n' ...
        'inter-animal distance (positive: moving apart), not locomotion speed.\n']);
    local_report_speed_coverage(psth_structure, animal_names, session_rows, ...
        'RelativeSpeed', high_z, ceiling_bin_width, min_nonplay_s, min_play_s, ...
        min_ceiling_nonplay_s, min_ceiling_play_s, ...
        curve_bin_width, min_bin_nonplay_s, min_bin_play_s);
end

function [r_distance, r_speed] = local_relative_speed_correlations(psth_structure, session_rows)
r_distance = nan(numel(session_rows), 1);
r_speed = nan(numel(session_rows), 1);
for s = 1:numel(session_rows)
    mdl = psth_structure(session_rows(s)).model_data;
    ok = ~isnan(mdl.RelativeSpeed) & ~isnan(mdl.Distance) & ~isnan(mdl.Speed);
    if nnz(ok) < 10
        continue
    end
    r_distance(s) = corr(mdl.RelativeSpeed(ok), mdl.Distance(ok));
    r_speed(s) = corr(mdl.RelativeSpeed(ok), mdl.Speed(ok));
end
end

function local_report_speed_coverage(psth_structure, animal_names, session_rows, ...
    x_name, high_z, ceiling_bin_width, min_nonplay_s, min_play_s, ...
    min_ceiling_nonplay_s, min_ceiling_play_s, ...
    curve_bin_width, min_bin_nonplay_s, min_bin_play_s)

n_sess = numel(session_rows);
animals = cell(n_sess, 1);
partners = nan(n_sess, 1);
nonplay_high_s = nan(n_sess, 1);
play_high_s = nan(n_sess, 1);
enough_high = false(n_sess, 1);
diff_high = nan(n_sess, 1);
diff_high_z = nan(n_sess, 1);
top_lo = nan(n_sess, 1);
top_hi = nan(n_sess, 1);
top_nonplay_s = nan(n_sess, 1);
top_play_s = nan(n_sess, 1);
diff_top = nan(n_sess, 1);
diff_top_z = nan(n_sess, 1);
max_nonplay = nan(n_sess, 1);
max_play = nan(n_sess, 1);
X_cell = cell(n_sess, 1);
Y_cell = cell(n_sess, 1);
play_cell = cell(n_sess, 1);
hop_s = nan(n_sess, 1);
wind_s = nan(n_sess, 1);

for s = 1:n_sess
    entry = psth_structure(session_rows(s));
    [X, Y, is_play, hop_s(s), wind_s(s)] = local_session_vectors(entry, x_name);
    animals{s} = char(string(animal_names{session_rows(s), 1}));
    partners(s) = animal_names{session_rows(s), 2};
    X_cell{s} = X;
    Y_cell{s} = Y;
    play_cell{s} = is_play;
    if isempty(X)
        nonplay_high_s(s) = 0;
        play_high_s(s) = 0;
        continue
    end

    session_sd = std(Y, 'omitnan');
    high = X >= high_z;
    nonplay_high_s(s) = nnz(high & ~is_play) * hop_s(s);
    play_high_s(s) = nnz(high & is_play) * hop_s(s);
    enough_high(s) = nonplay_high_s(s) >= min_nonplay_s && play_high_s(s) >= min_play_s;
    if any(high & is_play) && any(high & ~is_play)
        diff_high(s) = mean(Y(high & is_play), 'omitnan') - mean(Y(high & ~is_play), 'omitnan');
        if session_sd > 0
            diff_high_z(s) = diff_high(s) / session_sd;
        end
    end

    if any(~is_play)
        max_nonplay(s) = max(X(~is_play));
    end
    if any(is_play)
        max_play(s) = max(X(is_play));
    end

    [top_lo(s), top_hi(s), top_play_s(s), top_nonplay_s(s), diff_top(s)] = ...
        local_highest_supported_bin(X, Y, is_play, hop_s(s), ceiling_bin_width, ...
        min_ceiling_play_s, min_ceiling_nonplay_s);
    if session_sd > 0 && ~isnan(diff_top(s))
        diff_top_z(s) = diff_top(s) / session_sd;
    end
end

[~, order] = sort(nonplay_high_s, 'descend');
wind_ref = wind_s(find(~isnan(wind_s), 1, 'first'));
hop_ref = hop_s(find(~isnan(hop_s), 1, 'first'));
fprintf('\n%s, electrode 1. Frame hop %.4f s, window %.3f s.\n', x_name, hop_ref, wind_ref);
fprintf(['Durations are spans of overlapping spectrogram frames. %.1f s is about %.1f windows.\n'], ...
    min_nonplay_s, min_nonplay_s / wind_ref);
fprintf('A session has enough high %s when non-play >= %.1f s and play >= %.1f s at z >= %.1f.\n', ...
    x_name, min_nonplay_s, min_play_s, high_z);

for k = 1:n_sess
    s = order(k);
    fprintf('%s P%d\n', animals{s}, partners(s));
    fprintf('  z >= %.1f: non-play %.2f s, play %.2f s, max non-play z %.2f, max play z %.2f\n', ...
        high_z, nonplay_high_s(s), play_high_s(s), max_nonplay(s), max_play(s));
    fprintf('  difference at z >= %.1f: %s\n', high_z, local_diff_text(diff_high(s), diff_high_z(s), enough_high(s)));
    if isnan(top_lo(s))
        fprintf('  no %.2f-wide bin has >= %.1f s play and >= %.1f s non-play\n', ...
            ceiling_bin_width, min_ceiling_play_s, min_ceiling_nonplay_s);
    else
        fprintf('  highest such bin [%.2f, %.2f]: non-play %.2f s, play %.2f s, %s\n', ...
            top_lo(s), top_hi(s), top_nonplay_s(s), top_play_s(s), ...
            local_diff_text(diff_top(s), diff_top_z(s), top_lo(s) >= high_z));
    end
end

n_enough = nnz(enough_high);
fprintf('\nSessions with enough non-play at %s z >= %.1f: %d / %d\n', ...
    x_name, high_z, n_enough, n_sess);
local_print_difference_test(diff_high(enough_high), animals(enough_high), ...
    sprintf('%s z >= %.1f, sessions with enough non-play', x_name, high_z));

top_is_high = ~isnan(top_lo) & top_lo >= high_z;
fprintf('Sessions whose highest supported bin starts at z >= %.1f: %d / %d\n', ...
    high_z, nnz(top_is_high), n_sess);
local_print_difference_test(diff_top(top_is_high), animals(top_is_high), ...
    sprintf('highest supported %s bin at z >= %.1f', x_name, high_z));
fprintf('Sessions with any supported bin, including bins below z = %.1f: %d / %d\n', ...
    high_z, nnz(~isnan(diff_top)), n_sess);
local_print_difference_test(diff_top(~isnan(diff_top)), animals(~isnan(diff_top)), ...
    sprintf('highest supported %s bin of each session', x_name));

x_max = high_z + 1;
for s = 1:n_sess
    if ~isempty(X_cell{s})
        x_max = max(x_max, max(X_cell{s}));
    end
end
x_max = min(6, ceil(x_max / curve_bin_width) * curve_bin_width);
n_curve_edges = round((x_max - (-1)) / curve_bin_width) + 1;
curve_edges = linspace(-1, x_max, n_curve_edges);
[curve_x, diff_all, n_all] = local_diff_curves(X_cell, Y_cell, play_cell, hop_s, ...
    curve_edges, min_bin_play_s, min_bin_nonplay_s, true(n_sess, 1));
[~, diff_cov, n_cov] = local_diff_curves(X_cell, Y_cell, play_cell, hop_s, ...
    curve_edges, min_bin_play_s, min_bin_nonplay_s, enough_high);

y_label = 'Play − non-play (log10 amplitude)';
figure('Name', ['High-speed coverage: ' x_name])
subplot(2, 2, 1)
local_plot_diff_curve(curve_x, diff_all, x_name, y_label, ...
    sprintf('All electrode-1 sessions, %s', x_name))
subplot(2, 2, 2)
local_plot_diff_curve(curve_x, diff_cov, x_name, y_label, ...
    sprintf('Sessions with >= %.1f s non-play at z >= %.1f', min_nonplay_s, high_z))
subplot(2, 2, 3)
local_plot_coverage_count(curve_x, n_all, x_name, min_bin_play_s, min_bin_nonplay_s)
subplot(2, 2, 4)
local_plot_coverage_count(curve_x, n_cov, x_name, min_bin_play_s, min_bin_nonplay_s)

figure('Name', ['High-speed sessions: ' x_name])
labels = cell(n_sess, 1);
for s = 1:n_sess
    labels{s} = sprintf('%s P%d', animals{order(s)}, partners(order(s)));
end
subplot(1, 2, 1)
bar_colors = repmat([0.75 0.75 0.75], n_sess, 1);
bar_colors(enough_high(order), :) = repmat([0 0 0], nnz(enough_high), 1);
duration_bars = barh(nonplay_high_s(order));
duration_bars.FaceColor = 'flat';
duration_bars.CData = bar_colors;
xline(min_nonplay_s, 'k:')
set(gca, 'YTick', 1:n_sess, 'YTickLabel', labels, 'YDir', 'reverse')
xlabel(sprintf('Non-play frame span at %s z \\geq %.1f (s)', x_name, high_z))
title('Black: both play and non-play pass the minimum')

subplot(1, 2, 2)
hold on
has_top = ~isnan(top_lo) & ~isnan(diff_top);
below = has_top & ~top_is_high;
above = has_top & top_is_high;
scatter(top_lo(below), diff_top(below), 28, [0.6 0.6 0.6], 'filled')
scatter(top_lo(above), diff_top(above), 36, 'k', 'filled')
yline(0, 'k:')
xline(high_z, 'k:')
xlabel(sprintf('Start of highest supported %s bin (z)', x_name))
ylabel(y_label)
title('Each session at its highest supported speed')
if any(below) && any(above)
    legend({'Bin starts below the high-speed cut', 'Bin starts at the high-speed cut'}, ...
        'Location', 'best')
end
end

function [X, Y, is_play, hop_s, wind_s] = local_session_vectors(entry, x_name)
mdl = entry.model_data;
X = mdl.(x_name)(:);
Y = mdl.Power(:);
is_play = local_play_mask(mdl.Play);
if numel(X) ~= numel(Y) || numel(is_play) ~= numel(Y)
    error('Analyze_speed_freq_relation:LengthMismatch', ...
        '%s, Power, and Play have different lengths.', x_name);
end
keep = ~isnan(X) & ~isnan(Y);
X = X(keep);
Y = Y(keep);
is_play = is_play(keep);
hop_s = local_spectrogram_hop(entry);
wind_s = 0.25;
if isfield(entry, 'wind_length') && ~isempty(entry.wind_length)
    wind_s = entry.wind_length;
end
end

function is_play = local_play_mask(play_col)
play_col = play_col(:);
if islogical(play_col)
    is_play = play_col;
elseif isnumeric(play_col)
    is_play = play_col ~= 0;
else
    names = string(play_col);
    is_play = names == "true" | names == "1";
end
end

function hop_s = local_spectrogram_hop(entry)
fs = 2500;
if isfield(entry, 'wind_length') && isfield(entry, 'wind_overlap') ...
        && ~isempty(entry.wind_length) && ~isempty(entry.wind_overlap)
    hop_samples = entry.wind_length * fs - floor(entry.wind_overlap * fs);
    hop_s = hop_samples / fs;
else
    hop_s = 0.0012;
end
if ~isfinite(hop_s) || hop_s <= 0
    hop_s = 0.0012;
end
end

function [lo, hi, play_s, nonplay_s, diff_val] = local_highest_supported_bin( ...
    X, Y, is_play, hop_s, bin_width, min_play_s, min_nonplay_s)

lo = NaN;
hi = NaN;
play_s = NaN;
nonplay_s = NaN;
diff_val = NaN;
if isempty(X)
    return
end
edge_lo = floor(min(X) / bin_width) * bin_width;
edge_hi = ceil(max(X) / bin_width) * bin_width;
if edge_hi <= edge_lo
    edge_hi = edge_lo + bin_width;
end
n_edges = round((edge_hi - edge_lo) / bin_width) + 1;
edges = linspace(edge_lo, edge_hi, n_edges);
for b = numel(edges) - 1:-1:1
    if b == numel(edges) - 1
        in_bin = X >= edges(b) & X <= edges(b + 1);
    else
        in_bin = X >= edges(b) & X < edges(b + 1);
    end
    this_play_s = nnz(in_bin & is_play) * hop_s;
    this_nonplay_s = nnz(in_bin & ~is_play) * hop_s;
    if this_play_s >= min_play_s && this_nonplay_s >= min_nonplay_s
        lo = edges(b);
        hi = edges(b + 1);
        play_s = this_play_s;
        nonplay_s = this_nonplay_s;
        diff_val = mean(Y(in_bin & is_play), 'omitnan') - mean(Y(in_bin & ~is_play), 'omitnan');
        return
    end
end
end

function [centers, diff_mat, n_mat] = local_diff_curves( ...
    X_cell, Y_cell, play_cell, hop_s, edges, min_play_s, min_nonplay_s, use_session)

centers = edges(1:end - 1) + diff(edges) / 2;
n_sess = numel(X_cell);
diff_mat = nan(n_sess, numel(centers));
n_mat = zeros(1, numel(centers));
for s = 1:n_sess
    if ~use_session(s) || isempty(X_cell{s})
        continue
    end
    X = X_cell{s};
    Y = Y_cell{s};
    is_play = play_cell{s};
    for b = 1:numel(centers)
        if b == numel(centers)
            in_bin = X >= edges(b) & X <= edges(b + 1);
        else
            in_bin = X >= edges(b) & X < edges(b + 1);
        end
        play_s = nnz(in_bin & is_play) * hop_s(s);
        nonplay_s = nnz(in_bin & ~is_play) * hop_s(s);
        if play_s >= min_play_s && nonplay_s >= min_nonplay_s
            diff_mat(s, b) = mean(Y(in_bin & is_play), 'omitnan') - mean(Y(in_bin & ~is_play), 'omitnan');
            n_mat(b) = n_mat(b) + 1;
        end
    end
end
end

function text = local_diff_text(diff_log, diff_z, passes)
if isnan(diff_log)
    text = 'no overlap';
    return
end
amp_pct = 100 * (10 ^ diff_log - 1);
if passes
    flag = 'passes minimum';
else
    flag = 'below the minimum duration';
end
text = sprintf('%+.4f log10 (%+.1f%% amplitude, %+.3f session z), %s', ...
    diff_log, amp_pct, diff_z, flag);
end

function local_print_difference_test(vals, animals, label)
vals = vals(:);
animals = animals(:);
keep = ~isnan(vals);
vals = vals(keep);
animals = animals(keep);
fprintf('%s: n = %d\n', label, numel(vals));
if numel(vals) < 2
    return
end
[~, p_val, ci, stats] = ttest(vals);
fprintf('  one-sample t: mean %.4f log10, 95%% CI [%.4f, %.4f], t(%d) = %.2f, p = %.4f\n', ...
    mean(vals), ci(1), ci(2), stats.df, stats.tstat, p_val);
if numel(unique(animals)) < 2
    fprintf('  LME skipped: only one animal\n');
    return
end
lme = fitlme(table(vals, categorical(animals), 'VariableNames', {'Difference', 'Animal'}), ...
    'Difference ~ 1 + (1|Animal)');
coef = fixedEffects(lme);
ci_lme = coefCI(lme);
fprintf('  LME Difference ~ 1 + (1|Animal): mean %.4f, 95%% CI [%.4f, %.4f], p = %.4f\n', ...
    coef, ci_lme(1), ci_lme(2), lme.Coefficients.pValue(1));
end

function local_plot_diff_curve(x, diff_mat, x_name, y_label, plot_title)
hold on
mean_diff = mean(diff_mat, 1, 'omitnan');
if ~any(isfinite(mean_diff))
    xlim([x(1) x(end)])
    ylim([-1 1])
    yline(0, 'k:')
    title(plot_title)
    xlabel(sprintf('%s (session z)', x_name))
    ylabel(y_label)
    text(mean(xlim), 0, 'No bin has 3 sessions', 'HorizontalAlignment', 'center')
    set(gca, 'Tag', 'diffcurve')
    return
end
p_val = nan(1, numel(x));
ci = nan(2, numel(x));
for b = 1:numel(x)
    y = diff_mat(:, b);
    y = y(~isnan(y));
    if numel(y) < 3
        continue
    end
    [~, p_val(b), ci(:, b)] = ttest(y);
end
finite_ci = ~any(isnan(ci), 1);
if any(finite_ci)
    fill([x(finite_ci) fliplr(x(finite_ci))], ...
        [ci(1, finite_ci) fliplr(ci(2, finite_ci))], ...
        'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
end
plot(x, mean_diff, 'k', 'LineWidth', 1.5)
yline(0, 'k:')
axis tight
y_lim = ylim;
if all(isfinite(y_lim)) && y_lim(2) > y_lim(1)
    fill_lo = y_lim(1) + 0.9 * (y_lim(2) - y_lim(1));
    fill_hi = y_lim(2);
    sig = p_val < 0.05;
    sig(isnan(p_val)) = false;
    start_end = [find(diff([0 sig 0]) == 1)' (find(diff([0 sig 0]) == -1) - 1)'];
    for j = 1:size(start_end, 1)
        fill(x(start_end(j, [1 2 2 1])), [fill_lo fill_lo fill_hi fill_hi], ...
            'k', 'EdgeColor', 'none', 'FaceAlpha', .35)
    end
end
xlabel(sprintf('%s (session z)', x_name))
ylabel(y_label)
title(plot_title)
set(gca, 'Tag', 'diffcurve')
end

function local_plot_coverage_count(x, n_mat, x_name, min_play_s, min_nonplay_s)
plot(x, n_mat, 'k', 'LineWidth', 1.5)
xlabel(sprintf('%s (session z)', x_name))
ylabel('Sessions in the bin')
title(sprintf('Bins with >= %.1f s play and >= %.1f s non-play', min_play_s, min_nonplay_s))
axis tight
end
