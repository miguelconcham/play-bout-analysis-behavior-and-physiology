%% Analyze_speed_freq_relation
% Theta power vs relative speed during play vs non-play bouts.

this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Event aligned LFP\Analyze_speed_freq_relation.m';
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
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
