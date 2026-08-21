%% Analyze_speed_freq_relation
% Theta power vs relative speed during play vs non-play bouts.

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';

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
plot(variable_grid, mean_powers_play - mean_powers_noplay, ':k');
hold on
plot(variable_grid, mean(mean_powers_play - mean_powers_noplay), 'k');

%% Full-dataset correlation
selection_index = ~any(isnan(data_together(:, :)), 2) & data_together(:, 1) < 1 & abs(data_together(:, 2)) < 3;
[c, p] = corr(data_together(selection_index, 1), data_together(selection_index, 2));
