%% Analyze_psth_full_spectrogram
% Compare pre- vs during-play-bout full-spectrum power using LME.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

%% Load precomputed spectrogram PSTH
disp('loading')
load([saving_folder,'\psth_structure_delta_full_spectrogram.mat'],'psth_structure');
load([saving_folder,'\animal_names_delta_full_spectrogram.mat'],'animal_names');
f = psth_structure(1).f;
disp('ready')

%% Aggregate pre/during/post power across sessions
animal_index = [1 1 1 2 2 2 3 3 3 4 5 5 6 6];
all_pre = [];
all_during = [];
all_post = [];
all_animal_index = [];
all_pre_all = [];
all_during_all = [];

for fn = [1:11, 13]
    pre_pow = psth_structure(fn).play_bout_tw_this(:, :, 1:500);
    during_pow = psth_structure(fn).play_bout_tw_this(:, :, 501:1500);
    post_pow = psth_structure(fn).play_bout_tw_this(:, :, 1501:2000);

    all_pre = [all_pre; squeeze(mean(mean(pre_pow, 1), 3))];
    all_during = [all_during; squeeze(mean(mean(during_pow, 1), 3))];
    all_post = [all_post; squeeze(mean(mean(post_pow, 1), 3))];

    all_pre_all = [all_pre_all; squeeze(mean(pre_pow, 3))];
    all_during_all = [all_during_all; squeeze(mean(during_pow, 3))];

    all_animal_index = [all_animal_index; repmat(animal_index(fn), size(during_pow, 1), 1)];
end

%% Normalize session-level spectra
all_pre_norm = log10(all_pre);
all_during_norm = log10(all_during);

for j = 1:size(all_pre, 1)
    all_pre_norm(j, :) = (all_pre_norm(j, :) - mean(all_pre_norm(j, :))) / std(all_pre_norm(j, :));
    all_during_norm(j, :) = (all_during_norm(j, :) - mean(all_pre_norm(j, :))) / std(all_pre_norm(j, :));
end

%% Plot trial-averaged difference with t-test significance bands
width = -10000;
power_diff = all_during_all - all_pre_all;
[~, p, ci] = ttest(power_diff);
h = p < 0.001;
original_y = mean(power_diff);

figure
fill([f fliplr(f)], [ci(1, :) fliplr(ci(2, :))], 'm', 'EdgeColor', 'none', 'FaceAlpha', .2)
hold on
plot(f, original_y, 'm')
axis tight
y_lim = ylim;

start_end = [find(diff([0 h 0]) == 1) find(diff([0 h 0]) == -1) - 1];
for j = 1:size(start_end, 1)
    fill(f(start_end(j, [1 2 2 1])), y_lim(2) + [0 0 width width], 'm')
end

%% LME at each frequency (random intercept for animal)
A = power_diff;
nFreq = size(A, 2);

beta_est = nan(nFreq, 1);
SE = nan(nFreq, 1);
CI_low_95 = nan(nFreq, 1);
CI_high_95 = nan(nFreq, 1);
CI_low_99 = nan(nFreq, 1);
CI_high_99 = nan(nFreq, 1);
pValue = nan(nFreq, 1);

for freq_n = 1:nFreq
    tbl = table(A(:, freq_n), categorical(all_animal_index), 'VariableNames', {'Power', 'Animal'});
    lme = fitlme(tbl, 'Power ~ 1 + (1|Animal)');

    coef = lme.Coefficients;
    ci_95 = coefCI(lme, 'Alpha', 0.05);
    ci_99 = coefCI(lme, 'Alpha', 0.01);

    beta_est(freq_n) = coef.Estimate(1);
    SE(freq_n) = coef.SE(1);
    CI_low_95(freq_n) = ci_95(1, 1);
    CI_high_95(freq_n) = ci_95(1, 2);
    CI_low_99(freq_n) = ci_99(1, 1);
    CI_high_99(freq_n) = ci_99(1, 2);
    pValue(freq_n) = coef.pValue(1);
end

sigFreq = find(pValue < 0.01);
fprintf('Significant frequencies (p < 0.01):\n');
disp(sigFreq)

%% Plot LME estimate with significant frequencies marked
figure
hold on
plot(f, beta_est, 'k', 'LineWidth', 2);
yline(0, 'k--');
plot(f(sigFreq), beta_est(sigFreq), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
xlabel('Frequency (Hz)')
ylabel('Estimated mean power')
title('LME estimate ± 95% confidence interval')
box off

%% Plot LME estimate vs session-mean difference
figure
sigFreq = find(pValue < 0.01);
plot(f, beta_est, 'k', 'LineWidth', 2)
hold on
plot(f, mean(all_during - all_pre, 1), ':k')
plot(f(sigFreq), beta_est(sigFreq), 'ro', 'MarkerFaceColor', 'r')
xlabel('Frequency (Hz)')
ylabel('Estimated mean power')
title('LME estimate (intercept) for each frequency')

%% Plot per-animal mean difference with t-test bands
width = -10000;
[~, p, ci] = ttest(all_during, all_pre, 'Alpha', 0.01);
h = p < 0.01;
original_y = mean(all_during - all_pre);

figure
hold on
plot(f, f * 0, ':k')
fill([f fliplr(f)], [ci(1, :) fliplr(ci(2, :))], 'm', 'EdgeColor', 'none', 'FaceAlpha', .2)
plot(f, original_y, 'm')
axis tight
y_lim = ylim;

start_end = [find(diff([0 h 0]) == 1) find(diff([0 h 0]) == -1) - 1];
for j = 1:size(start_end, 1)
    fill(f(start_end(j, [1 2 2 1])), y_lim(2) + [0 0 width width], 'm')
end
