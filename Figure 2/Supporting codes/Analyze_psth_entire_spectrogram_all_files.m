

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
figure_3_new_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure 3 Updated';


 % for  the entire spectrum
psth_structure = [];
wind_length     = 1;
wind_overlap    = .990;
min_separation = .200;
f               = 0.1:.1:30;
freq_pow_range  = [35 90];


%% load if needed
disp('loading')
load([saving_folder,'\psth_structure_delta_full_spectrogram.mat'],'psth_structure');
load([saving_folder,'\animal_names_delta_full_spectrogram.mat'],'animal_names');
disp('ready')

%%
animal_index = [1 1 1 2 2 2 3 3 3 4 5 5 6 6];
all_pre     = [];
all_during  = [];
all_post    = [];
all_animal_index = [];

all_pre_all     = [];
all_during_all  = [];


for fn=[1:11,13]

pre_pow         = psth_structure(fn).play_bout_tw_this(:, :, 1:500);
during_pow      = psth_structure(fn).play_bout_tw_this(:, :, 501:1500);
post_pow        = psth_structure(fn).play_bout_tw_this(:, :, 1501:2000);

all_pre     = [all_pre;squeeze(mean(mean(pre_pow,1),3))];
all_during  = [all_during;squeeze(mean(mean(during_pow,1),3))];
all_post    = [all_post;squeeze(mean(mean(post_pow,1),3))];


all_pre_all     = [all_pre_all;squeeze(mean(pre_pow,3))];
all_during_all  = [all_during_all;squeeze(mean(during_pow,3))];

all_animal_index = [all_animal_index;repmat(animal_index(fn),size(during_pow,1),1)];
end

all_pre_norm        = log10(all_pre);
all_during_norm     =log10(all_during);

for j=1:size(all_pre,1)
all_pre_norm(j,:) = (all_pre_norm(j,:) - mean(all_pre_norm(j,:)))/std(all_pre_norm(j,:));
all_during_norm(j,:) = (all_during_norm(j,:) - mean(all_pre_norm(j,:)))/std(all_pre_norm(j,:));
end



%% ploting sinlge playbout


lme(all_animal_index)

width = -10000;

[~,p,ci]        = ttest(all_during_all-all_pre_all);
h               = p<0.001;
original_y      = mean(all_during_all-all_pre_all);


figure

fill([f fliplr(f)], [ci(1,:) fliplr(ci(2,:))], 'm', 'EdgeColor','none', 'FaceAlpha',.2)
hold on
plot(f,original_y, 'm' )

axis tight
y_lim = ylim;


for j=1:size(start_end,1)
fill(f(start_end(j,[1 2 2 1])), y_lim(2) + [0 0 width width], 'm')
end

%% same but ussing a lme

%% Test whether power at each frequency differs from zero using a Linear Mixed-Effects Model
% A: nObservations x nFrequencies
% all_animal_index: nObservations x 1
A = all_during_all-all_pre_all;
%% Test whether power at each frequency differs from zero using a Linear Mixed-Effects Model
% A: nObservations x nFrequencies
% all_animal_index: nObservations x 1
A = all_during_all-all_pre_all;
nFreq = size(A,2);

beta_est = nan(nFreq,1);
SE        = nan(nFreq,1);
CI_low_95    = nan(nFreq,1);
CI_high_95   = nan(nFreq,1);
CI_low_99    = nan(nFreq,1);
CI_high_99   = nan(nFreq,1);
tStat     = nan(nFreq,1);
pValue    = nan(nFreq,1);

for freq_n = 1:nFreq

    % Build table
    tbl = table( ...
        A(:,freq_n), ...
        categorical(all_animal_index), ...
        'VariableNames',{'Power','Animal'});

    % Random intercept for animal
    lme = fitlme(tbl,'Power ~ 1 + (1|Animal)');

    % Fixed-effect statistics (intercept)
    coef = lme.Coefficients;
    ci_95 = coefCI(lme,'Alpha',0.05);   % 95% confidence intervals
    ci_99 = coefCI(lme,'Alpha',0.01);   

    beta_est(freq_n) = coef.Estimate(1);
    SE(freq_n)       = coef.SE(1);
    CI_low_95(freq_n)   = ci_95(1,1);
    CI_high_95(freq_n)  = ci_95(1,2);

     CI_low_99(freq_n)   = ci_99(1,1);
    CI_high_99(freq_n)  = ci_99(1,2);
    tStat(freq_n)    = coef.tStat(1);
    pValue(freq_n)   = coef.pValue(1);

end

%% Significant frequencies


fprintf('Significant frequencies (p < 0.05):\n');
disp(sigFreq)

%% Plot with 95% confidence intervals
%% Plot mean estimate with 95% confidence band
sigFreq = find(pValue < 0.01);
figure;
hold on

% Confidence interval (shaded)
% fill([f(:); flipud(f(:))], ...
%      [CI_low_95(:); flipud(CI_high_95(:))], ...
%      [0.8 0.8 0.8], ...
%      'EdgeColor','none', ...
%      'FaceAlpha',0.4);
% fill([f(:); flipud(f(:))], ...
%      [CI_low_99(:); flipud(CI_high_99(:))], ...
%      [0.8 0.8 0.8], ...
%      'EdgeColor','none', ...
%      'FaceAlpha',0.4);

% Mean estimate
plot(f, beta_est, 'k', 'LineWidth',2);

% Zero line
yline(0,'k--');

% Mark significant frequencies
plot(f(sigFreq), beta_est(sigFreq), ...
    'ro', 'MarkerFaceColor','r', 'MarkerSize',5);

xlabel('Frequency (Hz)')
ylabel('Estimated mean power')
title('LME estimate ± 95% confidence interval')

box off

%% Display significant frequencies


fprintf('Significant frequencies (p < 0.05):\n');
disp(sigFreq)

%% Optional plot
figure;
sigFreq = find(pValue < 0.01);
plot(f,beta_est,'k','LineWidth',2); hold on
plot(f,all_during-all_pre,':k')
plot(f(sigFreq),beta_est(sigFreq),'ro','MarkerFaceColor','r')
xlabel('Frequency index')
ylabel('Estimated mean power')
title('LME estimate (intercept) for each frequency')

%% ploting by mean per animal

width = -10000;

[~,p,ci] = ttest(all_during,all_pre, 'Alpha',0.01)
h = p<0.01;
original_y      = mean(all_during-all_pre);


figure
hold on
plot(f,f*0, ':k' )

fill([f fliplr(f)], [ci(1,:) fliplr(ci(2,:))], 'm', 'EdgeColor','none', 'FaceAlpha',.2)
plot(f,original_y, 'm' )



axis tight
y_lim = ylim;


for j=1:size(start_end,1)
fill(f(start_end(j,[1 2 2 1])), y_lim(2) + [0 0 width width], 'm')
end


%%
figure
[~,~,ci] = ttest(log10(all_during));
mean_diff = mean(log10(all_during));
fill([f fliplr(f)], [ci(1,:) fliplr(ci(2,:))], 'm', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(f, mean_diff, 'm')
hold on

[~,~,ci] = ttest((all_pre));
mean_diff = mean((all_pre));
fill([f fliplr(f)], [ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(f, mean_diff, 'k')


%%

[h,p,ci] = ttest(log10(all_during)-log10(all_pre));
mean_diff = mean(log10(all_during)-log10(all_pre));
figure
fill([f fliplr(f)], [ci(1,:) fliplr(ci(2,:))], 'm', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(f, mean_diff, 'm')
%%



[h,p,ci] = ttest((all_during_norm)-(all_pre_norm));
mean_diff = mean((all_during_norm)-(all_pre_norm));
figure
fill([f fliplr(f)], [ci(1,:) fliplr(ci(2,:))], 'm', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(f, mean_diff, 'm')



