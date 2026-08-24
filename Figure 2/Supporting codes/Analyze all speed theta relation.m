
list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single','B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};
list_of_paertner = {[1 2],[1 2],[1 2],[1 2],[1 2],[1 2],[1 2 3]};


saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\psth power by frequency and behavior';


% Spectrogram parameters
% wind_length     = 1;    % delta
% wind_overlap    = 0.99;   % delta

wind_length     = .250;    % for theta
wind_overlap    = 0.249;   % for theta

% f          = .1:.1:5;      % frequency range for spectrogram delta
f          = 5:.1:16;       % frequency range for spectrogram thea
% freq_range = [.1 5];       % freq band delta
freq_range = [6 12];       % freq band theta



%% load 

load([saving_folder,'\psth_structure_speed_delta.mat'],'psth_structure')
load([saving_folder,'\animal_names_speed_delta.mat'],'animal_names')

%% some prelimiary plots

stasts_table = cell(numel(psth_structure),6);

for j=1:numel(psth_structure)

stasts_table(j,:) = { psth_structure(j).lm.Coefficients.Estimate(4), ...
                psth_structure(j).lm.Coefficients.pValue(4),...
                psth_structure(j).lm.Coefficients.tStat(4),...
                  animal_names{j,:}};

table = psth_structure(j).model_data;

end

stasts_table =cell2table(stasts_table);
stasts_table.Properties.VariableNames = {'Estimate','pValue','tStat','Animal','Partner','Electrode'};

stasts_table = stasts_table(stasts_table.Electrode==1,:);

lme = fitlme(stasts_table, 'Estimate ~ 1 + (1|Animal)');
coef = fixedEffects(lme);       % estimated mean
ci = coefCI(lme);               % confidence interval
pval = lme.Coefficients.pValue; % p-value for intercept

fprintf('Mean = %.3f, 95%% CI [%.3f, %.3f], p = %.4f\n', coef, ci(1), ci(2), pval);

%%
n_grid = 101;
variable_grid = linspace(-2, 2, n_grid);

%% 
data_together = [];
mean_powers_play    = nan(numel(psth_structure) ,n_grid);
mean_powers_noplay  = nan(numel(psth_structure) ,n_grid);
for j=1:numel(psth_structure) 


tbl = psth_structure(j).model_data;
tbl = tbl(~isnan(tbl.Speed),:);

Y = tbl.Power;
Y = (Y - mean(Y, 'omitmissing'))/std(Y, 'omitmissing');
X = tbl.RelativeSpeed;
Xa = X(tbl.Play=='true');
Ya = Y(tbl.Play=='true');


Xb = X(tbl.Play=='false');
Yb =  Y(tbl.Play=='false');


data_together = [data_together;[X Y]];

edges = [variable_grid Inf];

% Assign each sample to a bin
[~,~,binA] = histcounts(Xa, edges);
[~,~,binB] = histcounts(Xb, edges);

% Pre-allocate as NaN
meanA = nan(size(variable_grid));
meanB = nan(size(variable_grid));

% Compute mean per bin (only for bins with data)
for i = 1:numel(variable_grid)
    meanA(i) = mean(Ya(binA == i), 'omitnan');
    meanB(i) = mean(Yb(binB == i), 'omitnan');
end


mean_powers_play(j,:)    = meanA;
mean_powers_noplay(j,:)  = meanB;
end
%%
[bin_count, ~, samples2plot] = histcounts(data_together(:,1), -5:.1:10);

figure
bins_with_values = unique(samples2plot)'
max_count = 500;
min_count = 50;

indexes2plot = [];

for j=bins_with_values
    indexes= find(samples2plot==j);
   if numel(indexes)>=max_count

       indexes  = datasample(indexes,max_count );
   end
   if numel(indexes)>=min_count
       indexes2plot = [indexes2plot;indexes];
   end
end
selection_index = ~any(isnan(data_together(indexes2plot,:)),2)  & data_together(indexes2plot,1)<1 & abs(data_together(indexes2plot,2))<3;


subplot(1,2,1)
plot(data_together(indexes2plot(selection_index),1), data_together(indexes2plot(selection_index),2), '.k', 'MarkerSize',    .1)
[c,p] = corr(data_together(indexes2plot(selection_index),1),data_together(indexes2plot(selection_index),2));
title(num2str([c,p]))
axis xy

subplot(1,2,2)
plot(variable_grid, mean_powers_play-mean_powers_noplay, ':k');
hold on
plot(variable_grid, mean(mean_powers_play-mean_powers_noplay), 'k');


%%
selection_index = ~any(isnan(data_together(:,:)),2)  & data_together(:,1)<1 & abs(data_together(:,2))<3;

[c,p] = corr(data_together((selection_index),1),data_together((selection_index),2));
