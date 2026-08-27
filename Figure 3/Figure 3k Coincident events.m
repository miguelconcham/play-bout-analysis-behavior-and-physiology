%% Figure 3 — coincident events (empty-bar plot)
% Coincidence (surprise) change after vs before lag 0, per area pair and
% cell type (trough / peak / unlocked). Filled bars: 95% CI excludes 0;
% empty bars: not significant (same style as the supplementary figure).
% Former sections VI, IX, and XI+ in this file are omitted.

%% I define folders
% Coincidence structs, neuron tables, and NPX raw data are under Data\.

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = [data_root, '\Analysis results\Cross_correlogram'];
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure cross correlations';

bin_size = 0.01;
hist_range = [-.25 .5];
psth_edges = hist_range(1):bin_size:hist_range(2);
areas2analyse = {'DLPAG'	'DR'	'LPAG'	'SupCol'	'VLPAG'};

%% II load data
saving_folder = [data_root, '\Analysis results\Cross_correlogram'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PLAY   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([saving_folder,'\surpirse_stats_struct_structure_play_jittering_long_interval2.mat'],'synch_structure');
load([saving_folder,'\surpirse_stats_struct_animal_names_play_jittering_long_interval2.mat'],'animal_names');
animal_names(12,:) = [];
time_precision = 0.005;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  NON PLAY   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([saving_folder,'\surpirse_stats_struct_structure_NOPLAY_dittering_long_interval2_ALL.mat'],'synch_structure');
% load([saving_folder,'\surpirse_stats_struct_animal_names_NOPLAY_dittering_long_interval2_ALL.mat'],'animal_names');
%
% bin_size = 0.002;
% hist_range = [-1.5 1];
% psth_edges = hist_range(1):bin_size:hist_range(2);
% time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
% smoth_wind_sec  =0.05;

% load([saving_folder,'\surpirse_stats_struct_PLAY_2ms_dittering 15-Apr-2026 12_31_12.mat'],'synch_structure');
% load([saving_folder,'\animal_names_PLAY_2ms_dittering15-Apr-2026 12_31_12.mat'],'animal_names');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% for high precision case %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([saving_folder,'\surpirse_stats_struct_PLAY_dittering_long_interval_2ms_ALL.mat'],'synch_structure');
% load([saving_folder,'\surpirse_stats_animal_names__PLAY_dittering_long_interval_2ms_ALL.mat'],'animal_names');
% time_precision = 0.003;

psth_edges = synch_structure(1).psth_edges;
bin_size = mean(diff(psth_edges));
hist_range = psth_edges([1 end]);
smoth_wind_sec = 0.05;
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));

%% III merge synch structures
% Concatenate coincidence histograms, percentile CCGs, areas, and cluster IDs.
synch_spikes_histogram = [];
pctl_spikes_histogram = [];
all_area_comb = [];
session_id = [];
clusters_id = [];
cluster_fr = [];

for fn = 1:numel(synch_structure)
    synch_spikes_histogram          = cat(2,synch_spikes_histogram,synch_structure(fn).synch_spikes_histogram);
    snch_pctls_this_session         = synch_structure(fn).synch_spikes_pctl;
    smoothed_synch = synch_structure(fn).synch_spikes_histogram;

    for j = 1:size(smoothed_synch,2)
        for neuron = 1:2
            smoothed_synch(neuron,j,:) = movmean(smoothed_synch(neuron,j,:),smoth_wind_sec/bin_size);
        end
    end
    snch_pctls_this_session(smoothed_synch==0) = 1;
    snch_pctls_this_session(isnan(smoothed_synch)) = NaN;
    pctl_spikes_histogram = cat(2,pctl_spikes_histogram,snch_pctls_this_session);
    all_area_comb         = [all_area_comb;synch_structure(fn).cluster_info.area(synch_structure(fn).synch_comb)];
    session_id            = [session_id;repmat(animal_names(fn,1),size(synch_structure(fn).synch_comb,1),1)];
    clusters_id           = [clusters_id;synch_structure(fn).cluster_info.cluster_id(synch_structure(fn).synch_comb)];
    cluster_fr            = [cluster_fr;synch_structure(fn).cluster_info.fr(synch_structure(fn).synch_comb)];
end

%% IV obtain area combinations
% Unique undirected pairs among SupCol, DLPAG, LPAG, VLPAG, DR.
area_combinations = cell2table(all_area_comb);
area_combinations.Properties.VariableNames = {'Neuron1','Neuron2'};
unique_area_combinations = unique(area_combinations, 'rows');

relevant_areas = {'SupCol','DLPAG','LPAG','VLPAG','DR'};
unique_area_combinations = unique_area_combinations(ismember(unique_area_combinations.Neuron1,relevant_areas) & ismember(unique_area_combinations.Neuron2,relevant_areas),:);

data = string(table2array(unique_area_combinations));
sortedRows = sort(data, 2);
[~, ~, ic] = unique(sortedRows, 'rows');
rowGroups = accumarray(ic, (1:size(data,1))', [], @(x) {sort(x')});
mirroredPairs = rowGroups(cellfun(@(x) length(x) > 1, rowGroups));

for j = 1:size(mirroredPairs,1)
    unique_area_combinations{mirroredPairs{j},:} = repmat(unique_area_combinations{mirroredPairs{j}(1),:},2,1);
end

unique_area_combinations = unique(unique_area_combinations, 'rows');

%% V load delta and theta phase data
% Map each coincidence pair onto all_neurons_TD and label trough / peak / unlocked.
saving_folder = [data_root, '\Analysis results\phase locking data'];

load([saving_folder,'\theta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
all_neurons_TD = all_neurons;
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner1'}))      = {'ThetaPartner1'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner2'}))      = {'ThetaPartner2'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Play'}))          = {'ThetaPlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'PrePlay'}))       = {'ThetaPrePlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'EntireSession'})) = {'ThetaEntireSession'};

load([saving_folder,'\delta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) = {'isRt'};
all_neurons_TD.DeltaPartner1      = all_neurons.Partner1;
all_neurons_TD.DeltaPartner2      = all_neurons.Partner2;
all_neurons_TD.DeltaEntireSession = all_neurons.EntireSession;
all_neurons_TD.DeltaPlay          = all_neurons.Play;
all_neurons_TD.DeltaPrePlay       = all_neurons.PrePlay;
all_neurons_TD.Exited             = nan(size(all_neurons_TD,1),1);
all_neurons_TD.Inhibited          = nan(size(all_neurons_TD,1),1);

nRows = size(clusters_id,1);
idx_pairs = zeros(nRows,2);

for i = 1:nRows
    s = session_id{i};
    id_pair = clusters_id(i,:);
    for j = 1:2
        if any(strcmp(all_neurons_TD.session, s) & ...
                all_neurons_TD.cluster_id == id_pair(j), ...
                1)
            idx_pairs(i,j) = find( ...
                strcmp(all_neurons_TD.session, s) & ...
                all_neurons_TD.cluster_id == id_pair(j), ...
                1);
        else
            disp('Missing neuron')
            idx_pairs(i,j) = NaN;
        end
    end
end

%% VII Compute effect distributions
% Per pair: coincidence after lag 0 minus coincidence before lag 0.
% Rows go into table4mixed_effect as [delta, group, area-pair].
% Group 1 = trough–trough, 2 = peak–peak, 3 = unlocked–unlocked.

time_ranges_2_measure = [-1.4 0; 0 .3];
border_effect = [-1.4 .9];
baseline_range = [-1.4 0];
smooth_wind_sec = 0.05;
min_fr = 5;
fix_size = 50;
n_montecarlo = 5000;

trough        = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained = all_neurons_TD.DeltaEntireSession.PPCPval>.1;

entrained_group1     = trough;
entrained_group2     = trough;
comparison_group1    = peak;
comparison_group2    = peak;
non_entrained_group1 = non_entrained;
non_entrained_group2 = non_entrained;

table4mixed_effect = [];
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
time2consider = time_centers>border_effect(1) & time_centers<border_effect(2);
sub_time_centers = time_centers(time2consider);

n_counts_used = nan(size(unique_area_combinations,1),6);

for j = 1:size(unique_area_combinations,1)

    area_combination_indexes = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main     = find(entrained_group1(idx_pairs(:,1)) & entrained_group2(idx_pairs(:,2)) & area_combination_indexes & mean(cluster_fr,2)>min_fr);
    cell_type_comb_comp     = find(comparison_group1(idx_pairs(:,1)) & comparison_group2(idx_pairs(:,2)) & area_combination_indexes & mean(cluster_fr,2)>min_fr);
    cell_type_comb_nonentr  = find(non_entrained_group1(idx_pairs(:,1)) & non_entrained_group2(idx_pairs(:,2)) & area_combination_indexes & mean(cluster_fr,2)>min_fr);

    min_val = min(numel(cell_type_comb_main),numel(cell_type_comb_comp));

    if min_val>fix_size && nchoosek(min_val, fix_size)>n_montecarlo
        K = fix_size;
        n_counts_used(j,3) = min(n_montecarlo,K);
    elseif min_val>0 && min_val<fix_size
        K = round(.5*min_val);
        n_counts_used(j,3) = min(n_montecarlo,K);
    else
        K = 0;
    end
    n_counts_used(j,6) = K;
    disp(K)

    % Trough
    if K>0
        this_surpise = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,time_centers>border_effect(1) & time_centers<border_effect(2))));
        if K<=min_val
            n_counts_used(j,1) = min(n_montecarlo,nchoosek(size(this_surpise,1), K));
            n_counts_used(j,4) = size(this_surpise,1);
        end
        for row_n = 1:size(this_surpise,1)
            this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
        end
        time_slice_1 = sub_time_centers>=time_ranges_2_measure(1,1) & sub_time_centers<time_ranges_2_measure(1,2);
        time_slice_2 = sub_time_centers>=time_ranges_2_measure(2,1) & sub_time_centers<time_ranges_2_measure(2,2);
        table4mixed_effect = [table4mixed_effect;[mean(this_surpise(:,time_slice_2),2)-mean(this_surpise(:,time_slice_1),2) ones(size(this_surpise,1),2)*[1  0;0 j]]];
    end

    % Peak
    this_surpise = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_comp,time_centers>border_effect(1) & time_centers<border_effect(2))));
    for row_n = 1:size(this_surpise,1)
        this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
    end
    if K>0
        if K<min_val
            n_counts_used(j,2) = min(n_montecarlo,nchoosek(size(this_surpise,1), K));
            n_counts_used(j,5) = size(this_surpise,1);
        end
        time_slice_1 = sub_time_centers>=time_ranges_2_measure(1,1) & sub_time_centers<time_ranges_2_measure(1,2);
        time_slice_2 = sub_time_centers>=time_ranges_2_measure(2,1) & sub_time_centers<time_ranges_2_measure(2,2);
        if min(size(this_surpise))>1
            table4mixed_effect = [table4mixed_effect;[mean(this_surpise(:,time_slice_2),2)-mean(this_surpise(:,time_slice_1),2) ones(size(this_surpise,1),2)*[2  0;0 j]]];
        else
            table4mixed_effect = [table4mixed_effect;[mean(this_surpise(time_slice_2))-mean(this_surpise(time_slice_1)) [2 j]]];
        end
    end

    % Unlocked
    this_surpise = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_nonentr,time_centers>border_effect(1) & time_centers<border_effect(2))));
    for row_n = 1:size(this_surpise,1)
        this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
    end
    if K>0 && size(this_surpise,1)>K
        if K<min_val
            n_counts_used(j,2) = min(n_montecarlo,nchoosek(size(this_surpise,1), K));
            n_counts_used(j,5) = size(this_surpise,1);
        end
        time_slice_1 = sub_time_centers>=time_ranges_2_measure(1,1) & sub_time_centers<time_ranges_2_measure(1,2);
        time_slice_2 = sub_time_centers>=time_ranges_2_measure(2,1) & sub_time_centers<time_ranges_2_measure(2,2);
        if min(size(this_surpise))>1
            table4mixed_effect = [table4mixed_effect;[mean(this_surpise(:,time_slice_2),2)-mean(this_surpise(:,time_slice_1),2) ones(size(this_surpise,1),2)*[3  0;0 j]]];
        else
            table4mixed_effect = [table4mixed_effect;[mean(this_surpise(time_slice_2))-mean(this_surpise(time_slice_1)) [3 j]]];
        end
    end
end

%% VIII Mixed / Linear model analysis
% Estimated means + 95% CI for group × area-pair.

tbl = array2table(table4mixed_effect, 'VariableNames',{'value','group','condition'});
alpha_level = 0.005;

tbl.group = categorical(tbl.group);
tbl.condition = categorical(tbl.condition);

hasSubject = ismember('subject', tbl.Properties.VariableNames);
if hasSubject
    tbl.subject = categorical(tbl.subject);
end

if hasSubject
    lme = fitlme(tbl, 'value ~ group*condition + (1|subject)');
    modelObj = lme;
    fprintf('\n=============================\n');
    fprintf('MIXED EFFECTS MODEL FITTED\n');
    fprintf('=============================\n');
    disp(lme)
    fprintf('\nANOVA table:\n');
    anovaTbl = anova(lme);
    disp(anovaTbl)
else
    lm = fitlm(tbl, 'value ~ group*condition');
    modelObj = lm;
    fprintf('\n=============================\n');
    fprintf('LINEAR MODEL FITTED\n');
    fprintf('=============================\n');
    disp(lm)
    fprintf('\nANOVA table:\n');
    anovaTbl = anova(lm, 'summary');
    disp(anovaTbl)
end

groupCats = categories(tbl.group);
condCats  = categories(tbl.condition);
[groupGrid, condGrid] = ndgrid(groupCats, condCats);

predTbl = table;
predTbl.group = categorical(groupGrid(:), groupCats);
predTbl.condition = categorical(condGrid(:), condCats);

if hasSubject
    subjCats = categories(tbl.subject);
    predTbl.subject = categorical(repmat(subjCats(1), height(predTbl), 1), subjCats);
end

if hasSubject
    [yhat, yCI, ~] = predict(lme, predTbl);
else
    [yhat, yCI] = predict(lm, predTbl, 'Alpha', alpha_level);
end

predTbl.estimate = yhat;
predTbl.CI_low   = yCI(:,1);
predTbl.CI_high  = yCI(:,2);
predTbl.aboveZero = predTbl.CI_low > 0 | predTbl.CI_high < 0;

predTbl.group_num = str2double(string(predTbl.group));
predTbl.cond_num  = str2double(string(predTbl.condition));
predTbl = sortrows(predTbl, {'group_num','cond_num'});

fprintf('\n=============================\n');
fprintf('ESTIMATED MEANS + 95%% CI\n');
fprintf('=============================\n');
disp(predTbl(:, {'group','condition','estimate','CI_low','CI_high','aboveZero'}))
estimated_means_table = predTbl;

% Combinations with no observations should not be plotted
predTbl.Nobs = zeros(height(predTbl),1);
for r = 1:height(predTbl)
    idx = tbl.group == predTbl.group(r) & tbl.condition == predTbl.condition(r);
    predTbl.Nobs(r) = sum(idx);
end
missingIdx = predTbl.Nobs == 0;
predTbl.estimate(missingIdx)  = NaN;
predTbl.CI_low(missingIdx)    = NaN;
predTbl.CI_high(missingIdx)   = NaN;
predTbl.aboveZero(missingIdx) = false;

%% Fig 3k: X Plot (with empty bar plots, as on Supp Figure)
% Filled bar: 95% CI excludes 0. Empty bar: not significant.

figure('Color','w','Position',[100 100 1000 500]);
hold on

x_tick_lables = cell(size(unique_area_combinations,1),1);
for j = 1:size(unique_area_combinations,1)
    x_tick_lables{j} = [unique_area_combinations.Neuron1{j}, ' ',unique_area_combinations.Neuron2{j}];
end
groups = categories(predTbl.group);
offset = linspace(-0.18, 0.18, numel(groups));

color_code = {
    [0.2 0.4 0.8]
    [0.8 0.2 0.3]
    [0.2 0.8 0.3]
};

barWidth = 0.28;

for g = 1:numel(groups)

    idx = predTbl.group == groups{g};

    x   = predTbl.cond_num(idx);
    y   = predTbl.estimate(idx);
    lo  = predTbl.CI_low(idx);
    hi  = predTbl.CI_high(idx);
    sig = predTbl.aboveZero(idx);

    [x, sortIdx] = sort(x);
    y   = y(sortIdx);
    lo  = lo(sortIdx);
    hi  = hi(sortIdx);
    sig = sig(sortIdx);

    validIdx = ~isnan(y) & ~isnan(lo) & ~isnan(hi);

    xPlot   = x(validIdx) + offset(g);
    yPlot   = y(validIdx);
    loPlot  = lo(validIdx);
    hiPlot  = hi(validIdx);
    sigPlot = sig(validIdx);

    for mm = 1:numel(xPlot)
        if sigPlot(mm)
            bar(xPlot(mm), yPlot(mm), barWidth, ...
                'FaceColor', color_code{g}, ...
                'FaceAlpha', 0.2, ...
                'EdgeColor', color_code{g}, ...
                'LineWidth', 1.2);
        else
            bar(xPlot(mm), yPlot(mm), barWidth, ...
                'FaceColor', 'none', ...
                'EdgeColor', color_code{g}, ...
                'LineWidth', 1.2);
        end
    end

    errorbar(xPlot, yPlot, yPlot-loPlot, hiPlot-yPlot, ...
        'k.', ...
        'LineWidth', 1.5, ...
        'CapSize', 8);
end

yline(0, 'k--', 'LineWidth', 1.2);

xlabel('Condition', 'FontSize', 12)
ylabel('Estimated mean value', 'FontSize', 12)
title('Model-estimated means with 95% confidence intervals', 'FontSize', 14)

xticks(1:15)
xticklabels(x_tick_lables)

legend(arrayfun(@(x) sprintf('Group %s', x), string(groups), 'UniformOutput', false), ...
    'Location', 'best')

box off
set(gca, 'FontSize', 11, 'LineWidth', 1.2)

% print(gcf,'-vector','-dsvg',[figure_folder,'/ bar plot with significnat increase or decrease  synchrony.svg'])
