
%% I define folders
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = [data_root, '\Analysis results\Cross_correlogram'];
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Codes\Figure codes\Figure cross correlations';
% mkdir(saving_folder)


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
%%%% for high precission case%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% load([saving_folder,'\surpirse_stats_struct_PLAY_dittering_long_interval_2ms_ALL.mat'],'synch_structure');
% load([saving_folder,'\surpirse_stats_animal_names__PLAY_dittering_long_interval_2ms_ALL.mat'],'animal_names');
% time_precision = 0.003;

psth_edges =synch_structure(1).psth_edges;
bin_size = mean(diff(psth_edges));
% bin_size = 0.001;
% bin_size = 0.002;
 hist_range = psth_edges([1 end]);
% hist_range = [-1.5 1];
smoth_wind_sec = 0.05;
% psth_edges = hist_range(1):bin_size:hist_range(2);
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
%% III merge synch structures
synch_spikes_histogram = [];
pctl_spikes_histogram = [];

all_area_comb = [];
session_id = [];
clusters_id = [];
cluster_fr = [];


for fn = 1:numel(synch_structure)

    synch_spikes_histogram          = cat(2,synch_spikes_histogram,synch_structure(fn).synch_spikes_histogram);
    snch_pctls_this_session         = synch_structure(fn).synch_spikes_pctl ;
    smoothed_synch = synch_structure(fn).synch_spikes_histogram;

    for j=1:size(smoothed_synch,2)
        for neuron = 1:2
        smoothed_synch(neuron,j,:) = movmean( smoothed_synch(neuron,j,:),smoth_wind_sec/bin_size);
        end
    end
    snch_pctls_this_session(smoothed_synch==0 )=1;
    snch_pctls_this_session(isnan(smoothed_synch))= NaN;
    pctl_spikes_histogram           = cat(2,pctl_spikes_histogram,snch_pctls_this_session       );
    all_area_comb                   = [all_area_comb;synch_structure(fn).cluster_info.area(synch_structure(fn).synch_comb)];
    session_id                      = [session_id;repmat(animal_names(fn,1),size(synch_structure(fn).synch_comb,1),1)];
    clusters_id                     = [clusters_id;synch_structure(fn).cluster_info.cluster_id(synch_structure(fn).synch_comb)];
    cluster_fr                      = [cluster_fr;synch_structure(fn).cluster_info.fr(synch_structure(fn).synch_comb)];

end


%% IV obtain area combinations

area_combinations = cell2table(all_area_comb);
area_combinations.Properties.VariableNames = {'Neuron1','Neuron2'};
unique_area_combinations = unique(area_combinations, 'rows');


relevant_areas = {'SupCol','DLPAG','LPAG','VLPAG','DR'};
unique_area_combinations = unique_area_combinations(ismember(unique_area_combinations.Neuron1,relevant_areas) & ismember(unique_area_combinations.Neuron2,relevant_areas),:);


data = string(table2array(unique_area_combinations));
sortedRows = sort(data, 2);

% 2. Get a unique ID for every unique combination
[~, ~, ic] = unique(sortedRows, 'rows');

% 3. Group the original row indexes based on those unique IDs
% 'rowGroups' will be a cell array where each cell contains the indexes for a combination
rowGroups = accumarray(ic, (1:size(data,1))', [], @(x) {sort(x')});

% 4. Filter to keep only the groups with more than 1 index (the mirrored pairs)
mirroredPairs = rowGroups(cellfun(@(x) length(x) > 1, rowGroups));

for j=1:size(mirroredPairs,1)

    unique_area_combinations{mirroredPairs{j},:} = repmat(unique_area_combinations{mirroredPairs{j}(1),:},2,1);
end

unique_area_combinations = unique(unique_area_combinations, 'rows');
%% V load delta and theta phase data
saving_folder = [data_root, '\Analysis results\phase locking data'];

load([saving_folder,'\theta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
all_neurons_TD = all_neurons;
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner1'}))         = {'ThetaPartner1'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner2'}))         = {'ThetaPartner2'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Play'}))             = {'ThetaPlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'PrePlay'}))          = {'ThetaPrePlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'EntireSession'}))    = {'ThetaEntireSession'};

load([saving_folder,'\delta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'}))  =     {'isRt'  };
all_neurons_TD.DeltaPartner1                            = all_neurons.Partner1;
all_neurons_TD.DeltaPartner2                            = all_neurons.Partner2;
all_neurons_TD.DeltaEntireSession                       = all_neurons.EntireSession;
all_neurons_TD.DeltaPlay                                = all_neurons.Play;
all_neurons_TD.DeltaPrePlay                             = all_neurons.PrePlay;
all_neurons_TD.Exited                                   = nan(size(all_neurons_TD,1),1);
all_neurons_TD.Inhibited                                = nan(size(all_neurons_TD,1),1);

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
                1);  % only one match expected
        else
            disp('Missing neuron')
            idx_pairs(i,j)= NaN;
        end
    end
end

%% checlk general statistics

general_pval = squeeze(mean(pctl_spikes_histogram));

figure
subplot(2,2,1)
histogram(1-general_pval(:), 0:0.01:1)
ylabel('Count (log-scale)')
xlabel('Coincidence level')
yscale log
subplot(2,2,2)
swarmchart(round(mean(cluster_fr,2),2),1-mean(general_pval,2), '.')
ylabel('Coincidence level')
xlabel('Fr (Hz)')



subplot(2,2,3)

all_synch_count = synch_spikes_histogram;
for j=1:size(all_synch_count,2)
    for nn=1:2
    all_synch_count(nn,j,:) = movmean(all_synch_count(nn,j,:), smoth_wind_sec/bin_size);
    end
end
all_synch_count = squeeze(mean(all_synch_count));

nan_counts = sum(isnan(all_synch_count),2);
plot(mean(all_synch_count,2), 1-mean(general_pval,2), '.')
xlabel('mean coincidence count')
ylabel('mean coincidence level')


subplot(2,2,4)
histogram(mean(all_synch_count(mean(general_pval,2)<0.05*10),2), 0:0.0001:0.04, 'EdgeColor','none')
% 
figure

range_95 = prctile (all_synch_count, [5 95]);
range_95(2,:) = movmean(range_95(2,:),50);
plot(time_centers, range_95)
hold on
plot(time_centers,median(all_synch_count(mean(general_pval,2)<0.05*10,:)))
yscale log
% 





%% VI plot statistics
y_lim = [0 1];
% y_lim_AUC = [.25 .65];
y_lim_AUC = [0 .65];
alpha_level         = 0.05/10;

% dotted_center = 0;
x_lim = [-1.5 1];
% border_effect = [-.9 .5];
% baseline_range = [-.9 -.5];

 border_effect = [-1.4 .9];
baseline_range = [-1.4 0];

baseline_test = true;
smooth_wind_sec = 0.05;

selected_tail = 'right';
% border_effect = [-.25 .5];
% baseline_range = [-.2 -.1];
smooth_pctl         = true;
smooth_vis          = false;
sign_widht          = .2;
plot_AUC            = true;
plot_median         = true;
center_auc          = false;

if center_auc
    y_lim_AUC = y_lim_AUC-.5;
    dotted_center = 0;
end
min_fr              = 5;
fr_range            = [5 40];
general_baseline    = false;
pair_count          = false;
var_range           = .1;

trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% below you select what to compare %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entrained_group1 = trough;
entrained_group2 = trough;
comparison_group1 = peak;
comparison_group2 = peak;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
time2consider = time_centers>border_effect(1)  & time_centers<border_effect(2);

sub_time_centers = time_centers(time2consider);

plot_color = {'b', 'r'};



AUC_PER_CONDITION   = nan(2,size(unique_area_combinations,1),sum(time2consider));
PVAL_PER_CONDITION  = nan(2,size(unique_area_combinations,1),sum(time2consider));

figure
for j=1:size(unique_area_combinations,1)

    subplot(5,3,j)
    hold on

    area_combination_indexes      = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main          = find(entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes & mean(cluster_fr,2)>fr_range(1) & mean(cluster_fr,2)<fr_range(2));
     % cell_type_comb_main          = find(((entrained_group1(idx_pairs(:,1)) &  comparison_group1(idx_pairs(:,2))) | (comparison_group1(idx_pairs(:,1)) &  entrained_group1(idx_pairs(:,2)))) & area_combination_indexes & mean(cluster_fr,2)>min_fr);
    cell_type_comb_comp          = find(comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) & area_combination_indexes  & mean(cluster_fr,2)>fr_range(1) & mean(cluster_fr,2)<fr_range(2));


    if general_baseline
        baseline_surpirse = squeeze(mean(pctl_spikes_histogram(:, area_combination_indexes,time_centers<=baseline_range(2) & time_centers>=baseline_range(1))));
        baseline_surpirse = baseline_surpirse(:);
    end


    if pair_count
        min_count = min(numel(cell_type_comb_main), numel(cell_type_comb_comp));

        cell_type_comb_main = randsample(cell_type_comb_main,min(ceil((1+var_range)*min_count), numel(cell_type_comb_main)), false);
        cell_type_comb_comp = randsample(cell_type_comb_comp,min(ceil((1+var_range)*min_count), numel(cell_type_comb_comp)), false);
    end

    this_surpise =1- squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,:)));
    if smooth_pctl
        for row_n=1:size(this_surpise,1)
            this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
        end
    end
    AUC_over_time =  nan(1,sum(sum(time2consider)));
    is_above_chance = nan(1,sum(sum(time2consider)));
    if min(size(this_surpise))>1        
        surprise = this_surpise(:, time2consider);
        if ~general_baseline
            if min(size(surprise))>1
                baseline_surpirse = surprise(:, sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
            elseif min(size(surprise))>0
                baseline_surpirse = surprise(sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
            end
            baseline_surpirse = median(baseline_surpirse,2, 'omitmissing');
        end
        for t = 1:size(surprise,2)
            [ is_above_chance(t), ~, stats] =ranksum(surprise(:, t), baseline_surpirse, 'tail',selected_tail);

            if ~baseline_test
                [ is_above_chance(t),] =signrank(surprise(:, t), .5, 'tail','right');
            end
            if ~plot_median
                n_obs = length(surprise(:,t));
                n_base = length(baseline_surpirse);
                U = stats.ranksum - n_obs*(n_obs+1)/2;
                AUC = U / (n_obs * n_base);
                AUC_over_time(t) = AUC;
            end
        end
        if plot_median
            dotted_center = median(baseline_surpirse, 'omitmissing');
            AUC_over_time = median(surprise, 'omitmissing');
        end
        if plot_AUC
            if smooth_vis
                AUC_over_time = movmean(AUC_over_time, smooth_wind_sec/bin_size);
            end
            if center_auc
                AUC_over_time = AUC_over_time-mean(AUC_over_time(sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1)));
            end
            plot(sub_time_centers, AUC_over_time, plot_color{1})

            hold on
            ylim(y_lim_AUC)
            y_lim = y_lim_AUC;
            fill_widht = sign_widht*range(y_lim);
            hold on
            plot(border_effect, [1 1]*dotted_center, [':',plot_color{1}])

            start_end = [find(diff([0 (is_above_chance<alpha_level) 0])==1)' find(diff([0 is_above_chance<alpha_level 0])==-1)'-1];
            if ~isempty(start_end)
                for event_n =1:size(start_end,1)
                    fill(sub_time_centers(start_end(event_n,[1 2 2 1])), y_lim(2)*[1 1 1 1] + [0 0 -fill_widht -fill_widht], plot_color{1}, 'FaceAlpha',.5, 'EdgeColor','none')
                end
            end
        else
            semilogy(time_centers(time2consider), 1-is_above_chance, plot_color{1})
            hold on
            semilogy(border_effect, 1-[alpha_level alpha_level], 'K')
        end
    end

    AUC_PER_CONDITION(1,j,:)  = AUC_over_time;
    PVAL_PER_CONDITION(1,j,:) = is_above_chance;


    this_surpise = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_comp,:)));
     if smooth_pctl
        for row_n=1:size(this_surpise,1)
            this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
        end
    end
    AUC_over_time =  nan(1,sum(sum(time2consider)));
    is_above_chance = nan(1,sum(sum(time2consider)));
    if min(size(this_surpise))>1        
        surprise = this_surpise(:, time2consider);
        if ~general_baseline
            if min(size(surprise))>1
                baseline_surpirse = surprise(:, sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
            elseif min(size(surprise))>0
                baseline_surpirse = surprise(sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
            end
            baseline_surpirse = median(baseline_surpirse,2, 'omitmissing');
        end
        for t = 1:size(surprise,2)
            [ is_above_chance(t), ~, stats] =ranksum(surprise(:, t), baseline_surpirse, 'tail',selected_tail);

            if ~baseline_test
                [ is_above_chance(t),] =signrank(surprise(:, t), .5, 'tail','right');
            end
            if ~plot_median
                n_obs = length(surprise(:,t));
                n_base = length(baseline_surpirse);
                U = stats.ranksum - n_obs*(n_obs+1)/2;
                AUC = U / (n_obs * n_base);
                AUC_over_time(t) = AUC;
            end
        end
        if plot_median
            dotted_center = median(baseline_surpirse, 'omitmissing');
            AUC_over_time = median(surprise, 'omitmissing');
        end
        if plot_AUC
            if smooth_vis
                AUC_over_time = movmean(AUC_over_time, smooth_wind_sec/bin_size);
            end
            if center_auc
                AUC_over_time = AUC_over_time-mean(AUC_over_time(sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1)));
            end
            plot(sub_time_centers, AUC_over_time, plot_color{2})

            hold on
            ylim(y_lim_AUC)
            y_lim = y_lim_AUC;
            fill_widht = sign_widht*range(y_lim);
            hold on
            plot(border_effect, [1 1]*dotted_center, [':',plot_color{2}])

            start_end = [find(diff([0 (is_above_chance<alpha_level) 0])==1)' find(diff([0 is_above_chance<alpha_level 0])==-1)'-1];
            if ~isempty(start_end)
                for event_n =1:size(start_end,1)
                    fill(sub_time_centers(start_end(event_n,[1 2 2 1])), y_lim(1)*[1 1 1 1] + [0 0 fill_widht fill_widht], plot_color{2}, 'FaceAlpha',.5, 'EdgeColor','none')
                end
            end
        else
            semilogy(time_centers(time2consider), 1-is_above_chance, plot_color{2})
            hold on
            semilogy(border_effect, 1-[alpha_level alpha_level], 'K')
        end
    end
    AUC_PER_CONDITION(2,j,:)  = AUC_over_time;
    PVAL_PER_CONDITION(2,j,:) = is_above_chance;



    title([unique_area_combinations.Neuron1{j}, ' G1:', num2str(numel(cell_type_comb_main)), ' G2:', num2str(numel(cell_type_comb_comp))])
    ylabel(unique_area_combinations.Neuron2{j})
    xlim(x_lim)
    pause(.1)

end

%% save if needed
print(gcf,'-vector','-dsvg',[figure_folder,'/mean conincidence per area peak and trough short xlim NON PLAY.svg'])

%%  VII Compute effect distributions

y_lim = [0 1];
y_lim_AUC = [0 .8];
dotted_center = 0;
% border_effect = [-.9 .5];
% baseline_range = [-.9 -.'5];
time_ranges_2_measure = [-1.4 0;0 .3];
 border_effect = [-1.4 .9];
 baseline_range = [-1.4 0];
 do_montecarlo = false;

baseline_test = true;
smooth_wind_sec = 0.05;

selected_tail = 'right';
% border_effect = [-.25 .5];
% baseline_range = [-.2 -.1];
smooth_pctl         = true;
smooth_vis          = false;
sign_widht          = .2;
plot_AUC            = true;
plot_median         = true;
center_auc          = false;

if center_auc
    y_lim_AUC = y_lim_AUC-.5;
    dotted_center = 0;
end
min_fr              = 5;
general_baseline    = false;
alpha_level         = 0.057/15;
pair_count          = false;
var_range           = .1;

trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;


entrained_group1 = trough;
entrained_group2 = trough;
comparison_group1 = peak;
comparison_group2 = peak;
non_entrained_group1 =non_entrained;
non_entrained_group2 =non_entrained;

table4mixed_effect = [];
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
time2consider = time_centers>border_effect(1)  & time_centers<border_effect(2);

sub_time_centers = time_centers(time2consider);
fix_size    = 50;

n_montecarlo = 5000;


mean_base           = nan(2,n_montecarlo,size(unique_area_combinations,1));
mean_effect         = nan(2,numel(time_ranges_2_measure),n_montecarlo,size(unique_area_combinations,1));

n_counts_used = nan(size(unique_area_combinations,1),6);


for j=1:size(unique_area_combinations,1)

    

    area_combination_indexes      = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main          = find(entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes & mean(cluster_fr,2)>min_fr);

    cell_type_comb_comp          = find(comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) & area_combination_indexes  & mean(cluster_fr,2)>min_fr);


    cell_type_comb_nonentr          = find(non_entrained_group1(idx_pairs(:,1)) &  non_entrained_group2(idx_pairs(:,2)) & area_combination_indexes  & mean(cluster_fr,2)>min_fr);

    min_val = min(numel(cell_type_comb_main),numel(cell_type_comb_comp));


    if  min_val>fix_size && nchoosek(min_val, fix_size)>n_montecarlo
        K = fix_size;
        n_counts_used(j,3) =min(n_montecarlo,K);
    elseif min_val>0 && min_val<fix_size
        K = round(.5*min_val);
        n_counts_used(j,3) =min(n_montecarlo,K);
    else
        K=0;
    end
    n_counts_used(j,6) = K;
    disp(K)
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%% Trough %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if K>0      

        this_surpise =1- squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,time_centers>border_effect(1)  & time_centers<border_effect(2))));
         if K<=min_val
            n_counts_used(j,1) =min(n_montecarlo,nchoosek(size(this_surpise,1), K));
             n_counts_used(j,4) =size(this_surpise,1);
            
        end
        for row_n=1:size(this_surpise,1)
            this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
        end


        if min(size(this_surpise))>1
            baseline_surpirse = this_surpise(:, sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
        elseif min(size(this_surpise))>0
            baseline_surpirse = this_surpise(sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
        end
        baseline_surpirse = median(baseline_surpirse,2, 'omitmissing');
        time_slice_1 = sub_time_centers>=time_ranges_2_measure(1,1) &  sub_time_centers<time_ranges_2_measure(1,2);
        time_slice_2 = sub_time_centers>=time_ranges_2_measure(2,1) &  sub_time_centers<time_ranges_2_measure(2,2);

        table4mixed_effect = [table4mixed_effect;[mean(this_surpise(:,time_slice_2),2)-mean(this_surpise(:,time_slice_1),2) ones(size(this_surpise,1),2)*[1  0;0 j]]];
        if do_montecarlo
        for nm =1:n_counts_used(j,1)
            for tn = 1:size(time_ranges_2_measure,1)
                time_slice = sub_time_centers>=time_ranges_2_measure(tn,1) &  sub_time_centers<time_ranges_2_measure(tn,2);
                selection_slice = randsample(numel(baseline_surpirse),K);
                dotted_center = median(baseline_surpirse(selection_slice), 'omitmissing');
                mean_value = mean(this_surpise(selection_slice,time_slice) - repmat(baseline_surpirse(selection_slice), 1,sum(time_slice)), 'omitmissing');

                % mean_effect_pos (1,tn,nm,j) = mean(mean_value(mean_value>0));
                % mean_effect_neg (1,tn,nm,j) = mean(mean_value(mean_value<0));
                mean_effect(1,tn,nm,j)      = mean(mean_value);

            end
            mean_base(1,nm,j) = mean(mean(this_surpise(selection_slice,:)));
        end
        end
    end

   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%% peak %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    this_surpise = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_comp,time_centers>border_effect(1)  & time_centers<border_effect(2))));
    for row_n=1:size(this_surpise,1)
        this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
    end

    if K>0
        if K<min_val
            n_counts_used(j,2) =min(n_montecarlo,nchoosek(size(this_surpise,1), K));
             n_counts_used(j,5) =size(this_surpise,1);
        end
        if min(size(this_surpise))>1
            baseline_surpirse = this_surpise(:, sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
        elseif min(size(this_surpise))>0
            baseline_surpirse = this_surpise(sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
        end
        baseline_surpirse = median(baseline_surpirse,2, 'omitmissing');
        time_slice_1 = sub_time_centers>=time_ranges_2_measure(1,1) &  sub_time_centers<time_ranges_2_measure(1,2);
        time_slice_2 = sub_time_centers>=time_ranges_2_measure(2,1) &  sub_time_centers<time_ranges_2_measure(2,2);
        if min(size(this_surpise))>1
            table4mixed_effect = [table4mixed_effect;[mean(this_surpise(:,time_slice_2),2)-mean(this_surpise(:,time_slice_1),2) ones(size(this_surpise,1),2)*[2  0;0 j]]];
        else
             table4mixed_effect = [table4mixed_effect;[mean(this_surpise(time_slice_2))-mean(this_surpise(time_slice_1)) [2 j]]];
        end

        if do_montecarlo
        for nm =1:n_counts_used(j)
            for tn = 1:size(time_ranges_2_measure,1)
                time_slice = sub_time_centers>=time_ranges_2_measure(tn,1) &  sub_time_centers<time_ranges_2_measure(tn,2);
                selection_slice = randsample(numel(baseline_surpirse),K);
                dotted_center = median(baseline_surpirse(selection_slice), 'omitmissing');
                if min(size(this_surpise))>1
                    mean_value = mean(this_surpise(selection_slice,time_slice) - repmat(baseline_surpirse(selection_slice), 1,sum(time_slice)), 'omitmissing');
                else
                    mean_value = mean(this_surpise(time_slice) - repmat(baseline_surpirse(selection_slice), 1,sum(time_slice)), 'omitmissing');
                end
                    % mean_effect_pos (2,tn,nm,j) = mean(mean_value(mean_value>0));
                % mean_effect_neg (2,tn,nm,j) = mean(mean_value(mean_value<0));
                mean_effect(2,tn,nm,j)      = mean(mean_value);
            end
            mean_base(2,nm,j) = mean(mean(this_surpise(selection_slice,:)));
        end
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%% non entrained%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    this_surpise = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_nonentr,time_centers>border_effect(1)  & time_centers<border_effect(2))));
    for row_n=1:size(this_surpise,1)
        this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
    end

    if K>0 && size(this_surpise,1)>K
        if K<min_val
            n_counts_used(j,2) =min(n_montecarlo,nchoosek(size(this_surpise,1), K));
             n_counts_used(j,5) =size(this_surpise,1);
        end
        if min(size(this_surpise))>1
            baseline_surpirse = this_surpise(:, sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
        elseif min(size(this_surpise))>0
            baseline_surpirse = this_surpise(sub_time_centers<=baseline_range(2) & sub_time_centers>=baseline_range(1));
        end
        baseline_surpirse = median(baseline_surpirse,2, 'omitmissing');
        time_slice_1 = sub_time_centers>=time_ranges_2_measure(1,1) &  sub_time_centers<time_ranges_2_measure(1,2);
        time_slice_2 = sub_time_centers>=time_ranges_2_measure(2,1) &  sub_time_centers<time_ranges_2_measure(2,2);
        if min(size(this_surpise))>1
            table4mixed_effect = [table4mixed_effect;[mean(this_surpise(:,time_slice_2),2)-mean(this_surpise(:,time_slice_1),2) ones(size(this_surpise,1),2)*[3  0;0 j]]];
        else
             table4mixed_effect = [table4mixed_effect;[mean(this_surpise(time_slice_2))-mean(this_surpise(time_slice_1)) [3 j]]];
        end


    end



end



%% VIII Mixed / Linear model analysis
% =========================
%  Mixed / Linear model analysis
%  Estimated means + 95% CI + plot
%  =========================
%  1) Load table
%  -------------------------
tbl = array2table(table4mixed_effect, 'VariableNames',{'value','group','condition'});
alpha_level = 0.005;
% REQUIRED column names:
% tbl.value
% tbl.group
% tbl.condition
%
% OPTIONAL:
% tbl.subject

tile_side = 'right';
% -------------------------
%  2) Convert grouping vars
%  -------------------------
tbl.group = categorical(tbl.group);
tbl.condition = categorical(tbl.condition);

hasSubject = ismember('subject', tbl.Properties.VariableNames);

if hasSubject
    tbl.subject = categorical(tbl.subject);
end

% -------------------------
%  3) Fit model
%  -------------------------
if hasSubject
    % Mixed model if repeated observations exist
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
    % Standard linear model if all rows are independent
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

% -------------------------
%  4) Build prediction table
%  -------------------------
groupCats = categories(tbl.group);
condCats  = categories(tbl.condition);

[groupGrid, condGrid] = ndgrid(groupCats, condCats);

predTbl = table;
predTbl.group = categorical(groupGrid(:), groupCats);
predTbl.condition = categorical(condGrid(:), condCats);

% If mixed model, provide one valid subject category for prediction
if hasSubject
    subjCats = categories(tbl.subject);
    predTbl.subject = categorical(repmat(subjCats(1), height(predTbl), 1), subjCats);
end

% -------------------------
%  5) Get estimated means + 95% CI
%  -------------------------
if hasSubject
    [yhat, yCI, ~] = predict(lme, predTbl);
else
    [yhat, yCI] = predict(lm, predTbl, 'Alpha', alpha_level);
end

predTbl.estimate = yhat;
predTbl.CI_low   = yCI(:,1);
predTbl.CI_high  = yCI(:,2);
predTbl.aboveZero = predTbl.CI_low > 0 | predTbl.CI_high < 0;



% -------------------------
%  6) Add numeric versions for easy sorting/plotting
%  -------------------------
predTbl.group_num = str2double(string(predTbl.group));
predTbl.cond_num  = str2double(string(predTbl.condition));

predTbl = sortrows(predTbl, {'group_num','cond_num'});

% -------------------------
%  7) Display and save result table
%  -------------------------
fprintf('\n=============================\n');
fprintf('ESTIMATED MEANS + 95%% CI\n');
fprintf('=============================\n');

disp(predTbl(:, {'group','condition','estimate','CI_low','CI_high','aboveZero'}))

% Optional: save results
estimated_means_table = predTbl;
% writetable(predTbl, 'estimated_condition_means.csv');


%%  IX 8) Plot estimated means with 95% CI
%     (assign NaN to missing group-condition combinations)
%  -------------------------
% -------------------------
% Detect which group-condition combinations actually exist in raw data
% -------------------------
color_code = {'b','r','g'};
predTbl.Nobs = zeros(height(predTbl),1);

for r = 1:height(predTbl)
    idx = tbl.group == predTbl.group(r) & tbl.condition == predTbl.condition(r);
    predTbl.Nobs(r) = sum(idx);
end

% Any group-condition combo with no observations should not be plotted/interpreted
missingIdx = predTbl.Nobs == 0;

predTbl.estimate(missingIdx)  = NaN;
predTbl.CI_low(missingIdx)    = NaN;
predTbl.CI_high(missingIdx)   = NaN;
predTbl.aboveZero(missingIdx) = false;

% Plot


x_tick_lables = cell(size(unique_area_combinations,1),1);
for j=1:size(unique_area_combinations,1)
     x_tick_lables{j} = [unique_area_combinations.Neuron1{j}, ' ',unique_area_combinations.Neuron2{j}];
end
figure('Color','w','Position',[100 100 1000 500]); 
hold on

groups = categories(predTbl.group);
offset = linspace(-0.12, 0.12, numel(groups));   % works even if >2 groups

for g = 1:numel(groups)

    idx = predTbl.group == groups{g};

    x  = predTbl.cond_num(idx);
    y  = predTbl.estimate(idx);
    lo = predTbl.CI_low(idx);
    hi = predTbl.CI_high(idx);

    % sort by condition
    [x, sortIdx] = sort(x);
    y  = y(sortIdx);
    lo = lo(sortIdx);
    hi = hi(sortIdx);

    % Only keep valid (non-missing) entries
    validIdx = ~isnan(y) & ~isnan(lo) & ~isnan(hi);

    % Plot error bars only for valid points
bar(x(validIdx) + offset(g),y(validIdx), 'FaceColor',color_code{g}, 'FaceAlpha',.2, 'EdgeColor','none')
    errorbar(x(validIdx) + offset(g), ...
             y(validIdx), ...
             y(validIdx)-lo(validIdx), ...
             hi(validIdx)-y(validIdx), ...
             'o-', ...
             'LineWidth', 1.8, ...
             'MarkerSize', 7, ...
             'CapSize', 8,...
             'Color',color_code{g});
    
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

yline(0, 'k--', 'LineWidth', 1.2);

xlabel('Condition', 'FontSize', 12)
ylabel('Estimated mean value', 'FontSize', 12)
title('Model-estimated means with 95% confidence intervals', 'FontSize', 14)


legend(arrayfun(@(x) sprintf('Group %s', x), string(groups), 'UniformOutput', false), ...
    'Location', 'best')

box off
set(gca, 'FontSize', 11, 'LineWidth', 1.2)

%% X Plot (with empty bar plotS, as on Supp Figure)
figure('Color','w','Position',[100 100 1000 500]); 
hold on

x_tick_lables = cell(size(unique_area_combinations,1),1);
for j=1:size(unique_area_combinations,1)
     x_tick_lables{j} = [unique_area_combinations.Neuron1{j}, ' ',unique_area_combinations.Neuron2{j}];
end
groups = categories(predTbl.group);
offset = linspace(-0.18, 0.18, numel(groups));   % a bit wider for bars

% Example colors (edit as you want)
color_code = {
    [0.2 0.4 0.8]
    [0.8 0.2 0.3]
    [0.2 0.8 0.3]
};

barWidth = 0.28;

for g = 1:numel(groups)

    idx = predTbl.group == groups{g};

    x  = predTbl.cond_num(idx);
    y  = predTbl.estimate(idx);
    lo = predTbl.CI_low(idx);
    hi = predTbl.CI_high(idx);
    sig = predTbl.aboveZero(idx);

    % sort by condition
    [x, sortIdx] = sort(x);
    y   = y(sortIdx);
    lo  = lo(sortIdx);
    hi  = hi(sortIdx);
    sig = sig(sortIdx);

    % only keep valid (non-missing) entries
    validIdx = ~isnan(y) & ~isnan(lo) & ~isnan(hi);

    xPlot  = x(validIdx) + offset(g);
    yPlot  = y(validIdx);
    loPlot = lo(validIdx);
    hiPlot = hi(validIdx);
    sigPlot = sig(validIdx);

    % ---- draw bars one by one so each can have different fill ----
    for mm = 1:numel(xPlot)

        if sigPlot(mm)
            % Significant above zero -> filled transparent bar
            bar(xPlot(mm), yPlot(mm), barWidth, ...
                'FaceColor', color_code{g}, ...
                'FaceAlpha', 0.2, ...
                'EdgeColor', color_code{g}, ...
                'LineWidth', 1.2);
        else
            % Not significant above zero -> empty bar
            bar(xPlot(mm), yPlot(mm), barWidth, ...
                'FaceColor', 'none', ...
                'EdgeColor', color_code{g}, ...
                'LineWidth', 1.2);
        end
    end

    % ---- error bars ----
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
%% display mean values -------------------------
%  9) Optional: separate matrices for easy downstream use
%  -------------------------
nGroups = numel(groupCats);
nConds  = numel(condCats);

meanMat = nan(nGroups, nConds);
ciLowMat = nan(nGroups, nConds);
ciHighMat = nan(nGroups, nConds);
aboveZeroMat = false(nGroups, nConds);

for g = 1:nGroups
    for c = 1:nConds
        idx = predTbl.group == groupCats{g} & predTbl.condition == condCats{c};
        meanMat(g,c) = predTbl.estimate(idx);
        ciLowMat(g,c) = predTbl.CI_low(idx);
        ciHighMat(g,c) = predTbl.CI_high(idx);
        aboveZeroMat(g,c) = predTbl.aboveZero(idx);
    end
end

fprintf('\n=============================\n');
fprintf('MATRIX OUTPUTS\n');
fprintf('=============================\n');
fprintf('meanMat      = estimated means\n');
fprintf('ciLowMat     = lower 95%% CI\n');
fprintf('ciHighMat    = upper 95%% CI\n');
fprintf('aboveZeroMat = 1 if lower CI > 0\n\n');

disp('meanMat = ')
disp(meanMat)

disp('aboveZeroMat = ')
disp(aboveZeroMat)

datafor_thickness = meanMat;
datafor_thickness(~aboveZeroMat) = NaN;


%% plot only the bars
figure
barWidth = .8
for j=1:3
 subplot(3,1,j)
 hold on
 for mm=1:numel(condCats)
     k = str2double(condCats{mm});
     if aboveZeroMat(j,mm)
         % Significant above zero -> filled transparent bar
         bar(k, meanMat(j,mm), barWidth, ...
             'FaceColor', color_code{j}, ...
             'FaceAlpha', 0.2, ...
             'EdgeColor', color_code{g}, ...
             'LineWidth', 1.2);

     end
 end
    ylim([-.5 .5])
    yticks(-.5:.05:.5)
    
xlim([0 16])
xticks(1:15)
xticklabels(x_tick_lables)
xtickangle( gca , 90 )
set(gca, 'TickDir', 'out')
   
    % bar(meanMat(j,:))
end

%%
print(gcf,'-vector','-dsvg',[figure_folder,'/ bar plot with significnat increase or decrease  synchrony NON PLAY.svg'])



%% ploting imagesc
y_lim = [.3 .6];
mov_wind = 10;
c_lim = [.2 .6];
smooth_wind_sec = .1;
nan_pctls = any(isnan(pctl_spikes_histogram),3);
 nan_pctls = squeeze(any(nan_pctls));
log_scale = false;
modulated_pctg = Inf;
min_fr = 5;
trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
alpha_level     = 0.05/15;

possible_titles = {'Trough','Peak','Non-Entrained'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% below select 1 2 or 3 depending what to plot %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sg_title = possible_titles{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(sg_title, 'Trough')
    entrained_group1 = trough;
    entrained_group2 = trough;
elseif strcmp(sg_title, 'Peak')
    entrained_group1 = peak;
    entrained_group2 = peak;
else
    entrained_group1 = non_entrained;
    entrained_group2 = non_entrained;
end


plot_color = 'b';

z_scored_limit = Inf;
figure
smoth_wind = 10;
smoth_wind2 = 10;

parametric_cond = false;
x_lim = border_effect;

plot_ci = true;
plot_logit = ~plot_ci;

subplot_order =  reshape(1:30,5,6)';

imagesc_columns = subplot_order(1:2:end,:)';
imagesc_columns = imagesc_columns(:);


plot_columns = subplot_order(2:2:end,:)';
plot_columns = plot_columns(:);

for j=1:size(unique_area_combinations,1)

  

    area_combination_indexes      = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main          = entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes  & mean(cluster_fr,2)>min_fr & ~nan_pctls';

    % modulated = sum(mean())


    this_surpise = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,:)));
     for row_n=1:size(this_surpise,1)
        this_surpise(row_n, :) = movmean(this_surpise(row_n,:),smooth_wind_sec/bin_size);
     end
    surprise                    = this_surpise;
    if min(size(surprise))>1
  
    no_nan =all(~isnan(surprise),2);

    time_synch = sum(surprise,2);
    [~,order] = sort(time_synch);
    
    subplot(6,5,imagesc_columns(j))
    hold on
    imagesc(time_centers, 1:numel(order),surprise(order,:))
    clim(c_lim)
    xlim(x_lim) 
    plot([0 0], [1 numel(order)], 'w')  
     title([unique_area_combinations.Neuron1{j}, ' #', num2str(sum(cell_type_comb_main))])
     ylabel(unique_area_combinations.Neuron2{j})  

    xticks([x_lim(1) 0 x_lim(2)/2 x_lim(2)])
    set(gca, 'TickDir', 'out')

    subplot(6,5,plot_columns(j))
    plot(time_centers, median(surprise))
    xlim(border_effect)
    ylim(y_lim)

    end
    


end
sgtitle(sg_title)

%% save if needed
print(gcf,'-vector','-dsvg',[figure_folder,'/imagesc conincidence per area ',sg_title,'.svg'])
%% poting percetnag eof time synchronized
xtick_abel = cell(size(unique_area_combinations,1),1);
for j=1:size(unique_area_combinations,1)

     xtick_abel{j} = [unique_area_combinations.Neuron1{j}, ' to ', unique_area_combinations.Neuron2{j}];

end

figure
hold on
alpha_level = 0.01;
time2comapre = [-.2 .5];
index = time_centers(time2consider)>time2comapre(1) & time_centers(time2consider)<time2comapre(2);
A = squeeze(sum(PVAL_PER_CONDITION(1,:,index)<alpha_level,3))*bin_size/range(time2comapre);
B = squeeze(sum(PVAL_PER_CONDITION(2,:,index)<alpha_level,3))*bin_size/range(time2comapre);
PA = squeeze(min(PVAL_PER_CONDITION(1,:,index),[],3));
PB = squeeze(min(PVAL_PER_CONDITION(2,:,index),[],3));
for i = 1:numel(A)
    % Bar for A (upwards)
    if PA(i) < alpha_level
        bar(2*i-1, A(i), 'FaceColor', 'b', 'EdgeColor', 'b'); % blue if significant
    else
        bar(2*i -1, A(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end

    % Bar for B (downwards, using negative value)
    if PB(i) < alpha_level
        bar(2*i, B(i), 'FaceColor', 'r', 'EdgeColor', 'r'); % red if significant
    else
        bar(2*i, B(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end
end

xticks(1.5:2:(2*numel(A)))

xticklabels(xtick_abel)










%% save if needed
print(gcf,'-vector','-dsvg',[figure_folder,'/percetnag eof time synchronized Trough and non entrained, controled by count.svg'])


%% checking correlations for all observations


time2comapre = [-1 0];
index = time_centers(time2consider)>time2comapre(1) & time_centers(time2consider)<time2comapre(2);
figure
corrleation_r = nan(2,size(unique_area_combinations,1));
corrleation_p = nan(2,size(unique_area_combinations,1));
x_tick_lables = cell(size(unique_area_combinations,1),1);
for j=1:size(unique_area_combinations,1)

     subplot(5,3,j)
     plot(squeeze(AUC_PER_CONDITION(1,j,index)),squeeze(AUC_PER_CONDITION(2,j,index)), '.')
     [c,p] = corr(squeeze(AUC_PER_CONDITION(1,j,index)),squeeze(AUC_PER_CONDITION(2,j,index)), 'Type','Spearman');
     corrleation_r(1,j) = c;
     corrleation_p(1,j)  = p;

     legend([num2str(c), ' ', num2str(p)])

     title([unique_area_combinations.Neuron1{j}, ' to ', unique_area_combinations.Neuron2{j}])
     x_tick_lables{j} = [unique_area_combinations.Neuron1{j}, ' to ', unique_area_combinations.Neuron2{j}];

end



time2comapre = [0  .5];
index = time_centers(time2consider)>time2comapre(1) & time_centers(time2consider)<time2comapre(2);
figure
for j=1:size(unique_area_combinations,1)

     subplot(5,3,j)
     plot(squeeze(AUC_PER_CONDITION(1,j,index)),squeeze(AUC_PER_CONDITION(2,j,index)), '.')
     [c,p] = corr(squeeze(AUC_PER_CONDITION(1,j,index)),squeeze(AUC_PER_CONDITION(2,j,index)), 'Type','Spearman');
     legend([num2str(c), ' ', num2str(p)])
    corrleation_r(2,j) = c;
     corrleation_p(2,j)  = p;
     title([unique_area_combinations.Neuron1{j}, ' to ', unique_area_combinations.Neuron2{j}])

end

%% plot correlation fro surpirse values (all observation)
figure
hold on
alpha_level = 0.01;
A = corrleation_r(1,:);
B = corrleation_r(2,:)
PA = corrleation_p(1,:)
PB = corrleation_p(2,:);
for i = 1:numel(A)
    % Bar for A (upwards)
    if PA(i) < alpha_level
        bar(2*i-1, A(i), 'FaceColor', 'k', 'EdgeColor', 'k'); % blue if significant
    else
        bar(2*i -1, A(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end

    % Bar for B (downwards, using negative value)
    if PB(i) < alpha_level
        bar(2*i, B(i), 'FaceColor', 'm', 'EdgeColor', 'm'); % red if significant
    else
        bar(2*i, B(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end
end
xticks(1.5:2:(2*numel(A)))

xticklabels(x_tick_lables)
%% XI estimate correaltion surroage distribution

fix_size        =10 ;
n_montecarlo    = 5000;
time2comapre_cell = {[-1 0], [0 1]};




 border_effect = [-1.4 .9];
 time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
time2consider = time_centers>border_effect(1)  & time_centers<border_effect(2);
index = cell(numel(time2comapre_cell),1);
for j=1:numel(time2comapre_cell)
index{j} = time_centers(time2consider)>time2comapre_cell{j}(1) & time_centers(time2consider)<time2comapre_cell{j}(2);
end


smooth_wind_sec = 0.05;
min_fr              = 5;

trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
figure

entrained_group1 = trough;
entrained_group2 = trough;
comparison_group1 = peak;
comparison_group2 = peak;


sub_time_centers = time_centers(time2consider);

plot_color = {'b', 'r'};



correlations_c_montecarlo  = nan(numel(time2comapre_cell),n_montecarlo,size(unique_area_combinations,1));
correlations_p_montecarlo  = nan(numel(time2comapre_cell),n_montecarlo,size(unique_area_combinations,1));
n_counts_used = nan(size(unique_area_combinations,1),5);
figure

for j=1:size(unique_area_combinations,1)

    

    area_combination_indexes      = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main          = find(entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes & mean(cluster_fr,2)>min_fr);

    cell_type_comb_comp          = find(comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) & area_combination_indexes  & mean(cluster_fr,2)>min_fr);
    
    
    min_val = min(numel(cell_type_comb_main),numel(cell_type_comb_comp));


    if  min_val>fix_size && nchoosek(min_val, fix_size)>n_montecarlo
        K = fix_size;
        n_counts_used(j,4) =min(n_montecarlo,K);
    elseif min_val>0 && min_val<fix_size
        K = round(.75*min_val);
        n_counts_used(j,4) =min(n_montecarlo,K);
    else
        K=0;
    end
    disp(K)
    n_counts_used(j,5) = K;
  
    if K>0      

        this_surpise_main =1- squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,time_centers>border_effect(1)  & time_centers<border_effect(2))));
        
        for row_n=1:size(this_surpise_main,1)
            this_surpise_main(row_n, :) = movmean(this_surpise_main(row_n,:),smooth_wind_sec/bin_size);
        end
        this_surpise_comp = 1-squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_comp,:)));
        for row_n=1:size(this_surpise_comp,1)
            this_surpise_comp(row_n, :) = movmean(this_surpise_comp(row_n,:),smooth_wind_sec/bin_size);
        end

        if K<=min_val
            n_counts_used(j,1) =min(min(n_montecarlo,nchoosek(size(this_surpise_main,1), K)),nchoosek(size(this_surpise_comp,1), K));
            n_counts_used(j,2) =size(this_surpise_main,1);
            n_counts_used(j,3) =size(this_surpise_comp,1);

        end

        for mont_n =1:n_counts_used(j,1)

            for ttc = 1:numel(time2comapre_cell)
                median_coinc_over_time_main = median(this_surpise_main(randsample(numel(cell_type_comb_main), K),index{ttc}), 'omitmissing');
                median_coinc_over_time_comp = median(this_surpise_comp(randsample(numel(cell_type_comb_comp), K),index{ttc}), 'omitmissing');
                [c,p] = corr(median_coinc_over_time_main',median_coinc_over_time_comp', 'Type','Spearman');
                correlations_c_montecarlo(ttc,mont_n,j) = c;
                correlations_p_montecarlo(ttc,mont_n,j) = p;

            end
        end
    end

    

    subplot(5,3,j)
    histogram(squeeze(correlations_c_montecarlo(1,:,j)), -1:0.01:1, 'EdgeColor','none', 'Facecolor', 'k', 'FaceAlpha',  .2)
    hold on
    histogram(squeeze(correlations_c_montecarlo(2,:,j)), -1:0.01:1, 'EdgeColor','none', 'Facecolor', 'm', 'FaceAlpha',  .2)
    pause(.1)
end


%% plot and estimate cohen_d and p values

figure

Cohen_d     = nan(size(unique_area_combinations,1),2);
Mean_corr   = nan(size(unique_area_combinations,1),2); 
Cohen_p     = nan(size(unique_area_combinations,1),2);

for j=1:size(unique_area_combinations,1)
    subplot(5,3,j)
    histogram(squeeze(correlations_c_montecarlo(1,:,j)), -1:0.01:1, 'EdgeColor','none', 'Facecolor', 'k', 'FaceAlpha',  .2)
    Cohen_d(j,1) = mean(squeeze(correlations_c_montecarlo(1,:,j)), 'omitmissing')/std(squeeze(correlations_c_montecarlo(1,:,j)), 'omitmissing');
    Mean_corr(j,1) = mean(squeeze(correlations_c_montecarlo(1,:,j)), 'omitmissing');
    [~,Cohen_p(j,1)] = ttest(squeeze(correlations_c_montecarlo(1,:,j)));
    hold on
    histogram(squeeze(correlations_c_montecarlo(2,:,j)), -1:0.01:1, 'EdgeColor','none', 'Facecolor', 'm', 'FaceAlpha',  .2)
    Cohen_d(j,2) = mean(squeeze(correlations_c_montecarlo(2,:,j)), 'omitmissing')/std(squeeze(correlations_c_montecarlo(1,:,j)), 'omitmissing');
    Mean_corr(j,2) = mean(squeeze(correlations_c_montecarlo(2,:,j)), 'omitmissing');
    [~,Cohen_p(j,2)] = ttest(squeeze(correlations_c_montecarlo(2,:,j)));
    title([unique_area_combinations.Neuron1{j}, ' to ', unique_area_combinations.Neuron2{j}])
    pause(.1)
    ylim([0 200])

end


x_tick_lables = cell(size(unique_area_combinations,1),1);
for j=1:size(unique_area_combinations,1)

    
     x_tick_lables{j} = [unique_area_combinations.Neuron1{j}, ' to ', unique_area_combinations.Neuron2{j}];

end

%% save if needed

print(gcf,'-vector','-dsvg',[figure_folder,'/distribution of correlations for non play.svg'])



%% plot distribution of correlations
figure
hold on
alpha_level = 0.01;
A = Cohen_d(:,1);
B = Cohen_d(:,2);
% y_lim = [-2 2];

A = Mean_corr(:,1);
B = Mean_corr(:,2);
y_lim = [-.25 .25];
PA = Cohen_p(:,1);
PB = Cohen_p(:,2);
for i = 1:numel(A)
    % Bar for A (upwards)
    if PA(i) < alpha_level
        bar(2*i-1, A(i), 'FaceColor', 'k', 'EdgeColor', 'k'); % blue if significant
    else
        bar(2*i -1, A(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end

    % Bar for B (downwards, using negative value)
    if PB(i) < alpha_level
        bar(2*i, B(i), 'FaceColor', 'm', 'EdgeColor', 'm'); % red if significant
    else
        bar(2*i, B(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end
    plot([2*i 2*i]+.5, y_lim, ':k')
end
xticks(1.5:2:(2*numel(A)))
ylim(y_lim)

xticklabels(x_tick_lables)


%% save if needed
print(gcf,'-vector','-dsvg',[figure_folder,'/correaltion direction by structure correlation r NONPLAY behaviors.svg'])
