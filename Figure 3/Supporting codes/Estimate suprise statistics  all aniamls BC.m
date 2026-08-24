

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Cross_correlogram';
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Codes\Figure codes\Figure cross correlations';
% mkdir(saving_folder)
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];


animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};
animal2exclude = {''};
animal_list(ismember(animal_file_names,animal2exclude)) = [];
animal_names ={};
n_strctut = 1;


time_precision = 0.005;

bin_size = 0.01;
hist_range = [-.25 .5];

psth_edges = hist_range(1):bin_size:hist_range(2);
areas2analyse = {'DLPAG'	'DR'	'LPAG'	'SupCol'	'VLPAG'};
%%
tic


for fn = 1:numel(animal_list)

    if fn==1
        synch_structure = GENERATE_SURPIRSE_DYNAMICS_STATS_BC([npx_Raw_Data, '\', animal_list(fn).name],bin_size, hist_range,time_precision,areas2analyse);

        n_strctut = n_strctut+numel(synch_structure);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(synch_structure),1) num2cell(1:numel(synch_structure))']];
    else
        transt_psth  = GENERATE_SURPIRSE_DYNAMICS_STATS_BC([npx_Raw_Data, '\', animal_list(fn).name],bin_size, hist_range,time_precision,areas2analyse);

        for sub_j=1:numel(transt_psth)

            synch_structure(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({animal_list(fn).name},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]];

    end
    toc
    save([saving_folder,'\surpirse_stats_struct_structure_play_jitter2_BC.mat'],'synch_structure', '-v7.3');
    save([saving_folder,'\surpirse_stats_struct_animal_names_play_jitter2_BC.mat'],'animal_names');
    disp('saved')

end

%%

% disp('saving')
% save([saving_folder,'\surpirse_stats_struct_structure_no_play.mat'],'synch_structure', '-v7.3');
% save([saving_folder,'\surpirse_stats_struct_animal_names_no_play.mat'],'animal_names');

%%
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Cross_correlogram';

load([saving_folder,'\surpirse_stats_struct_structure_no_play.mat'],'synch_structure');
load([saving_folder,'\surpirse_stats_struct_animal_names_no_play.mat'],'animal_names');
% synch_structure(2) = [];

bin_size = 0.01;
hist_range = [-.25 .5];

psth_edges = hist_range(1):bin_size:hist_range(2);
% 
% load([saving_folder,'\surpirse_struct_structure_no_play.mat'],'synch_structure');
% load([saving_folder,'\surpirse_struct_animal_names_no_play.mat'],'animal_names');
%% correct structure
corrected_synch_structure = synch_structure;

%%
synch_spikes_histogram = [];
pctl_spikes_histogram = [];
shifted_spikes_histogram = [];

all_area_comb = [];
session_id = [];
clusters_id = [];

for fn = 1:numel(corrected_synch_structure)

    synch_spikes_histogram          = cat(2,synch_spikes_histogram,corrected_synch_structure(fn).synch_spikes_histogram);
    pctl_spikes_histogram           = cat(2,pctl_spikes_histogram,corrected_synch_structure(fn).synch_spikes_pctl        );
    all_area_comb                   = [all_area_comb;corrected_synch_structure(fn).these_neurons_areas(corrected_synch_structure(fn).synch_comb)];
    session_id                      = [session_id;repmat(animal_names(fn,1),size(corrected_synch_structure(fn).synch_comb,1),1)];
    clusters_id                     = [clusters_id;corrected_synch_structure(fn).cluster_info.cluster_id(corrected_synch_structure(fn).synch_comb)];

end


%%

area_combinations = cell2table(all_area_comb);
area_combinations.Properties.VariableNames = {'Neuron1','Neuron2'};
unique_area_combinations = unique(area_combinations, 'rows');


relevant_areas = {'SupCol','DLPAG','LPAG','VLPAG','DR'};
unique_area_combinations = unique_area_combinations(ismember(unique_area_combinations.Neuron1,relevant_areas) & ismember(unique_area_combinations.Neuron2,relevant_areas),:);


%%
%% load delta and theta phase data
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';

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

%%
y_lim = [0 1];
y_lim = [-1000 0]
border_effect = [-.24 .49];

plot_AUC = false;
modulated_pctg = Inf;
log_scale = true;

trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
alpha_level     = 0.05/15;
 figure

entrained_group1 = trough;
entrained_group2 = trough;
comparison_group1 = non_entrained;
comparison_group2 = non_entrained;
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
time2consider = time_centers>border_effect(1)  & time_centers<border_effect(2);

plot_color = 'r';
smoth_wind = 10;





for j=1:size(unique_area_combinations,1)

    subplot(5,5,j)
    hold on

    area_combination_indexes      = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main          = entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes ;

    cell_type_comb_comp          = comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) & area_combination_indexes ;



    this_cc = squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,time2consider)));
    if min(size(this_cc))>1
        baseline_surpirse = this_cc(:, time_centers(time2consider)<-.1);
    elseif min(size(this_cc))>0
        baseline_surpirse = this_cc(time_centers(time2consider)<-.1);
    end

    baseline_surpirse = baseline_surpirse(:);
    is_abobe_change = nan(1,size(this_cc,2));
    AUC_over_time =  nan(1,size(this_cc,2));
    if ~isempty(this_cc)
        surprise                    = this_cc;

        for raw_n=1:size(surprise,1)

            surprise(raw_n,:) = movmean(surprise(raw_n,:),smoth_wind);
        end

        for t = 1:size(surprise,2)

            [ is_abobe_change(t), ~, stats] =ranksum(surprise(:, t), baseline_surpirse);
            % [ is_abobe_change(t), ~, stats] =signrank(surprise(:, t), .5);

            n_obs = length(surprise(:,t));
            n_base = length(baseline_surpirse);
            U = stats.ranksum - n_obs*(n_obs+1)/2;
            AUC = U / (n_obs * n_base);
            AUC_over_time(t) = AUC;
        end

    end


    if plot_AUC
        plot(time_centers(time2consider), AUC_over_time)
        hold on
        plot([-.25 .5], [.5 .5], ':k')
    else
        semilogy(time_centers(time2consider), 1-is_abobe_change)
        hold on
        semilogy([-.25 .5], [0.95 0.95], 'K')
    end

     title([unique_area_combinations.Neuron1(j), ' ', num2str(size(surprise,1))])
     ylabel(unique_area_combinations.Neuron2(j))       

end
%% ploting imagesc
y_lim = [0 .5];
mov_wind = 10;
c_lim = [0 10];
log_scale = true;
modulated_pctg = Inf;
trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
alpha_level     = 0.05/15;

entrained_group1 = peak;
entrained_group2 = peak;
comparison_group1 = non_entrained;
comparison_group2 = non_entrained;
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));

plot_color = 'b';

z_scored_limit = Inf;
figure
smoth_wind = 10;
smoth_wind2 = 10;

parametric_cond = false;
x_lim = [-.25 .5];

plot_ci = true;
plot_logit = ~plot_ci;


for j=1:size(unique_area_combinations,1)

    subplot(5,5,j)
    hold on

    area_combination_indexes      = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main          = entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes ;



    this_cc = squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,:)));
   

    surprise                    = this_cc;
    if min(size(surprise))>1
  
    no_nan =all(~isnan(surprise),2);

    time_synch = sum(surprise,2);
    [~,order] = sort(time_synch);
    if log_scale
        surprise = log10(surprise);
        c_lim = log10([0.25 1] );
    end
  

    imagesc(time_centers, 1:numel(order),surprise(order,:))
    clim(c_lim)
    xlim(x_lim)
    % ylim(y_lim)

    hold on
    index = max(this_cc(:, time_centers>=x_lim(1) & time_centers<=x_lim(2)),[],2)>-Inf;


    hold on
    % plot([0 0], y_lim, ':k')
    xticks([x_lim(1) 0 x_lim(2)/2 x_lim(2)])
    set(gca, 'TickDir', 'out')
    % ylim(y_lim)

    end
    title([unique_area_combinations.Neuron1(j), ' ', num2str(sum(modulated))])
    ylabel(unique_area_combinations.Neuron2(j))


end

%% ploting distribnutruion of p values
% time_range2test = [0 .25];
time_range2test = [-.25 0];
mov_wind = 10;
c_lim = [2 10];
log_scale = true;


trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
alpha_level     = 0.05/15;

entrained_group1 = trough;
entrained_group2 = trough;
comparison_group1 = peak;
comparison_group2 = peak;
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));

plot_color = 'b';

z_scored_limit = Inf;
figure
smoth_wind = 10;
smoth_wind2 = 10;
mean_values_and_ACU = nan(size(unique_area_combinations,1),6);


x_lim = [-.25 .5];

plot_ci = true;
plot_logit = ~plot_ci;
combination_names = cell(size(unique_area_combinations,1),1);

for j=1:size(unique_area_combinations,1)

    subplot(5,5,j)
    hold on

    area_combination_indexes       = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main             = entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes ;
    cell_type_comb_comparisson       = comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) & area_combination_indexes ;



    this_cc = squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_main,:)));
    surprise = this_cc;
    surprise = surprise*2000/2001;
    index = time_centers>=time_range2test(1) & time_centers<time_range2test(2);
    prob_distr_main= mean(surprise(:,index),2,'omitmissing');
    if ~isempty(prob_distr_main)
        mean_values_and_ACU(j,1) = median(prob_distr_main);
         [p,~,stats] = signrank(prob_distr_main-.5);
        mean_values_and_ACU(j,5) = p;
    histogram(prob_distr_main, 0:0.005:1,'FaceColor','b',  'Normalization','percentage', 'FaceAlpha',.5, 'EdgeColor','none')
    end
    
    this_cc = squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_comparisson,:)));
    



    surprise = this_cc;
    surprise = ((surprise*2000)  +1)/2002;
    index = time_centers>=time_range2test(1) & time_centers<time_range2test(2);
    if size(surprise,2)>1
    prob_distr_comp= mean(surprise(:,index),2,'omitmissing');
    else
        prob_distr_comp= mean(surprise(index),'omitmissing');
    end

    

    hold on
    if ~isempty(prob_distr_comp)
    
        histogram(prob_distr_comp, 0:0.005:1,'FaceColor','r',  'Normalization','percentage', 'FaceAlpha',.5, 'EdgeColor','none')
        mean_values_and_ACU(j,2) = mean(prob_distr_comp);

        [p,~,stats] = signrank(prob_distr_comp-.5);
        mean_values_and_ACU(j,6) = p;
        if ~isempty(prob_distr_main)
             


            [p,~,stats] = ranksum(prob_distr_main,prob_distr_comp,'tail','right');

            mean_values_and_ACU(j,4) = p;
            % AUC effect size
            nA = length(prob_distr_main);
            nB = length(prob_distr_comp);

            R1 = stats.ranksum;
            U1 = R1 - nA*(nA+1)/2;
            AUC = U1 / (nA*nB);
             mean_values_and_ACU(j,3) = AUC;
            title([unique_area_combinations.Neuron1{j}, ' ', num2str(numel(prob_distr)), ' ', num2str(p), ' ', num2str(AUC)])
           
        end

    end
     combination_names{j} = [unique_area_combinations.Neuron1{j}, ' ',unique_area_combinations.Neuron2{j}];
    ylabel(unique_area_combinations.Neuron2(j))


end

%%

figure
subplot(1,2,1)
hold on
A = mean_values_and_ACU(:,2);
P = mean_values_and_ACU(:,6);
for i = 1:numel(A)
    % Bar for A (upwards)
    if P(i) < 0.01
        bar(i, A(i), 'FaceColor', 'b', 'EdgeColor', 'b'); % blue if significant
    else
        bar(i , A(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end

end

subplot(1,2,2)
hold on
A = mean_values_and_ACU(:,1);
P = mean_values_and_ACU(:,5);
for i = 1:numel(A)
    % Bar for A (upwards)
    if P(i) < 0.01
        bar(i, A(i), 'FaceColor', 'b', 'EdgeColor', 'b'); % blue if significant
    else
        bar(i , A(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end

end

%%
figure
hold on
A = mean_values_and_ACU(:,1);
B = mean_values_and_ACU(:,2);
P = mean_values_and_ACU(:,4);
for i = 1:numel(A)
    % Bar for A (upwards)
    if P(i) < 0.01
        bar(2*i-1, A(i), 'FaceColor', 'b', 'EdgeColor', 'b'); % blue if significant
    else
        bar(2*i -1, A(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end

    % Bar for B (downwards, using negative value)
    if P(i) < 0.01
        bar(2*i, B(i), 'FaceColor', 'r', 'EdgeColor', 'r'); % red if significant
    else
        bar(2*i, B(i), 'FaceColor', 'none', 'EdgeColor', 'k'); % black outline if not
    end
end

xticks(1.5:2:(2*numel(A)))

xticklabels(combination_names)


%%



trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
alpha_level     = 0.05/15;

entrained_group1 = trough;
entrained_group2 = non_entrained; 
comparison_group1 = non_entrained;
comparison_group2 = non_entrained;
time_centers = .5*(psth_edges(1:end-1)+psth_edges(2:end));
plot_color = 'r';
z_scored_limit = Inf;
figure
smoth_wind = 10;
smoth_wind2 = 10;
x_lim = [-.25 .5];
plot_ci = true;
plot_logit = ~plot_ci;

effect_size = cell(size(unique_area_combinations,1),3);
p_Val_per_cell = nan(size(unique_area_combinations,1),3);
p_Val_length_per_cell = nan(size(unique_area_combinations,1),3);

max_gap = 1;
for j=1:size(unique_area_combinations,1)

    for cell_type = 1:3

        if cell_type ==1
            entrained_group1 = trough;
            entrained_group2 = trough;
        elseif cell_type ==2
            entrained_group1 = peak;
            entrained_group2 = peak;
        else
            entrained_group1 = non_entrained;
            entrained_group2 = non_entrained;
        end

       

    area_combination_indexes       = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_main          = (entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) &  area_combination_indexes );
    cell_type_comb_non_entrained   = ( comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) &  area_combination_indexes);


    this_cc = squeeze(mean(synch_spikes_histogram(:,cell_type_comb_main,:)));
    for k = 1:size(this_cc,1)
        this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
    end
    if size(this_cc,2)==1
        this_cc = this_cc';
    end
    this_cc_control = squeeze(mean(shifted_spikes_histogram(:,cell_type_comb_main,:))); %desynch_spikes_histogram,shifted_spikes_histogram
    for k = 1:size(this_cc_control,1)
        this_cc_control(k,:) = movmean( this_cc_control(k,:),smoth_wind);
    end
    if size(this_cc_control,2)==1
        this_cc_control = this_cc_control';
    end
    surprise = this_cc- this_cc_control;

    for raw = 1:size(surprise,1)
        if std(surprise(raw,time_centers<0))>0
            surprise(raw,:) = ( surprise(raw,:)-mean(surprise(raw,time_centers<0)))/std(surprise(raw,time_centers<0));
        else
            surprise(raw,:) = zscore(surprise(raw,:));
        end
    end


    plot(time_centers, mean(surprise, 'omitmissing'),plot_color)


    index = max(this_cc(:, time_centers>x_lim(1) & time_centers<x_lim(2)),[],2)>-Inf;
    [h, p, ci] = ttest(surprise(index,:), 0,'alpha',alpha_level);
            
        h = h & time_centers>=x_lim(1) & time_centers<=x_lim(2);


        h = h(:)';           % make row vector
        h_out = h;

        d = diff([0 h 0]);   % detect transitions
        starts1 = find(d == 1);
        ends1   = find(d == -1) - 1;

        % connect consecutive 1-groups if gap <= max_gap
        for k = 1:length(ends1)-1
            gap_start = ends1(k) + 1;
            gap_end   = starts1(k+1) - 1;
            gap_len   = gap_end - gap_start + 1;

            if gap_len <= max_gap
                h_out(gap_start:gap_end) = 1;
            end
        end

        h = h_out;

    mu = mean(surprise(index,:), 1, 'omitnan');          % 1 x T
    sd = std(surprise(index,:), 0, 1, 'omitnan');        % 1 x T
    n  = sum(~isnan(surprise(index,:)), 1);              % 1 x T

    d = mu ./ sd;


     effect_size{j,cell_type} = d(time_centers>=x_lim(1) & time_centers<=x_lim(2));

        [L,n]= bwlabeln(h);
    
        if n>1
            pval_stats = nan(n,2);

            for nn=1:n
                pval_stats(nn,1) = median(p(L==nn));
                pval_stats(nn,2) =sum(L==nn);
            end
            % pval_stats = pval_stats(pval_stats(:,2)>2,:);

            if isempty(pval_stats)

                p_Val_per_cell(j,cell_type)        = 1;
                p_Val_length_per_cell(j,cell_type) = 0;
            else
                p_Val_per_cell(j,cell_type)        = max(pval_stats(:,1));
                p_Val_length_per_cell(j,cell_type) = max(pval_stats(:,2));
            end
        else

            p_Val_per_cell(j,cell_type)        = 1;
            p_Val_length_per_cell(j,cell_type) = 0;
        end



    end
end

%%


n = 5;

[I,J] = find(triu(true(n)));   % keep i <= j
idx = (I-1)*n + J;


sub_time = time_centers(time_centers>=x_lim(1) & time_centers<=x_lim(2));;
figure
colorspertype = {'b','r','g'}
legnedfor_xticks = {};
for j=1:numel(idx)
    for cell_type = 1:3
        subplot(3,numel(idx),j + (cell_type-1)*numel(idx))
        plot(sub_time,effect_size{idx(j),cell_type} , colorspertype{cell_type})
        title([unique_area_combinations{idx(j),1}{1}, ' to ', unique_area_combinations{idx(j),2}{1}])
        legnedfor_xticks{j} = [unique_area_combinations{idx(j),1}{1}, ' to ', unique_area_combinations{idx(j),2}{1}];
        ylim([-.6 .6])
    end
end

%%

figure

subplot(3,1,1)
bar(1- (p_Val_per_cell(idx,1)*15))
xticks(1:15)
xticklabels(legnedfor_xticks)
ylim([0.90 1])

subplot(3,1,2)
bar(1- (p_Val_per_cell(idx,2)*15))
ylim([0.90 1])


subplot(3,1,3)
bar(1- (p_Val_per_cell(idx,3)*15))'
xticks(1:15)
xticklabels(legnedfor_xticks)
ylim([0.90 1])




figure
bin_length = mean(diff(time_centers));
subplot(3,1,1)
bar(p_Val_length_per_cell(idx,1)*bin_length)
xticks(1:15)
xticklabels(legnedfor_xticks)
% ylim([0.90 1])

subplot(3,1,2)
bar(p_Val_length_per_cell(idx,2)*bin_length)
% ylim([0.90 1])


subplot(3,1,3)
bar(p_Val_length_per_cell(idx,3)*bin_length)
xticks(1:15)
xticklabels(legnedfor_xticks)
% ylim([0.90 1])

%%
cross_corr_figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Codes\Figure codes\Figure cross correlations';

print(gcf,'-vector','-dsvg',[cross_corr_figure_folder,'/Trough cells exeess coactvation.svg'])