

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


time_precision = 0.003;

bin_size = 0.001;
hist_range = [-1.5 1];

psth_edges = hist_range(1):bin_size:hist_range(2);
areas2analyse = {'DLPAG'	'DR'	'LPAG'	'SupCol'	'VLPAG'};
behaviors4playbout      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
% behaviors4playbout={'Grooming', 'PounceI','Rearing', 'Sniffing','Scratching', 'Bite'};
synch_structure = [];
%%
tic


for fn = n_strctut:numel(animal_list)

    if fn==1
        synch_structure = GENERATE_SURPIRSE_DYNAMICS_STATS_DITTERING([npx_Raw_Data, '\', animal_list(fn).name],bin_size, hist_range,time_precision,areas2analyse,behaviors4playbout);

        n_strctut = n_strctut+numel(synch_structure);
        animal_names = [animal_names;[repmat(animal_list(fn).name,numel(synch_structure),1) num2cell(1:numel(synch_structure))']];
    else
        transt_psth  = GENERATE_SURPIRSE_DYNAMICS_STATS_DITTERING([npx_Raw_Data, '\', animal_list(fn).name],bin_size, hist_range,time_precision,areas2analyse,behaviors4playbout);

        for sub_j=1:numel(transt_psth)

            synch_structure(n_strctut) = transt_psth(sub_j);
            n_strctut = n_strctut+1;
        end
        animal_names = [animal_names;[repmat({animal_list(fn).name},numel(transt_psth),1) num2cell(1:numel(transt_psth))' ]];

    end
    toc
    d = datestr(datetime('now'));
    d = strrep(d, ':', '_');
    save([saving_folder,'\surpirse_stats_struct_PLAY_2ms_dittering ',d, '.mat'],'synch_structure', '-v7.3');
    save([saving_folder,'\animal_names_PLAY_2ms_dittering',d, '.mat'],'animal_names');
    disp('saved')

end

%%
save([saving_folder,'\surpirse_stats_struct_PLAY_dittering_long_interval_2ms_ALL.mat'],'synch_structure', '-v7.3');
    save([saving_folder,'\surpirse_stats_animal_names__PLAY_dittering_long_interval_2ms_ALL.mat'],'animal_names');
% 
% disp('saving')
% save([saving_folder,'\surpirse_stats_struct_structure_no_play.mat'],'synch_structure', '-v7.3');
% save([saving_folder,'\surpirse_stats_struct_animal_names_no_play.mat'],'animal_names');

%%
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Cross_correlogram';

load([saving_folder,'\surpirse_stats_struct_structure_play_dittering.mat'],'synch_structure');
load([saving_folder,'\surpirse_stats_struct_animal_names_play_dittering.mat'],'animal_names');
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
y_lim = [0 .25];

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

plot_color = 'b';
z_scored_limit = Inf;
smoth_wind = 10;
smoth_wind2 = 10;

parametric_cond = false;
x_lim = [-.25 .5];
plot_ci = true;
plot_logit = ~plot_ci;


for j=1:size(unique_area_combinations,1)

    subplot(5,5,j)
    hold on

    area_combination_indexes       = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_trough          = entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) & area_combination_indexes ;
    % cell_type_comb_trough          = area_combination_indexes ;



    this_cc = squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_trough,:)));
    surprise = this_cc;
    % surprise = this_cc;
    % surprise(surprise<0.05) = 0;
    % surprise(surprise>0.05) = 1;
    % % surprise(isnan(this_cc)) = NaN;
    % surprise = 1-surprise;
    % modulated = (sum(surprise,2, 'omitmissing')>0) & sum(surprise,2, 'omitmissing')<sum(~isnan(surprise),2);
    %    surprise = this_cc(modulated,:);

    surprise(surprise>0.95) = 0;
    surprise(surprise<0.95) = 1;

      surprise = 1-surprise;
    % for raw = 1:size(surprise,1)
    %
    %     if std(surprise(raw,time_centers<0))>0
    %         surprise(raw,:) = ( surprise(raw,:)-mean(surprise(raw,time_centers<0)))/std(surprise(raw,time_centers<0));
    %     else
    %         surprise(raw,:) = zscore(surprise(raw,:));
    %     end
    % end

    if ~isempty(surprise)
    plot(time_centers, mean(surprise, 'omitmissing'),plot_color)

    xlim(x_lim)
    ylim(y_lim)

    hold on
    index = max(surprise(:, time_centers>x_lim(1) & time_centers<x_lim(2)),[],2)>-Inf;


    hold on
    plot([0 0], y_lim, ':k')
    hold on
    plot([x_lim(1)  x_lim(2)], [0.05 0.05], 'k')
    xticks([x_lim(1) 0 x_lim(2)/2 x_lim(2)])
    set(gca, 'TickDir', 'out')
    ylim(y_lim)
    end
    title([unique_area_combinations.Neuron1(j), ' ', num2str(sum(all(~isnan(surprise),2)))])
    ylabel(unique_area_combinations.Neuron2(j))


end
%%
y_lim = [0 .5];

trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
alpha_level     = 0.05/15;

entrained_group1 = trough;
entrained_group2 = trough;
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

    area_combination_indexes       = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)));
    cell_type_comb_trough          = comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) & area_combination_indexes ;



    this_cc = squeeze(mean(pctl_spikes_histogram(:,cell_type_comb_trough,:)));

    surprise = this_cc;
    surprise(surprise<0.05) = 0;
    surprise(surprise>0.05) = 1;
    % surprise(isnan(this_cc)) = NaN;


    surprise = 1-surprise;
    modulated = (sum(surprise,2, 'omitmissing')>0) & sum(surprise,2, 'omitmissing')<sum(~isnan(surprise),2);
   


    surprise = this_cc(modulated,:);
      surprise = 1-surprise;

    % surprise(surprise<0.05) = 0;
    % surprise(surprise>0.05) = 1;
    % surprise(isnan(this_cc)) = NaN;
    
    no_nan =all(~isnan(surprise),2);

    time_synch = sum(surprise,2);
    [~,order] = sort(time_synch);

    % for raw = 1:size(surprise,1)
    %
    %     if std(surprise(raw,time_centers<0))>0
    %         surprise(raw,:) = ( surprise(raw,:)-mean(surprise(raw,time_centers<0)))/std(surprise(raw,time_centers<0));
    %     else
    %         surprise(raw,:) = zscore(surprise(raw,:));
    %     end
    % end


    imagesc(time_centers, 1:size(surprise,2),surprise(order,:))

    xlim(x_lim)
    % ylim(y_lim)

    hold on
    index = max(this_cc(:, time_centers>x_lim(1) & time_centers<x_lim(2)),[],2)>-Inf;


    hold on
    % plot([0 0], y_lim, ':k')
    xticks([x_lim(1) 0 x_lim(2)/2 x_lim(2)])
    set(gca, 'TickDir', 'out')
    % ylim(y_lim)

    title([unique_area_combinations.Neuron1(j), ' ', num2str(sum(modulated))])
    ylabel(unique_area_combinations.Neuron2(j))


end


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
    cell_type_comb_trough          = (entrained_group1(idx_pairs(:,1)) &  entrained_group2(idx_pairs(:,2)) &  area_combination_indexes );
    cell_type_comb_non_entrained   = ( comparison_group1(idx_pairs(:,1)) &  comparison_group2(idx_pairs(:,2)) &  area_combination_indexes);


    this_cc = squeeze(mean(synch_spikes_histogram(:,cell_type_comb_trough,:)));
    for k = 1:size(this_cc,1)
        this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
    end
    if size(this_cc,2)==1
        this_cc = this_cc';
    end
    this_cc_control = squeeze(mean(shifted_spikes_histogram(:,cell_type_comb_trough,:))); %desynch_spikes_histogram,shifted_spikes_histogram
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