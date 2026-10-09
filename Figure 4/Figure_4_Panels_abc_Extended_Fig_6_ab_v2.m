%%  Based on NeuronalResponse_AcrossAreas_DRbased_CoModulation_Warped_v4.m
%%  Matias Mugnaini 

%% load table

% addpath(genpath('Y:\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis'))
% savepath
% addpath(genpath('Y:\PlayNeuralData\NPX-OPTO PLAY NMM\MATLAB CODES'))
% savepath
% addpath(genpath('Y:\PlayNeuralData\NPX-OPTO PLAY NMM'))
% savepath

saving_folder = 'Y:\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';

load([saving_folder,'\delta_all_neurons_v2.mat'],'all_neurons');
all_neurons_TD = all_neurons;
all_neurons.area(ismember(all_neurons.area, {'isRT'}))  =     {'isRt'  };
all_neurons_TD.DeltaPartner1                            = all_neurons.Partner1;
all_neurons_TD.DeltaPartner2                            = all_neurons.Partner2;
all_neurons_TD.DeltaEntireSession                       = all_neurons.EntireSession;
all_neurons_TD.DeltaPlay                                = all_neurons.Play;
all_neurons_TD.DeltaPrePlay                             = all_neurons.PrePlay;
all_neurons_TD.Exited                                   = nan(size(all_neurons_TD,1),1);
all_neurons_TD.Inhibited                                = nan(size(all_neurons_TD,1),1);

sum(ismember(all_neurons_TD.session,'B4S2 0825 Single') & all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2))


%%
alpha_level = 0.01;
non_entrained_level = 0.1;
non_entrained = all_neurons_TD.DeltaEntireSession.PPCPval>non_entrained_level;
entrained_cells = all_neurons_TD.DeltaEntireSession.PPCPval<alpha_level;
trough_cells = (all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2 |  all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2) & entrained_cells;
real_peak = ~(all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2 |  all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2) & entrained_cells;

[modulated_play,exited_play,inhibited_play,exited_non_play,inhibited_non_play,modulated_non_play,all_neurons_TD] = load_modulated_play();

session_and_id = all_neurons_TD (:,{'session', 'cluster_id', 'depth','ch'});

%% PARAMETERS

n_cat={'good','mua'};   % Cluster type  
ca=1;
sigma2=4;       % for shade plot smoothness
lim=5;
hist_ed = [-lim lim];
histogram_edges=hist_ed;
bin_size        = 0.2;
    bins=diff(histogram_edges)/bin_size;


post_onset=1; titlepostorpre='Post-Onset';% 1 if you want to evaluate the period posterior to the onset
% post_onset=0; titlepostorpre='Pre-Onset';

bins_to_eval=5; % use 10 for bin_size = 0.1; 5 for bin_size = 0.2

wrap_bins=10;

rat_to_plot=1;
pick_partner=2; % 2, 3 or 4

PlayBout=1; % Choose 1 if you want PlayBouts, 0 for individual behaviors
playbout_tittle='PlayBout';save_responses=1;
beh_to_plot='Pounce_A';
length_duration_threshold=0.25;

thebehaviors = 'PBP';
short_beh_list=0;


redundancies=1;


B1D1 = [1 2 3];
B1S3 = [4 5 6];
B2S2 = [7 8 9];
B3D2 = 10;

%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%CHECK THISSSSSS
B4D4 = [11 11] ; %   0 if medial 1 if lateral %% Dual 4rd Batch Probe 1
B4S2 = [12 12];                              %% Single 4rd Batch Probe 1

thisrat=[ 1 2 3 4 5 6 7 8 9 10 12 12 ];
MorL = {'_','_','_',  '_','_','_', '_','_','_',  '_','Med','Lat'};
probes=[              1 1 1 1 1  1 1 1 1 0  3 1  1 1 1  1 3  ]; %    3 is medial and 1 is lateral     for B4S2, oposite if it is B4D4


IsItAttack=[0 0];
pre_beh_event=1;

categories_play=[1 2 3 4 5 6 7 8 9  ]; % 
categories_non_play=[10 11 12 13 14 15 16  ];
categories = categories_play;

pval_th = 0.05/numel(categories);

behaviors_BL={'CC','CB','CD','Escape','Evasion','Pin','Pounce_B','Pounce_A','Boxing','Rearing','Grooming','Sniffing','Scratch','Pounce_Ai','Pounce_Bi','Bite'}; %% 

behaviors2check = behaviors_BL;

%%

exp_info = readtable('Y:\PlayNeuralData\NPX-OPTO PLAY NMM\Experiments_info.csv');

behavior_files = dir('*.txt');

mixed_combinations_cell = cell(5,4);
your_area={'SupCol','DLPAG','LPAG','VLPAG','DR'};
play_betas={};
            nonplay_betas={};

            mats =[];
            mat_sig_betas=[];
            indexes_areas=[];
            exc_or_inh_indx = [];
            play_sel = [];
            mats_NE = [];
            mats_Pk = [];
            mats_Th = [];
            mats_mod_play = [];
            mats_mod_non_play = [];

            full_play_specificity = [];

for areas_id = 1:numel (your_area)

    area_n = find(ismember({'SupCol','DLPAG','LPAG','VLPAG','DR'},your_area{areas_id}));

            full_zscore = [];
            full_zscore_off = [];
            full_zscore_warp = [];
            full_zscore_warp_sh = [];
            full_FR = [];
            full_BL = [];
            full_COV = [];
            full_WF = [];
            full_WF_long = [];
            full_dyn_res= [];
            full_pie_exc= [];
            full_pie_inh= [];
            full_clusterID = [];

            full_indexes_onset = [];
            full_indexes_offset = [];

            full_sig_betas=[];
                        
            full_mat_exc=[];
            full_mat_inh=[];
            full_mat_all=[];
            full_is_play = [];

            full_depth=[];

            full_non_entrained=[];
            full_peak=[];
            full_through=[];
            full_modulated_play=[];

full_modulated_non_play=[];

depth=[];

bins_this_area = cell(wrap_bins,1) ;

    
    bf_i=1;
    for bf = thisrat 
    
        aux_mydate=behavior_files(bf).name(6:9);
        ani_ID=behavior_files(bf).name(1:4);

        if str2double(behavior_files(bf).name(2))>3
            mydate=['2024' aux_mydate];
        else
            mydate=['2023' aux_mydate];
        end

        path='Y:\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Responses_Matrix\ModelCriterion\';


                load([path 'ResponsesMatrix_PPB_p1andp2_' num2str(length_duration_threshold) 's_' playbout_tittle '_' mydate '_' behavior_files(bf).name(1:4) '_' MorL{bf_i} '.mat'],"initial_cluster","AREAS","RM_areas","all_psth_zscore", ...
                    "pre_onset","post_onset","pre_offset","post_offset","all_psth_zscore_offset","all_features","all_this_psth","this_psth_shuffled","est_full","pval_full","mod_wrap", ...
                    "good_clusters","depth_or_Chn","all_psth_FR","all_psth_FR_offset","AvFR","av_FR_nonbehavior","width_ms","COV","width_ms_long","mean_wf","est","pval","Criterion","dev_test_p")

        
        depth = depth_or_Chn;

                        %% Select your area and its limits
    

                        areas={};    
                        if strcmp(ani_ID,'B1D1')
                            if bf_i > 12 & bf_i < 16
                                ratID={'Batch1Dual1'};probe=probes(bf_i);probe_area='mPFC';
                            else
                                ratID={'Batch1Dual1'};probe=probes(bf_i);probe_area='PAG';
                            end
                        elseif strcmp(ani_ID,'B1S3')
                            ratID={'Batch1Single3'};probe=probes(bf_i);probe_area='PAG';
                        elseif strcmp(ani_ID,'B2S2')
                            ratID={'Batch2Single2'};probe=probes(bf_i);probe_area='PAG';
                        elseif strcmp(ani_ID,'B3D2')
                            ratID={'Batch3Dual2'};probe=probes(bf_i);probe_area='PAG';
                        elseif strcmp(ani_ID,'B4S2')
                            if bf_i == 11
                                 ratID={'Batch4Single2'};probe=probes(bf_i);probe_area='PAG';
                            elseif bf_i == 12
                                ratID={'Batch4Single2'};probe=probes(bf_i);probe_area='PAG';
                            end
                        elseif strcmp(ani_ID,'B4D4')
                            if bf_i == 16
                                 ratID={'Batch4Dual4'};probe=probes(bf_i);probe_area='PAG';
                            elseif bf_i == 17
                                ratID={'Batch4Dual4'};probe=probes(bf_i);probe_area='PAG';
                            end
                        end

            [PAG_columns,available_AREAS,allareas] = take_PAG_area(ratID,probe,probe_area,[]);
        
        sel_area =[];

        sel_area=find(strcmp(available_AREAS,your_area{areas_id})) ;
    
            if sel_area==numel(available_AREAS)                
                init_area=PAG_columns(sel_area);
                end_area=max(depth);
                % end_area=3840;
            else
                init_area=PAG_columns(sel_area);
                end_area=PAG_columns(sel_area+1);
            end
        
            this_area_index = depth >= init_area & depth <= end_area;

        %% PEAK and THROUGH cells
        
        sel_session = strcmp(session_and_id{:,1}, behavior_files(bf).name(1:end-4));

        aux_goodClusters = [] ; aux_Depth = [] ; aux_session = [] ;
        aux_goodClusters = session_and_id{sel_session,2};
        aux_Depth = session_and_id{sel_session,3};
        aux_Ch = session_and_id{sel_session,4};

        path2='Y:\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Responses_Matrix\ModelCriterion_Onset_Miguel\';

        
        load([path2 'all_matched_tables.mat'])
        

        if bf == 12  %||  bf == 11 
            this_rat=all_matched_tables(:,6);

            idx = strfind(table2array(this_rat)', depth_or_Chn);  % works for numeric arrays too

            sel=all_matched_tables(idx:idx+numel(depth_or_Chn)-1,1);
            sel_session = table2array(sel);
            non_entrained_thisrat = non_entrained(sel_session); 
            all_neurons_TD.cluster_id(sel_session)
            trough_cells_thisrat  = trough_cells (sel_session);     
            real_peak_thisrat     = real_peak    (sel_session); 
             modulated_play_thisrat =  modulated_play(sel_session);
             modulated_non_play_thisrat = modulated_non_play(sel_session);

        else

            [new_aux_Depth,ord_by_depth] = sort(aux_Depth);
    
            aux_session = session_and_id{ord_by_depth,1};
            aux_goodClusters = aux_goodClusters(ord_by_depth);
            aux_Depth = aux_Depth(ord_by_depth);
    
            sel_goodClusters = aux_goodClusters==good_clusters';
            double_CHECK =     aux_Depth==depth_or_Chn';
    
            if ~all(double_CHECK) & ~all(sel_goodClusters)
                display('WARNING of mismatch')
            end
            non_entrained_thisrat = non_entrained(sel_session);             non_entrained_thisrat = non_entrained_thisrat(ord_by_depth);
            trough_cells_thisrat  = trough_cells (sel_session);             trough_cells_thisrat  = trough_cells_thisrat (ord_by_depth); 
            real_peak_thisrat     = real_peak    (sel_session);             real_peak_thisrat     = real_peak_thisrat    (ord_by_depth);       
             modulated_play_thisrat =  modulated_play(sel_session);         modulated_play_thisrat = modulated_play_thisrat(ord_by_depth);
             modulated_non_play_thisrat = modulated_non_play(sel_session);  modulated_non_play_thisrat = modulated_non_play_thisrat(ord_by_depth);

        end      
    
            %%

            heatmap_values_exc  = est_full > 0 & pval_full < pval_th;
            heatmap_values_inh  = est_full < 0 & pval_full < pval_th;
            heatmap_values_all  = pval_full < pval_th;

            significant_betas_sel = heatmap_values_exc | heatmap_values_inh;
            significant_betas     = est_full;
            significant_betas(~significant_betas_sel) = NaN;

            indx_exc = any (heatmap_values_exc,2);
            indx_inh = any (heatmap_values_inh,2);

            full_mat_exc=[full_mat_exc; heatmap_values_exc(this_area_index,:)];
            full_mat_inh=[full_mat_inh; heatmap_values_inh(this_area_index,:)];
            full_mat_all=[full_mat_all; heatmap_values_all(this_area_index,:)];

            full_depth  =[full_depth; depth(this_area_index)'];

            full_sig_betas=[full_sig_betas; significant_betas(this_area_index,:)];

            full_non_entrained=       [full_non_entrained ; non_entrained_thisrat(this_area_index)];
            full_peak=                [full_peak ; real_peak_thisrat(this_area_index)];
            full_through=             [full_through ; trough_cells_thisrat(this_area_index)];
            full_modulated_play     = [full_modulated_play ; modulated_play_thisrat(this_area_index)];
            full_modulated_non_play = [full_modulated_non_play ; modulated_non_play_thisrat(this_area_index)];


            %% Deviance
                                play_models    = [ 1 2 3 4 ];
                                nonplay_models = [ 1 5 6 7 ];
                                
                                deviance=[];
                                model_aic_mat = false(numel(Criterion),numel(Criterion{1}));
                                for i = 1:numel(Criterion)
                                    for j= 1:numel(Criterion{1})
                                
                                        if j ==1
                                            deviance(i,j)=Criterion{i}{j};
                                        else
                                            deviance(i,j)=Criterion{i}{j}.deviance;
                                        end
                                    end
                                    deviance_array = deviance(i,:);
                                    number_of_param = [0 1 2 sum(est{i}{4}~=0)-1];
                                    [p_val_matrix, best_model]= assign_lower_deviance(deviance_array(play_models), number_of_param, 0.05);
                                    model_aic_mat(i,best_model)= true;
                                    
                                    number_of_param2 = [0 1 2 sum(est{i}{7}~=0)-1];
                                    [p_val_matrix, best_model]= assign_lower_deviance_nonplay(deviance_array(nonplay_models), number_of_param2, 0.05);
                                    if best_model ~= 1
                                        model_aic_mat(i,best_model+3)= true;
                                    end
                                
                                end

                                is_play = any(model_aic_mat(:,play_models(2:end)),2);

                                full_is_play = [full_is_play;is_play(this_area_index)];
                               

            %% %%%%%%%%%%%%%%%%% Warped play specificity %%%%%%%%%%%%%%%%% 
            %% Filtering both the non-modulated and the at least play responsive
            
            aux1=heatmap_values_exc(this_area_index,:);
            aux2=heatmap_values_inh(this_area_index,:);
            pre_sel_exc=find(aux1==1);
            pre_sel_inh=find(aux2==1);

            pre_matrix = zeros(size(aux1,1),size(aux2,2));
            pre_matrix(pre_sel_exc) = 1;
            pre_matrix(pre_sel_inh) = -1;

                    X_thisarea = pre_matrix(:,categories);
                    aux_X=X_thisarea;
                    aux_X(find(X_thisarea==-1))=1;
                    sel_cell_type = real_peak_thisrat;
                    aux3=is_play(this_area_index) & sel_cell_type(this_area_index);
                    sel_mod_0 = sum(aux_X,2)>0 & aux3;
                    X_thisarea_sel = X_thisarea(sel_mod_0,:);
                    variance_measure_mean_thisarea = 1 - mean(abs(X_thisarea_sel),2);
                    aux_depth = depth(this_area_index);
                    aux_depth_sel = aux_depth(sel_mod_0);
          
            sel_area2=find(strcmp(available_AREAS,your_area{areas_id})) ;

            if ~isempty(sel_area2)

                if sel_area2==numel(available_AREAS)
                    initial_area=PAG_columns(sel_area2);
                        ending_area=3840;
                else
                    initial_area=PAG_columns(sel_area2);
                    ending_area=PAG_columns(sel_area2+1);
                end
    

                    if ~isempty(aux_depth_sel)
                        values_per_bin={};
                        [values_per_bin, index, selected_idx] = estimate_wrapped_column_play_specificity(aux_depth_sel,  initial_area, ending_area,  wrap_bins, variance_measure_mean_thisarea');
                        output_bin = unique(index);
                        
                        if ~isempty(index)
                            for i = 1:numel(output_bin)

                               ind = output_bin(i);
                                
                               bins_this_area{ind} = [bins_this_area{ind}; values_per_bin{ind}'];
                                
                            end
                        end
                        
                    end
            end
            %%
    
            limit_day(bf_i)=numel(this_area_index);bf_i=bf_i+1;
                
    end

    acum=0;
    for bf2 = 1:numel(thisrat)
        limit_days(bf2)=acum+limit_day(bf2);
        acum=sum(limit_day(1:bf2));
    end
        
fntsize = 16;
  
%%

sel_exc=find(full_mat_exc==1);
sel_inh=find(full_mat_inh==1);
matrix = zeros(size(full_mat_exc,1),size(full_mat_exc,2));
matrix(sel_exc) = 1;
matrix(sel_inh) = -1;
mat=matrix;

mat_sig_betas = [mat_sig_betas; full_sig_betas];


mats_NE = [mats_NE ;full_non_entrained];
mats_Pk = [mats_Pk ;full_peak];
mats_Th = [mats_Th ;full_through];

mats_mod_play = [mats_mod_play ;full_modulated_play];
mats_mod_non_play = [mats_mod_non_play ;full_modulated_non_play];


mats = [mats; mat] ;
indx_area = ones(size(mat,1),1)*areas_id;
indexes_areas = [indexes_areas; indx_area];
play_sel = [play_sel; full_is_play];

    
    aux_warpped_specificity=[];
    for k = 1:wrap_bins
        if ~isempty(bins_this_area{k})
            aux_warpped_specificity(k) = mean(bins_this_area{k},"omitmissing");
        else
            aux_warpped_specificity(k) = NaN;
        end
    end

    full_play_specificity = [full_play_specificity; aux_warpped_specificity'];

end

       %% Plot classification overlap

                SELF_name   = cellfun(@(X) strcat('Self-', X), behaviors2check, 'UniformOutput',false);    
                OTHER_name  = cellfun(@(X) strcat('Other-', X), behaviors2check, 'UniformOutput',false);
                all_names = [SELF_name,OTHER_name];
                behaviors2check4={};

        %%

        X = mats(:,categories);
        co_modulation_number2 = sum(~mats(:, 1:18)==0,2);

        X_Play_and_nonPlay = mats(:, [categories, categories_non_play]);
        X_other = mats(:,categories+16); %% For Other-driven behaviors
 
        % sel_mod = sum(X,2)>0; %%%%% Select MODULATED cells only
        aux_X=X;
        aux_X(find(X==-1))=1;
        sel_mod = sum(aux_X,2)>0 & play_sel; %%%%% Select MODULATED cells and play only
        co_modulation_number2 = co_modulation_number2(sel_mod);
        X_sel = X(sel_mod,:);
        indexes_areas_sel = indexes_areas(sel_mod);


        X_other_sel = X_other(sel_mod,:);

        mat_sig_betas_sel =mat_sig_betas(sel_mod,:);
        
        mats_NE_sel = mats_NE(sel_mod);
        mats_Pk_sel = mats_Pk(sel_mod);
        mats_Th_sel = mats_Th(sel_mod);

        mats_mod_play_sel = mats_mod_play(sel_mod);
        mats_mod_non_play_sel = mats_mod_non_play(sel_mod);
        X_sel_Play_and_nonPlay = X_Play_and_nonPlay(sel_mod,:);

        [coeff, score, latent, tsquared, explained] = pca(zscore(X_sel));

        scores_80prc=1:9;

        
        %%
        figure('units','normalized','outerposition',[0 0 0.3 .3])
        cse=cumsum(explained);
        plot(scores_80prc,cse(scores_80prc))
        axis square

        %%
        eva = evalclusters(score(:,scores_80prc), 'linkage', 'gap', ...
                   'KList', 1:size(X_sel,2)*2, 'Distance','Correlation');
        bestK = eva.OptimalK;
        fprintf('Best number of clusters (gap statistic) = %d\n', bestK);
        %%

        nClusters = bestK;
        D = pdist(score(:,scores_80prc), 'euclidean'); 
        Z = linkage(D, 'ward'); % cluster on first 2 PCs
        leafOrder = optimalleaforder(Z, D);
        cluster_labels_responsetype = cluster(Z, 'maxclust', nClusters);

        Xr = X_sel(leafOrder,:);
        indexes_areas_r = indexes_areas_sel(leafOrder);
        cluster_labels_responsetype_r = cluster_labels_responsetype(leafOrder);
        Xr_Play_and_nonPlay = X_sel_Play_and_nonPlay(leafOrder,:);

        Xr_other = X_other_sel(leafOrder,:);

        mat_r_sig_betas = mat_sig_betas_sel(leafOrder,:);

        clusters = unique(cluster_labels_responsetype_r);
         
        mats_NE_r = mats_NE_sel(leafOrder);
        mats_Pk_r = mats_Pk_sel(leafOrder);
        mats_Th_r = mats_Th_sel(leafOrder);

        mats_mod_play_r = mats_mod_play_sel(leafOrder);
        mats_mod_non_play_r =mats_mod_non_play_sel(leafOrder);

        %%

        Stop = 0;

        %% Control: STABILITY per cluster %%
        %% Resampling-based assignment stability %%

        original_matrix = score;          % nNeurons x 9 matrix
        labels0  =  cluster_labels_responsetype;            % original cluster labels (1..K)
        K = nClusters;

        % % %  Bootstrap / subsample % % %

        nIter = 500;
        n = size(original_matrix,1);        
        C = zeros(n,n);   % co-clustering counts
        N = zeros(n,n);   % number of times both neurons appear

        for it = 1:nIter
            idx = randperm(n, round(0.8*n));
            Xsub = original_matrix(idx,:);
            alt_D = pdist(Xsub, 'euclidean'); 
            
            alt_Z = linkage(alt_D, 'ward');        % or your method
            labels_sub = cluster(alt_Z, 'maxclust', K);
            
            for i = 1:length(idx)
                for j = i+1:length(idx)
                    ii = idx(i);
                    jj = idx(j);
                    N(ii,jj) = N(ii,jj) + 1;
                    if labels_sub(i) == labels_sub(j)
                        C(ii,jj) = C(ii,jj) + 1;
                    end
                end
            end
        end

        P = nan(n,n);
        mask = N > 0;
        P(mask) = C(mask) ./ N(mask);
        P = P + P';

        for i = 1:size(original_matrix,1)

            sameCluster = labels0 == labels0(i);
            sameCluster(i) = false;

            if sum(sameCluster) == 0
                stability(i) = NaN;
            else
                stability(i) = mean(P(i, sameCluster), 'omitnan');
            end

        end

        for k = 1:K

            idxk = labels0 == k;
            clusterStability(k) = mean(stability(idxk));   

        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 4 PANEL A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%  Plot GLM-based neurons modulation by behavior with Cluster-based dedrogram - Unchanged

        Xr_excinh = Xr;
        indexes_areas_r_excinh = indexes_areas_r;
        cluster_labels_responsetype_r_excinh = cluster_labels_responsetype_r;
        excinh ='Modulated';


        figure('units','normalized','outerposition',[0 0 0.35 1]);
        

        subplot(4,2,[1 3 5 7]);
                                            % --- you already have this from your dendrogram call:

                                [h,T,outperm] = dendrogram(Z,0,'Orientation','left','Reorder',leafOrder, 'ClusterIndices', cluster_labels_responsetype);

                                % number of leaves (should equal numel(outperm))
                                nLeaves = numel(outperm);

                                % Prepare leafColors: one RGB per leaf position (in plotted order 1..nLeaves)
                                leafColors = nan(nLeaves,3);

                                ax = gca;
                                axDir = get(ax,'YDir'); % 'normal' or 'reverse'

                                for hh = 1:numel(h)
                                    xd = get(h(hh),'XData');
                                    yd = get(h(hh),'YData');

                                    % For 'left' orientation, leaves are Y positions roughly 1..nLeaves
                                    % find elements where yd is (close to) an integer between 1 and nLeaves
                                    leaf_mask = round(yd)==yd & yd>=1 & yd<=nLeaves; % crude integer check
                                    if any(leaf_mask)
                                        leaf_positions = unique(yd(leaf_mask));
                                        for lp = leaf_positions
                                            leafColors(lp, :) = get(h(hh), 'Color');
                                        end
                                    end
                                end

                                % If some leafColors remain NaN (rare), fill them with the nearest non-NaN color:
                                missing = find(any(isnan(leafColors),2));
                                if ~isempty(missing)
                                    validIdx = find(~any(isnan(leafColors),2));
                                    for m = missing'
                                        % find nearest valid index
                                        [~, nearest] = min(abs(validIdx - m));
                                        leafColors(m,:) = leafColors(validIdx(nearest), :);
                                    end
                                end

                                % Now T gives cluster id for each leaf position (in the plotted order)
                                % Build clusterColorMap(clusterID, :) by picking the first leaf that has that cluster
                                clusterIDs = unique(T);
                                nClusters = max(clusterIDs);
                                clusterColorMap = zeros(nClusters,3);
                                for ci = clusterIDs(:)'
                                    leaf_pos_for_cluster = find(T == ci, 1, 'first'); % take first leaf of that cluster
                                    clusterColorMap(ci, :) = leafColors(leaf_pos_for_cluster, :);
                                end

        %%%%%

        set(h,'LineWidth',1.5);
        axis off
        title('Ward clustering');
        
        % Right: reordered matrix
        subplot(4,2,[2 4 6 8]);
        aux_Xr =Xr; aux_Xr(aux_Xr==1)=0.8; aux_Xr(aux_Xr==-1)=-0.8; %aux_Xr(aux_Xr==0)=nan;
        imagesc(aux_Xr);
        colormap([ 0 0.2 0.8; 1 1 1; 0.7 0.2 0]);  % white, blue, red
        caxis([-1 1]);
        set(gca, 'YDir','normal') 

        yticks(1:2:numel(cluster_labels_responsetype(leafOrder)))
        aux_lab = cluster_labels_responsetype(leafOrder);
        aux_lab = aux_lab(1:2:numel(aux_lab));
        yticklabels(aux_lab)
        xticks(1:numel(all_names(categories)))
        xticklabels([all_names(categories)])
        title('Reordered matrix');

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                    %% Calculate the dispertion of dots around the clusters centroids

                  unique_clusters = unique(aux_lab,'stable');
                  unique_colors = unique(clusterColorMap, 'rows','stable');
                    
                    %%%%%%%% SUM %%%%%%%%
                    % variance_measure_alternative_mean = sum(abs(Xr),2);  %%%% FLIP to match the plot
                    % variance_measure_alternative_mean_nonPlay = sum(abs(Xr_Play_and_nonPlay(:,categories_non_play)),2);                  
                    % variance_measure_alternative_mean_Other = sum(abs(Xr_other),2);

                    %%%%%%%% MEAN %%%%%%%%
                    variance_measure_alternative_mean = mean(abs(Xr),2);  %%%% FLIP to match the plot
                    variance_measure_alternative_mean_nonPlay = mean(abs(Xr_Play_and_nonPlay(:,categories_non_play)),2);                  
                    variance_measure_alternative_mean_Other = mean(abs(Xr_other),2);
                    
                    Xr_aux = Xr;Xr_aux(Xr_aux==1)= NaN;
                    variance_measure_alternative_mean_neg = mean(Xr_aux,2,'omitmissing');
                    Xr_aux = Xr;Xr_aux(Xr_aux==-1)= NaN;
                    variance_measure_alternative_mean_pos = mean(Xr_aux,2,'omitmissing');


                    %
                    K = numel(unique_clusters);
                    
                    
                    swarmX = []; % for plotting mean/std if needed
                    
                    order_p=[];g_order_p=[];
                    for k = 1:K
                        idx = unique_clusters(k);

                        yvals= variance_measure_alternative_mean(cluster_labels_responsetype_r == idx);

                        order_e(k) = mean(yvals);%/sum(yvals);
                        order_p = [order_p; yvals];
                        g_order_p = [g_order_p; yvals*0+k];

                        sig_from_0(k) = signrank(yvals,1);
                    end

                    
                    [ord,reOrdered_cluster]=sort(order_e)  ; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 4 PANEL B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                     %% Ordered in ascending dispertion

                    uniqueClusters = unique_clusters(reOrdered_cluster);
                    clusterColors   = unique_colors(reOrdered_cluster,:);
                    figure("Units","normalized","OuterPosition",[0 0 0.2 .8]); hold on

                    colors_turbo = flip(turbo(17));%colors_turbo = flipud(turbo(17));

                    newOrder = [];
                    excess_com = [];
                    g_excess_com = [];
                      
                    xlim([0 K+1]);
                    % ylim([-1 1])
                    ylim([0 1])
% 
                    for k = 1:K

                        idx = uniqueClusters(k);

                        yvals = variance_measure_alternative_mean(cluster_labels_responsetype_r == idx);

                        yvals_nonPlay = variance_measure_alternative_mean_nonPlay(cluster_labels_responsetype_r == idx);

                        yvals_Other = variance_measure_alternative_mean_Other(cluster_labels_responsetype_r == idx);

                        yvals_ind = 1 - yvals;  Wlcx_yvals= 1;
                        % yvals_ind = 1 - yvals_Other;  Wlcx_yvals= 1;                      
                        % yvals_ind =  ((1 - yvals) - (1 -yvals_nonPlay)) ; Wlcx_yvals= yvals_nonPlay;
                        % yvals_ind =  ((1 - yvals) - (1 -yvals_Other)) ; Wlcx_yvals= yvals_Other;
                        % yvals_ind =  ((1 - yvals_Other) - (1 -yvals)) ; Wlcx_yvals= yvals_Other;

                        newOrder = [newOrder; find(cluster_labels_responsetype_r == idx)];

                        xvals = yvals_ind*0+k;

                        jitter = 0.750;
                        scatter(xvals + (rand(size(xvals))-0.5)*jitter, yvals_ind, 50, ...
                                'MarkerFaceColor', colors_turbo(k,:), 'MarkerEdgeColor', 'k', 'LineWidth',0.5);

                        % Mean and std
                        meanVal = mean(yvals_ind);
                        stdVal = std(yvals_ind);
                        errorbar(k, meanVal, stdVal, 'k', 'LineWidth', 2, 'CapSize', 10);

                        yline(0,'--k')

                        hold on

                        excess_com = [excess_com; yvals_ind];
                        g_excess_com = [g_excess_com; yvals_ind*0+k];

                    end                    
                    xticks(1:K);
                    xticklabels(uniqueClusters)
                    xlabel('Cluster ID');
                    ylabel('Play specifity (self - other)');
                    set(gca,'XDir','reverse')
                    view([90 -90])
                     

                    [rho,pval] = corr(excess_com,g_excess_com,"Type","Spearman");%,"Rows","pairwise")
                    title(['r= ' num2str(round(rho,2)) '; p= ' num2str(round(pval,5))])

                    %% STATS             
                    
                    [~,~,stats]=kruskalwallis(excess_com,g_excess_com)
                    multcompare(stats)

%% 

                                            %% Reordered matrix and dendrogram (BLACK AND WHITE) for FIGURE 4 PANEL A
                                                    %% Play vs Non-play and Self vs Other

        figure('units','normalized','outerposition',[0 0 0.15 1]);

        aux_Xr_Self_Other = [Xr Xr_other]; aux_Xr_Self_Other(aux_Xr_Self_Other==1)=0.8; aux_Xr_Self_Other(aux_Xr_Self_Other==-1)=-0.8;
        imagesc(flip(aux_Xr_Self_Other(newOrder,:))); 

        colormap([ 0 0 0; 0.7 0.7 0.7; 1 1 1]);
        caxis([-1 1]);
        set(gca, 'YDir','normal') 

        yticks(1:2:numel(cluster_labels_responsetype_r(newOrder)))
        aux_lab2 = cluster_labels_responsetype_r(newOrder);
        aux_lab2 = aux_lab2(1:2:numel(aux_lab2));
        yticklabels(flip(aux_lab2))
        xticks(1:numel(all_names(categories)))
        xticklabels([all_names(categories)])


        %% PLOT Control: STABILITY per cluster %%

            figure('units','normalized','outerposition',[0 0 0.6 0.3]); hold on;
    
            for ii = 1:K
                k = uniqueClusters(ii);
                
                idxk = labels0 == k;
                y = stability(idxk);
                y = y(~isnan(y));   % remove singletons
                mean_y = mean(y,"omitmissing");
                
                x = ii + 0.15*(rand(size(y)) - 0.5);  % horizontal jitter
    
                plot([ii-.3 ii+.3],[mean_y mean_y],'k',LineWidth=2)
                
                scatter(x, y, 50, ...
                    'filled', ...
                    'MarkerFaceColor', 'b', ...
                    'MarkerFaceAlpha', 0.4, ...
                    'MarkerEdgeColor', 'none');
            end
    
            yline(0.8,'r')
            yline(0.5,'r--')
            
            xlim([0.5 K+0.5]);
            ylim([0 2.1]);
            xticks(1:K);
            xticklabels(uniqueClusters);
            xtickangle(45);
            
            ylabel('Assignment stability');
            xlabel('Cluster (ordered by mean specificity)');
            title('Stability of neurons within each cluster');
            box off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 4 PANEL C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Figure Specificity panel C

figure('units','normalized','outerposition',[0 0 0.2 1]);
ii=1;

numRegions = 5;
for r = 1:numRegions
    
        sel = (indexes_areas_r_excinh == r); % --- Subset indices for this structure
            CL=uniqueClusters;
        for k = 1:numel(CL)
            cl_i = CL(k);
            prc(k,1) = sum(mats_NE_r(cluster_labels_responsetype_r==(cl_i) & sel));
            prc(k,2) = sum(mats_Pk_r(cluster_labels_responsetype_r==(cl_i) & sel));
            prc(k,3) = sum(mats_Th_r(cluster_labels_responsetype_r==(cl_i) & sel));
            total_this_cluster = prc(k,1) + prc(k,2) + prc(k,3);
            prc(k,1) = prc(k,1)/total_this_cluster;
            prc(k,2) = prc(k,2)/total_this_cluster;
            prc(k,3) = prc(k,3)/total_this_cluster;

        end

        for k = 1:numel(CL)
            cl_i = CL(k);
            prc_cum(k,1) = sum(mats_NE_r(cluster_labels_responsetype_r==(cl_i) & sel));
            prc_cum(k,2) = sum(mats_Pk_r(cluster_labels_responsetype_r==(cl_i) & sel));
            prc_cum(k,3) = sum(mats_Th_r(cluster_labels_responsetype_r==(cl_i) & sel));
        end

        

        for cell_class = 1:4
            subplot(5,4,cell_class+(4*(r-1)))
            
            if cell_class~=4
                br = bar(flip(prc(:,cell_class))*100,"stacked")
                % br(1).BarWidth=1.1
                xticks(1:numel(CL));
                
                xlim([0.5 17.5])
                ylim([0 100])
                view([90 -90]);

            else
                cdf_h=[];
                for i = 1:size(prc_cum,2)
                    cdf_h(:,i) = cumsum(flip(prc_cum(:,i)))/sum(prc_cum(:,i));
                    plot(cdf_h(:,i),1:numel(CL),LineWidth=2)
        
                    hold on
                end

                [d,p] = kstest2(cdf_h(:,1),cdf_h(:,2))
                [d,p ] = kstest2(cdf_h(:,2),cdf_h(:,3))

                ylim([1 17])
            end

            if cell_class == 1
                title(your_area{r})
            end
        end

end


                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Extended Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 6 PANEL A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %% Peak vs Trough specificity

 
                
                        select_cell_types = logical([mats_NE_r, mats_Pk_r, mats_Th_r]);
                        select_cell_types = select_cell_types(newOrder,:);        
                        cluster_labels_ordered = cluster_labels_responsetype_r(newOrder);

                                        % NONENTRAINED 1 - PEAK 2 - TROUGH 3 ... in that order
                        select_cell_types_3 = zeros(size(select_cell_types,2),1);
                        for c = 1:3
                            select_cell_types_3(select_cell_types(:,c)) = c;
                        end

                        for bhs = 1:9
                            select_cell_types_4(:,bhs) =  select_cell_types_3;
                        end


                        sel_area = indexes_areas_r_excinh == 3;
                        
                        Xr_new       = Xr(sel_area,:);
                        Xr_other_new = Xr_other(sel_area,:);
                        types        = select_cell_types_4(:,1);
                        types_new    = types(sel_area,:);
                        cluster_labels_responsetype_r_new = cluster_labels_responsetype_r(sel_area);

                       
                        self_spec  = 1 - mean(abs(Xr_new(:,categories)),2,'omitnan');
                        other_spec = 1 - mean(abs(Xr_other_new(:,categories)),2,'omitnan');
                        specificity = self_spec - other_spec;
                        
                        peak   = specificity(types_new==2);
                        trough = specificity(types_new==3);
                        
                        figure; hold on
                        
                        jitter = 0.25;
                        
                        scatter(ones(size(peak))   +(rand(size(peak))-0.5)*jitter,...
                            peak,200,'MarkerEdgeColor',[1 .2 .2])
                        
                        scatter(2*ones(size(trough))+(rand(size(trough))-0.5)*jitter,...
                            trough,200,'MarkerEdgeColor',[.2 .2 1])
                        
                        errorbar(1,mean(peak),std(peak),'k','LineWidth',2)
                        errorbar(2,mean(trough),std(trough),'k','LineWidth',2)
                        
                        [p,~,stats] = ranksum(peak,trough);
                        
                        xlim([0.5 2.5])
                        xticks([1 2])
                        xticklabels({'Peak','Trough'})
                        ylabel('Specificity')
                        title(sprintf('Rank-sum p = %.4g',p))
                        box off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Extended Data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 6 PANEL B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %% Excited/Inhibited composition per delta locking type

            % for r = 1:numRegions

            sel_area = (indexes_areas_r_excinh == 3);

            Xr_self_other = [Xr Xr_other];
            types         = select_cell_types_4;

            % filter_by_area
            Xr_self_other = Xr_self_other(sel_area,:);
            types         = types        (sel_area,:);

            categories_full       = [categories categories*2 ];
            categories_full_index = [categories categories+16 ];

            figure
            
            delta_names = {'Unlocked','Peak','Trough'};

            
            for delta_type = 1:3
            
                subplot(1,3,delta_type)
                hold on
            
                % Select neurons belonging to this delta category
                delta_idx = types(:,1) == delta_type;
            
                betas = Xr_self_other(delta_idx,:);
            
                % Only behaviors
                betas = betas(:,categories_full);
            
                % Percentage excited and inhibited per behavior
            
                perc_exc = zeros(1,size(betas,2));
                perc_inh = zeros(1,size(betas,2));
                perc_ratio = zeros(1,size(betas,2));
            
            
                for bh = 1:size(betas,2)
            
                    b = betas(:,bh);
            
                    % remove NaNs
                    b = b(~isnan(b));
            
                    perc_exc(bh) =  0%sum(b>0)/numel(b)*100;
                    perc_inh(bh) =  0%sum(b<0)/numel(b)*100;
                    perc_ratio(bh) = (sum(b>0)/numel(b)*100) + (sum(b<0)/numel(b)*100);
            
                end
            
            
                % Matrix for stacked bar
                % First row = inhibited
                % Second row = excited
            
                perc_matrix = [perc_inh; perc_exc; perc_ratio]';
            
            
                b = bar(perc_matrix,'stacked');
            
            
                b(1).FaceColor = [0.3 0.3 1]; % inhibited
                b(2).FaceColor = [1 0.3 0.3]; % excited
                b(3).FaceColor = [0.3 1 0.3]; % excited
            
            
                ylim([0 50])
            
                xticks(1:length(categories_full_index))
                xticklabels([all_names(categories_full_index)])
                xtickangle(45)
            
                ylabel('% neurons')
            
                title(delta_names{delta_type})
            
                box off
            
            end
            
            
            legend({'Inhibited','Excited','Difference'})
            sgtitle('LPAG')

