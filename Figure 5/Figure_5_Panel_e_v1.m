%%  Based on Play_map_Dynamic_GLM_Cell_Res_Class_EarlyLateSust_Pounce_Entrained.m
%%  Matias Mugnaini 

%% Load ENTRAINED cells

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


%% Set to load individual behaviors PSTH

chose_beh = 1; chose_role = 1; play_cond = 0;

behavior_labels = {'POA', 'PWIA'};
base_conditions = {'POA', 'PWIA'};

behavior_labels_title = behavior_labels{chose_beh};
base_conditions_title = base_conditions{chose_beh};

choose_cell_type = 1; % 1 for Trough, 2 for Peak, 3 for NON-entrained

roles = {'self', 'other'};
roles_title = roles{chose_role};


% Add the common base variables
vars_to_load = [ "all_psth_self_warped", "all_psth_other_warped","all_psth_shuffled_self_warped","all_psth_shuffled_other_warped"];

% Add the per-behavior ones dynamically
% for b = 2:numel(behavior_labels) % start from 2 to skip 'play' (already added above)
    label = behavior_labels_title;
    vars_to_load = [vars_to_load, ...
        sprintf("all_psth_self_warped_%s", label), ...
        sprintf("all_psth_other_warped_%s", label)...
        sprintf("all_psth_shuffled_self_warped_%s", label), ...
        sprintf("all_psth_shuffled_other_warped_%s", label) ];
   
% end


%% Select rats

B1D1 = [1 2 3];
B1S3 = [4 5 6];
B2S2 = [7 8 9];
B3D2 = 10;


%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%CHECK THISSSSSS
B4D4 = [11 11] ; %   0 if medial 1 if lateral %% Dual 4rd Batch Probe 1
B4S2 = [12 12];                              %% Single 4rd Batch Probe 1

Npx2Probe={'Medial','Lateral'};
medial_lateral_probe=[0 0 0 0 0  0 0 0 0 0  1 0  0 0 0  0 1  ]; %   1 if medial and 0 if lateral
probes=[              1 1 1 1 1  1 1 1 1 0  3 1  1 1 1  1 3  ]; %    3 is medial and 1 is lateral     for B4S2, oposite if it is B4D4

MorL = {'_','_','_', '_','_', '_', '_','_','_', '_', 'Med','Lat', 'mPFC','mPFC','mPFC' , 'Med' , 'Lat'};
thisrat=[ B1D1 B1S3 B2S2 B3D2 B4S2 B1D1 B4D4];


PAG_mPFC=[            1 1 1 1 1  1 1 1 1 1  1 1  2 2 2  1 1  ]; % 1 for PAG , 2 for mPFC

length_duration_threshold=0.25;
wrap_bins=20;
alpha = 0.05%/5;

    ONOFFsets = 2;
%%


thisrat=[ B1D1 B1S3 B2S2 B3D2 B4S2  ];minus_areas=0;
bf_i = 1;


your_area={'DR','VLPAG','LPAG','DLPAG','SupCol'};

early_late=0; early_late_title = 'Early Onset'; % 0 for early and 5 for late ONset
% early_late=5; early_late_title = 'Late Onset'; % 0 for early and 5 for late ONset

%%

behaviors2check={'Pounce_A','CC','Pin','Boxing','Evasion','CB','Pounce_B','Escape','CD', ...
                'Rearing','Grooming','Scratch','Pounce_Ai','Pounce_Bi','Sniffing','Bite'};


    exc = 1; exc_title = 'Exc';mycol = [0.8500 0.3250 0.0980];
    % exc = 0; exc_title = 'Inh';  mycol = [0.3010 0.7450 0.9330]; 
    % exc = 2; exc_title = 'Mod';  mycol = [0.3010 0.7450 0.9330]; 
    lim_y=100;
    
    other = [ 0 numel(behaviors2check)-1]; other_title = {'Self','Other'}; 

        beh_to_choose = 1 + other; facealpha = 0.5; beh_title = 'Play Bout';      


    %%

AREAS=[];
depth = [];
full_areas = [];

all_pvalues=[];
all_coefficients= [];

full_dyn_res_warp=[];
full_pie_exc_warp=[];
full_pie_inh_warp=[];

psth_counts_exc={};
psth_counts_inh={};

behavior_files = dir('*.txt');
playbout_tittle='PlayBout';
% pval_th = 0.05;
 collected_psth_exc = cell(numel(your_area),5);
 collected_psth_inh = cell(numel(your_area),5);
 ids_exc = cell(numel(your_area),5);
 ids_inh = cell(numel(your_area),5);

 collected_psth = cell(numel(your_area),3);

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
       
        path_migue='Y:\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Responses_Matrix\ModelCriterion_Onset_Miguel\';

                 load([path_migue 'ResponsesMatrix_PPB_p1andp2_' num2str(length_duration_threshold) 's_' playbout_tittle '_' mydate '_' behavior_files(bf).name(1:4) '_' MorL{bf_i} '.mat'] , vars_to_load{5:8})


        depth = depth_or_Chn;
        ids_gcl = good_clusters;

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
            non_entrained_thisrat = non_entrained(sel_session);      non_entrained_thisrat = non_entrained_thisrat(ord_by_depth);
            trough_cells_thisrat  = trough_cells (sel_session);      trough_cells_thisrat  = trough_cells_thisrat (ord_by_depth); 
            real_peak_thisrat     = real_peak    (sel_session);      real_peak_thisrat     = real_peak_thisrat    (ord_by_depth);       
             modulated_play_thisrat =  modulated_play(sel_session); modulated_play_thisrat = modulated_play_thisrat(ord_by_depth);
             modulated_non_play_thisrat = modulated_non_play(sel_session);  modulated_non_play_thisrat = modulated_non_play_thisrat(ord_by_depth);
        end


        if choose_cell_type ==1 
            selected_cell_type = trough_cells_thisrat;
        elseif choose_cell_type == 2
            selected_cell_type = real_peak_thisrat;
        else
            selected_cell_type = non_entrained_thisrat;
        end

             %% ON SET and OFF SET modulation indexes WARP

              if  play_cond
                real_var     = sprintf('all_psth_%s_warped', roles{chose_role});
                shuffled_var = sprintf('all_psth_shuffled_%s_warped', roles{chose_role});
              else
                real_var     = sprintf('all_psth_%s_warped_%s', roles{chose_role}, label);
                shuffled_var = sprintf('all_psth_shuffled_%s_warped_%s', roles{chose_role}, label);
              end

            all_this_psth_indbeh = [];
            all_this_psth_indbeh = eval(real_var);

            full_zsc_warp_shuffle  = [];
            full_zsc_warp_shuffled = eval(shuffled_var);

            for raw_n=1:size(full_zsc_warp_shuffled,1)
                    mu = mean(all_this_psth_indbeh(raw_n,:), 'omitnan');
                    sd = std(all_this_psth_indbeh(raw_n,:), 'omitnan');
                    all_this_psth_indbeh(raw_n,:) = (all_this_psth_indbeh(raw_n,:) - mu) / sd;
                for raw_m = 1:size(full_zsc_warp_shuffled{raw_n},1)
                    full_zsc_warp_shuffled{raw_n}(raw_m,:) = (full_zsc_warp_shuffled{raw_n}(raw_m,:) - mu) / sd;
                    
                end
            end

        time_indxs=cell(5,1);
        time_indxs{1}= 16:18 ; time_indxs{2}= 22:27 ; time_indxs{3}= 33:38 ; time_indxs{4}= 21:38 ; time_indxs{5}= 42:45 ;

        FR_per_period=cell(5,1);
        FR_per_period{1}=mean(all_this_psth_indbeh(:,time_indxs{1}),2,"omitmissing");
        FR_per_period{2}=mean(all_this_psth_indbeh(:,time_indxs{2}),2,"omitmissing");
        FR_per_period{3}=mean(all_this_psth_indbeh(:,time_indxs{3}),2,"omitmissing");
        FR_per_period{4}=mean(all_this_psth_indbeh(:,time_indxs{4}),2,"omitmissing");
        FR_per_period{5}=mean(all_this_psth_indbeh(:,time_indxs{5}),2,"omitmissing");
  
        %% Inhibition and Excitation based on LOCAL difference with shuffle

    significant_indexes_exc=[];
    significant_indexes_inh=[];
    
    for j=1:size(full_zsc_warp_shuffled,1)

        % Early %%%%%%%%%%%%%%%%%%%%

        for t = 1:numel(time_indxs)
  
            mean_sh = mean(mean(full_zsc_warp_shuffled{j}, 'omitnan'));

            %%%% Compare against same period
            sh_FR_after_onset=mean(full_zsc_warp_shuffled{j}(:,time_indxs{t}),2,"omitmissing");
            
            %%%% Compare against all shuffle
            % sh_FR_after_onset=mean(full_zsc_warp_shuffled{j},2,"omitmissing");

            if FR_per_period{t}(j) > mean_sh 
                p_exc = mean(FR_per_period{t}(j) < sh_FR_after_onset);
                significant_indexes_exc(j,t) = (p_exc < alpha);
                significant_indexes_inh(j,t) = 0;
            elseif FR_per_period{t}(j) < mean_sh  
                p_inh = mean(FR_per_period{t}(j) > sh_FR_after_onset);
                significant_indexes_inh(j,t) = (p_inh < alpha);
                significant_indexes_exc(j,t) = 0;
            else
                significant_indexes_exc(j,t) = 0;
                significant_indexes_inh(j,t) = 0;
            end
        end
    end   

        significant_indexes_exc(:,4)=significant_indexes_exc(:,2) & significant_indexes_exc(:,3);
    significant_indexes_exc(:,2)=significant_indexes_exc(:,2);
    significant_indexes_exc(:,3)=significant_indexes_exc(:,3);

        significant_indexes_inh(:,4)=significant_indexes_inh(:,2) & significant_indexes_inh(:,3);
    significant_indexes_inh(:,2)=significant_indexes_inh(:,2);
    significant_indexes_inh(:,3)=significant_indexes_inh(:,3);

            significant_indexes_exc = logical(significant_indexes_exc) & selected_cell_type;
            significant_indexes_inh = logical(significant_indexes_inh) & selected_cell_type;
    
    %% Divide by Recorded Areas

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

       PAG_columns= [];
       start_out= [];
       sel_area = [];
       [PAG_columns,available_AREAS,allareas] = take_PAG_area(ratID,probe,probe_area,[]);
       initial_cluster=[1 initial_cluster];
      
        % bf
        % AREAS
        for j =1 :numel(your_area)

            correc_area_index = j; %find(ismember(your_area, AREAS{j }));

            sel_area=[];
            sel_area=find(strcmp(available_AREAS,your_area{j})) ;

            if ~isempty(sel_area)

                if sel_area==numel(available_AREAS)                
                    init_area=PAG_columns(sel_area);
                    end_area=3840;
                else
                    init_area=PAG_columns(sel_area);
                    end_area=PAG_columns(sel_area+1);
                end

                
    
                this_area_index = depth > init_area & depth < end_area;
                aux=all_this_psth_indbeh(this_area_index,:); this_area_index = [];
                collected_psth{correc_area_index}=[collected_psth{correc_area_index}; aux]; aux=[];
                
                for t = 1:numel(time_indxs)
                    this_area_index = depth > init_area & depth < end_area & significant_indexes_exc(:,t)';
                    if ~isempty(this_area_index)
                        aux=all_this_psth_indbeh(this_area_index,:); 
                        collected_psth_exc{correc_area_index,t}=[collected_psth_exc{correc_area_index,t}; aux]; aux=[];
                        aux_ids (:,1)= depth(this_area_index); aux_ids (:,2)= ids_gcl(this_area_index);
                        ids_exc{correc_area_index,t} = [ids_exc{correc_area_index,t}; aux_ids];aux_ids=[];
                    end

                    this_area_index = [];
                    this_area_index = depth > init_area & depth < end_area & significant_indexes_inh(:,t)';
                    if ~isempty(this_area_index)
                        aux=all_this_psth_indbeh(this_area_index,:); 
                        collected_psth_inh{correc_area_index,t}=[collected_psth_inh{correc_area_index,t}; aux]; aux=[];
                        aux_ids (:,1)= depth(this_area_index); aux_ids (:,2)= ids_gcl(this_area_index);
                        ids_inh{correc_area_index,t} = [ids_inh{correc_area_index,t}; aux_ids];aux_ids=[];
                    end
                    this_area_index = [];
                end
            end

        end
        bf_i = bf_i + 1;
    end

%%

        if exc==1
            counts = collected_psth_exc;
        elseif exc==0
            counts = collected_psth_inh;
        else
            counts = collected_psth;
        end
        
res_cat_num = 1:numel(your_area);
percentage=[];
for n_area=1:numel(your_area)
    for t_ind =  1:numel(time_indxs)

        ids_exc_sel = ids_exc{6-n_area,t_ind};
               
        other_res_cat = find(~ismember(res_cat_num,t_ind));
        for j =other_res_cat
           aux = ismember(ids_exc{6-n_area,j},ids_exc_sel); aux(:,1) = sum(aux,2)>1;
           percentage(6-n_area,t_ind,j)=mean(aux(:,1),"omitmissing")*100;
        end
    
    end
end


n_areas=5;


%% Save
    % save([path 'All_rats_' roles_title '_' behavior_labels_title '.mat'], ...
    %     'collected_psth_exc', 'collected_psth_inh', 'collected_psth', 'ids_exc', 'ids_inh')

%%

% load([path 'All_rats_' roles_title '_' behavior_labels_title '.mat'], ...
%     'ids_exc', 'ids_inh','collected_psth')

stop = 1;
%%

path='Y:\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Responses_Matrix\ModelCriterion\';
your_area={'DR','VLPAG','LPAG','DLPAG','SupCol'};
time_indxs=cell(5,1);
time_indxs{1}= 16:18 ; time_indxs{2}= 22:27 ; time_indxs{3}= 33:38 ; time_indxs{4}= 21:38 ; time_indxs{5}= 42:45 ;

y_lim = 60;

n_areas=5;

res_cat_num = 1:numel(your_area);
percentage={};

 

            overlap_index = [];
            overlap_ps = [];

for i = 1:2
    if i == 1 
        chose_beh_lab = 'POA'; chose_role_lab = 'self';
        load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '_' num2str(choose_cell_type) '.mat'], 'ids_exc','collected_psth')
        IDS = ids_exc;
        % N_behA = 

        % chose_beh_lab = 'POA'; chose_role_lab = 'other';
        % load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '.mat'], 'ids_exc')
        % IDS_comp = ids_exc;

        chose_beh_lab = 'PWIA';
        % chose_beh_lab = 'BI';
        load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '_' num2str(choose_cell_type) '.mat'], 'ids_exc')
        IDS_comp = ids_exc;
    else
        chose_beh_lab = 'POA'; chose_role_lab = 'self';
        load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '_' num2str(choose_cell_type) '.mat'], 'ids_inh','collected_psth')
        IDS = ids_inh;

        % chose_beh_lab = 'POA'; chose_role_lab = 'other';
        % load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '.mat'], 'ids_inh')
        % IDS_comp = ids_inh;

        chose_beh_lab = 'PWIA';
        % chose_beh_lab = 'BI';
        load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '_' num2str(choose_cell_type) '.mat'], 'ids_inh')
        IDS_comp = ids_inh;
    end

    for n_area=1:numel(your_area)
        for t_ind =  1:numel(time_indxs)
    
            ids_exc_inh_sel = IDS{6-n_area,t_ind};
            ids_exc_inh_sel_comp = IDS_comp{6-n_area,t_ind};

            if isempty(ids_exc_inh_sel) || isempty(ids_exc_inh_sel_comp)
                percentage{i}(n_area, t_ind, :) = NaN;
                continue
            end

                    overlap = ismember(ids_exc_inh_sel, ids_exc_inh_sel_comp, 'rows'); 

                    percentage{i}(n_area, t_ind) = sum(overlap)/(size(ids_exc_inh_sel,1) + size(ids_exc_inh_sel_comp,1)) * 100;


            %%%% OVERLAP indexes %%%%
            A               = ids_exc_inh_sel;
            B               = ids_exc_inh_sel_comp;        
            nA = size(A,1);
            nB = size(B,1);
            N  = size(collected_psth{6-n_area},1);
            nAB = sum(overlap);
            P = (nA/N) * (nB/N);
            expectedAB = N * P;
            overlapIndex = nAB / expectedAB;
            p_greater = 1 - binocdf(nAB - 1, N, P);   % more overlap than expected
            p_lower   = binocdf(nAB, N, P);           % less overlap than expected
            overlap_index (n_area,t_ind,i) = overlapIndex;
            overlap_ps (n_area,t_ind,i) = p_lower;

            %%%%%%%%%% MODULATED Cells %%%%%%%%%%%

            if i == 1
                chose_beh_lab = 'POA'; chose_role_lab = 'self';
                load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '_' num2str(choose_cell_type) '.mat'], 'ids_inh','collected_psth')
                IDS2 = ids_inh;
        
                % chose_beh_lab = 'POA'; chose_role_lab = 'other';
                % load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '.mat'], 'ids_inh')
                % IDS_comp2 = ids_inh;
        
                % chose_beh_lab = 'BI';
                chose_beh_lab = 'PWIA';
                load([path 'All_rats_' chose_role_lab '_' chose_beh_lab '_' num2str(choose_cell_type) '.mat'], 'ids_inh')
                IDS_comp2 = ids_inh;

                ids_exc_inh_sel2 = IDS2{6-n_area,t_ind}; %poa self inh
                ids_exc_inh_sel_comp2 = IDS_comp2{6-n_area,t_ind}; %PWIA self inh
 
                nA = size(ids_exc_inh_sel,1) + size(ids_exc_inh_sel2,1);
                nB = size(ids_exc_inh_sel_comp,1) + size(ids_exc_inh_sel_comp2,1);
                N  = size(collected_psth{6-n_area},1);

                overlap2 = ismember(ids_exc_inh_sel2, ids_exc_inh_sel_comp2, 'rows'); 

                nAB = sum([overlap;overlap2]);
                P = (nA/N) * (nB/N);
                expectedAB = N * P;
                overlapIndex = nAB / expectedAB;
                p_greater = 1 - binocdf(nAB - 1, N, P);   % more overlap than expected
                p_lower   = binocdf(nAB, N, P);           % less overlap than expected
                overlap_index_mod (n_area,t_ind,1) = overlapIndex;
                overlap_ps_mod (n_area,t_ind,1) = p_lower;
            end
        end
    end

    
end

%%
figure("Units","normalized",'outerposition',[0 0 0.2 0.85])

for t_ind =  1:numel(time_indxs)
    for n_area=1:numel(your_area)
    
            subplot(numel(your_area),1,n_area)
        
        counts_sel_exc = percentage{1}(n_area, t_ind);
        counts_sel_inh = percentage{2}(n_area, t_ind);

        % overlap_ps
        
        if counts_sel_exc > counts_sel_inh
            yl_txt = counts_sel_exc + 5;
        else
            yl_txt = counts_sel_inh + 5;
        end
        if overlap_ps(n_area,t_ind,1) < alpha && ~(overlap_ps(n_area,t_ind,1) == 0)
            text(t_ind-.2, yl_txt, '*', 'Color','k', 'FontSize',14, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold');  
        end
        if overlap_ps(n_area,t_ind,2) < alpha && ~(overlap_ps(n_area,t_ind,2) == 0)
            text(t_ind+.2, yl_txt, '*', 'Color','k', 'FontSize',14, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontWeight','bold');  
        end

        h=bar(t_ind,[counts_sel_exc counts_sel_inh]',"grouped");
        h(1).FaceColor=[.8 .5 .2];
        h(2).FaceColor=[.2 .5 .8];
        h(1).BarWidth=1;
        h(2).BarWidth=1;
        hold on
    
        xline(1.5,'g',LineWidth=3)
        xline(4.5,'r',LineWidth=3)
        ylim([0 y_lim])
        xticks(1:5)
        xlim([0 6])
        xticklabels({'Ant','Early','Late','Sust','Post'})
        title(your_area{6-n_area})
        axis square
                
    end
end


sgtitle(['Overlap percentages for Self POA and PWIA'  ])
% sgtitle(['Overlap percentages for Self POA and Other POA'  ])
