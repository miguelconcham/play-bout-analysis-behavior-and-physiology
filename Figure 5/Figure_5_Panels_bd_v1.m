%%  Based on Play_map_Dynamic_GLM_Wrap_Pounce_EntrainmentComparison.m
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

behavior_labels = {'POA', 'PWIA', 'CH', 'ES'};
base_conditions = {'POA', 'PWIA', 'CH', 'ES'};

behavior_labels_title = behavior_labels{chose_beh};
base_conditions_title = base_conditions{chose_beh};

choose_cell_type = 1; % 1 for Through, 2 for Peak, 3 for NON-entrained

    exc = 1; mycol = [0.8500 0.3250 0.0980];
    % exc = 0;  mycol = [0.3010 0.7450 0.9330]; 

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

MorL = {'_','_','_','_','_', '_','_','_','_','_', 'Med','Lat', 'mPFC','mPFC','mPFC' , 'Med' , 'Lat'};
% thisrat=[ B1D1 B1S3 B2S2 B3D2 B4S2 B1D1 B4D4];



PAG_mPFC=[            1 1 1 1 1  1 1 1 1 1  1 1  2 2 2  1 1  ]; % 1 for PAG , 2 for mPFC

length_duration_threshold=0.25;
wrap_bins=20;
alpha = 0.05;

    time_indxs=cell(4,1);
    time_indxs{1}= 1:5 ; time_indxs{2}= 7:11 ; time_indxs{3}= 21:25 ;time_indxs{4}= 26:30;
    ONOFFsets = 2;
%%


thisrat=[ B1D1 B1S3 B2S2 B3D2 B4S2  ];minus_areas=0;
bf_i = 1;

your_area={'DR','VLPAG','LPAG','DLPAG','SupCol'};


%% % Select Modulation Type %

    %%%%


    lim_y=45;
    

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
pval_th = 0.05;
 collected_psth_exc = cell(numel(your_area),1);
 collected_psth_inh = cell(numel(your_area),1);

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
        time_indxs{1}= 16:18 ; time_indxs{2}= 21:28 ; time_indxs{3}= 31:38 ; time_indxs{4}= 21:38 ; time_indxs{5}= 41:44 ;

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

            significant_indexes_exc = any(significant_indexes_exc(:,[2 3 ]),2) & selected_cell_type;
            significant_indexes_inh = any(significant_indexes_inh(:,[2 3 ]),2) & selected_cell_type;
    
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

        for j =1 :numel(your_area)
            
            correc_area_index = j; %find(ismember(your_area, your_area{j }));

            if isempty(collected_psth_exc{correc_area_index})
                collected_psth_exc{correc_area_index}=cell(wrap_bins,1);
            end
            if isempty(collected_psth_inh{correc_area_index})
                collected_psth_inh{correc_area_index}=cell(wrap_bins,1);
            end
            
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
    
                
                    if isempty(depth(significant_indexes_exc'))
                        psth_counts_exc{correc_area_index}(bf,:)=nan(1,wrap_bins);
                    else
                        [ aux_exc , numbers_per_area,index,neto_exc_selected] = estimate_wrapped_column(depth,  init_area, end_area,  wrap_bins, significant_indexes_exc');
                        output_bin = unique(index);

                        psth_counts_exc{correc_area_index}(bf_i,:)=(aux_exc./numbers_per_area)*100;
                    end
                
                    if isempty(depth(significant_indexes_inh))
                        psth_counts_inh{correc_area_index}(bf,:)=nan(1,wrap_bins);
                    else
                        [ aux_inh , numbers_per_area,index,neto_inh_selected] = estimate_wrapped_column(depth, init_area, end_area, wrap_bins, significant_indexes_inh');
                      
                        psth_counts_inh{correc_area_index}(bf_i,:)=(aux_inh./numbers_per_area)*100;
                    end

            end

        end
        bf_i = bf_i + 1;
    end

    %%
    %%%%

    figure('units','normalized','outerposition',[0 0 0.1 0.7])

    if exc==1
        counts = psth_counts_exc;
        
    else
        counts = psth_counts_inh;
    end
    
    counts_full=[];
    for stitch_i = 1 :numel(your_area)

        counts_full = [counts_full; nanmean(counts{stitch_i})'];

        % aux = sum(counts{stitch_i})/full_numbers_per_area{stitch_i}*100;
        % counts_full = [counts_full; aux'];

    end

    % percentages = counts_full/sum(counts_full)*100;

    counts_full(isnan(counts_full))=0;
    percentages = counts_full;
    
    full_areas = unique(full_areas);
    bins=1:20:max(depth);
    up_ylim=3840;
    pval_th = 0.05;
    sigma =12;
    yLabel='Percentage of modulated cells';

        res_area=boxFilter(percentages',5);

        ar=area(res_area); ar.FaceColor = mycol;
        ar.FaceAlpha = 0.25;
    
        xlim([0 numel(your_area)*wrap_bins])
        ylim([0 lim_y])

        
        ylabel(yLabel)
            xtic=xticks;
            camroll(+90)
            set(gca,'YDir','reverse') 

        for jj=1:numel(your_area)
    
            xline(jj*wrap_bins,'--','Color','r','LineWidth',3)
            txt = your_area(jj);
            text(jj*wrap_bins-1.5,lim_y-5,txt,'FontSize',14,'Color','red')
        end
    
        hold on


    sgtitle([roles_title ' - ' behavior_labels_title])
    


%%

psth_length = 60;
max_range = 10:40;
    if exc==1
        counts = psth_counts_exc;
        
    else
        counts = psth_counts_inh;
    end
    max_location = [];
    max_values = [];
   this_area_psth_full=[];
    for stitch_i = 1 :numel(your_area)
        this_area_psth = nan(numel(collected_psth_exc{stitch_i}), psth_length);;
        for bn = 1:numel(collected_psth_exc{stitch_i})

            if size(collected_psth_exc{stitch_i}{bn},1)==1
                this_area_psth(bn,:) =smooth(collected_psth_exc{stitch_i}{bn},5);
            elseif size(collected_psth_exc{stitch_i}{bn},1)>1
                this_area_psth(bn,:) =smooth(mean(collected_psth_exc{stitch_i}{bn}, 'omitmissing'),5);
            end            
        end

        for t=1:size(this_area_psth,2)
            this_area_psth(:,t) = smooth( this_area_psth(:,t),4);
        end
        [max_value_this_area,max_loc_this_area] = max(this_area_psth(:,max_range), [],2);
        max_loc_this_area = max_range(max_loc_this_area)';
        max_loc_this_area(all(isnan(this_area_psth),2)) = NaN;
        max_location = [max_location;max_loc_this_area];
        max_values = [max_values;max_value_this_area];

        this_area_psth_full = [this_area_psth_full; this_area_psth];
    end


%%
        figure
        colormap("jet");
        imagesc(this_area_psth_full)
        axis xy
        hold on
        plot([20 20], [1 100], 'g',LineWidth=4)
        plot([40 40], [1 100], 'r',LineWidth=4)
        plot(max_location, 1:numel(max_location), 'r',LineWidth=3)
        clim([-.5 1])

         for jj=1:numel(your_area)
    
            yline(jj*wrap_bins,'--','Color','k','LineWidth',3)
            txt = your_area(jj);
            text(psth_length-5,jj*wrap_bins-1.5,txt,'FontSize',14,'Color','red')
        end

       %%
