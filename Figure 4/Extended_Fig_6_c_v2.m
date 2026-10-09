%%  Based on Play_map_Dynamic_GLM_Wrap_Decoder_V2_EntrainmentComparison_stacked.m
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
PAG_mPFC=[            1 1 1 1 1  1 1 1 1 1  1 1  2 2 2  1 1  ]; % 1 for PAG , 2 for mPFC

length_duration_threshold=0.25;
wrap_bins=20;
alpha = 0.05;
bxfilterpar=8;
nperm = 1000;

lim_y=40;


%%

thisrat=[ B1D1 B1S3 B2S2 B3D2 B4S2  ];minus_areas=0;

your_area={'DR','VLPAG','LPAG','DLPAG','SupCol'};


%% % Select BEHAVIOR %

                sel_beh = 13;selforother = 1;

behaviors2check={'Pounce_A','CC','Pin','Boxing','Evasion','CB','Pounce_B','Escape','CD', ...
                'Rearing','Grooming','Scratch','Pounce_Ai','Pounce_Bi','Sniffing','Bite'};

    %%

AREAS=[];
depth = [];
full_areas = [];

all_pvalues=[];
all_coefficients= [];

full_dyn_res_warp=[];
full_pie_exc_warp=[];
full_pie_inh_warp=[];

behaviors_to_check_labels = {'Pounce neck', 'Pounce w/Imm','Discriminant'};
behaviors_to_check = numel(behaviors_to_check_labels); % Pounce, Pounce w/Imm and Separation
selfother_title = {'Self','Other','SelfvsOther'};sel_selfother = 3;

p_accuracy_all       = cell(numel(your_area),behaviors_to_check);
p_accuracy_all_shuff = cell(numel(your_area),behaviors_to_check);
all_aux_exc          = cell(numel(your_area),behaviors_to_check);
all_numbers_per_are  = cell(numel(your_area),behaviors_to_check);

for j=1:numel(your_area)
    for k=1:behaviors_to_check
        p_accuracy_all{j,k}      = nan(numel(thisrat), wrap_bins);
        all_aux_exc{j,k}         = zeros(1,wrap_bins);
        all_numbers_per_are{j,k} = zeros(1,wrap_bins);
    end
end

ids_exc = cell(numel(your_area),behaviors_to_check);

behavior_files = dir('*.txt');

playbout_tittle='PlayBout';

pval_th = 0.05;
collected_psth_exc = cell(numel(your_area),behaviors_to_check);

all_choose_cell_type = [1, 2, 3];

 for sel_behavior = 1:behaviors_to_check 

    all_percentage = nan(100,3);

    for cell_type_ind = 1:3

        choose_cell_type = all_choose_cell_type(cell_type_ind);

        bf_i = 1;
    
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
           
             path=['Y:\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Responses_Matrix\ModelCriterion\Decoder\' selfother_title{sel_selfother} '\'];
    
                    load([path 'ResponsesMatrix_PPB_p1andp2_' num2str(length_duration_threshold) 's_' playbout_tittle '_' mydate '_' behavior_files(bf).name(1:4) '_' MorL{bf_i} '.mat'], ...
                        "accuracy_distr","auc_distr","beta_distr","beta_p_distr","confusion_distr","p_accuracy_distr","p_auc_distr","p_confusion_distr")
    
    
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
            
                %% Set Decoder values
            
            sel_behavior
            p_acc = p_accuracy_distr{sel_behavior}(:,1) < alpha & selected_cell_type;
            % accuracy_distr{sel_behavior}(:,1);
            
            sum(p_acc) 
        
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
                if ~isempty(correc_area_index)
    
                    if isempty(collected_psth_exc{correc_area_index,sel_behavior})
                        collected_psth_exc{correc_area_index,sel_behavior}=cell(wrap_bins,1);
                    end
    
                    sel_area=[]; 
                    sel_area=find(strcmp(available_AREAS,your_area{j})) ;
                  
                    if ~isempty(sel_area)
    
                        if sel_area==numel(available_AREAS)
                            init_area=PAG_columns(sel_area);
                            end_area=max(depth);
                        else
                            init_area=PAG_columns(sel_area);
                            end_area=PAG_columns(sel_area+1);
                        end
    
                        this_area_index = depth > init_area & depth < end_area ;
    
                        if isempty(depth(this_area_index))
                            p_accuracy_all{correc_area_index,sel_behavior}(bf_i,:)=nan(1,wrap_bins);
                        else
                           [ aux_exc , numbers_per_area,index,neto_exc_selected] = estimate_wrapped_column(depth,  init_area, end_area,  wrap_bins, p_acc');
    
                            output_bin = unique(index);
    
                            p_accuracy_all{correc_area_index,sel_behavior}(bf_i,:)=(aux_exc./numbers_per_area)*100;
                            
                            for p = 1:nperm
                                r_ind = randi(numel(aux_exc),1,numel(aux_exc))  ;
                                p_accuracy_all_shuff{correc_area_index,sel_behavior}(bf_i,:,p)=(aux_exc(r_ind)./numbers_per_area(r_ind))*100;
                            end
    
                            all_aux_exc{correc_area_index,sel_behavior}         = all_aux_exc{correc_area_index,sel_behavior} + aux_exc ;
                            all_numbers_per_are{correc_area_index,sel_behavior} = all_numbers_per_are{correc_area_index,sel_behavior} + numbers_per_area ;
                            
                            aux_ids (:,1)= depth(this_area_index & p_acc'); aux_ids (:,2)= ids_gcl(this_area_index & p_acc');
                            ids_exc{correc_area_index,sel_behavior} = [ids_exc{correc_area_index,sel_behavior}; aux_ids];aux_ids=[];
                        end
    
    
                    else
                         p_accuracy_all{correc_area_index,sel_behavior}(bf_i,:)=nan(1,wrap_bins);
                    end
    
                end
            end
            bf_i = bf_i + 1;
        end
    
        %%
        
        perct = p_accuracy_all(:,sel_behavior);
        perct_shuff = p_accuracy_all_shuff(:,sel_behavior);
    
        pvals=[];
        for stitch_i = 1 :numel(your_area)
            pval=[];
            aux = nanmean(p_accuracy_all{stitch_i,sel_behavior});
            aux_shuff = nanmean(p_accuracy_all_shuff{stitch_i,sel_behavior});
    
            for along = 1:numel(aux)
    
                for ps=1:nperm
                    shuffles(ps)=aux_shuff(1,along,ps);
                end
                pval(along) = mean (aux(along) > shuffles);
            end
    
            pvals=[pvals pval];
    
        end
    
        %% Beaware that it's averaging percentages for each session. It already got percentages on each one of them
    
        
    
            counts = p_accuracy_all(:,sel_behavior);
    
        counts_full=[];
        counts_full_no_mean=[];
    
        for stitch_i = 1 :numel(your_area)
    
            counts_full = [counts_full; nanmean(counts{stitch_i})'];
            counts_full_no_mean= [ counts_full_no_mean;counts{stitch_i}];
        
        end
    
          
        counts_full(isnan(counts_full))=0;
        percentages = counts_full;

        all_percentage(:,choose_cell_type) = percentages;

    end

    figure('units','normalized','outerposition',[0 0 0.1 0.7])
    
    full_areas = unique(full_areas);

    yLabel='Percentage of modulated cells';
    
    res_area = [];
    for i = all_choose_cell_type
        res_area(:,i)=boxFilter(all_percentage(:,i)',bxfilterpar);        
    end   

        ar=area(res_area); %ar.FaceColor = [0.1 0.7 0.1];
        
        xlim([0 numel(your_area)*wrap_bins])
        ylim([0 lim_y])

        
        ylabel(yLabel)
            xtic=xticks;
            

        for jj=1:numel(your_area)
    
            xline(jj*wrap_bins,'--','Color','r','LineWidth',3)
            txt = your_area(jj);
            text(jj*wrap_bins-1.5,lim_y-5,txt,'FontSize',14,'Color','red')
        end
    
        hold on


    sgtitle([selfother_title{sel_selfother} ' - ' behaviors_to_check_labels{sel_behavior}])
% end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHUFFLE and plot on top


    camroll(+90)
            set(gca,'YDir','reverse') 
    
   

 end

