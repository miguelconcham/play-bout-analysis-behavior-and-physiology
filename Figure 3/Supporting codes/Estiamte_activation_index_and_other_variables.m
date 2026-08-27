%% 1 DATA LOADING and prealocting variables
cd('\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology')
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Peak trough dyanmics';
data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
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






% 
% behavior_labels = {'play', 'CH', 'POA','POB', 'PWIA','PWIB', 'EV', ...
%     'RE',  'ES', 'CD', 'SN', 'CB', 'GR'};
% 
% base_conditions = {'','play', 'CH', 'POA','POB', 'PWIA','PWIB', 'EV', ...
%     'RE',  'ES', 'CD', 'SN', 'CB', 'GR'};



behavior_labels = {'play', 'CH', 'POA','POB', 'PWIA','PWIB', 'EV', ...
    'RE',  'ES', 'CD', 'SN', 'CB', 'GR', 'PI','BO', 'SC','BI'};

base_conditions = {'','play', 'CH', 'POA','POB', 'PWIA','PWIB', 'EV', ...
    'RE',  'ES', 'CD', 'SN', 'CB', 'GR', 'PI','BO', 'SC','BI'};

roles = {'Self', 'Other'};

for b = 1:numel(base_conditions)
    for r = 1:numel(roles)
        if isempty(base_conditions{b})
            field_name = sprintf('PsthWarped%s', roles{r});
        else
            field_name = sprintf('PsthWarped%s%s', base_conditions{b}, roles{r});
        end
        all_neurons_TD.(field_name) = nan(size(all_neurons_TD,1), 60);
    end
end

animal_labels  = {'self','other'};
shuffle_label  = 'shuffled_';
for a = 1:numel(animal_labels)
    animal = animal_labels{a};

    for b = 1:numel(behavior_labels)
        behavior = behavior_labels{b};

        % Build field name Exited
        field_name = sprintf('Exited_%s_%s', animal, behavior);

        % Initialize the field with NaN
        all_neurons_TD.(field_name) = nan(height(all_neurons_TD), 1);

        % Build field name Inhibited
        field_name = sprintf('Inhibited_%s_%s', animal, behavior);

        % Initialize the field with NaN
        all_neurons_TD.(field_name) = nan(height(all_neurons_TD), 1);


    end
end


all_neurons_TD.PsthOnset  = nan(size(all_neurons_TD,1),400);
all_neurons_TD.PsthOffset  = nan(size(all_neurons_TD,1),400);
all_neurons_TD.PsthOnlyPB  = nan(size(all_neurons_TD,1),400);

psth_map = {};

% Combine your behavior and partner lists
n_partners = 2;
all_labels = [base_conditions, arrayfun(@(pn) sprintf('Partner%d', pn), 1:n_partners, 'UniformOutput', false)];

for b = 1:numel(all_labels)
    base = all_labels{b};

    % --- Handle the 'play' (empty) condition ---
    if isempty(base) || strcmp(base, 'play')
        var_self   = 'all_psth_self_warped';
        var_other  = 'all_psth_other_warped';
        field_self  = 'PsthWarpedSelf';
        field_other = 'PsthWarpedOther';
    else
        % --- Handle behavior or partner-specific cases ---
        var_self   = sprintf('all_psth_self_warped_%s', base);
        var_other  = sprintf('all_psth_other_warped_%s', base);
        field_self  = sprintf('PsthWarped%sSelf', base);
        field_other = sprintf('PsthWarped%sOther', base);
    end

    % Add mappings to the cell array
    psth_map = [psth_map; {field_self, var_self}; {field_other, var_other}];
end


vars_to_load = { ...
    "initial_cluster", "AREAS", "allareas", ...
    "all_psth_zscore", "all_psth_zscore_offset", "all_psth_FR_pblimit", ...
    "all_psth_PB_warped","all_psth_PB_warped_corrected", "good_clusters", "depth_or_Chn", ...
    "all_psth_FR", "all_psth_FR_offset",...
    "all_mean_std","all_psth_shuffled_PB_warped","all_psth_shuffled_PB_warped_corrected"};

% Add the common base variables
vars_to_load = [vars_to_load, ...
    "all_psth_self_warped", "all_psth_other_warped","all_psth_shuffled_self_warped","all_psth_shuffled_other_warped"];

% Add the per-behavior ones dynamically
for b = 2:numel(behavior_labels) % start from 2 to skip 'play' (already added above)
    label = behavior_labels{b};
    vars_to_load = [vars_to_load, ...
        sprintf("all_psth_self_warped_%s", label), ...
        sprintf("all_psth_other_warped_%s", label)...
        sprintf("all_psth_shuffled_self_warped_%s", label), ...
        sprintf("all_psth_shuffled_other_warped_%s", label)];
end
n_partners = 2;
for pn = 1:n_partners
    partner_label = sprintf('Partner%d', pn);

    vars_to_load = [vars_to_load, ...
        sprintf('all_psth_self_warped_%s', partner_label), ...
        sprintf('all_psth_other_warped_%s', partner_label), ...
        sprintf('all_psth_shuffled_self_warped_%s', partner_label), ...
        sprintf('all_psth_shuffled_other_warped_%s', partner_label), ...
        sprintf('all_psth_PB_warped_%s', partner_label), ...
        sprintf('all_psth_shuffled_PB_warped_%s', partner_label) ...
        ];
end




%% 2 DEFINE RAT/PROBE combination

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

time_indxs = 21:40;
ONOFFsets = 2;


%% 3 Select rat to estiamte
thisrat=[ B1D1, B1S3, B2S2, B3D2, B4S2  ];minus_areas=0;
bf_i = 1;

your_area={'DR','VLPAG','LPAG','DLPAG','SupCol'};

playbout_tittle='PlayBout';



behaviors2check={'Pounce_A','CC','Pin','Boxing','Evasion','CB','Pounce_B','Escape','CD', ...
    'Rearing','Grooming','Scratch','Pounce_Ai','Pounce_Bi','Sniffing','Bite'};


%% 4 ACTUALLY load behavior files and  load data
AREAS=[];
depth = [];
full_areas = [];

all_pvalues=[];
all_coefficients= [];

behavior_files = dir('*.txt');
for bf = thisrat

    aux_mydate=behavior_files(bf).name(6:9);
    ani_ID=behavior_files(bf).name(1:4);

    if str2double(behavior_files(bf).name(2))>3
        mydate=['2024' aux_mydate];
    else
        mydate=['2023' aux_mydate];
    end

    path_mati='\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Responses_Matrix\ModelCriterion\';

    load([path_mati 'ResponsesMatrix_PPB_p1andp2_' num2str(length_duration_threshold) 's_' playbout_tittle '_' mydate '_' behavior_files(bf).name(1:4) '_' MorL{bf_i} '.mat'],"initial_cluster","AREAS","RM_areas","all_psth_zscore", ...
        "pre_onset","post_onset","pre_offset","post_offset","all_psth_zscore_offset","all_this_psth","this_psth_shuffled","est_full","pval_full","mod_wrap", ...
        "good_clusters","depth_or_Chn","all_psth_FR","all_psth_FR_offset")

    path_migue='\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Responses_Matrix\ModelCriterion_Onset_Miguel\';
    load([path_migue 'ResponsesMatrix_PPB_p1andp2_' num2str(length_duration_threshold) 's_' playbout_tittle '_' mydate '_' behavior_files(bf).name(1:4) '_' MorL{bf_i} '.mat'] , vars_to_load{:})



    depth = depth_or_Chn;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Inhibition and Excitation based on LOCAL difference with shuffle %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    aux                     =        all_this_psth; %% PSTH warped

    full_zsc_warp_shuffled = this_psth_shuffled; %% PSTH warped SHUFFLED
    significant_indexes_exc =nan(size(full_zsc_warp_shuffled,2),1);
    significant_indexes_inh =nan(size(full_zsc_warp_shuffled,2),1);
    for j=1:size(full_zsc_warp_shuffled,2)


        sh_FR_after_onset=mean(full_zsc_warp_shuffled{j}(:,time_indxs),2,"omitmissing");

        FR_after_onset=mean(aux(:,time_indxs),2,"omitmissing");
        p_exc = mean(FR_after_onset(j) < sh_FR_after_onset); % p-value for being greater
        p_inh = mean(FR_after_onset(j) > sh_FR_after_onset); % p-value for being smaller

        significant_indexes_exc(j) = (p_exc < alpha);
        significant_indexes_inh(j) = (p_inh < alpha);

    end



    Significance = struct();
    for a = {'self','other'}
        for b = behavior_labels
            Significance.(a{1}).(b{1}).exc = nan(numel(good_clusters),1);
            Significance.(a{1}).(b{1}).inh = nan(numel(good_clusters),1);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% --- now cycle trough remeining psth and shuffles %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    animal_labels  = {'self','other'};
    shuffle_label  = 'shuffled_';  % only compare real vs shuffled


    animal_labels = {'self', 'other'};
    shuffle_label = 'shuffled_';
    partner_labels = arrayfun(@(pn) sprintf('Partner%d', pn), 1:n_partners, 'UniformOutput', false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Combine behaviors and partners into one list %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_labels = [behavior_labels, partner_labels];

    for a = 1:numel(animal_labels)
        for b = 1:numel(all_labels)
            label = all_labels{b};

            % --- Build PSTH variable names dynamically ---
            if strcmp(label, 'play')
                real_var     = sprintf('all_psth_%s_warped', animal_labels{a});
                shuffled_var = sprintf('all_psth_shuffled_%s_warped', animal_labels{a});
            else
                real_var     = sprintf('all_psth_%s_warped_%s', animal_labels{a}, label);
                shuffled_var = sprintf('all_psth_shuffled_%s_warped_%s', animal_labels{a}, label);
            end

            % --- Skip missing variables ---
            if ~exist(real_var, 'var') || ~exist(shuffled_var, 'var')
                fprintf('Skipping missing variable: %s or %s\n', real_var, shuffled_var);
                continue;
            end

            % --- Extract PSTH matrices ---
            aux = eval(real_var);
            psth_warp_shuffled = eval(shuffled_var);

            % --- Compute significance for each neuron/cluster ---
            if ~isempty(psth_warp_shuffled) && numel(psth_warp_shuffled) >= 1 && size(psth_warp_shuffled{1}, 2) > 1
                n_clusters = size(aux, 1);
                Significance.(animal_labels{a}).(label).exc = nan(n_clusters, 1);
                Significance.(animal_labels{a}).(label).inh = nan(n_clusters, 1);

                for j = 1:n_clusters
                    sh_FR_after_onset = mean(psth_warp_shuffled{j}(:, time_indxs), 2, "omitmissing");
                    FR_after_onset    = mean(aux(j, time_indxs), 2, "omitmissing");

                    p_exc = mean(FR_after_onset < sh_FR_after_onset); % Excitation test
                    p_inh = mean(FR_after_onset > sh_FR_after_onset); % Inhibition test

                    Significance.(animal_labels{a}).(label).exc(j) = (p_exc < alpha);
                    Significance.(animal_labels{a}).(label).inh(j) = (p_inh < alpha);
                end
            else
                fprintf('No shuffled data available for %s (%s)\n', animal_labels{a}, label);
            end
        end
    end


    session_name = behavior_files(bf).name;
    session_name =  strsplit(session_name, '.');
    session_name = session_name{1};
    pctg_cell_matching= 100*sum(ismember(all_neurons_TD.session,session_name) & ismember(all_neurons_TD.cluster_id, good_clusters))/sum(ismember(all_neurons_TD.session,session_name));
    disp([session_name , ' ',  num2str(pctg_cell_matching), '% match'])



    significant_indexes_exc = logical(significant_indexes_exc);
    significant_indexes_inh = logical(significant_indexes_inh);
    all_neurons_TD_MASK = (ismember(all_neurons_TD.session,session_name) & ismember(all_neurons_TD.cluster_id, good_clusters));
    loaded_data_mask = (ismember(good_clusters,all_neurons_TD.cluster_id(ismember(all_neurons_TD.session,session_name))));
    td_idx = find(all_neurons_TD_MASK);
    gc_idx = find(loaded_data_mask);
    td_ids = all_neurons_TD.cluster_id(td_idx);
    gc_ids = good_clusters(gc_idx);
    [tf, loc] = ismember(td_ids, gc_ids);     

     

   if all(all_neurons_TD.cluster_id(td_idx(tf)) == good_clusters(gc_idx(loc(tf)))')
       disp('Reassignment correct')
   else
       disp(['Warning, missasignment in ' session_name])
   end
    
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% NOW ADD DATA TO TABLE %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% 1 Exited %%%%
    all_neurons_TD.Exited( td_idx(tf)) =  significant_indexes_exc(gc_idx(loc(tf)));
    %%%% 2 Inhibited %%%%
    all_neurons_TD.Inhibited( td_idx(tf)) = ...
        significant_indexes_inh(gc_idx(loc(tf)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% all remaining inhibited an exited %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for a = 1:numel(animal_labels)
        animal = animal_labels{a};

        for b = 1:numel(behavior_labels)
            behavior = behavior_labels{b};

            % Build field name
            field_name = sprintf('Exited_%s_%s', animal, behavior);
            % Assign significance results
            all_neurons_TD.(field_name)(td_idx(tf)) = ...
                Significance.(animal).(behavior).exc(gc_idx(loc(tf)));

            % Build field name
            field_name = sprintf('Inhibited_%s_%s', animal, behavior);
            % Assign significance results
            all_neurons_TD.(field_name)(td_idx(tf)) = ...
                Significance.(animal).(behavior).inh(gc_idx(loc(tf)));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%
    %%%% 3 PsthWarped %%%%
    %%%%%%%%%%%%%%%%%%%%%%

    all_neurons_TD.PsthWarped( td_idx(tf),:) = ...
        all_psth_PB_warped(gc_idx(loc(tf)),:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% 4  ALL REMAINNG PSTH WARPED %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    all_neurons_TD.PsthOnset( td_idx(tf),:) = ...
        all_psth_zscore(gc_idx(loc(tf)),:);

    for i = 1:size(psth_map,1)

        field_name  = psth_map{i,1};
        source_name = psth_map{i,2};

        if exist(source_name, 'var')
            data = eval(source_name);
            all_neurons_TD.(field_name)(td_idx(tf),:) = ...
                data(gc_idx(loc(tf)),:);
        else
            warning('Variable %s not found in workspace, skipping %s.', source_name, field_name);
        end
    end

    bf_i = bf_i+1; %now next animal
end


%% 5 Estimating population respones (logit) for all behaviors and all areas
alpha_level         = 0.01;
non_entrained_lvl   = 0.1;
smooth_window       = 5;
warped_time         = (((1:60)/20)*5) - 5 ;
non_wraped_time     = (1:50)/5 - 5;
baseline            = [-Inf 0];

onset_time          = [0 0];
offset_time         = [5 5];
freq2use            = 'DeltaEntireSession'; %options:   ThetaEntireSession DeltaEntireSession
area_list           = {'SupCol' 'DLPAG'	'LPAG'	'VLPAG' 'DR' };
cell_type_names     = {'Peak-locked','Trough-locked','Unlocked'};
plot_angle          = true;
response_type       = 'All';

n_perm                      = 10000;
all_psth_together           = [];
all_mean_responses_together = [];
all_indexes                 = [];
psth_list                   = [{'PsthWarped'}';psth_map(:,1)];
psth_list([  4 5]) = [];
all_comp                    = cell(numel(psth_list),10);
psth_n = 1;
all_psth_cell               = cell(3,numel(area_list));
response_time               = [0 5];
all_activation_order        = cell(numel(psth_list),numel(area_list),4);


session_list                = unique(all_neurons_TD.session(~isnan(all_neurons_TD.Exited)));
session_index               =~ismember(all_neurons_TD.session,session_list([]) );


for psth2use_cell = psth_list'

    psth2use            =  psth2use_cell{1};

    if strcmp(psth2use, 'PsthOnlyPB')
        time2use = mig_edges_centers;
    elseif strcmp(psth2use, 'PsthOnset') || strcmp(psth2use, 'PsthOffset')
        time2use= non_wraped_time;
    else
        time2use = warped_time;
    end
    baseline_index      = time2use>=baseline(1) & time2use<=baseline(2);
    response_time_index = time2use>=response_time(1) & time2use<=response_time(2);
    all_areas_together  = [];
    indexes             = [];


    all_areas_mean_responses    = nan(numel(area_list),3);
    weighted_respones           = nan(numel(area_list),3);
    all_areas_p_val             = nan(numel(area_list),3);
    all_areas_p_val_pctg        = nan(numel(area_list),3);
    all_areas_multcomp          = nan(numel(area_list),3);
    all_areas_kw                = nan(numel(area_list),3);
    absolute_mdoulation         = nan(numel(area_list),3);
    activation_consistency      = nan(numel(area_list),3);
    activation_consistency_p    = nan(numel(area_list),3);
    pop_activation              = nan(numel(area_list),3);

    for an=1:numel(area_list)

        area_index      = ismember(all_neurons_TD.area, area_list{an});
        angletype_peak  = ~(all_neurons_TD.(freq2use).PreferedAngle>=pi/2 | all_neurons_TD.(freq2use).PreferedAngle<-pi/2);
        no_nan          = ~isnan(all_neurons_TD.Inhibited);
        if strcmp(response_type, 'All')
            response_index = true(size(all_neurons_TD, 1),1);
        else
            response_index  = all_neurons_TD.(response_type)==1;
        end

        entreined           = all_neurons_TD.(freq2use).PPCPval<=alpha_level ;
        not_entrained       = all_neurons_TD.(freq2use).PPCPval>non_entrained_lvl;
        peak_index          = session_index &  area_index & no_nan & response_index & entreined & angletype_peak  ;
        trough_index        = session_index &  area_index & no_nan & response_index & entreined & ~angletype_peak ;
        nonentrained_index  = session_index &  area_index & no_nan & response_index & not_entrained ;

        all_indexes                         = zeros(size(nonentrained_index));
        all_indexes(peak_index)             = 1;
        all_indexes(trough_index)           = 2;
        all_indexes(nonentrained_index)     = 3;
        all_psth_indexes                    = all_indexes(area_index & no_nan & response_index);
        all_activation_order{psth_n, an,1}  = all_psth_indexes;
        all_activation_order{psth_n, an,3}  = all_neurons_TD.session(area_index & no_nan & response_index);
        
        %Estimate peak poplation responses and activation index
        peak_psth                           = all_neurons_TD.(psth2use)(peak_index, :);
        for j=1:size(peak_psth,1)
            peak_psth(j,:) = smooth( peak_psth(j,:),smooth_window);
            if  std( peak_psth(j,baseline_index), 'omitmissing')>0.01
                peak_psth(j,:) = ( peak_psth(j,:) - mean( peak_psth(j,baseline_index), 'omitmissing'))/ std( peak_psth(j,baseline_index), 'omitmissing');
            else
                peak_psth(j,:) = ( peak_psth(j,:) - mean( peak_psth(j,:), 'omitmissing'))/ std( peak_psth(j,:), 'omitmissing');
            end
        end
        %Estimate trough poplation responses and activation index
        trough_psth  =  all_neurons_TD.(psth2use)(trough_index, :);
        for j=1:size(trough_psth,1)
            trough_psth(j,:) = smooth( trough_psth(j,:),smooth_window);
            if  std( trough_psth(j,baseline_index), 'omitmissing')>0.01
                trough_psth(j,:) = ( trough_psth(j,:) - mean( trough_psth(j,baseline_index), 'omitmissing'))/ std( trough_psth(j,baseline_index), 'omitmissing');
            else
                trough_psth(j,:) = ( trough_psth(j,:) - mean( trough_psth(j,:), 'omitmissing'))/ std( trough_psth(j,:), 'omitmissing');
            end
        end
        %Estimate non-entrained poplation responses and activation index
        nonentrained_psth  =  all_neurons_TD.(psth2use)(nonentrained_index, :);
        for j=1:size(nonentrained_psth,1)
            nonentrained_psth(j,:) = smooth( nonentrained_psth(j,:),smooth_window);
            if  std( nonentrained_psth(j,baseline_index), 'omitmissing')>0.01
                nonentrained_psth(j,:) = ( nonentrained_psth(j,:) - mean( nonentrained_psth(j,baseline_index), 'omitmissing'))/ std( nonentrained_psth(j,baseline_index), 'omitmissing');
            else
                nonentrained_psth(j,:) = ( nonentrained_psth(j,:) - mean( nonentrained_psth(j,:), 'omitmissing'))/ std( nonentrained_psth(j,:), 'omitmissing');
            end
        end
        all_psth = all_neurons_TD.(psth2use)(area_index & no_nan & response_index, :);
        for j=1:size(all_psth,1)
            all_psth(j,:) = smooth( all_psth(j,:),smooth_window);
            if  std( all_psth(j,baseline_index), 'omitmissing')>0.01
                all_psth(j,:) = ( all_psth(j,:) - mean( all_psth(j,baseline_index), 'omitmissing'))/ std( all_psth(j,baseline_index), 'omitmissing');
            else
                all_psth(j,:) = ( all_psth(j,:) - mean( all_psth(j,:), 'omitmissing'))/ std( all_psth(j,:), 'omitmissing');
            end
        end
      

        [~,all_activation_order{psth_n, an,2}] = sort(all_psth);
        ranks = 1:size(all_psth,1);
        for ti=1:size(all_psth,2)
            all_activation_order{psth_n, an,2}(:,ti) =size(all_psth,1)- ranks(all_activation_order{psth_n, an,2}(:,ti))+1;
        end
        all_activation_order{psth_n, an,4} = all_psth;

        null_psth_peak = nan(n_perm, size(peak_psth,2));
        for pn = 1:n_perm

            sub_selection = all_psth(randperm(size(all_psth,1),size(peak_psth,1)),:);

            null_psth_peak(pn,:) = median(sub_selection, 'omitmissing');
        end
        % peak_pctl_activation = 100*mean(null_psth>mean(peak_psth, 'omitmissing'));
        peak_pctl_activation = mean(null_psth_peak<median(peak_psth, 'omitmissing'));
        peak_pctl_activation(peak_pctl_activation==0) = 1/n_perm;
        peak_pctl_activation(peak_pctl_activation==1) = (n_perm-1)/n_perm;
        peak_pctl_activation = log10(peak_pctl_activation./(1-peak_pctl_activation));
        peak_pctl_inhibition = mean(null_psth_peak>mean(peak_psth, 'omitmissing'));

        null_psth_trough = nan(n_perm, size(trough_psth,2));
        for pn = 1:n_perm

            sub_selection = all_psth(randperm(size(all_psth,1),size(trough_psth,1)),:);

            null_psth_trough(pn,:) = median(sub_selection, 'omitmissing');
        end
        % trough_pctl_activation = 100*mean(null_psth>mean(trough_psth, 'omitmissing'));
        trough_pctl_activation = mean(null_psth_trough<median(trough_psth, 'omitmissing'));

        trough_pctl_activation(trough_pctl_activation==0) = 1/n_perm;
        trough_pctl_activation(trough_pctl_activation==1) = (n_perm-1)/n_perm;
        trough_pctl_activation = log10(trough_pctl_activation./(1-trough_pctl_activation));
        trough_pctl_inhibition = mean(null_psth_trough>median(trough_psth, 'omitmissing'));

        null_psth_nonentrained = nan(n_perm, size(nonentrained_psth,2));
        for pn = 1:n_perm

            sub_selection = all_psth(randperm(size(all_psth,1),size(nonentrained_psth,1)),:);

            null_psth_nonentrained(pn,:) = median(sub_selection, 'omitmissing');
        end
        % nonentrained_pctl_activation = 100*median(null_psth>median(nonentrained_psth, 'omitmissing'));
        nonentrained_pctl_activation = mean(null_psth_nonentrained<median(nonentrained_psth, 'omitmissing'));
        nonentrained_pctl_activation(nonentrained_pctl_activation==0) = 1/n_perm;
        nonentrained_pctl_activation(nonentrained_pctl_activation==1) = (n_perm-1)/n_perm;
        nonentrained_pctl_activation = log10(nonentrained_pctl_activation./(1-nonentrained_pctl_activation));
        nonentrained_pctl_inhibition = mean(null_psth_nonentrained>median(nonentrained_psth, 'omitmissing'));


        all_psth_cell{1,an}= [all_psth_cell{1,an};peak_pctl_activation];
        all_psth_cell{2,an}= [all_psth_cell{2,an};trough_pctl_activation];
        all_psth_cell{3,an}= [all_psth_cell{3,an};nonentrained_pctl_activation];


        peak_pop = median(peak_pctl_activation(response_time_index));
        trough_pop = median(trough_pctl_activation(response_time_index));
        nonentrained_pop = median(nonentrained_pctl_activation(response_time_index));

        pop_activation(an,:) = [peak_pop trough_pop nonentrained_pop];

        rpv_g1 = mean(peak_psth(:,response_time_index),2);
        if size(peak_psth,1)>1
            h = ttest(peak_psth);
        else
            h = true(size(peak_psth));
        end
        all_areas_p_val_pctg(an,1)=sum(h(response_time_index))/numel(h(response_time_index));
        weighted_respones(an,1) = mean(mean(peak_psth(:,h==1 & response_time_index ), 'omitmissing'), 'omitmissing');
        absolute_mdoulation(an,1) = mean(abs(mean(peak_psth(:,h==1 & response_time_index ), 'omitmissing')), 'omitmissing');

        rpv_g2 = mean(trough_psth(:,response_time_index),2);
       
        if size(trough_psth,1)>1
            h = ttest(trough_psth);
        else
            h = true(size(trough_psth));
        end
        all_areas_p_val_pctg(an,2)=sum(h(response_time_index))/numel(h(response_time_index));
        weighted_respones(an,2) = mean(mean(trough_psth(:,h==1 & response_time_index ), 'omitmissing'), 'omitmissing');
        absolute_mdoulation(an,2) = mean(abs(mean(trough_psth(:,h==1 & response_time_index ), 'omitmissing')), 'omitmissing');

        rpv_g3 = mean(nonentrained_psth(:,response_time_index),2);
        if size(nonentrained_psth,1)>1
            h = ttest(nonentrained_psth);
        else
            h = true(size(nonentrained_psth));
        end
        all_areas_p_val_pctg(an,3)=sum(h(response_time_index))/numel(h(response_time_index));
        weighted_respones(an,3) = mean( mean(nonentrained_psth(:,h==1 & response_time_index ), 'omitmissing'), 'omitmissing');
        absolute_mdoulation(an,3) = mean(abs(mean(nonentrained_psth(:,h==1 & response_time_index ), 'omitmissing')), 'omitmissing');


        all_areas_mean_responses(an,:)      = [mean(rpv_g1, 'omitmissing') mean(rpv_g2, 'omitmissing') mean(rpv_g3, 'omitmissing')];
        if sum(~isnan(rpv_g1)) >3
            p1=signrank(rpv_g1);
        else
            p1 = NaN;
        end

        if sum(~isnan(rpv_g2)) >3
            p2=signrank(rpv_g2);
        else
            p2 = NaN;
        end
        if sum(~isnan(rpv_g3)) >3
            p3=signrank(rpv_g3);
        else
            p3 = NaN;
        end

        all_areas_p_val(an,:)               = [p1 p2 p3];
        [ p, ~, stats] = kruskalwallis([rpv_g1;rpv_g2;rpv_g3],[rpv_g1*0;(rpv_g2*0 + 1);(rpv_g3*0 + 2)], 'off');
        if numel(stats.gnames) > 1
            comp_output = multcompare(stats, 'Display','off');
        else
            comp_output = NaN(1,3);
        end
        all_areas_kw(an,:)                  = p;
        all_areas_multcomp(an,:)             = comp_output(:,3);
      
    end

    all_comp(psth_n,:) = {all_areas_mean_responses,all_areas_p_val,all_areas_multcomp,all_areas_kw,all_areas_p_val_pctg,weighted_respones,absolute_mdoulation,activation_consistency,activation_consistency_p,pop_activation};
    psth_n = psth_n+1;


end
play_song([],[],[])

% 6 Saving — checkpoint for plotting 

checkpoint_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\activation index';
checkpoint_file = [checkpoint_folder, '\activation_index_section6.mat'];


disp('Saving section-6 checkpoint')
save(checkpoint_file, ...
    'all_activation_order', ...   % per-neuron PSTHs, cell-type index, session (sect. 7, IRI plots)
    'all_psth_cell', ...          % logit activation-index traces (sect. 8+)
    'psth_list', ...              % behavior/condition names for y-labels
    'area_list', ...
    'cell_type_names', ...
    'warped_time', ...
    '-v7.3');
disp(['Saved: ', checkpoint_file])



%% 7 SAVE DATA FOR COMPARISION of activation Index with delta power per behavior (former Supplementary figure)
saved_data = struct;

% play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25 ];
% non_play_behaviors = [10 11 12 13 16 17 22 23 26 27  32 33 34 35];

play_behaviors     = [4 6 8 14 18 20 24 5 7 9 15 19 21 25];
non_play_behaviors = [10 12 16 22 26 34 11 13 17 23 27 35];
time_analysis = (1:size(all_psth_cell{1,1},2))/4 - 5;

for an = 1:5
    for j = 1:3

        play_matrix     = all_psth_cell{j,an}(play_behaviors,:);
        non_play_matrix = all_psth_cell{j,an}(non_play_behaviors,:);

        % --- choose ONE consistent stacking order ---
        stacked_matrix = [play_matrix; non_play_matrix];

        % --- matching labels (must follow SAME order) ---
        stacked_labels = psth_list([play_behaviors, non_play_behaviors]);

        % --- store ---
        saved_data(an,j).stacked_matrix = stacked_matrix;
        saved_data(an,j).labels         = stacked_labels;
        saved_data(an,j).time           = time_analysis;

        saved_data(an,j).area           = area_list{an};
        saved_data(an,j).cell_type      = cell_type_names{j};

        saved_data(an,j).play_behaviors     = play_behaviors;
        saved_data(an,j).non_play_behaviors = non_play_behaviors;

    end
end

% Keep both locations for now: phase locking data folder and the repo checkpoint folder.
% save([saving_folder, '\stacked_psth_by_area_celltype.mat'], 'saved_data')
if ~exist(checkpoint_folder, 'dir'), mkdir(checkpoint_folder); end
save([checkpoint_folder, '\stacked_psth_by_area_celltype.mat'], 'saved_data')
disp(['Also saved: ', checkpoint_folder, '\stacked_psth_by_area_celltype.mat'])


%% 8 SAVE checkpoint for example plot (section 9)
% Needs the filled all_neurons_TD table from sections 1–4 (PSTHs, Exited /
% Inhibited, delta/theta locking). After this save, section 9 can be run
% on its own from a fresh MATLAB session.
checkpoint_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\activation index';

example_checkpoint = [checkpoint_folder, '\activation_index_example_plot.mat'];
vars_example = {'all_neurons_TD', 'psth_map'};
if exist('mig_edges_centers', 'var')
    vars_example{end+1} = 'mig_edges_centers';
end
disp('Saving example-plot checkpoint')
save(example_checkpoint, vars_example{:}, '-v7.3');
disp(['Saved: ', example_checkpoint])

%% Saving variables for ploting sections A1 to A7
% Checkpoint for Figure 3/Fig A1-A7 behavior modulation.m
checkpoint_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\activation index';
if ~exist(checkpoint_folder, 'dir')
    mkdir(checkpoint_folder);
end
a17_file = [checkpoint_folder, '\activation_index_A1_A7.mat'];

if ~exist('time2use', 'var') || isempty(time2use)
    time2use = warped_time;
end

% Short behavior names as a row (same as Play_map A1: labels(play_index_order)').
labels = cell(1, numel(psth_list));
labels{1} = 'PlayBout';
for i = 2:numel(psth_list)
    nm = psth_list{i};
    if numel(nm) >= 11
        labels{i} = nm(11:end);
    else
        labels{i} = nm;
    end
end

disp('Saving A1–A7 checkpoint')
save(a17_file, ...
    'all_activation_order', ...
    'area_list', ...
    'psth_list', ...
    'labels', ...
    'time2use', ...
    'warped_time', ...
    'n_perm', ...
    '-v7.3');
disp(['Saved: ', a17_file])


