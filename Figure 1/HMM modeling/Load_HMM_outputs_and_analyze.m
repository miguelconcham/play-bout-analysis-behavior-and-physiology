%% 1 SET PARAMS
function [] = Load_HMM_outputs_and_analyze(dir_output,animal_code)

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Codes repository\Data';

synch_folder = [data_root, '\Synch data'];
traking_folder = [data_root, '\Traking backups'];
behavior_folder = [data_root, '\Behavior backups'];
call_folder = [data_root, '\CallDetectionBackup'];
hmm_raw_data = [data_root, '\HMM data\HMM raw data'];
animal_code_params =strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
date                = animal_code_params{2};
repeated_animal     = animal_code_params{3}; %Single2, Dual Single
% animal_code         =  [animal_batch ' ' date ' ' repeated_animal]
file2save = animal_code;
play_state      = 0;
non_play_state  = 1;

load([synch_folder, '\',animal_code, '\synch_model_video2audio.mat'], 'synch_model_video2audio')


play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};
beh_bin = 0.01; %for hmm stimate
conv_length = 1; % convuolution lenbgth for smothing all data
list_of_trakings = dir([traking_folder, '\*',animal_code, '*']);

bin_size = 0.01;
psth_bin = [-20 20];
psth_limits = max(abs(psth_bin));


restrict2Partnerssession = true;
%% 2 Load data
disp('Loading Data')
Call_file = dir([call_folder, '\*',animal_code, '*.xlsx']);   %load call data
Behavior_file = dir([behavior_folder, '\*',animal_code, '*.txt']); %load behavior data


CallStats = readtable([call_folder, '\',Call_file.name]);
CallStats.Properties.VariableNames = cellfun(@(x) strrep(x, '_', ''),CallStats.Properties.VariableNames, 'UniformOutput',false );
Behavior =   readtable([behavior_folder,'\',Behavior_file.name]);
Behavior(:,2) = [];
Behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
Behavior(ismember(Behavior.Type,'Partners session'),:)




Behavior.Type2 = Behavior.Type;
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'}; %% Merging behaviors to Type2
Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
Behavior.Type2(ismember( Behavior.Type2,'')) = {'Other'};
Behavior(ismember(Behavior.Animal, 'Reversal'),:) = [];
animal_types = unique(Behavior.Animal);
animal_types(ismember(animal_types,'Session_structure'))=[];
Behavior.Start = predict(synch_model_video2audio, Behavior.Start);
Behavior.End = predict(synch_model_video2audio, Behavior.End);


min_time2analysis = Behavior.Start(ismember(Behavior.Type,'Partners session'))  ;
max_time2analysis = Behavior.End(ismember(Behavior.Type,'Partners session'))  ;
ListOfPartners = [min_time2analysis max_time2analysis];
[~, time_ordered] = sort(ListOfPartners(:,1));
ListOfPartners = ListOfPartners(time_ordered,:);

for pt = 1:size(ListOfPartners,1) %Create separate file sfore ach partner
    clear traking_structure
    load([traking_folder, '\',list_of_trakings(pt).name], 'traking_structure')
    traking_time = predict(synch_model_video2audio,(traking_structure.frames2stract/30)');
    if restrict2Partnerssession
        index2inlcude = traking_time>=ListOfPartners(pt,1) &  traking_time<=ListOfPartners(pt,2);
    else
        index2inlcude = true(size(traking_structure.time));
    end

    traking_time = traking_time(index2inlcude);



    partner_names = animal_types;
    partner_names(ismember(animal_types, repeated_animal)) = [];

    config.Behavior = Behavior;
    config.repeated_animal =repeated_animal;
    config.animal_types =animal_types        ;
    config.play_behaviors =play_behaviors      ;
    config.beh_bin =beh_bin             ;
    config.conv_length =conv_length;
    config.behavior_window = 0;


    [play_bouts_table] = play_bout(config);


    %% 3 EStimate behavioral varialbes
    disp('Estimating behavioral varialbes')
    %
    % triking_time =(triking_time-synch_model_video2audio.Coefficients.Estimate(1))/synch_model_video2audio.Coefficients.Estimate(2)
    % predict(synch_model_video2audio, triking_time')';


    % triking_time =predict(synch_model_video2audio, triking_time(index2inlcude)')';
    animal_pos      = traking_structure.animal_pos(index2inlcude,:);
    animal_pos(:,1) = smoothdata(animal_pos(:,1), 'loess',5);
    animal_pos(:,2) = smoothdata(animal_pos(:,2), 'loess',5);

    partner_pos = traking_structure.partner_pos(index2inlcude,:);
    partner_pos(:,1) = smoothdata(partner_pos(:,1), 'loess',5);
    partner_pos(:,2) = smoothdata(partner_pos(:,2), 'loess',5);

    inner_border = traking_structure.inner_border;
    all_external_border = traking_structure.outher_border;
    arena_center = mean(traking_structure.inner_border);
    all_inner_border = [inner_border;inner_border(1,:)] - arena_center;

    all_external_border = [all_external_border;all_external_border(1,:)] - arena_center;

    time_range  = [min(traking_time) max(traking_time)];
    binned_time = (beh_bin*floor(time_range(1)/beh_bin)):beh_bin:(beh_bin*ceil(time_range(2)/beh_bin));
    % binned_video_time = binned_time;



    animal_pos_x = interp1(traking_time,animal_pos(:,1), binned_time, 'cubic');
    animal_pos_y = interp1(traking_time,animal_pos(:,2), binned_time, 'cubic');

    animal_pos          = [animal_pos_x', animal_pos_y']-arena_center;
    animal_angle        = cart2pol(animal_pos(:,1), animal_pos(:,2));
    animal_angle_speed  = angdiff(animal_angle);
    animal_angle_acc    = angdiff(animal_angle_speed);

    animal_velocity     = diff(animal_pos);
    animal_acc          = diff(animal_velocity);
    animal_acc          = sqrt(sum(animal_acc.*animal_acc,2));
    animal_speed        = sqrt(sum(animal_velocity.*animal_velocity,2));


    partner_pos_x       = interp1(traking_time,partner_pos(:,1), binned_time);
    partner_pos_y       = interp1(traking_time,partner_pos(:,2), binned_time);
    partner_pos         = [partner_pos_x', partner_pos_y']-arena_center;
    partner_angle       = cart2pol(partner_pos(:,1), partner_pos(:,2));
    partner_angle_speed = angdiff(partner_angle);
    partner_angle_acc   = angdiff(partner_angle_speed);


    partner_velocity    = diff(partner_pos);
    partner_acc         =  diff(partner_velocity);
    partner_speed       = sqrt(sum(partner_velocity.*partner_velocity,2));
    partner_acc         = sqrt(sum(partner_acc.*partner_acc,2));

    relative_pos        = animal_pos-partner_pos;
    relative_velocity   = animal_velocity-partner_velocity;
    relative_distance   = sqrt(sum(relative_pos.*relative_pos,2));
    relative_angle      = cart2pol(relative_pos(:,1), relative_pos(:,2));
    relative_angle_speed=  cart2pol(relative_velocity(:,1), relative_velocity(:,2));
    relative_angle_acc  =  angdiff(relative_angle_speed);

    relative_speed      = diff(relative_distance);
    relative_acc        = diff(relative_speed);
    relative_speed      = sqrt(sum(relative_speed.*relative_speed,2));
    relative_acc        = sqrt(sum(relative_acc.*relative_acc,2));

    animal_spatial_properties   = [  animal_speed(1:end-1)  abs(animal_angle_speed(1:end-1)) abs(animal_angle_acc) animal_acc];
    partner_spatial_properties  = [  partner_speed(1:end-1)  abs(partner_angle_speed(1:end-1)) abs(partner_angle_acc) partner_acc];
    relative_spatial_properties = [relative_distance(1:end-2,:)   relative_speed(1:end-1)  abs(relative_angle_speed(1:end-1)) abs(relative_angle_acc) relative_acc];
    all_spatial_properties      = [animal_spatial_properties,partner_spatial_properties,relative_spatial_properties];
    spatial_property_names      = {  'AnimalSpeed','AnimalAngleSpeed','AnimalAngleAcc','AnimalAcc',...
        'PartnerSpeed','PartnerAngleSpeed','PartnerAngleAcc','PartnerAcc',...
        'RelativeDistance','RelativeSpeed','RelativeAngleSpeed','RelativeAngleAcc', 'RelativeAcc'};
    % call_prop_list = {'PrincipalFrequencykHz', 'SlopekHzs', 'Sinuosity', 'DeltaFreqkHz', 'FrequencyStandardDeviationkHz'};
    call_prop_list              = {};
    call_property_names         = ['NumCalls',call_prop_list];
    % call_property_names = [];
    ALL_VARIABLE_NAMES          = [spatial_property_names,call_property_names];
    %
    call_properties_shifted     = zeros(numel(binned_time)-1,numel(call_prop_list)+1);
    call_properties             = call_properties_shifted;
    % call_properties = [];
    for j=1:numel(binned_time)-1
        is_ther_call = CallStats.BeginTimes<=binned_time(j) &  CallStats.EndTimes>=binned_time(j+1);
        call_properties(j,1) = any(is_ther_call);
        call_properties_shifted(j,1) = any(is_ther_call);
        if sum(is_ther_call)==1
            call_properties(j,2:end) = (CallStats{is_ther_call,call_prop_list});
            call_properties_shifted(j,2:end) = (CallStats{is_ther_call,call_prop_list} -min( CallStats{:,call_prop_list}))/range( CallStats{:,call_prop_list});
        elseif sum(is_ther_call)>1
            call_properties(j,2:end) = (mean(CallStats{is_ther_call,call_prop_list})-min(CallStats{:,call_prop_list}))/range( CallStats{:,call_prop_list});
            call_properties_shifted(j,2:end) = mean(CallStats{is_ther_call,call_prop_list});
        end
    end

    call_properties_shifted(isnan(call_properties_shifted)) = 0;
    all_spatial_properties(isnan(all_spatial_properties))   = 0;


    all_properties      = [all_spatial_properties,call_properties(1:end-1,:)];
    conv_all_properties = all_properties;

    conv_length         = .5;
    conv_time           = linspace(-conv_length, conv_length, 1 +(conv_length/beh_bin));
    conv_fun            = normpdf(conv_time);
    conv_fun            = conv_fun/sum(conv_fun);
    conv_call_properties_non_shifted  = call_properties(1:end-1,:);

    for j=1:size(conv_all_properties,2)

        aux_conv = conv(conv_all_properties(:,j), conv_fun);
        conv_all_properties(:,j) = aux_conv(round(.5*(conv_length/beh_bin) +1):(end-round(.5*(conv_length/beh_bin))));

        if j<=size(conv_call_properties_non_shifted,2)
            aux_conv = conv(conv_call_properties_non_shifted(:,j), conv_fun);
            conv_call_properties_non_shifted(:,j) = aux_conv(round(.5*(conv_length/beh_bin) +1):(end-round(.5*(conv_length/beh_bin))));
        end
    end

    zscored_properties = conv_all_properties;
    for j=1:size(conv_all_properties,2)
        zscored_properties(:,j) = zscore(zscored_properties(:,j));
    end

    %% 4 estimate playbout time (iportant to obtain confussion matrix: true positive, tc
    disp('Estimating playbpout times')
    play_bout_time = any(binned_time>=play_bouts_table(:,1) & binned_time<=play_bouts_table(:,2)); %this variable is important
    %% 5 here i create the input array to run HMM in Pytohon. Classification matrix is importnat to obtain confusion matrtix
    disp('Creating unput array for HMM and saving')
    % prepar_coeff =call_score(:,[1 2]);


    classification_matrix  = [zscored_properties,play_bout_time(1:end-2)'];


    elemets_classified = ~any(isnan(classification_matrix),2);
    adjusted_binned_time = binned_time(~any(isnan(classification_matrix),2));
    classification_matrix(any(isnan(classification_matrix),2),:) = [];
    classification_matrix4HMM = classification_matrix(:, 1:end-1);



    % un-comment next line only to creatine matrix again.
    % writeNPY(classification_matrix4HMM, [file2save,'.npy'])

    %% 6 Here I load the HMM output from python
    disp('Loading HMM states from python output')
    variable_name2states =  [hmm_raw_data, '\',file2save, ' P', num2str(pt),'_states_K2.npy'];
    variable_name3states =  [hmm_raw_data, '\',file2save, ' P', num2str(pt),'_states_K3.npy'];
    time_limit = (adjusted_binned_time>= ListOfPartners(pt,1) & adjusted_binned_time<=  ListOfPartners(pt,2) );
    % time_limit = (adjusted_binned_time>= min(min(play_bouts_table)) & adjusted_binned_time<=   max(max(play_bouts_table)) );

    hmm_states =  readNPY(variable_name2states);
    hmm_3states =  readNPY(variable_name3states);
    % hmm_2states_es =  readNPY('2states_entiresession.npy');

    play_state = 0;
    non_play_state = 1;

    classificator_tp = sum(ismember(hmm_states(time_limit),play_state) & classification_matrix((time_limit),end) ==1);
    classificator_tn = sum(ismember(hmm_states(time_limit),non_play_state) & classification_matrix((time_limit),end) ==0);
    classificator_fn = sum(ismember(hmm_states(time_limit),non_play_state) & classification_matrix((time_limit),end) ==1);
    classificator_fp = sum(ismember(hmm_states(time_limit),play_state) & classification_matrix((time_limit),end) ==0);

    if classificator_tp<classificator_fn
        play_state=1;
        non_play_state =0;
        classificator_tp = sum(ismember(hmm_states(time_limit),play_state) & classification_matrix((time_limit),end) ==1);
        classificator_tn = sum(ismember(hmm_states(time_limit),non_play_state) & classification_matrix((time_limit),end) ==0);
        classificator_fn = sum(ismember(hmm_states(time_limit),non_play_state) & classification_matrix((time_limit),end) ==1);
        classificator_fp = sum(ismember(hmm_states(time_limit),play_state) & classification_matrix((time_limit),end) ==0);
    end



    figure
    bar([([classificator_tp classificator_fn])/(classificator_fn+classificator_tp) [classificator_tn classificator_fp]/(classificator_tn+ classificator_fp)])
    xticklabels({'TruePositives','FalseNegative','TrueNegative','FalsePositive'})

    title([file2save, ' P', num2str(pt)])
    pause(.1)
    
    confusion_matrix = [([classificator_tp classificator_fn])/(classificator_fn+classificator_tp) [classificator_tn classificator_fp]/(classificator_tn+ classificator_fp)];
    save([dir_output, '\',file2save, ' P', num2str(pt), ' confusion_matrix'], 'confusion_matrix')


    %% 7 (ESTIAMTE HMM ONSET OFFSET: beg_end_times) + plot hmm and play bouts together
    disp('ESTIAMTE HMM ONSET OFFSET')
    hmm_states = hmm_states==play_state;

    if hmm_states(end)==1
        end_state= [find(diff(hmm_states)==-1);numel(hmm_states)];
    else
        end_state= find(diff(hmm_states)==-1);
    end
    if hmm_states(1)==1
        start_state= [1;find(diff(hmm_states)==1)];
    else
        start_state= find(diff(hmm_states)==1);
    end
    beg_end_times = adjusted_binned_time([start_state end_state  ]);


    A = hmm_3states';

    % Find where the value changes
    changePoints = [1, find(diff(A) ~= 0) + 1, length(A) + 1];

    % Preallocate results
    results = [];

    % Loop over each segment
    for i = 1:length(changePoints) - 1
        startIdx = changePoints(i);
        endIdx = changePoints(i+1) - 1;
        value = A(startIdx);
        results = [results; value, startIdx, endIdx];
    end

    beg_end_times_3states =[double(results(:,1)),  adjusted_binned_time(results(:, [2 3]))];

    %% 8 detect behaviors within hmm states (define play and transition behaviors)
    disp('detecting behaviors within hmm states')
    Behavior2estimate = Behavior(Behavior.Start>=binned_time(1) & Behavior.End<=binned_time(end),:);

    behavior_types  = unique(Behavior2estimate.Type);
    proportions     = nan(numel(behavior_types),3);
    numbers         = nan(numel(behavior_types),1);

    for j=1:numel(behavior_types)

        beh_index = ismember(Behavior2estimate.Type, behavior_types{j});
        transition = any((Behavior2estimate.Start(beh_index)<=beg_end_times(:,1)' &  Behavior2estimate.End(beh_index)>=beg_end_times(:,1)') | ...
            (Behavior2estimate.Start(beh_index)<=beg_end_times(:,2)' &  Behavior2estimate.End(beh_index)>=beg_end_times(:,2)') ,2);

        within    = any(Behavior2estimate.Start(beh_index)>=beg_end_times(:,1)' &  Behavior2estimate.End(beh_index)<=beg_end_times(:,2)',2);
        outside   = ~transition & ~within;

        proportions(j,1) = sum(transition)/numel(transition);
        proportions(j,2) = sum(within)/numel(transition);
        proportions(j,3) = sum(outside)/numel(transition);
        numbers(j)      = numel(transition);
    end


    play_behavior_struct = [];
    play_behavior_struct.proportions = proportions;
    play_behavior_struct.behavior_types = behavior_types;
    play_behavior_struct.numbers = numbers;

    save([dir_output, '\',file2save, ' P', num2str(pt), ' play_behavior_struct_'], 'play_behavior_struct', '-v7.3')
    %% 9 estimate behavior onset to hmm states
    disp('estimating behavior onset to hmm states')
    play_behaviors_type1      = {'Pounce_A', 'Pounce_B', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD'};

    psth_edges = psth_bin(1):bin_size:psth_bin(2);
    filled_hmm_states = beg_end_times(beg_end_times(:,1)>=ListOfPartners(pt,1) & beg_end_times(:,1)<=ListOfPartners(pt,2),:);
    filled_hmm_3states = beg_end_times_3states(beg_end_times_3states(:,2)>=ListOfPartners(pt,1) & beg_end_times_3states(:,3)<=ListOfPartners(pt,2),:);
    % filled_hmm_states = beg_end_times;
    % filled_play_bouts=filled_play_bouts(filled_play_bouts(:,2)<250,:);
    behavior_tpes           = unique(Behavior.Type);
    behavior_onset          = zeros(numel(behavior_tpes),size(filled_hmm_states,1), numel(psth_edges));
    behavior_offset         = zeros(numel(behavior_tpes),size(filled_hmm_states,1), numel(psth_edges));
    behavior_onset_3states  = zeros(3,numel(behavior_tpes),size(filled_hmm_3states,1), numel(psth_edges));
    behavior_offset_3states = zeros(3,numel(behavior_tpes),size(filled_hmm_3states,1), numel(psth_edges));

    % properties2plot =
    Behavior2estimate = Behavior(Behavior.Start>=binned_time(1) & Behavior.End<=binned_time(end),:);
    is_there_play_beh   =  zeros(size(filled_hmm_states,1), numel(psth_edges));
    is_there_noplay_beh  =  zeros(size(filled_hmm_states,1), numel(psth_edges));

    for bn=1:numel(behavior_tpes)

        these_beh = find(ismember(Behavior2estimate.Type,behavior_tpes{bn}));
        for pb=1:size(filled_hmm_states,1)
            for sub_beh = these_beh'
                beh_start = Behavior2estimate.Start(sub_beh);
                beh_end = Behavior2estimate.End(sub_beh);
                time_index = psth_edges+filled_hmm_states(pb,1)>beh_start & psth_edges+filled_hmm_states(pb,1)<beh_end;
                behavior_onset(bn,pb,time_index) = behavior_onset(bn,pb,time_index)+1;
                if ismember(behavior_tpes{bn}, play_behaviors_type1)
                    is_there_play_beh(pb,time_index)=is_there_play_beh(pb,time_index)+1;
                else
                    is_there_noplay_beh(pb,time_index)=is_there_noplay_beh(pb,time_index)+1;
                end
                time_index = psth_edges+filled_hmm_states(pb,2)>beh_start & psth_edges+filled_hmm_states(pb,2)<beh_end;
                behavior_offset(bn,pb,time_index) = behavior_offset(bn,pb,time_index)+1;
            end
        end



        for pb=1:size(filled_hmm_3states,1)
            for sub_beh = these_beh'
                beh_start = Behavior2estimate.Start(sub_beh);
                beh_end = Behavior2estimate.End(sub_beh);

                time_index = psth_edges+filled_hmm_3states(pb,2)>beh_start & psth_edges+filled_hmm_3states(pb,2)<beh_end;
                behavior_onset_3states(filled_hmm_3states(pb,1)+1,bn,pb,time_index) = behavior_onset_3states(filled_hmm_3states(pb,1)+1,bn,pb,time_index)+1;

                time_index = psth_edges+filled_hmm_3states(pb,3)>beh_start & psth_edges+filled_hmm_3states(pb,3)<beh_end;
                behavior_offset_3states(filled_hmm_3states(pb,1)+1,bn,pb,time_index) = behavior_offset_3states(filled_hmm_3states(pb,1)+1,bn,pb,time_index)+1;
            end
        end

    end

    is_there_play_beh(is_there_play_beh>1)=1;
    behavior_onset_offset_struct =[];
    behavior_onset_offset_struct.behavior_onset             = behavior_onset;
    behavior_onset_offset_struct.behavior_offset            = behavior_offset;
    behavior_onset_offset_struct.psth_edges                 = psth_edges;
    behavior_onset_offset_struct.filled_play_bouts          = filled_hmm_states;
    behavior_onset_offset_struct.behavior_onset_3states     = behavior_onset_3states;
    behavior_onset_offset_struct.behavior_offset_3states    = behavior_offset_3states;
    behavior_onset_offset_struct.filled_hmm_3states         = filled_hmm_3states;
    behavior_onset_offset_struct.behavior_tpes              = behavior_tpes;
    save([dir_output, '\',file2save, ' P', num2str(pt),' behavior_onset_offset_struct'], 'behavior_onset_offset_struct', '-v7.3')
    %% 10 estimate variable onset to hmm states
    disp('estimate variable onset to hmm states')

    % filled_hmm_states = beg_end_times;

    psth_edges = psth_bin(1):beh_bin:psth_bin(2);
    range2check = [0 .250];


    beh_properties_onset = zeros(numel(ALL_VARIABLE_NAMES),size(filled_hmm_states,1), numel(psth_edges));
    beh_properties_offset = zeros(numel(ALL_VARIABLE_NAMES),size(filled_hmm_states,1), numel(psth_edges));


    beh_properties_onset_3states    = zeros(3,numel(ALL_VARIABLE_NAMES),size(filled_hmm_3states,1), numel(psth_edges));
    beh_properties_offset_3states   = zeros(3,numel(ALL_VARIABLE_NAMES),size(filled_hmm_3states,1), numel(psth_edges));

    mean_responses      = zeros(numel(ALL_VARIABLE_NAMES),size(filled_hmm_states,1));
    mean_responses2     = zeros(numel(ALL_VARIABLE_NAMES),size(filled_hmm_states,1));
    response_index      = zeros(numel(ALL_VARIABLE_NAMES),size(filled_hmm_states,1));
    is_there_play_bout  = zeros(size(filled_hmm_states,1), numel(psth_edges));
    is_this_hmm         = zeros(size(filled_hmm_states,1), numel(psth_edges));
    is_there_hmm        = zeros(size(filled_hmm_states,1), numel(psth_edges));
    what_3states_is     = zeros(size(filled_hmm_states,1), numel(psth_edges));

    for hmm_n=1:size(filled_hmm_states,1)
        hmm_start   = filled_hmm_states(hmm_n,1);
        hmm_end     = filled_hmm_states(hmm_n,2);
        [~,pb_time_index ]= min(abs(binned_time-hmm_start));
        for j=0:2
            triple_states_index =  any(psth_edges+hmm_start>=beg_end_times_3states(beg_end_times_3states(:,1)==j,2) & psth_edges+hmm_start<=beg_end_times_3states(beg_end_times_3states(:,1)==j,3));
            what_3states_is(hmm_n,triple_states_index) = j;
        end
        is_there_hmm_index =  any(psth_edges+hmm_start>=filled_hmm_states(:,1) & psth_edges+hmm_start<=filled_hmm_states(:,2));
        is_there_hmm(hmm_n,is_there_hmm_index) = 1;

        is_this_hmm_index =psth_edges>0 & psth_edges<hmm_end-hmm_start;
        is_this_hmm(hmm_n,is_this_hmm_index) = 1;

        is_there_pb_index =  any(psth_edges+hmm_start>=play_bouts_table(:,1) & psth_edges+hmm_start<=play_bouts_table(:,2));
        is_there_play_bout(hmm_n,is_there_pb_index)=1;

        for variable_n=1:numel(ALL_VARIABLE_NAMES)
            [~,pb_time_index ]= min(abs(binned_time-hmm_start));
            time_entire_index = (pb_time_index-psth_limits/beh_bin):(pb_time_index+psth_limits/beh_bin);
            fitting_index_conv      = ismember(1:size(conv_all_properties,1),time_entire_index);
            fittin_index_all_prop   = ismember(time_entire_index,1:size(conv_all_properties,1));
            beh_properties_onset(variable_n,hmm_n,fittin_index_all_prop) = conv_all_properties(fitting_index_conv, variable_n);

            [~,pb_time_index ]      = min(abs(binned_time-hmm_end));
            time_entire_index       = (pb_time_index-psth_limits/beh_bin):(pb_time_index+psth_limits/beh_bin);
            fitting_index_conv      = ismember(1:size(conv_all_properties,1),time_entire_index);
            fittin_index_all_prop   = ismember(time_entire_index,1:size(conv_all_properties,1));
            beh_properties_offset(variable_n,hmm_n,fittin_index_all_prop) = conv_all_properties(fitting_index_conv, variable_n);

            mean_responses(variable_n,hmm_n)  = mean(beh_properties_onset(variable_n,hmm_n,psth_edges>=0 & psth_edges<=hmm_end-hmm_start),3);
            mean_responses2(variable_n,hmm_n) = mean(beh_properties_onset(variable_n,hmm_n,psth_edges>=.250 & psth_edges<=.5),3);
            response_index(variable_n,hmm_n)  = mean(beh_properties_onset(variable_n,hmm_n,psth_edges>=0 & psth_edges<=hmm_end-hmm_start),3)/mean(beh_properties_onset(variable_n,hmm_n,psth_edges>=-5 & psth_edges<=0),3);

        end

    end



    for hmm_n=1:size(filled_hmm_3states,1)
        hmm_start   = filled_hmm_3states(hmm_n,2);
        hmm_end     = filled_hmm_3states(hmm_n,3);
        [~,pb_time_index ]= min(abs(binned_time-hmm_start));

        for variable_n=1:numel(ALL_VARIABLE_NAMES)
            [~,pb_time_index ]= min(abs(binned_time-hmm_start));
            time_entire_index = (pb_time_index-psth_limits/beh_bin):(pb_time_index+psth_limits/beh_bin);
            fitting_index_conv      = ismember(1:size(conv_all_properties,1),time_entire_index);
            fittin_index_all_prop   = ismember(time_entire_index,1:size(conv_all_properties,1));
            beh_properties_onset_3states(filled_hmm_3states(hmm_n,1)+1,variable_n,hmm_n,fittin_index_all_prop) = conv_all_properties(fitting_index_conv, variable_n);

            [~,pb_time_index ]      = min(abs(binned_time-hmm_end));
            time_entire_index       = (pb_time_index-psth_limits/beh_bin):(pb_time_index+psth_limits/beh_bin);
            fitting_index_conv      = ismember(1:size(conv_all_properties,1),time_entire_index);
            fittin_index_all_prop   = ismember(time_entire_index,1:size(conv_all_properties,1));
            beh_properties_offset_3states(filled_hmm_3states(hmm_n,1)+1,variable_n,hmm_n,fittin_index_all_prop) = conv_all_properties(fitting_index_conv, variable_n);

        end

    end



    variable_onset_struct =[];
    variable_onset_struct.beh_properties_onset              = beh_properties_onset;
    variable_onset_struct.beh_properties_offset             = beh_properties_offset;
    variable_onset_struct.mean_responses                    = mean_responses;
    variable_onset_struct.psth_edges                        = psth_edges;
    variable_onset_struct.variable_types                    = ALL_VARIABLE_NAMES;
    variable_onset_struct.filled_play_bouts                 = filled_hmm_states;
    variable_onset_struct.beh_properties_onset_3states      = beh_properties_onset_3states;
    variable_onset_struct.beh_properties_offset_3states     = beh_properties_offset_3states;
    variable_onset_struct.filled_play_bouts_3states          = filled_hmm_3states;
    save([dir_output, '\',file2save, ' P', num2str(pt), ' variable_onset_struct'], 'variable_onset_struct', '-v7.3')

    %% 11 plot hmm, pb, and  hmm no-pb

    % colormap(1-gray)

    [hmm_length_ordered, pb_order] = sort(diff(filled_hmm_states'));
    properly_labeled = filled_hmm_states(pb_order,1)<=Inf;

    figure('units','normalized','outerposition',[0 0 .25 1]);

    colormap(gray)
    imagesc(psth_edges, 1:sum(properly_labeled),  what_3states_is(pb_order(properly_labeled),:))
    hold on
    axis xy

    plot([0 0], [.5 sum(properly_labeled)+.5], 'k')
    plot(hmm_length_ordered,1:sum(properly_labeled),'k')
    title('Original Assignment')

    %% 11.B RE_ASIGN STATES NUMBERS
    % x = input('write re assignment order');
    % re_assignment = x;
    what_3states_is_corrected = what_3states_is;
    state_intercept = nan(3,1);
    for sn = 1:3
        state_intercept(sn) = sum(sum(is_there_play_beh & what_3states_is==sn-1));
    end


    [~, re_assignment] = sort(state_intercept, 'descend');
    original_order  = fliplr(1:3);
    re_assignment =  original_order(re_assignment);


    for j=0:2

        what_3states_is_corrected(what_3states_is==j)=re_assignment(j+1);
    end
    % what_3states_is_corrected;



    figure('units','normalized','outerposition',[.25 0 .25 1]);
    colormap(gray)
    imagesc(psth_edges, 1:sum(properly_labeled),  what_3states_is_corrected(pb_order(properly_labeled),:))
    hold on
    axis xy
    title('Current Re-assignment')
    pause(.1)
    disp('Current reasignment')
    disp(re_assignment)
    x = input('write re assignment order');
    re_assignment = x;

    for j=0:2
        what_3states_is_corrected(what_3states_is==j)=re_assignment(j+1);
    end

    what_3states_is = what_3states_is_corrected;



    figure('units','normalized','outerposition',[.5 0 .25 1]);
    colormap(gray)
    imagesc(psth_edges, 1:sum(properly_labeled),  what_3states_is(pb_order(properly_labeled),:))
    hold on
    axis xy
    title('Final Re-assignment')
    pause(.1)

    %% 12 plot var distribution between  hmm, pb, and  hmm no-pb
    disp('Estiamting variables for model')
    no_hmm_index        = is_there_hmm==0;
    % hmm_and_pb_index    = is_there_hmm==1 & is_there_play_beh==1;
    % hmm_and_NOpb_index  = is_there_hmm==1 & is_there_play_beh==0;
    hmm_and_pb_index    = is_this_hmm==1 & is_there_play_beh==1;
    hmm_and_NOpb_index  = is_this_hmm==1 & is_there_play_beh==0;


    variables_for_model = [];
    model_properties = zeros(size(squeeze(beh_properties_onset(1,:,:))));

    for variable_n=1:numel(ALL_VARIABLE_NAMES)
        this_matrix = zscore(squeeze(beh_properties_onset(variable_n,:,:)));
        no_hmm          =  this_matrix(no_hmm_index)';
        hmm_and_pb      = this_matrix(hmm_and_pb_index)';
        hmm_and_NOpb    = this_matrix(hmm_and_NOpb_index)';
        if variable_n<numel(ALL_VARIABLE_NAMES)
            variables_for_model = [variables_for_model,[hmm_and_pb,hmm_and_NOpb]'];
        else
            variables_for_model = [variables_for_model,[hmm_and_pb,hmm_and_NOpb]',[(hmm_and_pb*0)+1,hmm_and_NOpb*0]'];
        end


    end

    %% 13 creat GMM  hmm, pb, and  hmm no-pb
    disp('creat GMM  hmm, pb, and  hmm no-pb')
    % [yfit,scores] = SVM_model.predictFcn(variables_for_model(:, 1:end-1)) ;
    table4model = array2table(variables_for_model);
    table4model.Properties.VariableNames(1:end-1) = ALL_VARIABLE_NAMES;
    table4model.Properties.VariableNames(end) = {'PlayBout'};
    GLM_model = fitglm(table4model,'Distribution','binomial');
    all_Variables = nan(size(beh_properties_onset,2)*size(beh_properties_onset,3), size(beh_properties_onset,1));

    for j=1:14
        for k=1:size(beh_properties_onset,2)

            all_Variables((size(beh_properties_onset,3)*(k-1) +1):(size(beh_properties_onset,3)*k),j) = zscore(beh_properties_onset(j,k,:));
        end
    end

    ALL_varialbes_table = array2table(all_Variables);
    ALL_varialbes_table.Properties.VariableNames   =table4model.Properties.VariableNames(1:end-1);
    [yfit,scores]  = GLM_model.predict(all_Variables) ;
    % [yfit,scores]  = SVM_model.predictFcn(ALL_varialbes_table) ;
    entire_figure1 = nan(size(beh_properties_onset,3),size(beh_properties_onset,2))';
    entire_figure2 = nan(size(beh_properties_onset,3),size(beh_properties_onset,2))';
    %% 14 Estimate play probability
    disp('Estimate play probability')
    for k=1:size(beh_properties_onset,2)

        score2 = scores((size(beh_properties_onset,3)*(k-1) +1):(size(beh_properties_onset,3)*k),2);
        score1 =scores((size(beh_properties_onset,3)*(k-1) +1):(size(beh_properties_onset,3)*k),1);

        entire_figure1(k,:)=score1 ;
        entire_figure2(k,:)=score2 ;

    end
    %% 15 save prediction structure
    disp('Saving Structure data')
    prediction_struct = [];
    prediction_struct.entire_figure     = entire_figure1;
    prediction_struct.filled_play_bouts = filled_hmm_states;
    prediction_struct.psth_edges        = psth_edges;
    prediction_struct.GLM_model         = GLM_model;
    prediction_struct.is_there_play_bout= is_there_play_bout;
    prediction_struct.is_there_play_beh = is_there_play_beh;
    prediction_struct.is_there_hmm      = is_there_hmm;
    prediction_struct.is_this_hmm       = is_this_hmm;
    prediction_struct.what_3states_is   = what_3states_is;
    prediction_struct.re_assignment     = re_assignment;
    save([dir_output, '\',file2save, ' P', num2str(pt),' prediction_struct'], 'prediction_struct', '-v7.3')
    disp('Structure saved')
    %
    % model_properties(hmm_and_pb_index==1) = scores(variables_for_model(:,end)==0,1);
    % model_properties(hmm_and_NOpb_index==1) = scores(variables_for_model(:,end)==1,1);

    %% 16 estimate real call onsets offsets and call properties

    filled_hmm_states = beg_end_times;
    psth_edges = psth_bin(1):bin_size:psth_bin(2);

    behavior_tpes = unique(Behavior.Type);
    call_onset = zeros(size(filled_hmm_states,1), numel(psth_edges));
    call_offset = zeros(size(filled_hmm_states,1), numel(psth_edges));

    properties2estimate = {'PrincipalFrequencykHz','LowFreqkHz','HighFreqkHz','DeltaFreqkHz','FrequencyStandardDeviationkHz','SlopekHzs','Sinuosity','PeakFreqkHz'};
    properties2estimate = properties2estimate(ismember(properties2estimate, CallStats.Properties.VariableNames));
    properties_call_onset   = nan(numel(properties2estimate),size(filled_hmm_states,1), numel(psth_edges));
    properties_call_offset  = nan(numel(properties2estimate),size(filled_hmm_states,1), numel(psth_edges));


    % filled_play_bouts = beg_end_times
    n_calls_per_hmm         = nan(size(filled_hmm_states,1),1);

    for pb=1:size(filled_hmm_states,1)
        n_calls_per_hmm(pb) = nnz(CallStats.BeginTimes>=filled_hmm_states(pb,1) & CallStats.EndTimes<=filled_hmm_states(pb,2));
        call_list = find(CallStats.EndTimes>=psth_edges(1)+filled_hmm_states(pb,1) &  CallStats.BeginTimes<=psth_edges(end)+filled_hmm_states(pb,2))';
        for call = call_list
            call_start = CallStats.BeginTimes(call);
            call_end = CallStats.EndTimes(call);



            time_index = psth_edges+filled_hmm_states(pb,1)>call_start & psth_edges+filled_hmm_states(pb,1)<call_end;
            call_onset(pb,time_index) = call_onset(pb,time_index)+1;

            for p_n = 1:numel(properties2estimate)

                properties_call_onset(p_n,pb,time_index) = CallStats{call,properties2estimate{p_n}};

            end

            call_start = CallStats.BeginTimes(call);
            call_end = CallStats.EndTimes(call);



            time_index = psth_edges+beg_end_times(pb,2)>call_start & psth_edges+beg_end_times(pb,2)<call_end;
            call_offset(pb,time_index) = call_offset(pb,time_index)+1;


            for p_n = 1:numel(properties2estimate)

                properties_call_offset(p_n,pb,time_index) = CallStats{call,properties2estimate{p_n}};

            end

        end
    end


    call_struct = [];
    call_struct.call_onset              = call_onset;
    call_struct.call_offset             = call_offset;
    call_struct.filled_play_bouts       = filled_hmm_states;
    call_struct.psth_edges              = psth_edges;
    call_struct.properties_call_onset   = properties_call_onset;
    call_struct.properties_call_offset  = properties_call_offset;
    call_struct.properties2estimate     = properties2estimate;
    call_struct.n_calls_per_hmm         = n_calls_per_hmm;
    save([dir_output, '\',file2save, ' P', num2str(pt),' call_struct'], 'call_struct', '-v7.3')
end

end
