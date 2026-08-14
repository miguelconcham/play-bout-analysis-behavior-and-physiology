%% 1 SET PARAMS
function [] = Create_HMM_inputs(dir_output,animal_code)

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\Codes repository\Data';

synch_folder = [data_root, '\Synch data'];
traking_folder = [data_root, '\Traking backups'];
behavior_folder = [data_root, '\Behavior backups'];
call_folder = [data_root, '\CallDetectionBackup'];
animal_code_params =strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
date                = animal_code_params{2};
repeated_animal     = animal_code_params{3}; %Single2, Dual Single
% animal_code         =  [animal_batch ' ' date ' ' repeated_animal]
file2save = animal_code;
play_state      = 0;
non_play_state  = 1;

load([synch_folder, '\',animal_code, '\synch_model_video2audio.mat'], 'synch_model_video2audio')


play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB'};
beh_bin = 0.01; %for hmm stimate
conv_length = 1; % convuolution lenbgth for smothing all data
list_of_trakings = dir([traking_folder, '\*',animal_code, '*'])

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


    partner_pos_x       = interp1(traking_time,partner_pos(:,1), binned_time, 'cubic');
    partner_pos_y       = interp1(traking_time,partner_pos(:,2), binned_time, 'cubic');
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


    elemets_classified                                              = ~any(isnan(classification_matrix),2);
    adjusted_binned_time                                            = binned_time(~any(isnan(classification_matrix),2));
    classification_matrix(any(isnan(classification_matrix),2),:)    = [];
    classification_matrix4HMM                                       = classification_matrix(:, 1:end-1); % remove play states to train HMM!!!
    

    writeNPY(classification_matrix4HMM, [dir_output, '\',file2save, ' P', num2str(pt),'.npy'])
    save([dir_output, '\',file2save, ' P', num2str(pt),'_PropAndTime.mat'], 'adjusted_binned_time','all_properties','conv_all_properties' ,'classification_matrix4HMM','ALL_VARIABLE_NAMES')
end
end
