function [all_behavior_list,all_non_labeled_list,all_spatial_prop] =get_behavior_properties(animal_code, behaviors2check,call_prop_list, extended_time,get_non_labeled)

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';

synch_folder = [data_root, '\Synch data'];
traking_folder = [data_root, '\Traking backups'];
behavior_folder = [data_root, '\Behavior backups'];
call_folder = [data_root, '\CallDetectionBackup'];
animal_code_params =strsplit(animal_code, ' ');
repeated_animal     = animal_code_params{3}; %Single2, Dual Single
% animal_code         =  [animal_batch ' ' date ' ' repeated_animal]

load([synch_folder, '\',animal_code, '\synch_model_video2audio.mat'], 'synch_model_video2audio')
beh_bin = 0.01;

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
all_behavior_list    = {};
all_non_labeled_list = {};

for pt = 1:size(ListOfPartners,1) %Create separate file sfore ach partner
    clear traking_structure
    load([traking_folder, '\',list_of_trakings(pt).name], 'traking_structure')
    traking_time = predict(synch_model_video2audio,(traking_structure.frames2stract/30)');
    if restrict2Partnerssession
        index2inlcude = traking_time>=ListOfPartners(pt,1) &  traking_time<=ListOfPartners(pt,2);
    else
        index2inlcude = true(size(traking_structure.time));
    end

    if restrict2Partnerssession       
        behavior_index = ismember(Behavior.Type, behaviors2check) & Behavior.Start>=ListOfPartners(pt,1) &  Behavior.End<=ListOfPartners(pt,2);
        sum(behavior_index)
    end

    traking_time = traking_time(index2inlcude);


    %% 3 EStimate behavioral varialbes
    disp('Estimating behavioral variables')

    animal_pos      = traking_structure.animal_pos(index2inlcude,:);
    animal_pos(:,1) = smoothdata(animal_pos(:,1), 'loess',5);
    animal_pos(:,2) = smoothdata(animal_pos(:,2), 'loess',5);

    partner_pos = traking_structure.partner_pos(index2inlcude,:);
    partner_pos(:,1) = smoothdata(partner_pos(:,1), 'loess',5);
    partner_pos(:,2) = smoothdata(partner_pos(:,2), 'loess',5);

    inner_box = traking_structure.inner_box;
    arena_center = mean(traking_structure.inner_border);
    all_inner_box = [inner_box;inner_box(1,:)] - arena_center;

    all_inner_border = [traking_structure.inner_border;traking_structure.inner_border(1,:)]- arena_center;

    time_range  = [min(traking_time) max(traking_time)];
    binned_time = (beh_bin*floor(time_range(1)/beh_bin)):beh_bin:(beh_bin*ceil(time_range(2)/beh_bin));
    % binned_video_time = binned_time;



    animal_pos_x = interp1(traking_time,animal_pos(:,1), binned_time, 'cubic');
    animal_pos_y = interp1(traking_time,animal_pos(:,2), binned_time, 'cubic');
    animal_pos   = [animal_pos_x', animal_pos_y']-arena_center;
    L = sum(sqrt(diff(all_inner_border(:,1)).^2 + diff(all_inner_border(:,2)).^2));

    x_inner = all_inner_box(:,1);
    y_inner = all_inner_box(:,2);
    x_outer = all_inner_border(:,1);
    y_outer = all_inner_border(:,2);
    center = [0 0];
    animal_wall_pos_values = nan(size(animal_pos,1),1);

    for j=1:numel(animal_pos(:,1))
        x = animal_pos(j,1);
        y =animal_pos(j,2);
        animal_wall_pos_values(j) = normalizedRadialPosition(x, y, x_inner, y_inner, x_outer, y_outer, center);
    end

    animal_pos          = animal_pos/L;
    % [animal_pos_polar_th, animal_pos_polar_ro] = cart2pol(animal_pos(:,1),animal_pos(:,2));

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

    partner_wall_pos_values = nan(size(partner_pos,1),1);
    for j=1:numel(animal_pos(:,1))
        x = partner_pos(j,1);
        y =partner_pos(j,2);
        partner_wall_pos_values(j) = normalizedRadialPosition(x, y, x_inner, y_inner, x_outer, y_outer, center);
    end

    partner_pos         = partner_pos/L;
    % [partner_pos_polar_th, partner_pos_polar_ro] = cart2pol(partner_pos(:,1),partner_pos(:,2));
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
    relative_angle_speed= cart2pol(relative_velocity(:,1), relative_velocity(:,2));
    relative_angle_acc  = angdiff(relative_angle_speed);

    relative_speed      = diff(relative_distance);
    relative_acc        = diff(relative_speed);
    relative_speed      = sqrt(sum(relative_speed.*relative_speed,2));
    relative_acc        = sqrt(sum(relative_acc.*relative_acc,2));

    animal_spatial_properties   = [  animal_speed(1:end-1)  abs(animal_angle_speed(1:end-1))  abs(animal_angle_acc)  animal_acc  animal_wall_pos_values(1:end-2)];
    animal_spatial_properties   = (animal_spatial_properties - repmat(mean(animal_spatial_properties, 'omitmissing'), size(animal_spatial_properties,1),1))*diag(1./std(animal_spatial_properties, 'omitmissing'));
    animal_prop = {'Speed', 'AngularSpeed', 'AngulaAcc', 'Acc', 'WallPos'};
    partner_spatial_properties  = [  partner_speed(1:end-1) abs(partner_angle_speed(1:end-1)) abs(partner_angle_acc) partner_acc partner_wall_pos_values(1:end-2)];
    partner_spatial_properties   = (partner_spatial_properties - repmat(mean(partner_spatial_properties, 'omitmissing'), size(partner_spatial_properties,1),1))*diag(1./std(partner_spatial_properties, 'omitmissing'));
    shared_prop = {'RelativeDist', 'RelativeSpeed', 'RelativeAngularSpeed', 'RelativeAngularAcc', 'RelativeAcc'};
    relative_spatial_properties = [relative_distance(1:end-2,:)   relative_speed(1:end-1)  abs(relative_angle_speed(1:end-1)) abs(relative_angle_acc) relative_acc];
    relative_spatial_properties   = (relative_spatial_properties - repmat(mean(relative_spatial_properties, 'omitmissing'), size(relative_spatial_properties,1),1))*diag(1./std(relative_spatial_properties, 'omitmissing'));
    all_spatial_prop = [animal_prop,shared_prop];

    call_properties             = zeros(numel(binned_time)-1,numel(call_prop_list)+1);
    for j=1:numel(binned_time)-1
        is_ther_call = CallStats.BeginTimes<=binned_time(j) &  CallStats.EndTimes>=binned_time(j+1);
        call_properties(j,1) = any(is_ther_call);
        if sum(is_ther_call)==1
            call_properties(j,2:end) = (CallStats{is_ther_call,call_prop_list});
        elseif sum(is_ther_call)>1
            call_properties(j,2:end) = (mean(CallStats{is_ther_call,call_prop_list})-min(CallStats{:,call_prop_list}))/range( CallStats{:,call_prop_list});
        end
    end

    call_properties   = (call_properties - repmat(mean(call_properties, 'omitmissing'), size(call_properties,1),1))*diag(1./std(call_properties, 'omitmissing'));



    %%

    behavior_list = cell(sum(behavior_index),6);
    col_index = 1;
    for j=find(behavior_index)'
        beh_start =  Behavior.Start(j)+extended_time(1);
        beh_end   =  Behavior.End(j)+extended_time(2);
        [~, beg_loc] = min(abs(binned_time-beh_start));
        [~, end_loc] = min(abs(binned_time-beh_end));
        end_loc = min(end_loc,size(relative_spatial_properties,1));

        if ismember(Behavior.Animal(j), repeated_animal)
            behavior_list{col_index,1} = [animal_spatial_properties(beg_loc:end_loc,:) relative_spatial_properties(beg_loc:end_loc,:) call_properties(beg_loc:end_loc,:)];
            behavior_list{col_index,3} = 'Animal';
        else
            behavior_list{col_index,1} = [partner_spatial_properties(beg_loc:end_loc,:) relative_spatial_properties(beg_loc:end_loc,:) call_properties(beg_loc:end_loc,:)];
            behavior_list{col_index,3} = 'Partner';
        end
        behavior_list{col_index,2} = Behavior.Type{j};
        behavior_list{col_index,4} = animal_code;
        behavior_list{col_index,5} = j;
        behavior_list{col_index,6} = pt;
        col_index = col_index+1;
    end

    all_behavior_list = [all_behavior_list;behavior_list];

    if get_non_labeled
        disp('Estiamting not labeled behaviors')
        behavior_list_non_lebaled = {};
        behavior_lengths = Behavior.End(behavior_index)-Behavior.Start(behavior_index);

        n= 1;
        col_index = 1;
        start_time = ListOfPartners(pt,1);
        end_time = start_time+behavior_lengths(n);

        while end_time<=ListOfPartners(pt,2)
            beh_start       = start_time;
            beh_end         = end_time;
            [~, beg_loc]    = min(abs(binned_time-beh_start));
            [~, end_loc]    = min(abs(binned_time-beh_end));
            end_loc         = min(end_loc,size(animal_spatial_properties,1));

            behavior_list_non_lebaled{2*col_index -1,1}                 = [animal_spatial_properties(beg_loc:end_loc,:) relative_spatial_properties(beg_loc:end_loc,:) call_properties(beg_loc:end_loc,:)];
            behavior_list_non_lebaled{2*col_index -1,2}                 = 'Animal';
            behavior_list_non_lebaled{2*col_index -1,3}                 = [start_time end_time];
            behavior_list_non_lebaled{2*col_index -1,4}                 = pt;

            behavior_list_non_lebaled{2*col_index,1}                    = [partner_spatial_properties(beg_loc:end_loc,:) relative_spatial_properties(beg_loc:end_loc,:) call_properties(beg_loc:end_loc,:)];
            behavior_list_non_lebaled{2*col_index,2}                    = 'Partner';

            behavior_list_non_lebaled{2*col_index,3}                    = [start_time end_time];
            behavior_list_non_lebaled{2*col_index,4}                    = pt;


            n = n+1;
            col_index = col_index+1;
            if n>numel(behavior_lengths)
                n=1;
            end
            start_time  = end_time;
            end_time    = start_time+behavior_lengths(n);


        end
        all_non_labeled_list = [all_non_labeled_list;behavior_list_non_lebaled];
    end
end

end






