function [play_bouts_table, play_bout_behaviors,interplay_bout_behaviors] = play_bout(config)

field_expected_names = {'Behavior','repeated_animal','animal_types','play_behaviors','beh_bin','conv_length','behavior_window'};
fields = fieldnames(config);

if any(~ismember(field_expected_names,fields))
    disp('The following fields are missing')
    disp(field_expected_names(~ismember(field_expected_names,fields)))
    play_bouts_table    = [];
    play_bout_behaviors = [];
else


    Behavior            = config.Behavior;
    repeated_animal     = config.repeated_animal;
    animal_types        = config.animal_types;
    play_behaviors      = config.play_behaviors;
    beh_bin             = config.beh_bin;
    conv_length         = config.conv_length;
    behavior_window     = config.behavior_window;



    extra_time = 10;

    partner_names = animal_types;
    partner_names(ismember(animal_types, repeated_animal)) = [];
    animal_sessions = nan(numel(partner_names),2);
    for j=1:numel(partner_names)
        animal_sessions(j,:) = [min(Behavior.Start(ismember(Behavior.Animal,partner_names{j}  ))) max(Behavior.End(ismember(Behavior.Animal,partner_names{j} )))];
    end
    [~,session_order] = sort(animal_sessions(:,1));
    partner_names = partner_names(session_order);
    animal_sessions = animal_sessions(session_order,:);
    animal_types = [repeated_animal;partner_names];
    play_behavior_events = Behavior{ismember(Behavior.Type2, play_behaviors),{'Start', 'End'}};
    play_behavior_animal = Behavior{ismember(Behavior.Type2, play_behaviors),{'Animal'}};

    beh_time = 0:beh_bin:(beh_bin*round(max(max(play_behavior_events+extra_time))/beh_bin));
    behavior_freq = zeros(numel(animal_types), numel(beh_time) );
    for an =1:numel(animal_types)

        this_animal_behavior = play_behavior_events(ismember(play_behavior_animal, animal_types{an}),:);

        for bn=1:size(this_animal_behavior,1)
            beh_start = this_animal_behavior(bn,1);
            beh_end = this_animal_behavior(bn,2);
            time_index = beh_time>=beh_start & beh_time<=beh_end;
            behavior_freq(an,time_index) = behavior_freq(an,time_index)+1;
        end
    end
    conv_time = linspace(-conv_length, conv_length, 1 +(conv_length/beh_bin));
    conv_fun =normpdf(conv_time);
    conv_fun = conv_fun/sum(conv_fun);

    conv_behavior = behavior_freq;
    for an =1:numel(animal_types)
        conv_this_beh = conv(behavior_freq(an,:), conv_fun);
        conv_behavior(an,:) = conv_this_beh(round(.5*(conv_length/beh_bin) +1):(end-round(.5*(conv_length/beh_bin))));
    end

    data2segment = mean(conv_behavior);
    data2segment(data2segment>0)=1;
    L               = bwlabeln(data2segment);
    binned_beh      = L>0;

    events_start    = find(diff(binned_beh)==1);
    event_end       = find(diff(binned_beh)==-1);
    event_end       = [event_end max(find(binned_beh>0))];
    event_end       = unique(event_end);

    bout_lengths    = beh_time(event_end)-beh_time(events_start);
    sorted_behavior_start = Behavior;
    [~,beh_order]   = sort(Behavior.Start);
    sorted_behavior_start = sorted_behavior_start(beh_order,:);

    sorted_behavior_end    = Behavior;
    [~,beh_order]   = sort(Behavior.End);
    sorted_behavior_end = sorted_behavior_end(beh_order,:);

    original_start_end = beh_time([events_start' event_end']);
    play_bouts_table    = nan(numel(bout_lengths),2);
    play_bout_behaviors = cell(numel(bout_lengths),1);

    bout_start=[];bout_end=[];
    for bn=1:numel(bout_lengths)
        bout_start = original_start_end(bn,1);
        bout_end = original_start_end(bn,2);

        first_play_behavior         = min(find(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<=bout_end & ismember(sorted_behavior_start.Type2, play_behaviors)));
        last_play_behavior          = max(find(sorted_behavior_end.End<=bout_end & sorted_behavior_end.End>=bout_start & ismember(sorted_behavior_end.Type2, play_behaviors)));
        bout_start_2  (bn)                  = sorted_behavior_start.Start(first_play_behavior);
        bout_end_2  (bn)                  =  sorted_behavior_end.End(last_play_behavior);

        behaviors2include = Behavior.Start>=bout_start_2(bn)-behavior_window & Behavior.Start<=bout_end_2(bn)+behavior_window;
        play_bout_behaviors{bn} = Behavior(behaviors2include,:);
        play_bouts_table(bn,:) = [bout_start_2(bn) bout_end_2(bn)];
    end

    for bn=1:numel(bout_lengths)-1
        % bout_start = original_start_end(bn+1,1);
        % bout_end = original_start_end(bn,2);
        % 
        % first_play_behavior         = min(find(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<=bout_end & ismember(sorted_behavior_start.Type2, play_behaviors)));
        % last_play_behavior          = max(find(sorted_behavior_end.End<=bout_end & sorted_behavior_end.End>=bout_start & ismember(sorted_behavior_end.Type2, play_behaviors)));
        % bout_start                  = sorted_behavior_start.Start(first_play_behavior);
        % bout_end                    =  sorted_behavior_end.End(last_play_behavior);

        behaviors2include = Behavior.Start>=bout_end_2(bn)+behavior_window & Behavior.Start<=bout_start_2(bn+1)-behavior_window;
        interplay_bout_behaviors{bn} = Behavior(behaviors2include,:);
        % play_bouts_table(bn,:) = [bout_start bout_end];
    end
end
end
