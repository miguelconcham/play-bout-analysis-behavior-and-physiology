%% defining loading folders and fix parameters
repo_root     = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology';
data_root     = [repo_root, '\Data'];
saving_folder = [repo_root, '\Figure 1\outputs'];
behavior_data = [data_root, '\Behavior backups'];
%if set to true, it will save the figures from next sextion in saving
%folder:
save_bool = false;
current_folder = cd;
cd(behavior_data)
behavior_files = dir('*.txt');

beh_meta_data = cell(numel(behavior_files),3);

for bf = 1:numel(behavior_files)
    this_file_metadata = strsplit(behavior_files(bf).name, ' ');
    this_file_metadata{3}= strsplit(this_file_metadata{3}, '.');
    this_file_metadata{3} = this_file_metadata{3}{1};
    this_file_metadata{2} = str2double(this_file_metadata{2});

    beh_meta_data(bf,:) =this_file_metadata;
end

beh_meta_data = cell2table(beh_meta_data);
beh_meta_data.Properties.VariableNames = {'Animal','Session','Implanted'};
lag_time = 20;
beh_bin = 0.01;
conv_length = 1;

hist_range = [0 20];
hist_bins = .1;
hist_edges = hist_range(1):hist_bins:hist_range(2);
hist_edges_centers = .5*(hist_edges(1:end-1)+hist_edges(2:end));


cc_matrix           = nan(numel(behavior_files),2*(lag_time/beh_bin) +1 );
inter_bout_matrix   = nan(numel(behavior_files),numel(hist_edges_centers));
bout_length_matrix  = nan(numel(behavior_files),numel(hist_edges_centers));
partner_specific_comb_bout_length    = [];
partner_specific_comb_bout_interval = [];
animal_comb = [];

play_bout_structure = [];
renage_around = 1.25;

this_session_first_last_table = []
inter_play_bout_table = [];

play_bout_corr_pre = [];
play_bout_corr_post = [];

play_bout_corr_pre_zs = [];
play_bout_corr_post_zs = [];
%% Fig 1G load behavior data, and estimate play bout data (it will also generate play bout plots likje in fig 1G)
% a lot of the estimates are not used in the article

restrain2morethan1beh = false;
ALL_behavior_matrices   = [];
play_behaviors      = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'CB', 'CD', 'Escape'};
beh_matrix_range    = [-20 20];
beh_matrix_bin      = 0.01;
beh_matrix_indexes  = beh_matrix_range(1):beh_matrix_bin:beh_matrix_range(2);

all_behavior_timeplots = cell(numel(behavior_files),1);
all_behavior_timeplots_conv  = cell(numel(behavior_files),1);

for bf = 1:numel(behavior_files)
    Behavior =   readtable(behavior_files(bf).name);
    Behavior(:,2) = [];
    Behavior.Properties.VariableNames = {'Animal', 'Start', 'End', 'Length', 'Type'};
    Behavior.Type2 = Behavior.Type;
    Behavior.Type2(ismember(Behavior.Type2, {'Pounce_A', 'Pounce_B'})) = {'Pounce'}; %% Merging behaviors to Type2
    Behavior.Type2(ismember(Behavior.Type2, {'Pounce_Ai', 'Pounce_Bi'})) = {'PounceI'};
    Behavior.Type2(ismember( Behavior.Type2,'')) = {'Other'};
    Behavior(ismember(Behavior.Animal, 'Reversal'),:) = [];

    new_behavior = [];

    sessions_structure =  Behavior(ismember(Behavior.Animal, 'Session_structure'),:)
    sessions_structure(ismember(sessions_structure.Type, 'Tickling'),:) = [];

    for sn=1:size(sessions_structure,1)
        new_behavior = [new_behavior;Behavior(Behavior.Start>=sessions_structure.Start(sn) & Behavior.End<=sessions_structure.End(sn),:)];
    end

    Behavior = new_behavior;
    Behavior(ismember(Behavior.Animal, 'Session_structure'),:) = [];
    animal_types = unique(Behavior.Animal);
    repeated_animal = beh_meta_data.Implanted(bf); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    partner_names = animal_types;
    partner_names(ismember(animal_types, repeated_animal)) = [];

    config.Behavior = Behavior;
    config.repeated_animal =repeated_animal;
    config.animal_types =animal_types        ;
    config.play_behaviors =play_behaviors      ;
    config.beh_bin =beh_bin             ;
    config.conv_length =conv_length         ;

    [play_bouts_table] = play_bout(config) ;

    animal_sessions =sessions_structure{:, {'Start','End'}};

    main_partner = nan(size(animal_sessions,1),numel(partner_names));

    for  sn=1:size(sessions_structure,1)
        for pt=1:numel(partner_names)

            main_partner(sn,pt) = sum(ismember(Behavior.Animal(Behavior.Start>=sessions_structure.Start(sn) & Behavior.End<=sessions_structure.End(sn),:),partner_names{pt}));
        end
    end



    [~, loc]  =max(main_partner,[],2)
    partner_per_session = partner_names(loc);

    for pt=1:numel(partner_names)
        Behavior(Behavior.Start>=sessions_structure.Start(pt) & Behavior.End<=sessions_structure.End(pt) & ~ismember(Behavior.Animal, partner_per_session{pt}),:) = [];
    end


    new_behavior = [];


    for sn=1:size(sessions_structure,1)
        new_behavior = [new_behavior;Behavior(Behavior.Start>=sessions_structure.Start(sn) & Behavior.End<=sessions_structure.End(sn),:)];
    end
    Behavior = new_behavior;

    for j=1:numel(partner_names)
        animal_sessions(j,:) = [min(Behavior.Start(ismember(Behavior.Animal,partner_names{j}  ))) max(Behavior.End(ismember(Behavior.Animal,partner_names{j} )))];
    end

    [~,session_order] = sort(animal_sessions(:,1));
    partner_names = partner_names(session_order);
    animal_sessions = animal_sessions(session_order,:);
    animal_types = [repeated_animal;partner_names];
    play_behavior_events = Behavior{ismember(Behavior.Type2, play_behaviors),{'Start', 'End'}};
    play_behavior_animal = Behavior{ismember(Behavior.Type2, play_behaviors),{'Animal'}};

    beh_time = 0:beh_bin:(beh_bin*round(max(max(play_behavior_events+10))/beh_bin));
    behavior_freq = zeros(numel(animal_types), numel(beh_time) );
    for an =1:numel(animal_types)

        this_animal_behavior = play_behavior_events(ismember(play_behavior_animal, animal_types{an}),:);

        for bn=1:size(this_animal_behavior,1)
            beh_start = this_animal_behavior(bn,1);
            beh_end = this_animal_behavior(bn,2);
            time_index = beh_time>=beh_start & beh_time<=beh_end;
            behavior_freq(an,time_index) = 1; %% if the time bin includes the behavior, is equal to 1
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

    all_behavior_timeplots_conv{bf} = conv_behavior;
    all_behavior_timeplots{bf}      = behavior_freq;
    data2segment = mean(conv_behavior);
    % data2segment(data2segment>0)=1;
    L               = bwlabeln(data2segment); %% start and end of the event
    binned_beh      = L>0;

    events_start    = find(diff(binned_beh)==1);
    event_end       = find(diff(binned_beh)==-1);
    event_end       = [event_end max(find(binned_beh>0))];
    event_end       = unique(event_end);
    bout_lengths    = beh_time(event_end)-beh_time(events_start);

    sorted_behavior_start = Behavior;
    [~,beh_order]   = sort(Behavior.Start);
    sorted_behavior_start = sorted_behavior_start(beh_order,:);
    n_behavior = nan(numel(bout_lengths),1);
    sorted_behavior_end    = Behavior;
    [~,beh_order]   = sort(Behavior.End);
    sorted_behavior_end = sorted_behavior_end(beh_order,:);
    original_start_end = beh_time([events_start' event_end']);
    bout_start_end  = nan(numel(bout_lengths),2);

    for bn=1:numel(bout_lengths)
        bout_start = original_start_end(bn,1);
        bout_end = original_start_end(bn,2);
        first_play_behavior         = min(find(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<=bout_end & ismember(sorted_behavior_start.Type2, play_behaviors)));
        last_play_behavior          = max(find(sorted_behavior_end.End<=bout_end & sorted_behavior_end.End>=bout_start & ismember(sorted_behavior_end.Type2, play_behaviors)));
        first_play_behavior_type    = sorted_behavior_start.Type2(first_play_behavior);
        last_play_behavior_type     = sorted_behavior_end.Type2(last_play_behavior);
        bout_start                  = sorted_behavior_start.Start(first_play_behavior);
        bout_end                    = sorted_behavior_end.End(last_play_behavior);

        bout_start_end(bn,:)=[bout_start bout_end];


        nbehaviors_within = sorted_behavior_start.Type2(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<bout_end);
        n_behavior(bn) = numel(nbehaviors_within );
        players_within  = sorted_behavior_start.Animal(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<bout_end);
        mutual = numel(unique(players_within))>1;
        nbehaviors_within_SE =  [sorted_behavior_start.Start((sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<bout_end)) ...
            sorted_behavior_start.End(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<bout_end)];
        nbehaviors_within_SE(nbehaviors_within_SE(:,2)>=bout_end,2) = bout_end;
        nbehaviors_within_SE(nbehaviors_within_SE(:,1)<=bout_start,1) = bout_start;
        nbehaviors_within_length =  nbehaviors_within_SE(:,2) - nbehaviors_within_SE(:,1);
        this_partner = [];
        start_search = 0;
        while isempty(this_partner)
            this_session = find(animal_sessions(:,1)-start_search<=bout_start & animal_sessions(:,2)+start_search>=bout_end);
            this_partner = partner_names(this_session);
            start_search = start_search+1;
        end
        merged_animal = players_within{1};
        marged_behaviors = nbehaviors_within{1};
        merged_behaviors_length = num2str(nbehaviors_within_length(1));
        for mb = 2:numel(nbehaviors_within)
            merged_animal           = [merged_animal,'/',players_within{mb}];
            marged_behaviors        = [marged_behaviors,'/', nbehaviors_within{mb}];
            merged_behaviors_length = [merged_behaviors_length, '/', num2str(nbehaviors_within_length(mb))];
        end
        n_beh_bout = sum(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<=bout_end);
        n_play_beh_bout = sum(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<=bout_end & ismember(sorted_behavior_start.Type2, play_behaviors));

        first_behavior = min(find(sorted_behavior_start.Start>=bout_start & sorted_behavior_start.Start<bout_end));
        first_behavior_type = sorted_behavior_start.Type2(first_behavior);
        last_behavior  = max(find(sorted_behavior_end.Start<bout_end & sorted_behavior_end.Start>=bout_start ));
        last_behavior_type = sorted_behavior_end.Type2(last_behavior);

        raw2add = [num2cell([bout_start bout_end])...
            [first_behavior_type last_behavior_type] ...
            num2cell([strcmp(sorted_behavior_start.Animal(first_behavior),repeated_animal) strcmp(sorted_behavior_end.Animal(last_behavior),repeated_animal)]) ...
            num2cell([sorted_behavior_start.Start(first_behavior) sorted_behavior_end.Start( last_behavior)]-bout_start) ...
            num2cell(bout_end - [sorted_behavior_start.End(first_behavior) sorted_behavior_end.End(last_behavior)]) ...
            [first_play_behavior_type  last_play_behavior_type] ...
            num2cell([strcmp(sorted_behavior_start.Animal(first_play_behavior),repeated_animal) strcmp(sorted_behavior_end.Animal(last_play_behavior),repeated_animal)]) ...
            num2cell([sorted_behavior_start.Start(first_play_behavior) sorted_behavior_end.Start( last_play_behavior)]-bout_start) ...
            num2cell(bout_end - [sorted_behavior_start.End(first_play_behavior) sorted_behavior_end.End(last_play_behavior)]) ,...
            num2cell([n_beh_bout n_play_beh_bout bout_end-bout_start bn bf numel(bout_lengths) mutual]) ...
            marged_behaviors,merged_behaviors_length,merged_animal,...
            table2cell(beh_meta_data(bf,:)), this_partner];
        this_session_first_last_table = [this_session_first_last_table;raw2add];

    end


    raw2interplay = [];


    this_inter_play_bout_table = cell(size(bout_start_end,1)-1,5);
    for pb_n =1:size(bout_start_end,1)-1

        bout_end = bout_start_end(pb_n,2);
        next_bout_start =  bout_start_end(pb_n+1,1);

        behavior_list = Behavior.Start>=bout_end &  Behavior.End <=next_bout_start;
        this_inter_play_bout_table{pb_n,1} = Behavior.Type2(behavior_list);
        this_inter_play_bout_table{pb_n,2} = Behavior.Length(behavior_list);
        this_inter_play_bout_table{pb_n,3} = next_bout_start-bout_end;
        this_inter_play_bout_table{pb_n,4} = numel(Behavior.Type2(behavior_list));
        this_inter_play_bout_table{pb_n,5}= bf;
    end
    inter_play_bout_table = [inter_play_bout_table;this_inter_play_bout_table];
    if restrain2morethan1beh
        index = n_behavior>1;
    else
        index = true(size(bout_start_end,1),1);
    end
    bout_lengths = bout_start_end(index,2) - bout_start_end(index,1);
    bout_intervals = bout_start_end(2:end,1)-bout_start_end(1:end-1,2);

    play_bout_corr_pre = [play_bout_corr_pre;[bout_lengths(1:end-1),bout_intervals]];
    play_bout_corr_post = [play_bout_corr_post;[bout_lengths(2:end),bout_intervals]];
    play_bout_corr_pre_zs = [play_bout_corr_pre_zs;zscore([bout_lengths(1:end-1),bout_intervals])];
    play_bout_corr_post_zs = [play_bout_corr_post_zs;zscore([bout_lengths(2:end),bout_intervals])];


    bout_length_matrix(bf,:) = histcounts(bout_lengths,hist_edges);
    inter_bout_matrix(bf,:) = histcounts(bout_intervals,hist_edges);
    beh_type_list           = unique(Behavior.Type2)';

    range2extract =    [-1 1]*round(max(bout_lengths)*1.25);
    beh_time2extract = range2extract(1):beh_bin:range2extract(2);
    behavior_matrixes = zeros(2,numel(beh_type_list), size(bout_start_end,1), numel(beh_matrix_indexes));
    for bout_n = 1:size(bout_start_end)

        bout_start  = bout_start_end(bout_n,1);
        bout_end    = bout_start_end(bout_n,2);


        for bt = 1:numel(beh_type_list)
            this_beh_index = ismember(Behavior.Type2,beh_type_list(bt));
            bool_index_start = any(bout_start+beh_matrix_indexes>=Behavior.Start(this_beh_index) & bout_start+beh_matrix_indexes<=Behavior.End(this_beh_index));
            behavior_matrixes(1,bt,bout_n,bool_index_start) = 1;

            bool_index_end = any(bout_end+beh_matrix_indexes>=Behavior.Start(this_beh_index) & bout_end+beh_matrix_indexes<=Behavior.End(this_beh_index));
            behavior_matrixes(2,bt,bout_n,bool_index_end) = 1;
        end
    end

    if bf ==1
        ALL_behavior_matrices.matrix = behavior_matrixes;
        ALL_behavior_matrices.behavior_types = beh_type_list;
    else

        ALL_behavior_matrices(bf).matrix = behavior_matrixes;
        ALL_behavior_matrices(bf).behavior_types = beh_type_list;
    end


    for pn =1:numel(partner_names)
        figure('units','normalized','outerposition',[0 .25 1 .75])
        t = tiledlayout(1,1);
        ax1 = axes(t);

        start_index = bout_start_end(:,1)>=animal_sessions(pn,1) & ...
            bout_start_end(:,1)<=animal_sessions(pn,2);
        end_index = bout_start_end(:,2)>=animal_sessions(pn,1) & ...
            bout_start_end(:,2)<=animal_sessions(pn,2);
        this_session_index = start_index |end_index ;
        bout_lengths = bout_start_end(this_session_index & index,2)-bout_start_end(this_session_index & index,1);
        this_session_index= find(this_session_index);
        bout_intervals = bout_start_end(this_session_index(2:end),1)-bout_start_end(this_session_index(1:end-1),2);
        plot(beh_time,mean([behavior_freq(1,:);behavior_freq(pn+1,:)]), 'k' )
        hold on
        mean2plot = mean([conv_behavior(1,:);conv_behavior(pn+1,:)]);
        mean2plot(mean2plot>0)=1; %% Conv equals 1
        plot(beh_time,mean2plot, 'b' )
        y_lim = ylim;
        if ~isempty(this_session_index)
            plot([bout_start_end(this_session_index,1)';bout_start_end(this_session_index,1)'], y_lim, 'g')
            plot([bout_start_end(this_session_index,2)';bout_start_end(this_session_index,2)'], y_lim, 'r')

        end
        xlim(animal_sessions(pn,:))
        ylabel([partner_names{pn},' ',beh_meta_data.Animal{bf}])

        xticks((bout_start_end(:,1) +bout_start_end(:,2))/2)
        xticklabels(strsplit(num2str(bout_lengths'), ' '));
        xtickangle(90)
        xlabel('Play-bout duration (s)')
        set(gca,'box','off')
        set(gca,'TickDir','out')
        ax1.XAxisLocation = "bottom";

        ax2= axes(t);
        ax2.Color = 'none';
        ax2.XAxisLocation = 'top';
        ax2.YTick = [];
        xlim(animal_sessions(pn,:))
        xticks((bout_start_end(this_session_index(2:end),1)+bout_start_end(this_session_index(1:end-1),2))/2)
        xticklabels(strsplit(num2str(bout_intervals'), ' '));
        xtickangle(90)
        xlabel('Inter play-bout interval duration (s)')
        set(gca,'TickDir','out')
        set(gca,'box','off')

        partner_specific_comb_bout_interval    = [partner_specific_comb_bout_interval;histcounts(bout_intervals,hist_edges)];
        partner_specific_comb_bout_length = [partner_specific_comb_bout_length;histcounts(bout_lengths,hist_edges)];
        animal_comb = [animal_comb;{partner_names{pn},beh_meta_data.Animal{bf}}];
        pause(.01)
        if save_bool
            saveas(gcf, [saving_folder, '\Play-bout distribution ', partner_names{pn},' ',beh_meta_data.Animal{bf}, ' ', num2str(beh_meta_data.Session(bf)), '.jpg'])
            saveas(gcf, [saving_folder, '\Play-bout distribution ', partner_names{pn},' ',beh_meta_data.Animal{bf}, ' ', num2str(beh_meta_data.Session(bf)), '.svg'])
        end

        close gcf
    end



    figure('units','normalized','outerposition',[0 .25 1 .75])
    t = tiledlayout(1,1);
    ax1 = axes(t);
    hold on
    for bn = 1:size(Behavior,1)

        beh_start   = Behavior.Start(bn);
        beh_end     =  Behavior.End(bn);
        animal_id   = find(ismember(animal_types,Behavior.Animal(bn)));
        if ismember(Behavior.Type2(bn),play_behaviors)
            fill([beh_start beh_end beh_end beh_start ], [-.5 -.5 .5 .5] + animal_id, 'r', 'EdgeColor',[.5 0 0])
        else
            fill([beh_start beh_end beh_end beh_start ], [-.5 -.5 .5 .5] + animal_id, 'k', 'EdgeColor',[.2 .2 .2])
        end
    end
    for bn =1:size(bout_start_end,1)

        plot(bout_start_end(bn,1)*[1 1], [.5 numel(animal_types)+.5], 'b')
        plot(bout_start_end(bn,2)*[1 1], [.5 numel(animal_types)+.5], 'b')
    end
    repeated_animal_index = ismember(Behavior.Animal,animal_types(1));

    repeated_animal_beh= [Behavior.Start(repeated_animal_index) Behavior.End(repeated_animal_index)];
    xticks((repeated_animal_beh(:,2)+repeated_animal_beh(:,1))/2)
    xticklabels(Behavior.Type2(repeated_animal_index))
    xtickangle(90)
    x_lim = xlim;

    ax2= axes(t);
    ax2.Color = 'none';
    ax2.XAxisLocation = 'top';
    ax2.YTick = [];
    xlim(x_lim)
    xtickangle(90)


    partner_animal_beh= [Behavior.Start(~repeated_animal_index) Behavior.End(~repeated_animal_index)];
    xticks((partner_animal_beh(:,2)+partner_animal_beh(:,1))/2)
    xticklabels(Behavior.Type2(~repeated_animal_index))

    for pn=1:size(animal_sessions,1)
        axes(ax1)
        xlim(animal_sessions(pn,:)+[-10 10])
        axes(ax2)
        xlim(animal_sessions(pn,:)+[-10 10])
        if save_bool
            print(gcf,'-vector','-dsvg',[saving_folder, '\Play-bout sheme ', partner_names{pn},' ',beh_meta_data.Animal{bf}, ' ', num2str(beh_meta_data.Session(bf)), '.svg'])
            saveas(gcf, [saving_folder, '\Play-bout sheme ', partner_names{pn},' ',beh_meta_data.Animal{bf}, ' ', num2str(beh_meta_data.Session(bf)), '.jpg'])
        end
    end

    close gcf

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% transform data in tables %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animal_comb = cell2table(animal_comb);
animal_comb.Properties.VariableNames = {'Player1','Player2'};


play_bout_table =  cell2table(this_session_first_last_table);

play_bout_table.Properties.VariableNames = {'BoutStart','BoutEnd',...
    'AnyBehTypeFirst','AnyBehTypeLast','AnyBehAnimalFirst','AnyBehAnimalLast',...
    'AnyBehBegLatencyFirst','AnyBehBegLatencyLast','AnyBehEndLatencyFirst','AnyBehEndLatencyLast',...
    'PlayTypeFirst','PlayTypeLast','PlayAnimalFirst','PlayAnimalLast',...
    'PlayBegLatencyFirst','PlayBegLatencyLast','PlayEndLatencyFirst','PlayEndLatencyLast',...
    '#Beh','#PlayBeh','BoutLength','Bout#','Sess#','NBouts','IsMutual', ...
    'BehaviorList', 'BehaviorListLength','BehaviorListAnimalId', 'Animal','Session','Implanted','Partner'};
%  A lot of the data here is not used anymore


inter_play_bout_table = cell2table(inter_play_bout_table);
inter_play_bout_table.Properties.VariableNames = {'BehaviorList','BehaviirLength','IPBLength','#Beh','SessionN'};

%% normalize bout lengh distribution and smooth


normalzed_bout_length_matrix = partner_specific_comb_bout_length;
normalzed_inter_bout_matrix = partner_specific_comb_bout_interval;


for j=1:size(normalzed_bout_length_matrix,1)
    normalzed_bout_length_matrix(j,:) = normalzed_bout_length_matrix(j,:)/sum(normalzed_bout_length_matrix(j,:));
    normalzed_bout_length_matrix(j,:) = movmean( normalzed_bout_length_matrix(j,:), .5/beh_bin);
end


for j=1:size(normalzed_inter_bout_matrix,1)
    normalzed_inter_bout_matrix(j,:) = normalzed_inter_bout_matrix(j,:)/sum(normalzed_inter_bout_matrix(j,:));
    normalzed_inter_bout_matrix(j,:) = movmean( normalzed_inter_bout_matrix(j,:), .5/beh_bin);
end

%% plot distribution of play bout length and inter play bout interval together with double exponential fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% IMPORTANT: the fiting may not converge the first time %%%%
%%%% YOu may need to run it multiple times                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_lim = [0 10];

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)

hold on
plot(hist_edges_centers,normalzed_bout_length_matrix, ':k' )
plot(hist_edges_centers,mean(normalzed_bout_length_matrix,"omitnan"), 'r', 'LineWidth',2 )
title('Play-bout duration')
ylabel('Prob')
xlabel('Time (s)')
xlim(x_lim)

subplot(2,2,2)
plot(hist_edges_centers,normalzed_inter_bout_matrix,':k' )
hold on
plot(hist_edges_centers,mean(normalzed_inter_bout_matrix,"omitnan") , 'r', 'LineWidth',2 )
title('Inter play-bout interval')
ylabel('Prob')
xlabel('Time (s)')
xlim(x_lim)

playing_animals2 =unique(animal_comb.Player2)';
playing_animals1 =unique(animal_comb.Player1)';

all_animal = unique([animal_comb.Player1;animal_comb.Player2])';

mean_animal_inter_bout = nan(numel(all_animal),size(normalzed_inter_bout_matrix,2));
mean_animal_bout_length = nan(numel(all_animal),size(normalzed_bout_length_matrix,2));

for j=1:numel(all_animal)
    this_animal = ismember(animal_comb.Player2, all_animal{j}) | ...
        ismember(animal_comb.Player1, all_animal{j});

    mean_animal_inter_bout(j,:) = mean(normalzed_inter_bout_matrix(this_animal,:));
    mean_animal_bout_length(j,:) = mean(normalzed_bout_length_matrix(this_animal,:));
end

animals2omit = round(std(mean_animal_bout_length,[],2),10)==0;

subplot(2,2,3)
hold on
plot(hist_edges_centers,mean_animal_bout_length(~animals2omit,:), 'k:' , 'HandleVisibility','off')
plot(hist_edges_centers,mean(mean_animal_bout_length(~animals2omit,:),"omitnan"), 'r', 'LineWidth',2 , 'HandleVisibility','off')
ylabel('Prob')
xlabel('Time (s)')

rep_hist_edges_centers = repmat(hist_edges_centers,sum(~animals2omit));
animals2include  =mean_animal_bout_length(~animals2omit,:);
no_nan = ~isnan(animals2include);
starting_points = rand(1,5);
starting_points(1) = 2.8;

fit_hs_exp = fit_exp_heaviside3(rep_hist_edges_centers(no_nan),animals2include(no_nan),starting_points );
x = fit_hs_exp.fited_param;
t = 0:0.001:20;

hs_function = heaviside(t-x(1));

F =  (x(4)*exp(-t.*(hs_function*(x(3)-x(2)) +x(2)))  + x(5) + hs_function*(x(4)*(exp(-x(1)*x(2))-exp(-x(1)*x(3)))));

hold on
plot(t, F, ':b', 'LineWidth',2 )
legend(['lambda 1 =', num2str(round(x(2),3)),  ' lambda 2 =', num2str(round(x(3),3))])
y_lim = ylim;

plot([x(1) x(1)], y_lim, 'b', 'HandleVisibility','off')
x_ticks = xticks;
x_ticks = [x_ticks x(1)];
x_ticks = sort(x_ticks);
xticks(x_ticks)
ylabel('Prob')
xlabel('Time (s)')
xlim(x_lim)

subplot(2,2,1)
hold on
plot(t, F, ':b', 'LineWidth',2 )
legend(['lambda 1 =', num2str(round(x(2),3)),  ' lambda 2 =', num2str(round(x(3),3))])
y_lim = ylim;
plot([x(1) x(1)], y_lim, 'b', 'HandleVisibility','off')
x_ticks = xticks;
x_ticks = [x_ticks x(1)];
x_ticks = sort(x_ticks);
xticks(x_ticks)
ylabel('Prob')
xlabel('Time (s)')
xlim(x_lim)

subplot(2,2,4)
plot(hist_edges_centers,mean_animal_inter_bout(~animals2omit,:),'k:' , 'HandleVisibility','off')
hold on
plot(hist_edges_centers,mean(mean_animal_inter_bout(~animals2omit,:),"omitnan") , 'r', 'LineWidth',2 , 'HandleVisibility','off')
y_lim = ylim;
plot([3.55 3.55], y_lim, 'b')
x_ticks = xticks;
x_ticks = [x_ticks 3.55];
x_ticks = sort(x_ticks);
xticks(x_ticks)
xlim(x_lim)

