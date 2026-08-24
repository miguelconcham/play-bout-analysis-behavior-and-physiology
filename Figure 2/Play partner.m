saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\psth power by frequency and behavior';

%%
disp('loading')


load([saving_folder,'\psth_structure_delta_updated.mat'],'psth_structure');
load([saving_folder,'\animal_names_delta_updated.mat'],'animal_names');
disp('ready')
%% merging_psth
smooth_wind = 20;
baseline_range = [-2 0]
animal_label = {'B1D1','B1S3','B2S2','B3D2', 'B4S2', 'B4D4'};
electorde_numner = [1 2];
bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges = psth_structure(1).hist_range;
wrap_range = psth_structure(1).range_time_wrap;
time = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index = time<baseline_range(2) & time>baseline_range(1);
baseline_index_time_wrap = 1:round((abs(wrap_range(1))/bin_size));
pre_time = [-.37 -.07]; ('defined as on the ')
pre_time_index = time<pre_time(2) & time>pre_time(1);
all_psth_onset                  = [];
all_behavior_onset_implanted        = [];
all_behavior_onset_non_implanted    = [];
all_session_onset                   = [];
all_playbout_onset                  = [];



all_psth_onset_behavior         = [];
all_psth_onset_only_playobut    = [];
all_psth_offset         = [];
all_psth_tw             = [];
all_psth_tw_3points     = [];
all_play_bouts          = [];
time_wrap_time          = [(baseline_index_time_wrap*bin_size) + wrap_range(1),linspace(0,1,psth_structure(1).n_bins_time_wrap),1 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];
time_wrap_3_points      = [(baseline_index_time_wrap*bin_size) + wrap_range(1),linspace(0,1-1/psth_structure(1).n_bins_time_wrap,psth_structure(1).n_bins_time_wrap), ...
   linspace(1,2-1/psth_structure(1).n_bins_time_wrap,psth_structure(1).n_bins_time_wrap),2 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];
partner_number = [];
session_index = [];
animal_index = [];
baseline_values = [];
play_per_partner = [];
play_per_partner_per_session = [];
baseline_values_per_session = [];
all_electrodes = [];
animal_index_per_session = [];

for j=1:numel(psth_structure)
    if contains(animal_names{j,1},animal_label)

        this_session_beahvior =  psth_structure(j).Behavior;
        animal_num      = find(cell2mat(cellfun(@(x) contains(animal_names{j,1},x), animal_label, 'UniformOutput',false)));
        electrode_num   =animal_names{j,2};
        animal_name_id = animal_names{j,1};
        animal_name_id = strsplit(animal_name_id, ' ');
        animal_name_id = animal_name_id{end};
        this_animal_playbouts = psth_structure(j).play_bouts_table;
        this_animal_lengths = diff(this_animal_playbouts');
        all_play_bouts = [all_play_bouts;this_animal_playbouts];
        partner_sessions = psth_structure(j).Behavior(ismember(psth_structure(j).Behavior.Animal,'Session_structure'),:);
        partner_sessions(ismember(partner_sessions.Type, 'Tickling'),:) = [];
        [~, loc]      = max(this_animal_playbouts(:,1)'>=partner_sessions.Start & this_animal_playbouts(:,2)'<partner_sessions.End,[],1);
        partner_number = [partner_number;loc'];


        animal_index_per_session= [animal_index_per_session;animal_num];




        this_psth_onset         = psth_structure(j).play_bout_onset;
        this_psth_onset_onlypb  = this_psth_onset;

        this_psth_behavior_implanted     = zeros(size(this_psth_onset));
        this_psth_behavior_non_implanted = zeros(size(this_psth_onset));
        this_playbout_onset              = zeros(size(this_psth_onset));
        this_psth_session                = zeros(size(this_psth_onset))+j;

        animal_index = [animal_index;repmat(animal_num,size(this_psth_onset,1),1)];
        session_index = [session_index;repmat(j,size(this_psth_onset,1),1)];

        for trial=1:size(this_psth_onset,1)
            this_playbout_onset(trial,:) = time>0 & time<this_animal_playbouts(trial,2)-this_animal_playbouts(trial,1);
            this_psth_behavior_implanted(trial,:) = any((time+this_animal_playbouts(trial,1))>this_session_beahvior.Start(ismember(this_session_beahvior.Animal,animal_name_id )) & ...
                (time+this_animal_playbouts(trial,1))<this_session_beahvior.End(ismember(this_session_beahvior.Animal,animal_name_id )));


            this_psth_behavior_non_implanted(trial,:) = any((time+this_animal_playbouts(trial,1))>this_session_beahvior.Start(~ismember(this_session_beahvior.Animal,{animal_name_id, 'Session_structure'} )) & ...
                (time+this_animal_playbouts(trial,1))<this_session_beahvior.End(~ismember(this_session_beahvior.Animal,{animal_name_id, 'Session_structure'} )));
            this_psth_onset(trial,:) = ( this_psth_onset(trial,:) - mean( this_psth_onset(trial,baseline_index), 'Omitmissing'))/std( this_psth_onset(trial,baseline_index), 'Omitmissing');
            this_psth_onset(trial,:) = movmean(this_psth_onset(trial,:), smooth_wind);
            % this_psth_onset_onlypb(trial,:) = this_psth_onset(trial,:);
            this_psth_onset(trial,time> this_animal_lengths(trial)) = NaN;
        end
        baseline_this_session=mean(this_psth_onset(:,pre_time_index),2);
        baseline_values = [baseline_values;baseline_this_session];

        play_per_partner_this_animal = nan(1,3);
        baseline_values_this_session = nan(1,3);
        for pn=1:size(partner_sessions,1)
            session_length = partner_sessions.End(pn) - partner_sessions.Start(pn);
            play_per_partner_this_animal(pn) = sum(this_animal_playbouts(loc==pn,2)-this_animal_playbouts(loc==pn,1))/session_length;

            baseline_values_this_session(pn)  = mean(baseline_this_session(loc==pn));
        end
        play_per_partner = [play_per_partner;repmat(play_per_partner_this_animal,size(this_psth_onset,1),1)];
        play_per_partner_per_session = [play_per_partner_per_session;play_per_partner_this_animal];
        baseline_values_per_session =[baseline_values_per_session;baseline_values_this_session] ;

        all_psth_onset                       = [all_psth_onset; this_psth_onset];
        all_behavior_onset_implanted         = [all_behavior_onset_implanted;this_psth_behavior_implanted];
        all_behavior_onset_non_implanted     = [all_behavior_onset_non_implanted;this_psth_behavior_non_implanted];
        all_psth_onset_only_playobut         = [all_psth_onset_only_playobut; this_psth_onset_onlypb];
        all_session_onset                    = [all_session_onset;this_psth_session];
        all_playbout_onset                   = [all_playbout_onset;this_playbout_onset];


        this_psth_onset     = psth_structure(j).play_bout_onset;
        this_psth_offset    = psth_structure(j).play_bout_offset;
        for trial=1:size(this_psth_offset,1)
            this_psth_offset(trial,:) = ( this_psth_offset(trial,:) - mean( this_psth_onset(trial,baseline_index)))/std( this_psth_onset(trial,baseline_index));
            this_psth_offset(trial,:) = movmean(this_psth_offset(trial,:), smooth_wind);
        end
        all_psth_offset = [all_psth_offset; this_psth_offset];

        this_psth_tw = psth_structure(j).play_bout_tw_this;
        for trial=1:size(this_psth_tw,1)
            this_psth_tw(trial,:) = ( this_psth_tw(trial,:) - mean( this_psth_tw(trial,baseline_index_time_wrap)))/std( this_psth_tw(trial,baseline_index_time_wrap));
        end
        all_psth_tw = [all_psth_tw; this_psth_tw];


        this_psth_tw = psth_structure(j).three_point_tw;
        for trial=1:size(this_psth_tw,1)
            this_psth_tw(trial,:) = ( this_psth_tw(trial,:) - mean( this_psth_tw(trial,baseline_index_time_wrap)))/std( this_psth_tw(trial,baseline_index_time_wrap));
        end
        all_psth_tw_3points = [all_psth_tw_3points; this_psth_tw];



        this_psth_ab = psth_structure(j).animal_behavior_onset;
        for trial=1:size(this_psth_ab,1)
            this_psth_ab(trial,:) = ( this_psth_ab(trial,:) - mean( this_psth_ab(trial,baseline_index_time_wrap)))/std( this_psth_ab(trial,baseline_index_time_wrap));
            this_psth_ab(trial,:) = movmean(this_psth_ab(trial,:), smooth_wind);
            this_psth_ab(trial,time> this_animal_lengths(trial)) = NaN;
        end
        all_psth_onset_behavior = [all_psth_onset_behavior; this_psth_ab];
        all_electrodes = [all_electrodes;electrode_num];
    end
end

play_bout_length = diff(all_play_bouts')';


[sorted_play_bout_length, order] = sort(play_bout_length);

%% plot animal responses

X_lim = [-1 2];
y_lim = [-1 3];
alpha= 0.01;

figure
staked_means = [];
time2use = time;
data_for_paired = [];
for an= 1:numel(animal_label)

    [sorted_play_bout_length, order] = sort(play_bout_length(animal_index==an,:));
    subplot(5,numel(animal_label),(1:numel(animal_label):2*numel(animal_label)) + an-1)
    array = all_psth_onset(animal_index==an,:);
    parthers = partner_number(animal_index==an);
    play_per_partner_this_session  = play_per_partner(animal_index==an,:);
    this_baseline_values = baseline_values(animal_index==an);
    imagesc(time2use,1:numel(sorted_play_bout_length),array(order,:) )
    xlim(X_lim)
    clim([-2 2])
    axis xy
    hold on
    imagesc_y_lim = ylim;
    plot([0 0], imagesc_y_lim, 'w')
    plot(sorted_play_bout_length, 1:numel(sorted_play_bout_length), 'w')


    subplot(5,numel(animal_label),((2*numel(animal_label) + 1):numel(animal_label):4*numel(animal_label)) + an-1)

    [~, ~, ci]  = ttest(array(:,:));
    no_nan = ~any(isnan(ci));
    fill([time2use(no_nan) fliplr(time2use(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(time2use,mean(array(:,:), 'omitmissing'), 'k')

    [~, ~, ci]  = ttest(array(parthers==2,:));
    no_nan = ~any(isnan(ci));
    fill([time2use(no_nan) fliplr(time2use(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'r', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(time2use,mean(array(parthers==2,:), 'omitmissing'), 'r')
    xlim(X_lim)
    ylim(y_lim)

      subplot(5,numel(animal_label),(5*(numel(animal_label)-1)) + an-1)
    expected_orders = {[1 2],[2 1]};
    this_baseline_values = zscore(this_baseline_values);
      swarmchart(cell2mat(expected_orders((play_per_partner_this_session(:,2)>play_per_partner_this_session(:,1) )+1)'),this_baseline_values, 'k.')
      hold on
      play_lengths= mean(play_per_partner_this_session);
      % play_lengths = play_lengths/max(play_lengths);
      this_staked_means = nan(numel(play_lengths),1);
      for j=1:3
          index = partner_number==j;
          if sum(index)>0
              this_staked_means(j) = mean(this_baseline_values(parthers==j));
              plot([.5 1.5] + j-1,[1 1]*mean(this_baseline_values(parthers==j)), 'r')
          end
      end

      staked_means =[staked_means; [play_lengths' this_staked_means,ones(numel(this_staked_means),1)*an (1:3)']];
      data_for_paired = [data_for_paired;play_lengths];
end

%%
figure

subplot(1,3,1)
x_rand = (rand(sum(all_electrodes==1),2)-.5)/2 + ones(sum(all_electrodes==1),2)*diag([1 2]);
plot(x_rand', play_per_partner_per_session(all_electrodes==1, 1:2)', ':k')
hold on
    plot(x_rand', play_per_partner_per_session(all_electrodes==1, 1:2)', '.k', 'MarkerSize', 15)

 p = signrank(play_per_partner_per_session(all_electrodes==1, 1), play_per_partner_per_session(all_electrodes==1, 2));
title(p)

reduced = play_per_partner_per_session(:,2)<play_per_partner_per_session(:,1);
axis square
ylabel('Proportion of play during session')
xticks([1 2])
xticklabels({'Partner 1','Partner 2'})


subplot(1,3,2)
table_data = [play_per_partner_per_session(all_electrodes==1,1)-play_per_partner_per_session(all_electrodes==1,2),baseline_values_per_session(all_electrodes==1,1)-baseline_values_per_session(all_electrodes==1,2), animal_index_per_session(all_electrodes==1)];
table_data = array2table(table_data);
table_data.Properties.VariableNames = {'x','y','Subject'};
lme = fitlme(table_data, 'y ~ 1 + x+ (1|Subject)');

% Fit model
lme = fitlme(table_data, 'y ~ 1 + x + (1|Subject)');

% Create x-grid for smooth prediction
xvals = linspace(min(table_data.x), max(table_data.x), 100)';
tbl_pred = table(xvals, 'VariableNames', {'x'});

% Get predicted mean and confidence intervals
tbl_pred = table(xvals, repmat(table_data.Subject(1), numel(xvals), 1), ...
    'VariableNames', {'x','Subject'});
[yhat, yCI] = predict(lme, tbl_pred, 'Conditional', false);% 'Conditional', false => marginal over random effects (population-level line)

% Plot data
 hold on
scatter(table_data.x, table_data.y, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.3)

% Plot confidence band (as shaded area)
fill([xvals; flipud(xvals)], ...
     [yCI(:,1); flipud(yCI(:,2))], ...
     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Plot fitted line
plot(xvals, yhat, 'b', 'LineWidth', 2);
xlabel('Play change (Partner 1- Partner 2)')
ylabel('Delta power change (Partner 1- Partner 2)')
title( num2str([lme.anova.FStat lme.anova.pValue]))
axis square

subplot(1,3,3)
no_nan = ~any(isnan(staked_means(:,1:2)),2);
plot(staked_means(:,1), staked_means(:,2), '.k')
hold on
sesions = unique(staked_means(:,3))
for ns=1:numel(sesions)
    x = staked_means(staked_means(:,3) ==sesions(ns),1);
    y = staked_means(staked_means(:,3) ==sesions(ns),2);
    parter_number = 1:3;

    [~, correct_order] =sort(x);
    parter_number = parter_number(correct_order);
    if ns==5
        plot( x(correct_order),y(correct_order) ,'.r')
        plot( x(correct_order),y(correct_order) ,'r')
        for pn=1:numel(parter_number)
            text(x(correct_order(pn)),y(correct_order(pn)), ['Pt #',num2str(parter_number(pn))], 'Color', 'r')
        end
    elseif ns==6
         plot( x(correct_order),y(correct_order) ,'.b')
        plot( x(correct_order),y(correct_order) ,'b')
        for pn=1:numel(parter_number)
            text(x(correct_order(pn)),y(correct_order(pn)), ['Pt #',num2str(parter_number(pn))], 'Color', 'b')
        end

    else

      plot( x(correct_order),y(correct_order) ,'.k')
       plot( x(correct_order),y(correct_order) ,'k')
       for pn=1:numel(parter_number)
       text(x(correct_order(pn)),y(correct_order(pn)), ['Pt #',num2str(parter_number(pn))])
       end
    end
end

[c,p] = corr(staked_means(no_nan,1), staked_means(no_nan,2), 'Type','Spearman')
title([c p])


% p = ranksum(baseline_values(partner_number==1), baseline_values(partner_number==2))
title(p)

axis square
xlabel('Proportion of play during session')
ylabel('Change in delta power (zscored)')
%%
staked_mean_responses = cell(3,1);
min_length = .0;
pn=1;
for partner_number_list = {1, 2, [1 2 3]};
    figure
    staked_mean_responses{pn} = [];
    for an= 1:numel(animal_label)
        animal_bool = animal_index==an & ismember(partner_number, partner_number_list{1});
        length_bool = play_bout_length>min_length;
        [sorted_play_bout_length, order] = sort(play_bout_length(animal_bool & length_bool,:));
        subplot(5,numel(animal_label),(1:numel(animal_label):2*numel(animal_label)) + an-1)

        array = all_psth_onset(animal_bool & length_bool,:);
        imagesc(time,1:numel(sorted_play_bout_length),array(order,:) )
        xlim(X_lim)
        clim([-2 2])
        axis xy
        hold on
        plot([0 0],[1 numel(sorted_play_bout_length)], 'w')
        plot(sorted_play_bout_length,1:numel(sorted_play_bout_length), 'w')
        title(animal_label{an})

        subplot(5,numel(animal_label),((2*numel(animal_label) + 1):numel(animal_label):5*numel(animal_label)) + an-1)

        [~, ~, ci]  = ttest(array);
        no_nan = ~any(isnan(ci));
        fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
        hold on
        plot(time,mean(array, 'omitmissing'), 'k')
        staked_mean_responses{pn} = [ staked_mean_responses{pn};mean(array, 'omitmissing')];
        xlim(X_lim)
    end
    pn =pn+1;
end

%%

load([saving_folder,'\results_play_bout_PBonly_zscore4_updated.mat'],'results');   %     here we use all play bouts

x_lim = [-1 2];
alpha = 0.05;


limited_time = results.time;


% 
est = results.partner_est;
ci = results.partner_ci ;      
d = results.partner_d  ;       
pvals_fdr = results.partner_pvals      ;

figure;
subplot(2,1,1); hold on;

% Shaded CI
what_to_plot = (staked_mean_responses{2}-staked_mean_responses{1})';
plot(time,what_to_plot, 'k')
no_nan = ~any(isnan(ci),2)
fill([limited_time(no_nan); flipud(limited_time(no_nan))], [ci(no_nan,1); flipud(ci(no_nan,2))], ...
    [0.8 0.8 1], 'EdgeColor','none','FaceAlpha',0.4);
plot(limited_time, est, 'b','LineWidth',2);

% Mark significant bins
sig_idx = pvals_fdr < alpha;
plot(limited_time(sig_idx), est(sig_idx), 'r*','MarkerSize',6);

ylabel('Mean Power');
title('Mixed-Effects Power (CI + FDR-corrected sig)');
grid on;
xlim(x_lim)

subplot(2,1,2);
plot(limited_time, d, 'k','LineWidth',2);
ylabel('Effect size (Cohen''s d)');
xlabel('Time (s)');
grid on;
xlim(x_lim)