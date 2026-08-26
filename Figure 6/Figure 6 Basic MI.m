
two_animals_data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\TWO ANIMALS';

%% LOAD DATA
% load the corresponding data to plot play or non-play related mutual
% information plots.

figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure Mutual informatin inputs';

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Theta psth';

disp('loading')
load([saving_folder,'\behavior_structure_20bins_delta.mat'],'behavior_structure');
load([saving_folder,'\mi_structure_20bins_delta.mat'],'mi_structure');
load([saving_folder,'\animal_names_mi_structure_20bins_delta.mat'],'animal_names');

% for all non playful behavior
% load([saving_folder,'\mi_structure_20bins_delta_no_play.mat'],'mi_structure');
% load([saving_folder,'\animal_names_mi_structure_20bins_delta_no_play.mat'],'animal_names');

% for aggresive behaviors
% load([saving_folder,'\mi_structure_20bins_delta_aggression.mat'],'mi_structure');
% load([saving_folder,'\animal_names_mi_structure_20bins_delta_aggression.mat'],'animal_names');
%%   mergind data

corr_win_size_sec         = .5;
baseline4corr             = [-1 0];
baseline4mi               = [-2 0];
overal_correlation_time     = [-1 2];

overal_correaltions = [];

sr = mean(diff(mi_structure(1).psth_time));

corr_wind_size = round(corr_win_size_sec/sr);

all_psth_onset          = [];
all_pb_length           = [];
mean_delta_this_animal  = [];
all_mi_t                = [];
all_mi_bc               = [];
all_mi_pct              = [];
all_mi_pct2             = [];
all_mi_zscored          = [];
global_mi_distr_zscored = [];
global_mi               = [];
global_mi_distr         = [];
global_mi_pctl          = [];
power_correlations      = [];
bc_power_correlations   = [];
animal_pair_index       = [];
all_animals_index       =  [];
session_index           = [];


animal_number = [1 2 3 4 4 5 5];

for ns=1:numel(mi_structure)

    disp(ns)
    psth_time = mi_structure(ns).psth_time;

    psth_zscored = mi_structure(ns).psth_zscored;

    power4cor = psth_zscored(:,:,psth_time>=overal_correlation_time(1) & psth_time<=overal_correlation_time(2));
    delta_1 = (power4cor(1,:,:));
    delta_1 = delta_1(:);   

    delta_2 = (power4cor(2,:,:));
    delta_2 = delta_2(:);

    no_nan = ~isnan(delta_1) & ~isnan(delta_2);
    [c,p] = corr(delta_1(no_nan),delta_2(no_nan), 'Type','Spearman');

    
    overal_correaltions = [overal_correaltions; [c,p]];

    corr_matrix = nan(size(squeeze(psth_zscored(1,:,:))));

    for trial=1:size(corr_matrix,1)
        for t=corr_wind_size:size(corr_matrix,2)

            corr_matrix(trial,t-round(corr_wind_size/2)) = corr(squeeze(psth_zscored(1,trial,(t-corr_wind_size+1):t)),squeeze(psth_zscored(2,trial,(t-corr_wind_size+1):t)), 'Type','Spearman');
        end
    end

    power_correlations = [power_correlations;corr_matrix];

    animal_pair_index = [animal_pair_index;ones(size(corr_matrix,1),1)*animal_number(ns)];

    all_animals_index = [all_animals_index;[ones(size(corr_matrix,1),1)*animal_number(ns);-ones(size(corr_matrix,1),1)*animal_number(ns)]];

    session_index = [session_index;[ones(size(corr_matrix,1),1)*ns]];


    global_mi               = [global_mi;mi_structure(ns).mi_global];
    global_mi_distr         = [global_mi_distr;mi_structure(ns).MI_rand_global'];

    all_mi_zscored          = [all_mi_zscored;(mi_structure(ns).mi_global-mean(mi_structure(ns).MI_rand_global))/std(mi_structure(ns).MI_rand_global)];
    global_mi_distr_zscored = [global_mi_distr_zscored;zscore(mi_structure(ns).MI_rand_global)'];

    this_pctl               = sum(mi_structure(ns).MI_rand_global>mi_structure(ns).mi_global)/numel(mi_structure(ns).MI_rand_global);
    global_mi_pctl          = [global_mi_pctl;this_pctl];
    
    PlayBout_table = mi_structure(ns).play_bouts_table;
    pb_length = diff(PlayBout_table')';

    bc_corr_matrix = corr_matrix;
    for j=1:numel(pb_length)
        bc_corr_matrix(j,:) =bc_corr_matrix(j,:) - mean(bc_corr_matrix(j,psth_time>=baseline4corr(1) & psth_time<baseline4corr(2)));
        psth_zscored(1,j, psth_time>pb_length(j)) = NaN;
        psth_zscored(2,j, psth_time>pb_length(j)) = NaN;
    end

        bc_power_correlations = [bc_power_correlations;bc_corr_matrix];

    mean_delta_this_animal = [mean_delta_this_animal;[mean(squeeze(psth_zscored(1,:,:)),'omitmissing');mean(squeeze(psth_zscored(2,:,:)),'omitmissing')]];
    psth_zscored = [squeeze(psth_zscored(1,:,:));squeeze(psth_zscored(2,:,:))];
    all_psth_onset = [all_psth_onset;psth_zscored];

    all_pb_length = [all_pb_length;[pb_length;pb_length]];

    mi_time = mi_structure(ns).mi_time;
    this_mi_t = mi_structure(ns).mi_t';
    all_mi_t = [all_mi_t;mi_structure(ns).mi_t'];
    baseline_range_index = mi_time>=baseline4mi(1) & mi_time<=baseline4mi(2);
    this_mi_t_bc = (this_mi_t - mean(this_mi_t(baseline_range_index)))/std(this_mi_t(baseline_range_index));
    all_mi_bc = [all_mi_bc;this_mi_t_bc];
    all_mi_pct = [all_mi_pct;mi_structure(ns).mi_t_pctl'];
    mi_t_pctl_2 =this_mi_t;
    MI_t_rand_play = mi_structure(ns).MI_t_rand_play;
    for t=1:numel(this_mi_t)
    mi_t_pctl_2(t) = sum(MI_t_rand_play>this_mi_t(t))/numel(MI_t_rand_play);
    end

     all_mi_pct2 = [all_mi_pct2;mi_t_pctl_2];

end

 
%% Fig 6d: distribution of correlation values

figure
x = ones(size(overal_correaltions,1),1) + (rand(size(overal_correaltions,1),1)-.5)*.5;
significan_index = overal_correaltions(:,2)<=0.05;
plot(x,overal_correaltions(:,1), '.k')


plot(x(significan_index),overal_correaltions(significan_index,1), '.r')
hold on
plot(x(~significan_index),overal_correaltions(~significan_index,1), '.k')

xlim([0 2])

title('Power correlation')
%% plot mean power psth together with imagesc (not used in figures)
y_lim = [-1 2.5];
 c_lim = [-2 2]
x_lim = mi_structure(ns).global_time_range;
min_length = 0.0;
baseline = [x_lim(1) 0];
% smoth_wind = round((1/3.5)/mean(diff(mi_structure(1).psth_time))); %halfe delta cycle
smoth_wind = 1; %no smoothing
stacked_mean_power = [];
for ns =1:6

   

    this_play_bout = mi_structure(ns).play_bouts_table;
  

    these_lengths = diff(this_play_bout')';
    [sorted_lengths, order] = sort(these_lengths);

    order = order(sorted_lengths>min_length);
    sorted_lengths = sorted_lengths(sorted_lengths>min_length);
    psth_zscored    = mi_structure(ns).psth_zscored;
     psth_time = mi_structure(ns).psth_time;
    for j=1:numel(these_lengths)
        psth_zscored(1,j, :) = (psth_zscored(1,j, :) - mean(psth_zscored(1,j, psth_time<baseline(2) & psth_time>baseline(1))))/std(psth_zscored(1,j, psth_time<baseline(2) & psth_time>baseline(1)));
        psth_zscored(1,j, :)  = movmean( psth_zscored(1,j, :) ,smoth_wind);
        % psth_zscored(1,j, psth_time>these_lengths(j)) = NaN;
        psth_zscored(2,j, :) = (psth_zscored(2,j, :) - mean(psth_zscored(2,j,psth_time<baseline(2) & psth_time>baseline(1))))/std(psth_zscored(2,j, psth_time<baseline(2) & psth_time>baseline(1)));
        psth_zscored(2,j, :)  = movmean( psth_zscored(2,j, :) ,smoth_wind);

        % psth_zscored(2,j, psth_time>these_lengths(j)) = NaN;
    end
    


    psth_time       = mi_structure(ns).psth_time;


    figure('units','normalized','outerposition',[0 0 .4 .5]);
    subplot(3,2,[1 3])
    imagesc(psth_time, 1:size(psth_zscored,2), squeeze(psth_zscored(1,order,:)))
    axis xy
    xlim(x_lim)
    hold on
    clim(c_lim)
    plot(sorted_lengths, 1:size(sorted_lengths,1), 'w')
    plot(sorted_lengths*0, 1:size(sorted_lengths,1), 'w')

    subplot(3,2,5)
    [h, ~, ci] = ttest(squeeze(psth_zscored(1,order,:)));
    h(isnan(h)) = 0;
    no_nan = ~any(isnan(ci));
    fill([psth_time(no_nan) fliplr(psth_time(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'k', 'FaceAlpha',.1, 'EdgeColor','none')
    hold on
    plot(psth_time, mean(squeeze(psth_zscored(1,order,:)), 'omitmissing'), 'k')
    stacked_mean_power = [stacked_mean_power;mean(squeeze(psth_zscored(1,order,:)))];
    plot([0 0],y_lim, 'r')
    plot(x_lim,[0 0], ':k')
    xlim(x_lim)
    ylim(y_lim)
    significan_points = psth_time([find(diff([0,h,0])==1)' (find(diff([0,h,0])==-1)-1)']);






    for k=1:size(significan_points,1)
        fill([significan_points(k,:) fliplr(significan_points(k,:))],[.9 .9 1 1]*y_lim(2), 'r', 'FaceAlpha',.2, 'EdgeColor','none')
    end
    

    subplot(3,2,[1 3 ]+ 1)
    imagesc(psth_time, 1:size(psth_zscored,2), squeeze(psth_zscored(2,order,:)))
    axis xy
    xlim(x_lim)
    hold on
       plot(sorted_lengths, 1:size(sorted_lengths,1), 'w')
    plot(sorted_lengths*0, 1:size(sorted_lengths,1), 'w')
    clim(c_lim)

    subplot(3,2,6)
    [h, ~, ci] = ttest(squeeze(psth_zscored(2,order,:)));
    h(isnan(h)) = 0;
    no_nan = ~any(isnan(ci));


    fill([psth_time(no_nan) fliplr(psth_time(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'k', 'FaceAlpha',.1, 'EdgeColor','none')
    hold on
    plot(psth_time, mean(squeeze(psth_zscored(2,order,:)), 'omitmissing'), 'k')
    plot([0 0],y_lim, 'r')
    plot(x_lim,[0 0], ':k')
       significan_points = psth_time([find(diff([0,h,0])==1)' (find(diff([0,h,0])==-1)-1)']);

    for k=1:size(significan_points,1)
        fill([significan_points(k,:) fliplr(significan_points(k,:))],[.9 .9 1 1]*y_lim(2), 'r', 'FaceAlpha',.2, 'EdgeColor','none')
    end
    xlim(x_lim)
    ylim(y_lim)


    figure('units','normalized','outerposition',[0 .4 1 .5]);
    subplot(3,6,[1 7])
 
    imagesc(psth_time, 1:size(psth_zscored,2), squeeze(psth_zscored(1,:,:)))
    axis xy
    xlim(x_lim)
    hold on
    clim(c_lim)
    plot(sorted_lengths, 1:size(sorted_lengths,1), 'w')
    plot(sorted_lengths*0, 1:size(sorted_lengths,1), 'w')
    for comb = 1:5
        subplot(3,6,[1 7 ]+ comb)
        imagesc(psth_time, 1:size(psth_zscored,2), squeeze(psth_zscored(2,randsample(size(psth_zscored,2),size(psth_zscored,2)),:)))
        axis xy
        xlim(x_lim)
        hold on
        plot(sorted_lengths, 1:size(sorted_lengths,1), 'w')
        plot(sorted_lengths*0, 1:size(sorted_lengths,1), 'w')
        clim(c_lim)

        pause(.1)
    end

end



%% FIg 6b: plot a psth together mi and shufled trials
y_lim = [-1 2.5];
 c_lim = [-2 2]
x_lim = mi_structure(ns).global_time_range;
min_length = 0.0;
baseline = [x_lim(1) -1];
% smoth_wind = round((1/3.5)/mean(diff(mi_structure(1).psth_time))); %halfe delta cycle
smoth_wind = 1; %no smoothing
stacked_mean_power = [];

for ns =1:7

    this_play_bout = mi_structure(ns).play_bouts_table;

    these_lengths = diff(this_play_bout')';
    [sorted_lengths, order] = sort(these_lengths);

    order = order(sorted_lengths>min_length);
    sorted_lengths = sorted_lengths(sorted_lengths>min_length);
    psth_zscored    = mi_structure(ns).psth_zscored;
    psth_time = mi_structure(ns).psth_time;
    for j=1:numel(these_lengths)
        psth_zscored(1,j, :) = (psth_zscored(1,j, :) - mean(psth_zscored(1,j, psth_time<baseline(2))))/std(psth_zscored(1,j, psth_time<baseline(2)));
        % psth_zscored(1,j, psth_time>these_lengths(j)) = NaN;
        psth_zscored(2,j, :) = (psth_zscored(2,j, :) - mean(psth_zscored(2,j,psth_time<baseline(2) )))/std(psth_zscored(2,j, psth_time<baseline(2) ));

        % psth_zscored(2,j, psth_time>these_lengths(j)) = NaN;
    end
    psth_time       = mi_structure(ns).psth_time;


    % [power_value,new_order] = sort(mean(mean(psth_zscored([1 2],:,psth_time>0 & psth_time<1),3)));
    [power_value2,new_order] = sort((mean(psth_zscored(2,:,psth_time>0 & psth_time<1),3)));
    power_value1 = squeeze(mean(psth_zscored(1,new_order,psth_time>0 & psth_time<1),3));


  

    power_figure = figure('units','normalized','outerposition',[0 0 1 .5]);

    subplot(1,7,1:2)
    plot(power_value2, 1:numel(power_value2), 'm')
    hold on
    plot(power_value1, 1:numel(power_value1), 'g')   
    legend({'Power anima 1','Power Animal 2'})
    title({'Sorted Power'})
    xlabel('Delta Power')



    mi_figure = figure('units','normalized','outerposition',[0 .5 1 .25]);
    subplot(1,7,1)
    sorted_lengths = these_lengths(new_order);
    imagesc(psth_time, 1:size(psth_zscored,2), squeeze(psth_zscored(1,new_order,:)))
    axis xy
    xlim(x_lim)
    hold on
    clim(c_lim)
    plot(sorted_lengths*0, 1:size(sorted_lengths,1), 'w')



    global_time_range = [-2 4];
    time4global = psth_time>global_time_range(1)  & psth_time<global_time_range(2);

    I1 = squeeze(psth_zscored(1,:,time4global));
    I2 = squeeze(psth_zscored(2,:,time4global));

    mi_global = (image_mutual_info(I1, I2, 50)-mean(mi_structure(ns).MI_rand_global))/std(mi_structure(ns).MI_rand_global);




    subplot(1,7,2)
    imagesc(psth_time, 1:size(psth_zscored,2), squeeze(psth_zscored(2,new_order,:)))
    axis xy
    xlim(x_lim)
    hold on
    plot(sorted_lengths*0, 1:size(sorted_lengths,1), 'w')
    clim(c_lim)

    title({'Sorted ', num2str(mi_global)})




    for comb = 1:5
        figure(mi_figure)
        subplot(1,7,2 +comb)
        rand_order = randsample(size(psth_zscored,2),size(psth_zscored,2));
        imagesc(psth_time, 1:size(psth_zscored,2), squeeze(psth_zscored(1,rand_order,:)))
        axis xy
        xlim(x_lim)
        hold on
        plot(sorted_lengths*0, 1:size(sorted_lengths,1), 'w')


        I1 = squeeze(psth_zscored(1,:,time4global));
        I2 = squeeze(psth_zscored(2,rand_order,time4global));

        mi_global = (image_mutual_info(I1, I2, 50)-mean(mi_structure(ns).MI_rand_global))/std(mi_structure(ns).MI_rand_global);
        xlabel('Time (s)')
        title({'Shuffled ', num2str(mi_global)})

        clim(c_lim)


        figure(power_figure)
          power_value1 = squeeze(mean(psth_zscored(1,rand_order,psth_time>0 & psth_time<1),3));
         subplot(1,7,2 +comb)
         plot(power_value2, 1:numel(power_value2), 'm')
         hold on
         plot(power_value1, 1:numel(power_value1), 'g')
          title({'Shuffled ', num2str(mi_global)})
           xlabel('Delta Power')


        pause(.1)
    end

        figure(mi_figure)
    sgtitle(animal_list(ns).name)

        figure(power_figure)
        sgtitle(animal_list(ns).name)

end

%% save needed figure
% 
% print(gcf,'-vector','-dsvg',[figure_folder,'/mi example and shufled versions example 17_6.svg'])
% print(gcf,'-vector','-dsvg',[figure_folder,'/power  example and shufled versions example 17_6.svg'])

%% Fig 6d: plot global mi distribution and observations


figure
hist_ranges = round([min(global_mi_distr_zscored(:)) max(global_mi_distr_zscored(:))],2)
hold on
for j=1:size(global_mi_distr_zscored,1)
histogram(global_mi_distr_zscored(j,:),hist_ranges(1):0.05:hist_ranges(2),'FaceColor', 'k', 'FaceAlpha',.25, 'EdgeColor', 'none')
end

plot([all_mi_zscored';all_mi_zscored'], [0 50], 'r')
plot(all_mi_zscored, all_mi_zscored*0, '.r', 'MarkerSize',12)
text( all_mi_zscored, all_mi_zscored*0 + 60, strsplit(num2str(round(global_mi_pctl,3)'), ' '),'Color',[1 0 0])

%% (not used any more) plot mean delta power across animals and confident intervals (using ttest)
x_lim = [-1 2];
figure
what2plot  = find(psth_time>=x_lim(1) & psth_time<=x_lim(2));
subplot(5,1,3:5)

plot(psth_time(what2plot),mean_delta_this_animal(:,what2plot), ':k','LineWidth',.05)

hold on
[~,~, ci] = ttest(mean_delta_this_animal(:,what2plot));
fill([psth_time(what2plot) fliplr(psth_time(what2plot))],[ci(1,:), fliplr(ci(2,:))],'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(psth_time(what2plot),mean(mean_delta_this_animal(:,what2plot), 'omitmissing'), 'r', 'LineWidth',2)
plot([psth_time(what2plot(1)) psth_time(what2plot(end))], [0 0], 'k')
xlim(x_lim)
ylim([-.5 1])
%% Fig 6g,h: ploting tiem resolve mi
y_lim = [-4 8];
x_lim = [-1 2];
figure
subplot(5,1,1:3)
colormap(jet(256))
[~, order] = sort(mean(all_mi_bc(:,mi_time>1& mi_time<3),2));

all_mi_bc_smoothed = all_mi_bc;
for j=1:size(all_mi_bc,1)
    all_mi_bc_smoothed(j,:) = movmean(   all_mi_bc_smoothed(j,:),20);
end
imagesc(mi_time,1:size(all_mi_bc,1),all_mi_bc_smoothed(order,:))
hold on
plot([0 0], [1 size(all_mi_bc,1)], 'w')
axis xy
clim([- 1 3])
xlim(x_lim)
yticks(1:size(animal_names,1))
yticklabels(animal_names(:,1))


subplot(5,1,4:5)
plot(mi_time,all_mi_bc_smoothed(order,:), ':k')
hold on
plot(mi_time,mean(all_mi_bc_smoothed(order,:)), 'k', 'LineWidth',2)
[h, ~, ci] = ttest(all_mi_bc_smoothed);

     no_nan = ~any(isnan(ci));
     fill([mi_time(no_nan) fliplr(mi_time(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'k', 'FaceAlpha',.25, 'EdgeColor','none')
     plot(mi_time,mean(all_mi_bc_smoothed, 'omitmissing'), 'k')

 starts = find(diff([0 h]) == 1);
     ends   = find(diff([h 0]) == -1)-1;
     ends(ends==0) = 1;
     % y-coordinates for rectangle height
     y_bottom = y_lim(2)*0.8;
     y_top    = y_lim(2)*.9;

     for k = 1:length(starts)
         xs = [mi_time(starts(k))  mi_time(ends(k)) ...
             mi_time(ends(k))   mi_time(starts(k))];
         ys = [y_bottom y_bottom y_top y_top];

         fill(xs, ys, 'b', 'FaceAlpha', .25, 'EdgeColor', 'none');
     end
     plot([0 0], y_lim, 'b')
     axis tight
     xlim(x_lim)

     ylim(y_lim)
    




%% (not used anymore)  plot mi significance plot (logit)

figure
% plot(mi_time,all_mi_bc, ':k')
subplot(1,2,1)
hold on
[~,~, ci] = ttest(all_mi_bc);
fill([mi_time fliplr(mi_time)],[ci(1,:), fliplr(ci(2,:))],'k', 'FaceAlpha',.2, 'EdgeColor','none')
plot(mi_time,mean(all_mi_bc, 'omitmissing'), 'k', 'LineWidth',2)
ylim([-2 4])
yyaxis right
[~,~, ci] = ttest(mean_delta_this_animal);
fill([psth_time fliplr(psth_time)],[ci(1,:), fliplr(ci(2,:))],'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(psth_time,mean(mean_delta_this_animal, 'omitmissing'), 'r', 'LineWidth',2)
hold on
plot([psth_time(1) psth_time(end)], [0 0], 'k')
xlim([-3 4])
ylim([-1 2])




subplot(1,2,2)
% plot(mi_time,1-all_mi_pct, ':k')
y = (1-all_mi_pct+0.001)/1.0002;
y = log10(y./(1-y));
hold on
plot(mi_time,mean(y, 'omitmissing'), 'k')
[~,~, ci] = ttest(y);
no_nan =all(~isnan(ci));
fill([mi_time(no_nan) fliplr(mi_time(no_nan))],[ci(1,no_nan), fliplr(ci(2,no_nan))],'k', 'FaceAlpha',.2, 'EdgeColor','none')
% ylim([.5 1])
% plot([psth_time(1) psth_time(end)], [0.95 0.95], 'r')
plot([psth_time(1) psth_time(end)], -log10([0.05 0.05]), 'r')

hold on
ylim([-2 4])

xlim([-1 2])

