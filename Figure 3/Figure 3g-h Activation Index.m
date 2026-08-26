%% Figure 3 — Activation index plots
% Checkpoint from Supporting codes/Estiamte_activation_index_and_other_variables.m.
% Sections 2–4: primary activation-index views (how extreme a cell type’s
% activity is relative to other neurons in the same structure).
% Sections 5–9: alternative summaries of that same within-structure ranking.

%% 1 LOADING
% Load precomputed activation-index and per-neuron PSTHs
% (all_psth_cell, all_activation_order, psth_list, area_list, cell_type_names, warped_time).

checkpoint_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\activation index';
checkpoint_file = [checkpoint_folder, '\activation_index_section6.mat'];

disp('Loading activation-index checkpoint')
load(checkpoint_file);
disp('Checkpoint loaded')


%% 2 Per-neuron PSTHs by cell type and area (play vs non-play). NOT SHOWN IN PAPAER
% Within each structure, neurons are split into peak-locked, trough-locked, and
% unlocked. Each row is one neuron’s z-scored, time-warped PSTH, median across
% play (top) vs non-play (middle) behaviors. Bottom: population mean ± 95% CI
% (play red, non-play black). This is the activity that later sections rank
% against the rest of the same structure. Trough cells are typically suppressed
% during play.
psth_n=1;
% play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25 28 29 30 31];
play_behaviors     = [ 4  6  8  14  18  20  24  5 7 9 15 19 21 25 ];
non_play_behaviors = [10  12  16  22  26   34  11 13 17 23 27 35];

for an = 1:numel(area_list)


    figure('units','normalized','outerposition',[0 0 .33 1]);
    peak_indexes            = find(all_activation_order{psth_n,an,1}==1);
    trough_indexes          = find(all_activation_order{psth_n,an,1}==2);
    non_entrained_indexes   = find(all_activation_order{psth_n,an,1}==3);
    indexes = {peak_indexes,trough_indexes,non_entrained_indexes};
    indexes_names = {'Peak','Trough','NonModulated'}
    for index_n = 1:3


        staked_play_medians      = [];
        staked_non_play_medias  = [];
        for nn = 1:numel(indexes{index_n})
            all_responses = nan(size(all_activation_order,2),size(all_activation_order{1,1,4},2));

            for psth_n = 1:size(all_activation_order,1)
                all_responses(psth_n,:) = all_activation_order{psth_n,an,4}(indexes{index_n}(nn),:);
            end



            staked_play_medians = [staked_play_medians;median(all_responses(play_behaviors,:), 'omitmissing')];
            staked_non_play_medias = [staked_non_play_medias;median(all_responses(non_play_behaviors,:), 'omitmissing')];


        end

        subplot(4,3,1+index_n-1)
        imagesc(warped_time,1:size(staked_play_medians,1),staked_play_medians)
        title(indexes_names{index_n})
        clim([-2 2])
        axis xy

        subplot(4,3,4+index_n-1)
        imagesc(warped_time,1:size(staked_non_play_medias,1),staked_non_play_medias)
        axis xy
        clim([-2 2])

        subplot(4,3,[7 10]+index_n-1)

        hold on
        matrix2plot = staked_non_play_medias;
        [~, ~, ci] =ttest(matrix2plot);
        no_nan = ~any(isnan(ci));

        fill([warped_time(no_nan) fliplr(warped_time(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'k','FaceAlpha', .25, 'EdgeColor', 'none')
        plot(warped_time, mean(matrix2plot, 'omitmissing'), 'k')

        matrix2plot = staked_play_medians;
        [~, ~, ci] =ttest(matrix2plot);
        no_nan = ~any(isnan(ci));

        fill([warped_time(no_nan) fliplr(warped_time(no_nan))],[ci(1,no_nan) ci(2,no_nan)],'r','FaceAlpha', .225, 'EdgeColor', 'none')
        plot(warped_time, mean(matrix2plot, 'omitmissing'), 'r')
        ylim([-1 1])

    end
    sgtitle(area_list(an))
end



%% 3: Fig 3h Activation index, mean ± CI 
% Logit activation index over warped time: how extreme the cell-type median is
% relative to a within-structure null (other neurons in the same area).
% Rows = structures, columns = peak / trough / unlocked. Play red, non-play black.
% Positive = that cell type ranks above the local population; negative = below.
figure
play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25]; % bite and pining excluded due to lower sample size 28 29 30 31
non_play_behaviors = [10 11 12 13 16 17 22 23 26 27  32 33 34 35];
y_lim = [-2 2];
time_analysis = 1:size(all_psth_cell{1,1},2);
for an=1:5
    for j=1:3
        subplot(5,3,3*(an-1) + j )
     
        play_matrix = all_psth_cell{j,an}(play_behaviors,:);
        non_play_matrix = all_psth_cell{j,an}(non_play_behaviors,:);

        for n=1:size(play_matrix,1)
            play_matrix(n,:) = movmean(play_matrix(n,:),1);
        end

        for n=1:size(non_play_matrix,1)
            non_play_matrix(n,:) = movmean(non_play_matrix(n,:),1);
        end
        hold on
        [~, ~, ci] = ttest(non_play_matrix);
        fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.2, 'EdgeColor','none')
         plot(mean(non_play_matrix), 'k')

          [~, ~, ci] = ttest(play_matrix);
        fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
        plot(mean(play_matrix), 'r')
        plot([20 20], y_lim, 'k')
        plot([40 40], y_lim, 'k')
        plot([0 60], [0 0], 'k:')

        hold on
       
        ylim(y_lim)
           if j==1
            ylabel({area_list{an}, 'ActivationIndex'})
           end

           if an==1
               title(cell_type_names{j})
           end

    end
end
%% print if needed
% print(gcf,'-vector','-dsvg',[figure_folder,'/mean and ci activation index.svg'])
%% 4 Fig 3g Activation index per behavior 
% Same logit ranking as section 3, shown as one row per behavior (play then
% non-play). Color is how strongly that cell type outranks (or is outranked by)
% the rest of the structure over time. White line separates play from non-play.

figure
% colormap(cividis)
% play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25 ];
% non_play_behaviors = [10 11 12 13 16 17 22 23 26 27  32 33 34 35];

play_behaviors     = [4 6 8 14 18 20 24 5 7 9 15 19 21 25];
non_play_behaviors = [10 12 16 22 26 34 11 13 17 23 27 35];

y_lim = [-2 2];
time_analysis =( 1:size(all_psth_cell{1,1},2))/4 -5;
for an=1:5
    for j=1:3
        subplot(5,3,3*(an-1) + j )
     
        play_matrix = all_psth_cell{j,an}(play_behaviors,:);
        non_play_matrix = all_psth_cell{j,an}(non_play_behaviors,:);      

       
        hold on
        stacked_matrix = [play_matrix;non_play_matrix];
        imagesc(time_analysis, 1:size(stacked_matrix,1),stacked_matrix)

        if j==1
            yticks(1:numel([non_play_behaviors,play_behaviors]))
            yticklabels(psth_list([play_behaviors,non_play_behaviors]))
        end

        if an==1
            title(cell_type_names{j})
        end

        axis xy
     
        plot([20 20]/4 - 5, [1 size(stacked_matrix,1)], 'k')
        plot([40 40]/4 - 5, [1 size(stacked_matrix,1)], 'k')
        plot([0 60]/4 - 5, [0 0]+numel(play_behaviors), 'w', 'LineWidth',2)

        plot([0 60]/4 - 5, [0 0]+numel(play_behaviors)/2, 'w:', 'LineWidth',2)


        plot([0 60]/4 - 5, [0 0]+numel(play_behaviors) +numel(non_play_behaviors)/2, 'w:', 'LineWidth',2)

        clim([-4 4])
        axis tight
       
    end
end


%% 5 (Not shown in paper) Self–other consistency of the within-structure ranking
% Alternative to the time-resolved index: correlation of activation-index
% waveforms across behaviors. For each cell type × structure, compare
% (1) matched self vs partner of the same behavior vs (2) unmatched behaviors,
% separately for play and non-play. Tests whether a cell type’s ranking inside
% the structure is behavior-specific (higher matched correlation) or a generic
% play-related signature. Title = rank-sum p (matched vs unmatched).

figure
 play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25 ];
 non_play_behaviors = [10 11 12 13 16 17 22 23 26 27  32 33 34 35];
for an=1:5
    for j=1:3
      
     
        play_matrix = all_psth_cell{j,an}(play_behaviors,:);
        non_play_matrix = all_psth_cell{j,an}(non_play_behaviors,:);

        for n=1:size(play_matrix,1)
            play_matrix(n,:) = movmean(play_matrix(n,:),1);
        end

        for n=1:size(non_play_matrix,1)
            non_play_matrix(n,:) = movmean(non_play_matrix(n,:),1);
        end

        C_play = corr(play_matrix');        % size: n_behaviors x n_behaviors
        C_nonplay = corr(non_play_matrix');

        matched_corrs_play = [];
        nonmatched_corrs_play = [];
        n_play = size(play_matrix,1);
        for k = 1:2:n_play
            % matched pair (same behavior: self vs other)
            matched_corrs_play(end+1) = C_play(k, k+1);

            % compare to all OTHER behaviors
            other_idx = setdiff(1:n_play, [k k+1]);
            nonmatched_corrs_play = [nonmatched_corrs_play, C_play(k, other_idx), C_play(k+1, other_idx)];
        end

        n_non = size(non_play_matrix,1);
        matched_corrs_non = [];
        nonmatched_corrs_non = [];
        for k = 1:2:n_non
            matched_corrs_non(end+1) = C_nonplay(k, k+1);
            other_idx = setdiff(1:n_non, [k k+1]);
            nonmatched_corrs_non = [nonmatched_corrs_non, C_nonplay(k, other_idx), C_nonplay(k+1, other_idx)];
        end
        all_var = {matched_corrs_play,nonmatched_corrs_play,matched_corrs_non,nonmatched_corrs_non};
          subplot(5,3,3*(an-1) + j )
          hold on
          for cc_n = 1:numel(all_var)
              x_rand = (rand(size(all_var{cc_n})) - .5)*.5;
              plot(x_rand + cc_n,all_var{cc_n}, '.k', 'MarkerSize',5)
              plot([-.25 .25]+cc_n,[0 0]+ mean(all_var{cc_n}), 'r' )
          end

          p1 = ranksum(all_var{1},all_var{2});
          p2 = ranksum(all_var{3},all_var{4});

          title(num2str([p1 p2]))


        
    end
end




%% 6 (Not shown in paper)  Bout-averaged activation index per behavior
% Alternative summary of the same ranking: mean logit during the bout window
% [0 5] s, one bar per behavior (play magenta, non-play black). Collapses the
% time course from sections 3–4 into a single within-structure rank per
% behavior (positive = cell type above the local population).

figure
colormap(jet)
play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25 28 29 30 31 ];
non_play_behaviors = [10 11 12 13 16 17 22 23 26 27  32 33 34 35];
y_lim = [-2 2];
time_analysis =( 1:size(all_psth_cell{1,1},2))/4 -5;
value_during_behavior = [0 5];
selection_index = time_analysis>=value_during_behavior(1) & time_analysis<=value_during_behavior(2);

for an=1:5

    
    for j=1:3
        subplot(5,3,3*(an-1) + j )
     
        play_matrix = all_psth_cell{j,an}(play_behaviors,:);
        non_play_matrix = all_psth_cell{j,an}(non_play_behaviors,:);


        hold on
        stacked_matrix = [play_matrix;non_play_matrix];
        logit_mean = mean(stacked_matrix(:,selection_index ),2);
        pctg_time_significant_positive = mean(stacked_matrix(:,selection_index )>-log10(0.05),2);
        pctg_time_significant_negative = mean(stacked_matrix(:,selection_index )<log10(0.05),2);
        modulation_effect = pctg_time_significant_positive;

        summary_p = 10.^logit_mean ./ (1 + 10.^logit_mean);
        
        % 
        % barh(1:numel(play_behaviors), summary_p(1:numel(play_behaviors)), 'FaceColor','m', 'EdgeColor','none')
        % barh((1:numel(non_play_behaviors)) + numel(play_behaviors), summary_p((1:numel(non_play_behaviors)) + numel(play_behaviors)), 'FaceColor','k', 'EdgeColor','none')
             
        barh(1:numel(play_behaviors), logit_mean(1:numel(play_behaviors)), 'FaceColor','m', 'EdgeColor','none')
        barh((1:numel(non_play_behaviors)) + numel(play_behaviors), logit_mean((1:numel(non_play_behaviors)) + numel(play_behaviors)), 'FaceColor','k', 'EdgeColor','none')
       
        % 
        % barh(1:numel(play_behaviors), modulation_effect(1:numel(play_behaviors)), 'FaceColor','m', 'EdgeColor','none')
        % barh((1:numel(non_play_behaviors)) + numel(play_behaviors), modulation_effect((1:numel(non_play_behaviors)) + numel(play_behaviors)), 'FaceColor','k', 'EdgeColor','none')
       xlim([-4 4])


    end
end

%% 7 (Not shown in paper)  Fraction of behaviors with a significant ranking
% Alternative to the mean index: how often the cell type’s rank in the structure
% crosses a significance threshold (|logit| > −log10(0.05)) at any time.
% Bars: play up, play down, non-play up, non-play down. Asks whether play
% increases the chance that a cell type sits at the extreme of its area.

figure

alpha_level = 0.05;
figure
behavior_range = 20:40;
for an=1:5
    for j=1:3
        subplot(5,3,3*(an-1) + j )
        play_matrix = all_psth_cell{j,an}(play_behaviors,:);
        non_play_matrix = all_psth_cell{j,an}(non_play_behaviors,:);

        for n=1:size(play_matrix,1)
            play_matrix(n,:) = movmean(play_matrix(n,:),2);
        end

        for n=1:size(non_play_matrix,1)
            non_play_matrix(n,:) = movmean(non_play_matrix(n,:),2);
        end
        pctg_play_001_p = sum(any((play_matrix)>abs(log10(alpha_level)),2))/size(play_matrix,1);
        pctg_non_play_001_p = sum(any((non_play_matrix)>abs(log10(alpha_level)),2))/size(non_play_matrix,1);
         pctg_play_001_n = sum(any((play_matrix)<-abs(log10(alpha_level)),2))/size(play_matrix,1);
        pctg_non_play_001_n = sum(any((non_play_matrix)<-abs(log10(alpha_level)),2))/size(non_play_matrix,1);
        bar([pctg_play_001_p pctg_play_001_n pctg_non_play_001_p pctg_non_play_001_n])
        ylim([0 1])
    end
end


%% 8 (Not shown in paper)  Within-structure relative activity (IRI heatmaps)
% Neuron-level alternative to the population activation index. For each
% structure, subtract the area-wide median PSTH from every neuron (IRI), then
% average those residuals by cell type. Heatmaps: behaviors × time. This is
% how peak / trough / unlocked neurons rank relative to the rest of the same
% structure, without the permutation/logit transform used in sections 3–7.
% The IRI value also ranks the activation of a neuron within a strucutre

% c_lim = [-6 6]
c_lim = [-2 2];

play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25 28 29 30 31 ];
non_play_behaviors = [10 11 12 13 16 17 22 23 26 27  34 35];
figure
colormap(jet)
nAreas = numel(area_list);
for an = 1:nAreas
nBeh = size(all_activation_order,1);
nNeurons = size(all_activation_order{1, an, 4},1);   % or whatever list you're looping over
nT = size(all_activation_order{1,an,4},2);

all_resp = nan(nBeh, nNeurons, nT);

for psth_n = 1:nBeh
    for i = 1:nNeurons
        
        neuron_id = (i);

        all_resp(psth_n, i, :) = all_activation_order{psth_n, an, 4}(neuron_id, :);

    end
end



peak_idx = find(all_activation_order{1,an,1} == 1);
trough_idx = find(all_activation_order{1,an,1} == 2);
non_idx = find(all_activation_order{1,an,1} == 3);
 

IRI= nan(size(all_resp));
pop_median = (median(all_resp,2,'omitnan'));

for psth_n = 1:nBeh
     for i = 1:nNeurons
         IRI(psth_n, i, :) = all_resp(psth_n, i, :) -pop_median(psth_n,1, :);
     end
end


% matrix2plot = all_resp;
matrix2plot = IRI;

trough_IRI =  squeeze(mean(matrix2plot(:,trough_idx,:),2, 'omitmissing'));
peak_IRI =  squeeze(mean(matrix2plot(:,peak_idx,:),2, 'omitmissing'));
unlocked_IRI =  squeeze(mean(matrix2plot(:,non_idx,:),2, 'omitmissing'));


cell_array = {peak_IRI,trough_IRI,unlocked_IRI};

non_play_behaviors = [10 11 12 13 16 17 22 23 26 27 32 33 34 35];
% non_play_behaviors = [10 11]
for type_n=1:3
subplot(nAreas,3,3*(an-1) + type_n)
imagesc(cell_array{type_n}([non_play_behaviors,play_behaviors],:))
clim(c_lim)
axis xy
hold on

plot([1 60], [1 1]*numel(non_play_behaviors), 'w', 'LineWidth',2)
 if an==1
            title(cell_type_names{type_n})
        end
end

subplot(nAreas,3,3*(an-1) + 1)
yticks(1:numel([non_play_behaviors,play_behaviors]))
yticklabels(psth_list([non_play_behaviors,play_behaviors]))

       
end

%% 9 (Not shown in paper) Within-structure relative activity, play vs non-play
% Same IRI as section 8 (neuron minus structure median, then mean by cell type),
% shown as mean ± 95% CI over time. Play red, non-play black. Asks whether a
% cell type’s rank inside the area shifts during play relative to other
% behaviors, using a signed residual instead of the logit activation index.

c_lim = [-2.5 2.5];
play_behaviors     = [ 4 5 6 7 8 9 14 15 18 19 20 21 24 25 28 29 30 31 ];
non_play_behaviors = [10 11 12 13 16 17 22 23 26 27  34 35];
figure
colormap(jet)

for an = 1:nAreas

    nBeh = size(all_activation_order,1);
    nNeurons = size(all_activation_order{1, an, 4},1);
    nT = size(all_activation_order{1,an,4},2);

    % -----------------------------
    % build response tensor
    % -----------------------------
    all_resp = nan(nBeh, nNeurons, nT);

    for psth_n = 1:nBeh
        for i = 1:nNeurons
            all_resp(psth_n, i, :) = all_activation_order{psth_n, an, 4}(i, :);
        end
    end

    % -----------------------------
    % IRI
    % -----------------------------
    pop_median = squeeze(median(all_resp,2,'omitnan'));

    IRI = all_resp - pop_median(:, ones(1,nNeurons), :);

    % -----------------------------
    % cell types
    % -----------------------------
    peak_idx   = find(all_activation_order{1,an,1} == 1);
    trough_idx = find(all_activation_order{1,an,1} == 2);
    non_idx    = find(all_activation_order{1,an,1} == 3);

    matrix2plot = IRI;

    trough_IRI   = squeeze(mean(matrix2plot(:,trough_idx,:),2,'omitnan'));
    peak_IRI     = squeeze(mean(matrix2plot(:,peak_idx,:),2,'omitnan'));
    unlocked_IRI = squeeze(mean(matrix2plot(:,non_idx,:),2,'omitnan'));

    cell_array = {peak_IRI,trough_IRI, unlocked_IRI};

    % -----------------------------
    % GROUPS
    % -----------------------------
 

    for type_n = 1:3

        subplot(nAreas,3,3*(an-1) + type_n)
        hold on

        data = cell_array{type_n};  % behavior × time

        % -----------------------------
        % PLAY
        % -----------------------------
        play_data = data(play_behaviors, :);
        [~,~,ci_play] = ttest(play_data);

        mu_play = mean(play_data,1,'omitnan');

        fill([1:nT fliplr(1:nT)], ...
             [ci_play(1,:) fliplr(ci_play(2,:))], ...
             'r','FaceAlpha',0.25,'EdgeColor','none');

        plot(mu_play,'r','LineWidth',2)

        % -----------------------------
        % NON-PLAY
        % -----------------------------
        non_data = data(non_play_behaviors, :);
        [~,~,ci_non] = ttest(non_data);

        mu_non = mean(non_data,1,'omitnan');

        fill([1:nT fliplr(1:nT)], ...
             [ci_non(1,:) fliplr(ci_non(2,:))], ...
             'k','FaceAlpha',0.25,'EdgeColor','none');

        plot(mu_non,'k','LineWidth',2)

        % -----------------------------
        % formatting
        % -----------------------------
        ylim(c_lim)
        axis xy
        % title(['Type ' num2str(type_n)])

        if an == 1 && type_n == 1
            legend({'play CI','play mean','non-play CI','non-play mean'})
        end

        if type_n==1
            ylabel({area_list{an}, 'ActivationIndex'})
        end

        if an==1
            title(cell_type_names{type_n})
        end

    end

end