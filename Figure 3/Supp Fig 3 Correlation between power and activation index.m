%% Supp Fig 3 — Correlation between delta power and activation index
% Relates two behavior-resolved measures: (i) the cell-type activation
% index (how that cell type ranks inside each structure) and (ii) PAG
% delta-power PSTHs. Behaviors are matched by name and self/other role,
% not by list order, then correlated in time and as bout averages.

%% 1 Load activation-index and delta-power results
% stacked_psth_by_area_celltype.mat : logit activation index per behavior,
%   area, and cell type (Estimate_activation_index...).
% results_full_lfp_per_behavior_time_warpped.mat : time-warped LME beta
%   for delta power, play vs non-play (Analyze_psth_all_behaviors).

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
figure_2_new_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Codes\Figure codes\FIgure 2 delta modulation per behavior';
saving_folder = [data_root, '\Analysis results\phase locking data'];

load([saving_folder,'\results_full_lfp_per_behavior_time_warpped.mat'])
load([saving_folder,'\stacked_psth_by_area_celltype.mat'])   
stacked_labels = saved_data(1,1).labels;

%% 2 Preview the two measures
% Quick look at both datasets before matching. Rows/columns should already
% resemble Fig. 3g (activation index) and the per-behavior delta PSTH.

%%% Activation-index heatmaps (area × cell type)
% Color = logit within-structure rank over warped time. Rows = behaviors
% (play then non-play). Columns in the grid = peak / trough / unlocked.


figure

y_lim = [-2 2];

for an = 1:5
    for j = 1:3

        subplot(5,3,3*(an-1) + j)

        % -----------------------------
        % load precomputed data
        % -----------------------------
        stacked_matrix = saved_data(an,j).stacked_matrix;
        time_analysis_activaton_index  = saved_data(an,j).time;
        labels         = saved_data(an,j).labels;

        play_behaviors     = saved_data(an,j).play_behaviors;
        non_play_behaviors = saved_data(an,j).non_play_behaviors;

        % -----------------------------
        % plot
        % -----------------------------
        imagesc(time_analysis_activaton_index, 1:size(stacked_matrix,1), stacked_matrix)
        axis xy
        hold on

        % -----------------------------
        % titles
        % -----------------------------
        if j == 1
            yticks(1:numel(labels))
            yticklabels(labels)
        end

        if an == 1
            title(saved_data(an,j).cell_type)
        end

        % -----------------------------
        % reference lines (same logic as before)
        % -----------------------------
        plot([20 20]/4 - 5, [1 size(stacked_matrix,1)], 'k')
        plot([40 40]/4 - 5, [1 size(stacked_matrix,1)], 'k')

        plot([0 60]/4 - 5, [0 0] + numel(play_behaviors), ...
            'w', 'LineWidth', 2)

        plot([0 60]/4 - 5, [0 0] + numel(play_behaviors)/2, ...
            'w:', 'LineWidth', 2)

        plot([0 60]/4 - 5, [0 0] + numel(play_behaviors) + ...
            numel(non_play_behaviors)/2, 'w:', 'LineWidth', 2)

        clim([-4 4])
        axis tight
    end
end

%%% Delta-power PSTHs (play vs non-play)
% Time-warped LME beta for delta power. Heatmap: all behaviors (self/other);
% traces: mean ± CI (play red, non-play black). This is the LFP measure
% later correlated with the activation index.
load([saving_folder,'\results_full_lfp_per_behavior_time_warpped.mat'])
all_behaviors = {'Pounce','CC','Boxing','Evasion','Pin','Escape','CB','CD', ...
                 'Grooming','Rearing','Sniffing','Scratching','PounceI','Bite'};

alpha_tresh =Inf;

matrix_play = results_time_wrapped(1).beta_std;
pval_play   = results_time_wrapped(1).p_beh;
pval_play_ai =- log10(pval_play).*matrix_play;


matrix_play(~any(pval_play<alpha_tresh,2),:) = NaN;
play_indexes = [group1,group1+numel(all_behaviors)];


matrix_noplay = results_time_wrapped(2).beta_std;
pval_noplay  = results_time_wrapped(2).p_beh;
pval_noplay_ai = -log10(pval_noplay).*matrix_noplay;

matrix_noplay(~any(pval_noplay<alpha_tresh,2),:) = NaN;
noplay_indexes = [group2,group2+numel(all_behaviors)];

figure
colormap(jet(256))


ax2 = subplot(3,1,2:3);
matrix2plot = matrix_noplay;
mean2plot = mean(matrix2plot);
[~, ~, ci] = ttest(matrix2plot);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)],[ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis_time_wrapped, mean2plot, 'k')

matrix2plot = matrix_play;
mean2plot = mean(matrix2plot);
[~, ~, ci] = ttest(matrix2plot);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)],[ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis_time_wrapped, mean2plot, 'r')
plot(time_analysis_time_wrapped([1 end]),[0 0], 'k:')
plot([0 0],[-4 4], 'k')
plot([5 5],[-4 4], 'k')

axis tight

set(gca,'TickDir','out')

ax1 = subplot(3,1,1);




matrix2plot = [matrix_play;matrix_noplay  ];
imagesc(time_analysis_time_wrapped,1:size(matrix2plot,1),matrix2plot)
axis xy

clim([-4 4])
hold on
plot(time_analysis_time_wrapped([1 end]), [1 1]*size(matrix_play,1), 'w', 'LineWidth',2)
plot([0 0], [0 1]*size(matrix2plot,1), 'w', 'LineWidth',2)
plot([5 5], [0 1]*size(matrix2plot,1), 'w', 'LineWidth',2)

yticks(1:size(matrix2plot,1))
yticklabels(all_behaviors_2x([play_indexes,noplay_indexes]))
set(gca,'TickDir','out')
colorbar


pos1 = ax1.Position;
pos2 = ax2.Position;

ax2.Position(1) = pos1(1);   % align left edges
ax2.Position(3) = pos1(3);   % optional: match width

%% 3 Bout-averaged delta power
% Mean LME beta during the bout window [0 5] s. Filled bars: any time bin
% with p < 0.01. Swarm compares play vs non-play (rank-sum). Companion
% summary of the delta-power side of the correlation.

alpha_level = 0.01;

figure
subplot(1,2,1)
y_lim = [-5 5];
time2compare = [0 5];
play_bars = mean(matrix_play(:,time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2)),2);
non_play_bars = mean(matrix_noplay(:,time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2)),2);

non_play_p_val_window   = any(pval_noplay(:,time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2))<alpha_level,2);
[non_play_p_val_window, non_play_order] = sort(non_play_p_val_window);


play_p_val_window        = any(pval_play(:,time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2))<alpha_level,2);
[play_p_val_window, play_order] = sort(play_p_val_window);
bar_p_val = [non_play_p_val_window' play_p_val_window'];



bar_values = [non_play_bars(non_play_order)' play_bars(play_order)'];

           
% bar(bar_values)
hold on;

for i = 1:length(bar_values)

    if bar_p_val(i) == 1
        % significant → filled
        bar(i, bar_values(i), 'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'k');
    else
        % non-significant → empty
        bar(i, bar_values(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    end

end

xlim([0 length(bar_values)+1]);
xticks(1:numel([non_play_bars' play_bars']))
xticklabels(all_behaviors_2x([noplay_indexes(non_play_order),play_indexes(play_order)]))
ylim(y_lim)
plot([1 1]*numel(non_play_bars)+ .5,y_lim, 'k')



subplot(1,2,2)

rand_x = (rand(size(bar_values)) - .5)*.5 +[non_play_bars'*0 play_bars'*0+1];
plot(rand_x,bar_values, 'k.', 'MarkerSize',5)
hold on
plot(rand_x(bar_p_val),bar_values(bar_p_val), 'r.', 'MarkerSize',10)


% s = swarmchart([non_play_bars'*0 play_bars'*0+1]+rand_x,[non_play_bars' play_bars'], 'k.');




hold on
plot([-.25 .25], mean(non_play_bars)*[1 1], 'r', 'MarkerSize',8)
hold on
plot([-.25 .25]+1, mean(play_bars)*[1 1], 'r', 'MarkerSize',8)
[p, h, stats]=ranksum(non_play_bars,play_bars);
title([num2str(stats.zval), ' ', num2str(stats.ranksum),' ',num2str(p)])

ylim(y_lim)

%% 4 Match behaviors and correlate activation index with delta power
% Decode short PSTH codes (GR, PO, …) and LFP names onto the same
% behavior × self/other keys, then correlate the two measures without
% assuming identical list order. Also stores time-course r and a
% cross-behavior control (mismatched pairs).

% Behavioral ontology:
% GR = Grooming, RE = Rearing, SN = Sniffing, BI = Bite,
% CB = Approach, CD = Darting, CC/CH = Chasing, EV = Evasion,
% ES = Escape, PO/POA/POB = Pounce, PW* = PWI (pounce immobility).

%%% 4.1 Decode activation-index (PSTH) labels

n = numel(stacked_labels);
psth_decoded = cell(n,2); % {behavior, target}

for i = 1:n
    s = stacked_labels{i};

    % ---- target ----
    if contains(s,'Self')
        target = 'self';
    else
        target = 'other';
    end

    % ---- behavior ----
    if contains(s,'GR')
        beh = 'Grooming';

    elseif contains(s,'RE')
        beh = 'Rearing';

    elseif contains(s,'SN')
        beh = 'Sniffing';

    elseif contains(s,'BI')
        beh = 'Bite';

    elseif contains(s,'CB')
        beh = 'Approach';

    elseif contains(s,'CD')
        beh = 'Darting';

    elseif contains(s,'CC') || contains(s,'CH')
        beh = 'Chasing';

    elseif contains(s,'EV')
        beh = 'Evasion';

    elseif contains(s,'ES')
        beh = 'Escape';

    elseif contains(s,'PW')
        beh = 'PWI';

    elseif contains(s,'PO')
        beh = 'Pounce';

    else
        beh = 'Unknown';
    end

    psth_decoded{i,1} = beh;
    psth_decoded{i,2} = target;
end


%%% 4.2 Decode delta-power (LFP) labels

ref = all_behaviors_2x([play_indexes, noplay_indexes]);

ref_decoded = cell(numel(ref),2);

for i = 1:numel(ref)
    s = ref{i};

    % ---- target ----
    if contains(s,'self')
        target = 'self';
    else
        target = 'other';
    end

    % ---- behavior ----
    if contains(s,'Grooming')
        beh = 'Grooming';

    elseif contains(s,'Rearing')
        beh = 'Rearing';

    elseif contains(s,'Sniffing')
        beh = 'Sniffing';

    elseif contains(s,'Bite')
        beh = 'Bite';

    elseif contains(s,'Approach') || contains(s,'CB')
        beh = 'Approach';

    elseif contains(s,'Darting') || contains(s,'CD')
        beh = 'Darting';

    elseif contains(s,'Chasing') || contains(s,'CC') || contains(s,'CH')
        beh = 'Chasing';

    elseif contains(s,'Evasion')
        beh = 'Evasion';

    elseif contains(s,'Escape')
        beh = 'Escape';

    elseif contains(s,'PounceI')
        beh = 'PWI';

    elseif contains(s,'Pounce')
        beh = 'Pounce';

    else
        beh = 'Unknown';
    end

    ref_decoded{i,1} = beh;
    ref_decoded{i,2} = target;
end


%%% 4.3 Comparison keys (order-independent)

psth_key = strcat(psth_decoded(:,1), "_", psth_decoded(:,2));
ref_key  = strcat(ref_decoded(:,1), "_", ref_decoded(:,2));


%%% 4.4 Membership mapping (PSTH ↔ LFP)

% PSTH → where in reference?
[present_in_ref, ref_idx] = ismember(psth_key, ref_key);

% Reference → where in PSTH?
[present_in_psth, psth_idx] = ismember(ref_key, psth_key);


%%% 4.5 Match rates

fprintf('PSTH entries found in reference: %.2f%%\n', 100*mean(present_in_ref));
fprintf('Reference entries found in PSTH: %.2f%%\n', 100*mean(present_in_psth));


%%% 4.6 Mismatches (debug)

missing_in_ref = find(~present_in_ref);
missing_in_psth = find(~present_in_psth);

disp('--- PSTH missing in reference (first 10) ---');
for i = 1:min(10,numel(missing_in_ref))
    k = missing_in_ref(i);
    fprintf('%s\n', psth_key{k});
end

disp('--- Reference missing in PSTH (first 10) ---');
for i = 1:min(10,numel(missing_in_psth))
    k = missing_in_psth(i);
    fprintf('%s\n', ref_key{k});
end


%%% 4.7 Mapping table (for later example plots)

mapping = table;
mapping.psth_key = psth_key;
mapping.present_in_ref = present_in_ref;
mapping.ref_index = ref_idx;
mapping.psth_label =psth_decoded(:,1);
mapping.psth_so =psth_decoded(:,2);



%%% 4.8 Bout window [0 5] s
% Time axes of the two measures differ; restrict both to the bout.

time_psth = saved_data(1,1).time;
time_lfp  = time_analysis_time_wrapped;

idx_psth = time_psth >= 0 & time_psth <= 5;
idx_lfp  = time_lfp  >= 0 & time_lfp  <= 5;


%%% 4.9 Delta-power activity in that window

matrix_play   = results_time_wrapped(1).beta_std;
matrix_noplay = results_time_wrapped(2).beta_std;

lfp_all = [matrix_play; matrix_noplay];

lfp_activity = mean(lfp_all(:, idx_lfp), 2);


%%% 4.10 Preallocate correlation outputs

R = nan(5,3);
P = nan(5,3);

time_correalion_p =nan(5,3,26);
time_correalion_r =nan(5,3,26);
crossed_correlations = cell(5,3);

figure

%%% 4.11 Main loop: align, correlate, scatter
% For each area × cell type: bout-averaged AI vs delta power (scatter),
% time-course r per behavior, and crossed r (wrong behavior pairs) as control.

for an = 1:5

    for j = 1:3
        subplot(5,3,(an-1)*3 + j)

        %%% PSTH (activation-index) activity in the bout window

        stacked_matrix = saved_data(an,j).stacked_matrix;
       

        psth_activity = mean(stacked_matrix(:, idx_psth), 2);


        %%% Matched AI and delta-power values (aligned by key)

        psth_vals = [];
        re_arranged_psth = [];
        re_arranged_lfp = [];
        lfp_vals  = [];

        for i = 1:numel(psth_key)

            % skip behaviors not found in reference
            if ~present_in_ref(i)
                continue
            end

            % corresponding LFP row
            lfp_i = ref_idx(i);

            psth_vals(end+1,1) = psth_activity(i);
            re_arranged_psth(end+1,:) = stacked_matrix(i,:);

            lfp_vals(end+1,1)  = lfp_activity(lfp_i);
             aux_lfp= lfp_all(lfp_i,:);
            lfp_in_psth_time = interp1(time_lfp,aux_lfp,time_psth);
            re_arranged_lfp(end+1,:) = lfp_in_psth_time(:);

            no_nan = ~isnan(re_arranged_lfp(end,:)) &  ~isnan(re_arranged_psth(end,:));
            [r,p] =  corr(re_arranged_lfp(end,no_nan)',re_arranged_psth(end,no_nan)');
            time_correalion_p(an,j,i) =p;
            time_correalion_r(an,j,i)    =r;

            
        end

        for play_behavior1 = 1:7
            for play_behavior2 = play_behavior1+1:7
                            no_nan = re_arranged_psth(play_behavior1,:) &  ~isnan(re_arranged_lfp(play_behavior2,:));

            r = corr(re_arranged_psth(play_behavior1,no_nan )',re_arranged_lfp(play_behavior2,no_nan  )');
            crossed_correlations{an,j} = [crossed_correlations{an,j} ;r];
            end
        end


        for play_behavior1 = 8:14
            for play_behavior2 = play_behavior1+1:14
                            no_nan = re_arranged_psth(play_behavior1,:) &  ~isnan(re_arranged_lfp(play_behavior2,:));

            r = corr(re_arranged_psth(play_behavior1,no_nan )',re_arranged_lfp(play_behavior2,no_nan  )');
            crossed_correlations{an,j} = [crossed_correlations{an,j} ;r];
            end
        end

         for no_play_behavior1 = 15:20
            for no_play_behavior2 = no_play_behavior1+1:20
                no_nan = re_arranged_psth(no_play_behavior1,:) &  ~isnan(re_arranged_lfp(no_play_behavior2,:));

                r = corr(re_arranged_psth(no_play_behavior1,no_nan )',re_arranged_lfp(no_play_behavior2,no_nan )');
                crossed_correlations{an,j} = [crossed_correlations{an,j} ;r];
            end
         end


         for no_play_behavior1 = 21:26
            for no_play_behavior2 = no_play_behavior1+1:26
                no_nan = re_arranged_psth(no_play_behavior1,:) &  ~isnan(re_arranged_lfp(no_play_behavior2,:));

                r = corr(re_arranged_psth(no_play_behavior1,no_nan )',re_arranged_lfp(no_play_behavior2,no_nan )');
                crossed_correlations{an,j} = [crossed_correlations{an,j} ;r];
            end
        end




        %%% Drop unmatched / NaN pairs

        valid = ~isnan(psth_vals) & ~isnan(lfp_vals);

        psth_vals = psth_vals(valid);
        lfp_vals  = lfp_vals(valid);

        %%% Bout-averaged correlation (scatter of behaviors)

        if numel(psth_vals) >= 3
            [r,p] = corr(psth_vals, lfp_vals);

            plot(psth_vals,lfp_vals, 'k.')
            title(num2str([p r]))

            R(an,j) = r;
            P(an,j) = p;

        end

    end
end


%% 5 Area × cell-type correlation matrix
% Bout-averaged r from section 4 (activation index vs delta power, 0–5 s).
% Rows = structures, columns = peak / trough / unlocked.

figure
imagesc(R)

clim([-1 1])
colorbar

xlabel('Cell type')
ylabel('Area')

title('PSTH vs LFP correlation (0–5 s)')

set(gca,'XTick',1:3)
set(gca,'XTickLabel',{'Peak','Trough','Unlocked'})

set(gca,'YTick',1:5)

%% 6 Per-behavior time-course correlations
% r of the full warped traces (activation index vs interpolated delta power),
% one point per matched behavior. Red: p ≤ 0.01; black: n.s.

alpha_level = 0.01;

figure
for an=1:5
    for j=1:3
        subplot(5,3,(an-1)*3 + j)

        p_this_condition = squeeze(time_correalion_p(an,j,:));
        r_this_cond = squeeze(time_correalion_r(an,j,:));
        rand_x = (rand(size(r_this_cond))-.5)/2;
        plot(rand_x(p_this_condition>alpha_level),r_this_cond(p_this_condition>alpha_level), 'k.', 'MarkerSize',4)
        hold on
        plot(rand_x(p_this_condition<=alpha_level),r_this_cond(p_this_condition<=alpha_level), 'r.', 'MarkerSize',8)
        xlim([-1 1])
        ylim([-1 1])

    end
end

%% 7 Distribution of time-course r
% Histogram of the per-behavior r values from section 6, for each structure
% and cell type. Title = t-test p (mean r vs 0). Crossed (mismatch) r is
% computed in section 4 as a control (histogram currently commented).

alpha_level = 0.01;
mean_r_per_area = nan(5,3);
p_val_per_area = nan(5,3);
figure
for an=1:5
    for j=1:3
        subplot(5,3,(an-1)*3 + j)
        r_this_cond = squeeze(time_correalion_r(an,j,:));
        control_r = crossed_correlations{an,j};
               % histogram(control_r,-1:.1:1, 'FaceColor','k', 'FaceAlpha',.2, 'Normalization','percentage')
        hold on
       histogram(r_this_cond,-1:.1:1, 'FaceColor','r', 'FaceAlpha',.8, 'Normalization','percentage')
        mean_r_per_area(an,j)  = mean(r_this_cond, 'omitmissing');
      [h,p]= ttest(r_this_cond);
      p_val_per_area(an,j)=p;
      title(p)

    end
end

%% 8 Mean r per structure
% Bars of mean time-course r, grouped peak / trough / unlocked within each
% area. Filled = p ≤ 0.05 from section 7; hollow = n.s.


color_code = {'r','b','g','r','b','g','r','b','g','r','b','g','r','b','g'};
area_names = {};

for j=1:5

area_names{j} = saved_data(j,1).area;
end

mean_r_per_area_t = mean_r_per_area';
p_val_t = p_val_per_area';
mean_r_per_area_t = mean_r_per_area_t(:);
p_val_t = p_val_t(:);
figure
hold on
for j=1:numel(mean_r_per_area_t)

    if p_val_t(j)>0.05
        bar(j,mean_r_per_area_t(j), 'FaceColor','none', 'EdgeColor',color_code{j});

    else
         bar(j,mean_r_per_area_t(j), 'FaceColor',color_code{j}, 'EdgeColor',color_code{j});

    end
end
for j=1:4
plot([3 3]*j+.5, [-.25 .25], 'k:')
end

xticks(2:3:16)
xticklabels(area_names)



%% 9 Example traces: trough cells in LPAG (as on former supplementary figure)
% Overlay, for each matched behavior: trough-cell activation index in LPAG
% (blue), delta power (red), and trough-cell activation index in SC (black).
% Title = behavior key and r (AI vs delta).

an=3;
j=2;
an_comp = 1;
stacked_matrix = saved_data(an,j).stacked_matrix;
comparison_stacked = saved_data(an_comp,j).stacked_matrix;
figure
for i=1:26

lfp_i = ref_idx(i);

example_AI = stacked_matrix(i,:);

example_lfp  = lfp_all(lfp_i,:);
example_lfp = interp1(time_lfp,example_lfp,time_psth);

subplot(4,7,i)

plot(time_psth, movmean(example_AI,1), 'b')
hold on
plot(time_psth, movmean(example_lfp,1), 'r')
[r,p] = corr(example_AI(1:end-1)',example_lfp(1:end-1)');
plot(time_psth, movmean(comparison_stacked(i,:),1), 'k')

title([strrep(mapping.psth_key(i), '_',' '),  ' ', num2str(r)])
end
legend({'Activation Index Trough/LPAG','Delta power','ACtivation Index Trough/SC'})


%% 10 Example traces: trough cells in VLPAG (as on former supplementary figure)
% Same overlay as section 9, for trough cells in VLPAG vs SC.
figure


an=4;
j=2;
an_comp = 1;
stacked_matrix = saved_data(an,j).stacked_matrix;
comparison_stacked = saved_data(an_comp,j).stacked_matrix;
figure
for i=1:26

lfp_i = ref_idx(i);

example_AI = stacked_matrix(i,:);

example_lfp  = lfp_all(lfp_i,:);
example_lfp = interp1(time_lfp,example_lfp,time_psth);

subplot(4,7,i)

plot(time_psth, movmean(example_AI,1), 'b')
hold on
plot(time_psth, movmean(example_lfp,1), 'r')
[r,p] = corr(example_AI(1:end-1)',example_lfp(1:end-1)');
plot(time_psth, movmean(comparison_stacked(i,:),1), 'k')

title([mapping.psth_key(i),  ' ', num2str(r)])
end
legend({'Activation Index Trough/VLPAG','Delta power','ACtivation Index Trough/SC'})

