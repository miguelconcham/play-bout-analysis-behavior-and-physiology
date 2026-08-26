

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
figure_2_new_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\FIgure 2 delta modulation per behavior';


wind_length     = 1;
wind_overlap    = .990;
min_separation = .200;
f               = .1:.5:6;
freq_pow_range  = [1 5];


%% load if needed
disp('loading')
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';

load([saving_folder,'\psth_structure_delta_all_behaviors_self_other_partner.mat'],'psth_structure');
load([saving_folder,'\animal_names_delta_all_behaviors_self_other_partner.mat'],'animal_names');
disp('ready')
%% merging_psth
smooth_wind         = 20;
baseline_range      = [-2 0];
animal_label        = {'B1D1','B1S3','B2S2','B3D2', 'B4S2', 'B4D4'};
electorde_numner    = [1 2];
bin_size            = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
n_bins_timewrap     = psth_structure(1).n_bins_time_wrap;
psth_ranges         = psth_structure(1).hist_range;
wrap_range          = psth_structure(1).range_time_wrap;
time                = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index      = time<baseline_range(2) & time>baseline_range(1);
baseline_index_time_wrap = 1:round((abs(wrap_range(1))/bin_size));

all_psth_onset_self          = cell(size(psth_structure(1).all_behavior_psth_self));
all_session_index_self       = cell(1,size(psth_structure(1).all_behavior_psth_self,2));
all_electrode_index_self     = cell(1,size(psth_structure(1).all_behavior_psth_self,2));
all_psth_onset_other         = cell(size(psth_structure(1).all_behavior_psth_other));
all_session_index_other      = cell(1,size(psth_structure(1).all_behavior_psth_other,2));
all_electrode_index_other    = cell(1,size(psth_structure(1).all_behavior_psth_other,2));
all_partner_order            = cell(1,size(psth_structure(1).all_behavior_psth_other,2));

time_wrap_time          = [(baseline_index_time_wrap*bin_size) + wrap_range(1),linspace(bin_size,5,psth_structure(1).n_bins_time_wrap),5 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];


for j=1:numel(psth_structure)
    if contains(animal_names{j,1},animal_label)

        this_session_beahvior =  psth_structure(j).Behavior;
        animal_num           = find(cell2mat(cellfun(@(x) contains(animal_names{j,1},x), animal_label, 'UniformOutput',false)));
        electrode_num        =animal_names{j,2};
        animal_name_id = animal_names{j,1};
        animal_name_id = strsplit(animal_name_id, ' ');      


        for bn=1:size(all_psth_onset_self,2)
            n_obs = size(psth_structure(j).all_behavior_psth_self {1,bn},1);
            for tp=1:3
                if tp<3 && n_obs>0
                    all_psth_onset_self{tp,bn} = [all_psth_onset_self{tp,bn};psth_structure(j).all_behavior_psth_self{tp,bn}(:, 1:4000)];
                else
                    all_psth_onset_self{tp,bn} = [all_psth_onset_self{tp,bn};psth_structure(j).all_behavior_psth_self{tp,bn}(:, :)];
                end

            end
            all_session_index_self{bn}= [all_session_index_self{bn};animal_num*ones(n_obs,1)];
            all_electrode_index_self{bn}= [all_electrode_index_self{bn};electrode_num*ones(n_obs,1)];
        end


        for bn=1:size(all_psth_onset_other,2)
            n_obs = size(psth_structure(j).all_behavior_psth_other {1,bn},1);
            for tp=1:3
                if tp<3 && n_obs>0
                    all_psth_onset_other{tp,bn} = [all_psth_onset_other{tp,bn};psth_structure(j).all_behavior_psth_other{tp,bn}(:, 1:4000)];
                else
                    all_psth_onset_other{tp,bn} = [all_psth_onset_other{tp,bn};psth_structure(j).all_behavior_psth_other{tp,bn}(:, :)];
                end

            end
            all_session_index_other{bn}= [all_session_index_other{bn};animal_num*ones(n_obs,1)];
            all_electrode_index_other{bn}= [all_electrode_index_other{bn};electrode_num*ones(n_obs,1)];
            all_partner_order{bn}= [all_partner_order{bn};psth_structure(j).partner_order_pb{bn}];
        end
    end

end





%% time_wrap_power

all_behaviors           =  {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD','Grooming','Rearing', 'Sniffing','Scratching','PounceI', 'Bite'};

figure
for b=1:14
    subplot(4,4,b)


    plot(time_wrap_time, mean(all_psth_onset_self{3,b}, 'omitmissing'))
    hold on
    yyaxis right
    plot(time_wrap_time,mean(all_psth_onset_other{3,b}, 'omitmissing'))
    xlim([-5 5])
    title(all_behaviors{b})
end

%% osnet pwoer
all_behaviors           =  {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD','Grooming','Rearing', 'Sniffing','Scratching','PounceI', 'Bite'};

figure
for b=1:14
    subplot(4,4,b)


    plot(time(2:end), mean(all_psth_onset_self{1,b}, 'omitmissing'))
    hold on
    yyaxis right
    plot(time(2:end),mean(all_psth_onset_other{1,b}, 'omitmissing'))
    xlim([-5 5])
    title(all_behaviors{b})
end


%%



n = numel(all_behaviors);

all_behaviors_2x = cell(1, 2*n);

% self half
for i = 1:n
    all_behaviors_2x{i} = [all_behaviors{i} ' self'];
end

% other half
for i = 1:n
    all_behaviors_2x{n+i} = [all_behaviors{i} ' other'];
end
%% creat lfp mixed model

% -------------------------
% TIME WINDOWS
% -------------------------
baseline_idx = time >= -2 & time <= 0;
analysis_idx = time >= -1 & time <= 2;

time_analysis = time(analysis_idx);

% -------------------------
% BEHAVIOR GROUPS
% -------------------------
group1 = [1:8];
group2 = [ 9 10 11 13 14];
selected_behaviors_all = [group1,group2];

beh_groups = {group1, group2};

% -------------------------
% OUTPUT STRUCTURE
% -------------------------
results = struct;


%%
% 
% for g = 1:2
% 
%     beh_list = beh_groups{g};
% 
%     all_y = [];
%     all_session = [];
%     all_behavior = [];
% 
%     % =========================
%     % STACK DATA
%     % =========================
%     for b = beh_list
% 
%         psth = all_psth_onset_self{1,b};
%         sess = all_session_index_self{b};
%         elec = all_electrode_index_self{b};
% 
%         % electrode filter
%         idx = elec == 1;
%         psth = psth(idx,:);
%         sess = sess(idx);
% 
%         % baseline z-score
%         mu = mean(psth(:,baseline_idx),2);
%         sd = std(psth(:,baseline_idx),0,2);
%         sd(sd==0) = NaN;
% 
%         psth = (psth - mu) ./ sd;
% 
%         valid = ~isnan(sd);
%         psth = psth(valid,:);
%         sess = sess(valid);
% 
%         % crop time
%         psth = psth(:,analysis_idx);
% 
%         nObs = size(psth,1);
% 
%         all_y = [all_y; psth];
%         all_session = [all_session; sess];
%         all_behavior = [all_behavior; repmat(b, nObs,1)];
% 
%     end
% 
% 
%     for b = beh_list
% 
%         psth = all_psth_onset_other{1,b};
%         sess = all_session_index_other{b};
%         elec = all_electrode_index_other{b};
% 
%         % electrode filter
%         idx = elec == 1;
%         psth = psth(idx,:);
%         sess = sess(idx);
% 
%         % baseline z-score
%         mu = mean(psth(:,baseline_idx),2);
%         sd = std(psth(:,baseline_idx),0,2);
%         sd(sd==0) = NaN;
% 
%         psth = (psth - mu) ./ sd;
% 
%         valid = ~isnan(sd);
%         psth = psth(valid,:);
%         sess = sess(valid);
% 
%         % crop time
%         psth = psth(:,analysis_idx);
% 
%         nObs = size(psth,1);
% 
%         all_y = [all_y; psth];
%         all_session = [all_session; sess];
%         all_behavior = [all_behavior; repmat(b+numel(all_behaviors), nObs,1)];
% 
%     end
% 
% 
%     self_other_beh_list = unique(all_behavior);
%     nBeh = length(beh_list);
%     nTime = length(time_analysis);
% 
%     beta = nan(2*nBeh, nTime);
%     p_beh = nan(2*nBeh, nTime);
%     p_global = nan(1, nTime);
% 
%     % =========================
%     % LME PER TIMEPOINT
%     % =========================
%     for t = 1:nTime
%         disp(num2str(time_analysis(t)))
% 
%         y = all_y(:,t);
%         valid = ~isnan(y);
% 
%         tbl = table( ...
%             y(valid), ...
%             categorical(all_session(valid)), ...
%             categorical(all_behavior(valid)), ...
%             'VariableNames', {'y','session','behavior'});
% 
%         y = y(valid);
%         sess = all_session(valid);
% 
%         v = var(y);
%         if numel(unique(tbl.session)) < 2
%             continue;
%         end
% 
%         % -------------------------
%         % GLOBAL EFFECT
%         % -------------------------
%         lme_g = fitlme(tbl, 'y ~ 1 + (1|session)');
%         p_global(t) = lme_g.Coefficients.pValue(1);
% 
%         % -------------------------
%         % PER BEHAVIOR EFFECT
%         % -------------------------
%         for i = 1:2*nBeh
%             lme_b = fitlme(tbl(tbl.behavior == categorical(self_other_beh_list(i)),:) , 'y ~ 1 + (1|session)');
% 
%             % LOW VARIANCE REGIME
%             % fallback: session-mean t-test
%             p_beh(i,t) = lme_b.Coefficients.pValue(1);
%             beta(i,t)  = lme_b.Coefficients.Estimate(1);
% 
%         end
%     end
% 
%     % -------------------------
%     % STORE RESULTS CLEANLY
%     % -------------------------
%     results(g).beta = beta;
%     results(g).p_beh = p_beh;
%     results(g).p_global = p_global;
%     results(g).behaviors = self_other_beh_list;
% 
% end
%%
all_behaviors = {'Pounce','CC','Boxing','Evasion','Pin','Escape','CB','CD', ...
                 'Grooming','Rearing','Sniffing','Scratching','PounceI','Bite'};

%% -------------------------
% TIME WINDOWS
% -------------------------
baseline_idx = time >= -2 & time <= 0;
analysis_idx = time >= -1 & time <= 2;

time_analysis = time(analysis_idx);

% -------------------------
% BEHAVIOR GROUPS
% -------------------------
group1 = 1:8;
group2 = [9 10 11 13 14];

beh_groups = {group1, group2};

%% -------------------------
% OUTPUT STRUCTURE
% -------------------------
results = struct;

%========================
% MAIN LOOP
%=========================
for g = 1:2

    beh_list = beh_groups{g};
    nBeh = length(beh_list);

    all_y = [];
    all_session = [];
    all_behavior = [];

    % =========================
    % STACK SELF DATA
    % =========================
    for b = beh_list

        psth = all_psth_onset_self{1,b};
        sess = all_session_index_self{b};
        elec = all_electrode_index_self{b};

        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);

        % baseline z-score
        mu = mean(psth(:,baseline_idx),2);
        sd = std(psth(:,baseline_idx),0,2);
        sd(sd==0) = NaN;

        psth = (psth - mu) ./ sd;

        valid = ~isnan(sd);
        psth = psth(valid,:);
        sess = sess(valid);

        psth = psth(:,analysis_idx);
        nObs = size(psth,1);

        all_y = [all_y; psth];
        all_session = [all_session; sess];
        all_behavior = [all_behavior; repmat(b, nObs,1)];

    end

    % =========================
    % STACK OTHER DATA
    % =========================
    for b = beh_list

        psth = all_psth_onset_other{1,b};
        sess = all_session_index_other{b};
        elec = all_electrode_index_other{b};

        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);

        % baseline z-score
        mu = mean(psth(:,baseline_idx),2);
        sd = std(psth(:,baseline_idx),0,2);
        sd(sd==0) = NaN;

        psth = (psth - mu) ./ sd;

        valid = ~isnan(sd);
        psth = psth(valid,:);
        sess = sess(valid);

        psth = psth(:,analysis_idx);
        nObs = size(psth,1);

        all_y = [all_y; psth];
        all_session = [all_session; sess];
        all_behavior = [all_behavior; repmat(b + numel(all_behaviors), nObs,1)];

    end

    % -------------------------
    % DEFINE BEHAVIOR LIST (ORDERED!)
    % -------------------------
    self_other_beh_list = [beh_list, beh_list + numel(all_behaviors)];

    nTime = length(time_analysis);

    % outputs
    beta = nan(2*nBeh, nTime);
    beta_std = nan(2*nBeh, nTime);   % <-- NEW (standardized beta)
    p_beh = nan(2*nBeh, nTime);
    p_global = nan(1, nTime);

    % =========================
    % LME PER TIMEPOINT
    % =========================
    for t = 1:nTime
        disp(num2str(time_analysis(t)))

        y = all_y(:,t);
        valid = ~isnan(y);

        tbl = table( ...
            y(valid), ...
            categorical(all_session(valid)), ...
            categorical(all_behavior(valid)), ...
            'VariableNames', {'y','session','behavior'});

        if numel(unique(tbl.session)) < 2
            continue;
        end

        % -------------------------
        % GLOBAL EFFECT
        % -------------------------
        lme_g = fitlme(tbl, 'y ~ 1 + (1|session)');
        p_global(t) = lme_g.Coefficients.pValue;

        % -------------------------
        % PER BEHAVIOR EFFECT
        % -------------------------
        for i = 1:2*nBeh

            beh_id = self_other_beh_list(i);

            tbl_i = tbl(tbl.behavior == categorical(beh_id), :);
            if isempty(tbl_i)
                continue;
            end

            lme_b = fitlme(tbl_i, 'y ~ 1 + (1|session)');

            % raw beta (mean)
            beta(i,t) = lme_b.Coefficients.Estimate;

            % p-value
            p_beh(i,t) = lme_b.Coefficients.pValue;

            % -------------------------
            % STANDARDIZED BETA (KEY ADDITION)
            % -------------------------
            se = lme_b.Coefficients.SE;

            if se > 0
                beta_std(i,t) = beta(i,t) / se;
            end

        end
    end

    % -------------------------
    % STORE RESULTS
    % -------------------------
    results(g).beta         = beta;
    results(g).beta_std     = beta_std;   % <-- NEW
    results(g).p_beh        = p_beh;
    results(g).p_global     = p_global;
    results(g).behaviors    = self_other_beh_list;

end
%% save result to avoid estimating model again
save([saving_folder,'\results_full_lfp_per_behavior.mat'], ...
     'results', ...
     'all_behaviors_2x', ...
     'group1', ...
     'group2', ...
     'time_analysis','all_behaviors', ...
     '-v7.3')
%% Load results if needed

load([saving_folder,'\results_full_lfp_per_behavior.mat'], ...
     'results', ...
     'all_behaviors_2x', ...
     'group1', ...
     'group2', ...
     'time_analysis','all_behaviors')
%% plot images and mean respones per behavior type
all_behaviors = {'Pounce','CC','Boxing','Evasion','Pin','Escape','CB','CD', ...
                 'Grooming','Rearing','Sniffing','Scratching','PounceI','Bite'};

alpha_tresh =Inf;

matrix_play = results(1).beta_std;
pval_play   = results(1).p_beh;
pval_play_ai =- log10(pval_play).*matrix_play;


matrix_play(~any(pval_play<alpha_tresh,2),:) = NaN;
play_indexes = [group1,group1+numel(all_behaviors)];


matrix_noplay = results(2).beta_std;
pval_noplay  = results(2).p_beh;
pval_noplay_ai = -log10(pval_noplay).*matrix_noplay;

matrix_noplay(~any(pval_noplay<alpha_tresh,2),:) = NaN;
noplay_indexes = [group2,group2+numel(all_behaviors)];

figure
colormap(jet(256))


ax2 = subplot(3,1,2:3);
matrix2plot = matrix_noplay;
mean2plot = mean(matrix2plot);
[~, ~, ci] = ttest(matrix2plot);
fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean2plot, 'k')

matrix2plot = matrix_play;
mean2plot = mean(matrix2plot);
[~, ~, ci] = ttest(matrix2plot);
fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean2plot, 'r')
plot(time_analysis([1 end]),[0 0], 'k:')
plot([0 0],[-2.5 2.5], 'k')
axis tight

set(gca,'TickDir','out')

ax1 = subplot(3,1,1);




matrix2plot = [matrix_play;matrix_noplay  ];
imagesc(time_analysis,1:size(matrix2plot,1),matrix2plot)
axis xy

clim([-4 4])
hold on
plot(time_analysis([1 end]), [1 1]*size(matrix_play,1), 'w', 'LineWidth',2)
plot([0 0], [0 1]*size(matrix2plot,1), 'w', 'LineWidth',2)
yticks(1:size(matrix2plot,1))
yticklabels(all_behaviors_2x([play_indexes,noplay_indexes]))
set(gca,'TickDir','out')
colorbar


pos1 = ax1.Position;
pos2 = ax2.Position;

ax2.Position(1) = pos1(1);   % align left edges
ax2.Position(3) = pos1(3);   % optional: match width


%% for saving figures used this code line
% print(gcf,'-vector','-dsvg',[figure_2_new_folder,'\images dn mean response colormap jet.svg'])
%%

figure
subplot(1,2,1)
y_lim = [-5 5];
time2compare = [-.2 .5];
play_bars = mean(matrix_play(:,time_analysis>time2compare(1) & time_analysis<time2compare(2)),2);
non_play_bars = mean(matrix_noplay(:,time_analysis>time2compare(1) & time_analysis<time2compare(2)),2);

non_play_p_val_window   = any(pval_noplay(:,time_analysis>time2compare(1) & time_analysis<time2compare(2))<0.05,2);
[non_play_p_val_window, non_play_order] = sort(non_play_p_val_window);


play_p_val_window        = any(pval_play(:,time_analysis>time2compare(1) & time_analysis<time2compare(2))<0.05,2);
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

ylim([-2.5 3.5])
%% for saving figures used this code line
% print(gcf,'-vector','-dsvg',[figure_2_new_folder,'\bar plots and scatter plot mean response.svg'])


%% self other play




figure
matrix_play = results(1).beta_std;
pval_play   = results(1).p_beh;


matrix_play(~any(pval_play<alpha_tresh,2),:) = NaN;
play_indexes = [group1,group1+numel(all_behaviors)];

index_self              = group1;
index_other             = 1:numel(play_indexes);
index_other(index_self) = [];


matrix2plot         = matrix_play(index_self,:);
mean2plot           = mean(matrix2plot);
[~, ~, ci]          = ttest(matrix2plot);

fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean2plot, 'b')

matrix2plot         = matrix_play(index_other,:);
mean2plot           = mean(matrix2plot);
[~, ~, ci]          = ttest(matrix2plot);
fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean2plot, 'r')
plot(time_analysis([1 end]),[0 0], 'k:')
plot([0 0],[-2 4], 'k')
axis tight
set(gca,'TickDir','out')
%%

print(gcf,'-vector','-dsvg',[figure_2_new_folder,'\mean response play self and other.svg'])

%% self other no play

figure

alpha_tresh = Inf;

matrix_noplay = results(2).beta_std;
pval_noplay  = results(2).p_beh;


matrix_noplay(~any(pval_noplay<alpha_tresh,2),:) = NaN;
no_play_indexes = [group2,group2+numel(all_behaviors)];

index_self              = 1:numel(group2);
index_other             = 1:numel(no_play_indexes);
index_other(index_self) = [];


matrix2plot         = matrix_noplay(index_self,:);
mean2plot           = mean(matrix2plot);
[~, ~, ci]          = ttest(matrix2plot);

fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean2plot, 'b')

matrix2plot         = matrix_noplay(index_other,:);
mean2plot           = mean(matrix2plot);
[~, ~, ci]          = ttest(matrix2plot);
fill([time_analysis fliplr(time_analysis)],[ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean2plot, 'r')
plot(time_analysis([1 end]),[0 0], 'k:')
plot([0 0],[-4 2], 'k')
axis tight
set(gca,'TickDir','out')

%% creating lfp respones timewrapped

%% -------------------------
% TIME WINDOWS
% -------------------------
baseline_idx = time_wrap_time >= -5 & time_wrap_time <= 0;
analysis_idx = time_wrap_time >= -5 & time_wrap_time <= 10;

time_analysis_time_wrapped = time_wrap_time(analysis_idx);

% -------------------------
% BEHAVIOR GROUPS
% -------------------------
group1 = 1:8;
group2 = [9 10 11 13 14];

beh_groups = {group1, group2};

%% -------------------------
% OUTPUT STRUCTURE
% -------------------------
results_time_wrapped = struct;

%========================
% MAIN LOOP
%=========================
for g = 1:2

    beh_list        = beh_groups{g};
    nBeh            = length(beh_list);

    all_y           = [];
    all_session     = [];
    all_behavior    = [];

    % =========================
    % STACK SELF DATA
    % =========================
    for b = beh_list

        psth = all_psth_onset_self{3,b};
        sess = all_session_index_self{b};
        elec = all_electrode_index_self{b};

        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);

        % baseline z-score
        mu = mean(psth(:,baseline_idx),2);
        sd = std(psth(:,baseline_idx),0,2);
        sd(sd==0) = NaN;

        psth = (psth - mu) ./ sd;

        valid = ~isnan(sd);
        psth = psth(valid,:);
        sess = sess(valid);

        psth = psth(:,analysis_idx);
        nObs = size(psth,1);

        all_y = [all_y; psth];
        all_session = [all_session; sess];
        all_behavior = [all_behavior; repmat(b, nObs,1)];

    end

    % =========================
    % STACK OTHER DATA
    % =========================
    for b = beh_list

        psth = all_psth_onset_other{3,b};
        sess = all_session_index_other{b};
        elec = all_electrode_index_other{b};

        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);

        % baseline z-score
        mu = mean(psth(:,baseline_idx),2);
        sd = std(psth(:,baseline_idx),0,2);
        sd(sd==0) = NaN;

        psth = (psth - mu) ./ sd;

        valid = ~isnan(sd);
        psth = psth(valid,:);
        sess = sess(valid);

        psth = psth(:,analysis_idx);
        nObs = size(psth,1);

        all_y = [all_y; psth];
        all_session = [all_session; sess];
        all_behavior = [all_behavior; repmat(b + numel(all_behaviors), nObs,1)];

    end

    % -------------------------
    % DEFINE BEHAVIOR LIST (ORDERED!)
    % -------------------------
    self_other_beh_list = [beh_list, beh_list + numel(all_behaviors)];

    nTime = length(time_analysis_time_wrapped);

    % outputs
    beta = nan(2*nBeh, nTime);
    beta_std = nan(2*nBeh, nTime);   % <-- NEW (standardized beta)
    p_beh = nan(2*nBeh, nTime);
    p_global = nan(1, nTime);

    % =========================
    % LME PER TIMEPOINT
    % =========================
    for t = 1:nTime
        disp(num2str(time_analysis_time_wrapped(t)))

        y = all_y(:,t);
        valid = ~isnan(y);

        tbl = table( ...
            y(valid), ...
            categorical(all_session(valid)), ...
            categorical(all_behavior(valid)), ...
            'VariableNames', {'y','session','behavior'});

        if numel(unique(tbl.session)) < 2
            continue;
        end

        % -------------------------
        % GLOBAL EFFECT
        % -------------------------
        lme_g = fitlme(tbl, 'y ~ 1 + (1|session)');
        p_global(t) = lme_g.Coefficients.pValue;

        % -------------------------
        % PER BEHAVIOR EFFECT
        % -------------------------
        for i = 1:2*nBeh

            beh_id = self_other_beh_list(i);

            tbl_i = tbl(tbl.behavior == categorical(beh_id), :);
            if isempty(tbl_i)
                continue;
            end

            lme_b = fitlme(tbl_i, 'y ~ 1 + (1|session)');

            % raw beta (mean)
            beta(i,t) = lme_b.Coefficients.Estimate;

            % p-value
            p_beh(i,t) = lme_b.Coefficients.pValue;

            % -------------------------
            % STANDARDIZED BETA (KEY ADDITION)
            % -------------------------
            se = lme_b.Coefficients.SE;

            if se > 0
                beta_std(i,t) = beta(i,t) / se;
            end

        end
    end

    % -------------------------
    % STORE RESULTS
    % -------------------------
    results_time_wrapped(g).beta         = beta;
    results_time_wrapped(g).beta_std     = beta_std;   % <-- NEW
    results_time_wrapped(g).p_beh        = p_beh;
    results_time_wrapped(g).p_global     = p_global;
    results_time_wrapped(g).behaviors    = self_other_beh_list;

end
%% save if needed

save([saving_folder,'\results_full_lfp_per_behavior_time_warpped.mat'], ...
     'results_time_wrapped', ...
     'all_behaviors_2x', ...
     'group1', ...
     'group2', ...
     'time_analysis_time_wrapped','all_behaviors', ...
     '-v7.3')
%% load if needed

load([saving_folder,'\results_full_lfp_per_behavior_time_warpped.mat'], ...
     'results_time_wrapped', ...
     'all_behaviors_2x', ...
     'group1', ...
     'group2', ...
     'time_analysis_time_wrapped','all_behaviors')
%%

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
plot([0 0],[-5 5], 'k')
plot([5 5],[-5 5], 'k')

axis tight

set(gca,'TickDir','out')

ax1 = subplot(3,1,1);




matrix2plot = [matrix_play;matrix_noplay  ];
imagesc(time_analysis_time_wrapped,1:size(matrix2plot,1),matrix2plot)
axis xy

clim([-10 10])
hold on
plot(time_analysis_time_wrapped([1 end]), [1 1]*size(matrix_play,1), 'w', 'LineWidth',2)
plot([0 0], [0 1]*size(matrix2plot,1), 'w', 'LineWidth',2)
yticks(1:size(matrix2plot,1))
yticklabels(all_behaviors_2x([play_indexes,noplay_indexes]))
set(gca,'TickDir','out')
colorbar


pos1 = ax1.Position;
pos2 = ax2.Position;

ax2.Position(1) = pos1(1);   % align left edges
ax2.Position(3) = pos1(3);   % optional: match width


%%
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

