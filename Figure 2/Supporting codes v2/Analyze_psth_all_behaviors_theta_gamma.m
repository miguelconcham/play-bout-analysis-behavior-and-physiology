%% Analyze_psth_all_behaviors_theta_gamma
% Per-behavior PSTH analysis (self vs partner) with LME and plotting.
% Onset-aligned and time-wrapped estimates, play and non-play.
% Separate from Analyze_psth_all_behaviors.m, which stays the delta script.
%
% band_list selects which saved Estimate_psth_all_behaviors_theta_gamma.m
% files to fit. 'delta' keeps the original result filenames. 'theta' and
% 'gamma' write results_full_lfp_per_behavior_<band>.mat and
% results_full_lfp_per_behavior_time_warpped_<band>.mat.
% Each LME is skipped when that file is already in saving_folder.
% The %% Save sections only write the workspace, so a finished fit can be
% saved without fitting again.
% Figure window titles are psth_<band>_<plot>, the same plot name for
% delta, theta, and gamma, and safe to use as a Windows file name.
% Self vs other plots mark times when that mean is significantly above 0
% (darker bar) or below 0 (lighter bar). Blue is self, red is other.
% Run Estimate_psth_all_behaviors_theta_gamma.m first so those PSTH files exist.

this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_psth_all_behaviors_theta_gamma.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');
figure_2_new_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\FIgure 2 delta modulation per behavior';

all_behaviors = {'Pounce', 'CC','Boxing', 'Evasion','Pin', 'Escape', 'CB', 'CD','Grooming','Rearing', 'Sniffing','Scratching','PounceI', 'Bite'};
band_list = {'delta', 'theta', 'gamma'};

for band_n = 1:numel(band_list)
band_name = band_list{band_n};
structure_file = ['psth_structure_' band_name '_all_behaviors_self_other_partner.mat'];
names_file = ['animal_names_' band_name '_all_behaviors_self_other_partner.mat'];
if strcmp(band_name, 'delta')
    results_file = 'results_full_lfp_per_behavior.mat';
    results_tw_file = 'results_full_lfp_per_behavior_time_warpped.mat';
    svg_name = 'mean response play self and other.svg';
else
    results_file = ['results_full_lfp_per_behavior_' band_name '.mat'];
    results_tw_file = ['results_full_lfp_per_behavior_time_warpped_' band_name '.mat'];
    svg_name = ['mean response play self and other ' band_name '.svg'];
end

%% Load precomputed behavior PSTH structures
disp(['loading ' band_name])
load(fullfile(saving_folder, structure_file), 'psth_structure');
load(fullfile(saving_folder, names_file), 'animal_names');
disp('ready')

%% Merge PSTH across sessions (self and partner)
animal_label = {'B1D1','B1S3','B2S2','B3D2', 'B4S2', 'B4D4'};
bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges = psth_structure(1).hist_range;
wrap_range = psth_structure(1).range_time_wrap;
time = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index_time_wrap = 1:round((abs(wrap_range(1))/bin_size));

all_psth_onset_self = cell(size(psth_structure(1).all_behavior_psth_self));
all_session_index_self = cell(1,size(psth_structure(1).all_behavior_psth_self,2));
all_electrode_index_self = cell(1,size(psth_structure(1).all_behavior_psth_self,2));
all_psth_onset_other = cell(size(psth_structure(1).all_behavior_psth_other));
all_session_index_other = cell(1,size(psth_structure(1).all_behavior_psth_other,2));
all_electrode_index_other = cell(1,size(psth_structure(1).all_behavior_psth_other,2));
all_partner_order = cell(1,size(psth_structure(1).all_behavior_psth_other,2));

time_wrap_time = [(baseline_index_time_wrap*bin_size) + wrap_range(1), linspace(bin_size,5,psth_structure(1).n_bins_time_wrap), 5 + (1:round((abs(wrap_range(2))/bin_size)))*bin_size];

for j = 1:numel(psth_structure)
    if contains(animal_names{j,1}, animal_label)
        animal_num = find(cell2mat(cellfun(@(x) contains(animal_names{j,1}, x), animal_label, 'UniformOutput', false)));
        electrode_num = animal_names{j,2};

        for bn = 1:size(all_psth_onset_self, 2)
            n_obs = size(psth_structure(j).all_behavior_psth_self{1,bn}, 1);
            for tp = 1:3
                if tp < 3 && n_obs > 0
                    all_psth_onset_self{tp,bn} = [all_psth_onset_self{tp,bn}; psth_structure(j).all_behavior_psth_self{tp,bn}(:, 1:4000)];
                else
                    all_psth_onset_self{tp,bn} = [all_psth_onset_self{tp,bn}; psth_structure(j).all_behavior_psth_self{tp,bn}(:, :)];
                end
            end
            all_session_index_self{bn} = [all_session_index_self{bn}; animal_num*ones(n_obs,1)];
            all_electrode_index_self{bn} = [all_electrode_index_self{bn}; electrode_num*ones(n_obs,1)];
        end

        for bn = 1:size(all_psth_onset_other, 2)
            n_obs = size(psth_structure(j).all_behavior_psth_other{1,bn}, 1);
            for tp = 1:3
                if tp < 3 && n_obs > 0
                    all_psth_onset_other{tp,bn} = [all_psth_onset_other{tp,bn}; psth_structure(j).all_behavior_psth_other{tp,bn}(:, 1:4000)];
                else
                    all_psth_onset_other{tp,bn} = [all_psth_onset_other{tp,bn}; psth_structure(j).all_behavior_psth_other{tp,bn}(:, :)];
                end
            end
            all_session_index_other{bn} = [all_session_index_other{bn}; animal_num*ones(n_obs,1)];
            all_electrode_index_other{bn} = [all_electrode_index_other{bn}; electrode_num*ones(n_obs,1)];
            all_partner_order{bn} = [all_partner_order{bn}; psth_structure(j).partner_order_pb{bn}];
        end
    end
end

%% Preview — time-wrapped mean PSTH per behavior
figure
name_psth_figure(band_name, 'timewrap_preview')
for b = 1:14
    subplot(4,4,b)
    plot(time_wrap_time, mean(all_psth_onset_self{3,b}, 'omitmissing'))
    hold on
    yyaxis right
    plot(time_wrap_time, mean(all_psth_onset_other{3,b}, 'omitmissing'))
    xlim([-5 5])
    title(all_behaviors{b})
end

%% Preview — onset-aligned mean PSTH per behavior
figure
name_psth_figure(band_name, 'onset_preview')
for b = 1:14
    subplot(4,4,b)
    plot(time(2:end), mean(all_psth_onset_self{1,b}, 'omitmissing'))
    hold on
    yyaxis right
    plot(time(2:end), mean(all_psth_onset_other{1,b}, 'omitmissing'))
    xlim([-5 5])
    title(all_behaviors{b})
end

%% Build self/other behavior labels
n = numel(all_behaviors);
all_behaviors_2x = cell(1, 2*n);
for i = 1:n
    all_behaviors_2x{i} = [all_behaviors{i} ' self'];
    all_behaviors_2x{n+i} = [all_behaviors{i} ' other'];
end

%% LME analysis — onset-aligned PSTH
% Skipped when results_full_lfp_per_behavior_<band>.mat is already saved.
group1 = 1:8;
group2 = [9 10 11 13 14];
beh_groups = {group1, group2};
if exist(fullfile(saving_folder, results_file), 'file') == 2
    disp(['onset LME already saved, loading ' results_file])
    load(fullfile(saving_folder, results_file), ...
        'results', 'all_behaviors_2x', 'group1', 'group2', 'time_analysis', 'all_behaviors')
else
baseline_idx = time >= -2 & time <= 0;
analysis_idx = time >= -1 & time <= 2;
time_analysis = time(analysis_idx);
results = struct;

for g = 1:2
    beh_list = beh_groups{g};
    nBeh = length(beh_list);
    all_y = [];
    all_session = [];
    all_behavior = [];

    for b = beh_list
        psth = all_psth_onset_self{1,b};
        sess = all_session_index_self{b};
        elec = all_electrode_index_self{b};
        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);
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

    for b = beh_list
        psth = all_psth_onset_other{1,b};
        sess = all_session_index_other{b};
        elec = all_electrode_index_other{b};
        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);
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

    self_other_beh_list = [beh_list, beh_list + numel(all_behaviors)];
    nTime = length(time_analysis);
    beta = nan(2*nBeh, nTime);
    beta_std = nan(2*nBeh, nTime);
    p_beh = nan(2*nBeh, nTime);
    p_global = nan(1, nTime);

    for t = 1:nTime
        disp(num2str(time_analysis(t)))
        y = all_y(:,t);
        valid = ~isnan(y);
        tbl = table(y(valid), categorical(all_session(valid)), categorical(all_behavior(valid)), 'VariableNames', {'y','session','behavior'});
        if numel(unique(tbl.session)) < 2
            continue;
        end
        lme_g = fitlme(tbl, 'y ~ 1 + (1|session)');
        p_global(t) = lme_g.Coefficients.pValue;
        for i = 1:2*nBeh
            beh_id = self_other_beh_list(i);
            tbl_i = tbl(tbl.behavior == categorical(beh_id), :);
            if isempty(tbl_i)
                continue;
            end
            lme_b = fitlme(tbl_i, 'y ~ 1 + (1|session)');
            beta(i,t) = lme_b.Coefficients.Estimate;
            p_beh(i,t) = lme_b.Coefficients.pValue;
            se = lme_b.Coefficients.SE;
            if se > 0
                beta_std(i,t) = beta(i,t) / se;
            end
        end
    end

    results(g).beta = beta;
    results(g).beta_std = beta_std;
    results(g).p_beh = p_beh;
    results(g).p_global = p_global;
    results(g).behaviors = self_other_beh_list;
end
end

%% Save onset-aligned LME
% This section only writes the workspace. Run it alone after a fit.
disp(['saving ' results_file])
save(fullfile(saving_folder, results_file), ...
    'results', 'all_behaviors_2x', 'group1', 'group2', 'time_analysis', 'all_behaviors', 'band_name', '-v7.3')

%% Plot colormap and mean responses (play vs non-play)
alpha_tresh = Inf;
matrix_play = results(1).beta_std;
pval_play = results(1).p_beh;
matrix_play(~any(pval_play<alpha_tresh,2),:) = NaN;
play_indexes = [group1, group1+numel(all_behaviors)];
matrix_noplay = results(2).beta_std;
pval_noplay = results(2).p_beh;
matrix_noplay(~any(pval_noplay<alpha_tresh,2),:) = NaN;
noplay_indexes = [group2, group2+numel(all_behaviors)];

figure
name_psth_figure(band_name, 'onset_play_nonplay')
colormap(jet(256))
ax2 = subplot(3,1,2:3);
matrix2plot = matrix_noplay;
[~, ~, ci] = ttest(matrix2plot);
fill([time_analysis fliplr(time_analysis)], [ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean(matrix2plot), 'k')
matrix2plot = matrix_play;
[~, ~, ci] = ttest(matrix2plot);
fill([time_analysis fliplr(time_analysis)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(time_analysis, mean(matrix2plot), 'r')
plot(time_analysis([1 end]), [0 0], 'k:')
plot([0 0], [-2.5 2.5], 'k')
axis tight
set(gca,'TickDir','out')
ax1 = subplot(3,1,1);
matrix2plot = [matrix_play; matrix_noplay];
imagesc(time_analysis, 1:size(matrix2plot,1), matrix2plot)
axis xy
clim([-4 4])
hold on
plot(time_analysis([1 end]), [1 1]*size(matrix_play,1), 'w', 'LineWidth',2)
plot([0 0], [0 1]*size(matrix2plot,1), 'w', 'LineWidth',2)
yticks(1:size(matrix2plot,1))
yticklabels(all_behaviors_2x([play_indexes, noplay_indexes]))
set(gca,'TickDir','out')
colorbar
ax2.Position(1) = ax1.Position(1);
ax2.Position(3) = ax1.Position(3);

%% Plot bar and scatter comparison (onset-aligned)
figure
name_psth_figure(band_name, 'onset_bars')
subplot(1,2,1)
y_lim = [-5 5];
time2compare = [-.2 .5];
play_bars = mean(matrix_play(:, time_analysis>time2compare(1) & time_analysis<time2compare(2)), 2);
non_play_bars = mean(matrix_noplay(:, time_analysis>time2compare(1) & time_analysis<time2compare(2)), 2);
non_play_p_val_window = any(pval_noplay(:, time_analysis>time2compare(1) & time_analysis<time2compare(2))<0.05, 2);
[non_play_p_val_window, non_play_order] = sort(non_play_p_val_window);
play_p_val_window = any(pval_play(:, time_analysis>time2compare(1) & time_analysis<time2compare(2))<0.05, 2);
[play_p_val_window, play_order] = sort(play_p_val_window);
bar_p_val = [non_play_p_val_window' play_p_val_window'];
bar_values = [non_play_bars(non_play_order)' play_bars(play_order)'];
hold on
for i = 1:length(bar_values)
    if bar_p_val(i) == 1
        bar(i, bar_values(i), 'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'k');
    else
        bar(i, bar_values(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
end
xlim([0 length(bar_values)+1]);
xticks(1:numel([non_play_bars' play_bars']))
xticklabels(all_behaviors_2x([noplay_indexes(non_play_order), play_indexes(play_order)]))
ylim(y_lim)
plot([1 1]*numel(non_play_bars)+.5, y_lim, 'k')

subplot(1,2,2)
rand_x = (rand(size(bar_values)) - .5)*.5 + [non_play_bars'*0 play_bars'*0+1];
plot(rand_x, bar_values, 'k.', 'MarkerSize', 5)
hold on
plot(rand_x(bar_p_val), bar_values(bar_p_val), 'r.', 'MarkerSize', 10)
plot([-.25 .25], mean(non_play_bars)*[1 1], 'r', 'MarkerSize', 8)
plot([-.25 .25]+1, mean(play_bars)*[1 1], 'r', 'MarkerSize', 8)
[p, ~, stats] = ranksum(non_play_bars, play_bars);
title([num2str(stats.zval), ' ', num2str(stats.ranksum), ' ', num2str(p)])
ylim([-2.5 3.5])

%% Plot self vs other — play behaviors
figure
name_psth_figure(band_name, 'onset_self_other_play')
index_self = group1;
index_other = 1:numel(play_indexes);
index_other(index_self) = [];
self_matrix = matrix_play(index_self,:);
other_matrix = matrix_play(index_other,:);
[~, ~, ci] = ttest(self_matrix);
fill([time_analysis fliplr(time_analysis)], [ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean(self_matrix), 'b')
[~, ~, ci] = ttest(other_matrix);
fill([time_analysis fliplr(time_analysis)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(time_analysis, mean(other_matrix), 'r')
plot(time_analysis([1 end]), [0 0], 'k:')
plot([0 0], [-2 4], 'k')
axis tight
set(gca,'TickDir','out')
mark_significant_psth(time_analysis, self_matrix, other_matrix)
print(gcf,'-vector','-dsvg',fullfile(figure_2_new_folder, [get(gcf, 'Name') '.svg']))

%% Plot self vs other — non-play behaviors
figure
name_psth_figure(band_name, 'onset_self_other_nonplay')
no_play_indexes = [group2, group2+numel(all_behaviors)];
index_self = 1:numel(group2);
index_other = 1:numel(no_play_indexes);
index_other(index_self) = [];
self_matrix = matrix_noplay(index_self,:);
other_matrix = matrix_noplay(index_other,:);
[~, ~, ci] = ttest(self_matrix);
fill([time_analysis fliplr(time_analysis)], [ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis, mean(self_matrix), 'b')
[~, ~, ci] = ttest(other_matrix);
fill([time_analysis fliplr(time_analysis)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(time_analysis, mean(other_matrix), 'r')
plot(time_analysis([1 end]), [0 0], 'k:')
plot([0 0], [-4 2], 'k')
axis tight
set(gca,'TickDir','out')
mark_significant_psth(time_analysis, self_matrix, other_matrix)

%% LME analysis — time-wrapped PSTH
% Skipped when results_full_lfp_per_behavior_time_warpped_<band>.mat is already saved.
if exist(fullfile(saving_folder, results_tw_file), 'file') == 2
    disp(['time-wrapped LME already saved, loading ' results_tw_file])
    load(fullfile(saving_folder, results_tw_file), ...
        'results_time_wrapped', 'all_behaviors_2x', 'group1', 'group2', 'time_analysis_time_wrapped', 'all_behaviors')
else
baseline_idx = time_wrap_time >= -5 & time_wrap_time <= 0;
analysis_idx = time_wrap_time >= -5 & time_wrap_time <= 10;
time_analysis_time_wrapped = time_wrap_time(analysis_idx);
results_time_wrapped = struct;

for g = 1:2
    beh_list = beh_groups{g};
    nBeh = length(beh_list);
    all_y = [];
    all_session = [];
    all_behavior = [];

    for b = beh_list
        psth = all_psth_onset_self{3,b};
        sess = all_session_index_self{b};
        elec = all_electrode_index_self{b};
        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);
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

    for b = beh_list
        psth = all_psth_onset_other{3,b};
        sess = all_session_index_other{b};
        elec = all_electrode_index_other{b};
        idx = elec == 1;
        psth = psth(idx,:);
        sess = sess(idx);
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

    self_other_beh_list = [beh_list, beh_list + numel(all_behaviors)];
    nTime = length(time_analysis_time_wrapped);
    beta = nan(2*nBeh, nTime);
    beta_std = nan(2*nBeh, nTime);
    p_beh = nan(2*nBeh, nTime);
    p_global = nan(1, nTime);

    for t = 1:nTime
        disp(num2str(time_analysis_time_wrapped(t)))
        y = all_y(:,t);
        valid = ~isnan(y);
        tbl = table(y(valid), categorical(all_session(valid)), categorical(all_behavior(valid)), 'VariableNames', {'y','session','behavior'});
        if numel(unique(tbl.session)) < 2
            continue;
        end
        lme_g = fitlme(tbl, 'y ~ 1 + (1|session)');
        p_global(t) = lme_g.Coefficients.pValue;
        for i = 1:2*nBeh
            beh_id = self_other_beh_list(i);
            tbl_i = tbl(tbl.behavior == categorical(beh_id), :);
            if isempty(tbl_i)
                continue;
            end
            lme_b = fitlme(tbl_i, 'y ~ 1 + (1|session)');
            beta(i,t) = lme_b.Coefficients.Estimate;
            p_beh(i,t) = lme_b.Coefficients.pValue;
            se = lme_b.Coefficients.SE;
            if se > 0
                beta_std(i,t) = beta(i,t) / se;
            end
        end
    end

    results_time_wrapped(g).beta = beta;
    results_time_wrapped(g).beta_std = beta_std;
    results_time_wrapped(g).p_beh = p_beh;
    results_time_wrapped(g).p_global = p_global;
    results_time_wrapped(g).behaviors = self_other_beh_list;
end
end

%% Save time-wrapped LME
% This section only writes the workspace. Run it alone after a fit.
disp(['saving ' results_tw_file])
save(fullfile(saving_folder, results_tw_file), ...
    'results_time_wrapped', 'all_behaviors_2x', 'group1', 'group2', 'time_analysis_time_wrapped', 'all_behaviors', 'band_name', '-v7.3')

%% Plot colormap and mean responses (time-wrapped)
matrix_play = results_time_wrapped(1).beta_std;
pval_play = results_time_wrapped(1).p_beh;
matrix_play(~any(pval_play<alpha_tresh,2),:) = NaN;
matrix_noplay = results_time_wrapped(2).beta_std;
pval_noplay = results_time_wrapped(2).p_beh;
matrix_noplay(~any(pval_noplay<alpha_tresh,2),:) = NaN;

figure
name_psth_figure(band_name, 'timewrap_play_nonplay')
colormap(jet(256))
ax2 = subplot(3,1,2:3);
[~, ~, ci] = ttest(matrix_noplay);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)], [ci(1,:) fliplr(ci(2,:))], 'k', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis_time_wrapped, mean(matrix_noplay), 'k')
[~, ~, ci] = ttest(matrix_play);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(time_analysis_time_wrapped, mean(matrix_play), 'r')
plot(time_analysis_time_wrapped([1 end]), [0 0], 'k:')
plot([0 0], [-5 5], 'k')
plot([5 5], [-5 5], 'k')
axis tight
set(gca,'TickDir','out')
ax1 = subplot(3,1,1);
matrix2plot = [matrix_play; matrix_noplay];
imagesc(time_analysis_time_wrapped, 1:size(matrix2plot,1), matrix2plot)
axis xy
clim([-10 10])
hold on
plot(time_analysis_time_wrapped([1 end]), [1 1]*size(matrix_play,1), 'w', 'LineWidth',2)
plot([0 0], [0 1]*size(matrix2plot,1), 'w', 'LineWidth',2)
yticks(1:size(matrix2plot,1))
yticklabels(all_behaviors_2x([play_indexes, noplay_indexes]))
set(gca,'TickDir','out')
colorbar
ax2.Position(1) = ax1.Position(1);
ax2.Position(3) = ax1.Position(3);

%% Plot bar and scatter comparison (time-wrapped)
alpha_level = 0.01;
figure
name_psth_figure(band_name, 'timewrap_bars')
subplot(1,2,1)
y_lim = [-5 5];
time2compare = [0 5];
play_bars = mean(matrix_play(:, time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2)), 2);
non_play_bars = mean(matrix_noplay(:, time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2)), 2);
non_play_p_val_window = any(pval_noplay(:, time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2))<alpha_level, 2);
[non_play_p_val_window, non_play_order] = sort(non_play_p_val_window);
play_p_val_window = any(pval_play(:, time_analysis_time_wrapped>time2compare(1) & time_analysis_time_wrapped<time2compare(2))<alpha_level, 2);
[play_p_val_window, play_order] = sort(play_p_val_window);
bar_p_val = [non_play_p_val_window' play_p_val_window'];
bar_values = [non_play_bars(non_play_order)' play_bars(play_order)'];
hold on
for i = 1:length(bar_values)
    if bar_p_val(i) == 1
        bar(i, bar_values(i), 'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'k');
    else
        bar(i, bar_values(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
end
xlim([0 length(bar_values)+1]);
xticks(1:numel([non_play_bars' play_bars']))
xticklabels(all_behaviors_2x([noplay_indexes(non_play_order), play_indexes(play_order)]))
ylim(y_lim)
plot([1 1]*numel(non_play_bars)+.5, y_lim, 'k')

subplot(1,2,2)
rand_x = (rand(size(bar_values)) - .5)*.5 + [non_play_bars'*0 play_bars'*0+1];
plot(rand_x, bar_values, 'k.', 'MarkerSize', 5)
hold on
plot(rand_x(bar_p_val), bar_values(bar_p_val), 'r.', 'MarkerSize', 10)
plot([-.25 .25], mean(non_play_bars)*[1 1], 'r', 'MarkerSize', 8)
plot([-.25 .25]+1, mean(play_bars)*[1 1], 'r', 'MarkerSize', 8)
[p, ~, stats] = ranksum(non_play_bars, play_bars);
title([num2str(stats.zval), ' ', num2str(stats.ranksum), ' ', num2str(p)])
ylim(y_lim)

%% Plot self vs other — play behaviors, time-wrapped
figure
name_psth_figure(band_name, 'timewrap_self_other_play')
index_self = group1;
index_other = 1:numel(play_indexes);
index_other(index_self) = [];
self_matrix = matrix_play(index_self,:);
other_matrix = matrix_play(index_other,:);
[~, ~, ci] = ttest(self_matrix);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)], [ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis_time_wrapped, mean(self_matrix), 'b')
[~, ~, ci] = ttest(other_matrix);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(time_analysis_time_wrapped, mean(other_matrix), 'r')
plot(time_analysis_time_wrapped([1 end]), [0 0], 'k:')
plot([0 0], [-5 5], 'k')
plot([5 5], [-5 5], 'k')
axis tight
set(gca,'TickDir','out')
mark_significant_psth(time_analysis_time_wrapped, self_matrix, other_matrix)

%% Plot self vs other — non-play behaviors, time-wrapped
figure
name_psth_figure(band_name, 'timewrap_self_other_nonplay')
index_self = 1:numel(group2);
index_other = 1:numel(noplay_indexes);
index_other(index_self) = [];
self_matrix = matrix_noplay(index_self,:);
other_matrix = matrix_noplay(index_other,:);
[~, ~, ci] = ttest(self_matrix);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)], [ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.2, 'EdgeColor','none')
hold on
plot(time_analysis_time_wrapped, mean(self_matrix), 'b')
[~, ~, ci] = ttest(other_matrix);
fill([time_analysis_time_wrapped fliplr(time_analysis_time_wrapped)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.2, 'EdgeColor','none')
plot(time_analysis_time_wrapped, mean(other_matrix), 'r')
plot(time_analysis_time_wrapped([1 end]), [0 0], 'k:')
plot([0 0], [-5 5], 'k')
plot([5 5], [-5 5], 'k')
axis tight
set(gca,'TickDir','out')
mark_significant_psth(time_analysis_time_wrapped, self_matrix, other_matrix)
end

function name_psth_figure(band_name, plot_id)
% Window title is the save name. Same plot_id for delta, theta, and gamma.
set(gcf, 'NumberTitle', 'off', 'Name', ['psth_' band_name '_' plot_id]);
end

function mark_significant_psth(time_axis, self_matrix, other_matrix)
% Strips above the axes. Blue is self, red is other.
% Darker rectangles: mean significantly above 0. Lighter: significantly below 0.
yl = ylim;
pad = 0.045 * diff(yl);
if ~(pad > 0)
    pad = 0.1;
end
ylim([yl(1), yl(2) + 2.6 * pad]);
shade_one_trace(time_axis, self_matrix, [0.10 0.25 0.80], yl(2) + 1.35 * pad, pad);
shade_one_trace(time_axis, other_matrix, [0.80 0.15 0.15], yl(2) + 0.15 * pad, pad);
end

function shade_one_trace(time_axis, values, color, y0, height)
[above, below] = trace_vs_zero(values);
paint_mask(time_axis, above, y0 + 0.55 * height, 0.45 * height, color, 0.75);
paint_mask(time_axis, below, y0, 0.45 * height, color, 0.30);
end

function [above, below] = trace_vs_zero(values)
n_time = size(values, 2);
above = false(1, n_time);
below = false(1, n_time);
mean_trace = mean(values, 1, 'omitnan');
for t = 1:n_time
    sample = values(:, t);
    sample = sample(~isnan(sample));
    if numel(sample) < 2 || std(sample) == 0
        continue
    end
    [~, p] = ttest(sample, 0);
    if p < 0.05 && mean_trace(t) > 0
        above(t) = true;
    elseif p < 0.05 && mean_trace(t) < 0
        below(t) = true;
    end
end
end

function paint_mask(time_axis, mask, y0, height, color, face_alpha)
if ~any(mask) || numel(time_axis) < 2
    return
end
dt = median(diff(time_axis));
edges = [false, mask(:)', false];
starts = find(diff(edges) == 1);
stops = find(diff(edges) == -1) - 1;
for k = 1:numel(starts)
    x0 = time_axis(starts(k)) - dt / 2;
    x1 = time_axis(stops(k)) + dt / 2;
    patch([x0 x1 x1 x0], [y0 y0 y0 + height y0 + height], color, ...
        'FaceAlpha', face_alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end
