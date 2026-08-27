%% Fig A1–A7 — behavior modulation by cell type and area
% Sections A1–A7 from Play_map_Dynamic_GLM_Wrap_IndexRes_PhaseComp.m
% Run independently after Estiamte_activation_index_and_other_variables.m
% has saved activation_index_A1_A7.mat.

checkpoint_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\activation index';
a17_file = [checkpoint_folder, '\activation_index_A1_A7.mat'];
figure_folder_peakandtrough = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Peak trough dyanmics';
if ~exist(figure_folder_peakandtrough, 'dir')
    mkdir(figure_folder_peakandtrough);
end

disp('Loading A1–A7 checkpoint')
load(a17_file);
disp('Checkpoint loaded')

%% A1) ploing neural response per behavior and estiamting significnative resposne per behavior and cell type
% Heatmaps of median z-scored rate (0–5 s) per neuron × behavior, split by
% peak / trough / unlocked. Also runs the permutation tests used later:
%   p_val_behavior_shuffle  — shuffle behavior order within each neuron
%   p_val_labelshuffle_*    — shuffle peak/trough/unlocked labels
%   p_val_comodulation      — whether the count of significant behaviors
%                             is itself unusual
% bar_responses is the median across neurons (one value per type/area/behavior).

response_time = [0 5];
response_time_index = time2use >= response_time(1) & time2use <= response_time(2);

selected_behaviors = [1 4 5 6 7 8 9 14 15 18 19 20 21 24 25];
play_index_order   = [1 4 14 18 20 24 6 8 5 15 19 21 25 7 9 10 12 16 22 26 11 13 17 23 27];
play_border        = numel(selected_behaviors);
self_play_order    = 8;
sel_other_no_play  = 20;
all_borders        = [1 self_play_order play_border sel_other_no_play];
c_lim              = [-5 5];
matrix_type        = 4; % 2 or 4

bar_responses                   = nan(3, numel(area_list), numel(play_index_order));
p_val_behavior_shuffle         = nan(3, numel(area_list), numel(play_index_order));
activation_by_area_by_type      = {};
n_behaviors                     = numel(play_index_order);


p_val_labelshuffle_one_sided = nan(3, numel(area_list), n_behaviors);
p_val_labelshuffle_two_sided = nan(3, numel(area_list), n_behaviors);
logit_labelshuffle_one_sided  = nan(3, numel(area_list), n_behaviors);
logit_labelshuffle_two_sided  = nan(3, numel(area_list), n_behaviors);
p_val_comodulation            = nan(3, numel(area_list));


% ---- main loop ----
promedio_per_area = [];
for k = 1:numel(area_list)
    all_neurons = size(all_activation_order{1,k,2}, 1);
    all_activation = nan(all_neurons, size(all_activation_order,1), ...
                         size(all_activation_order{1,k,matrix_type}, 2));

    for bn = 1:size(all_activation_order,1)
        all_activation(:,bn,:) = all_activation_order{bn,k,matrix_type};
    end

    % logical mask for behaviors
    play_index = play_index_order;
    % play_index(selected_behaviors) = true;

    % define cell groups
    cell_groups = {
        all_activation_order{1,k,1}==1, 'Peak cells';
        all_activation_order{1,k,1}==2, 'Trough cells';
        all_activation_order{1,k,1}==3, 'Non-entrained cells';
    };

    figure('units','normalized','outerposition',[0 0 .5 1]);
    session_list = unique(all_activation_order{1,k,3});
    promedio_per_type =[];
   stacked_labels = {};
   staked_session_labels = {};
    for ci = 1:3
        subplot(3,1,ci)
        cells_mask = cell_groups{ci,1};
        y_labels   = all_activation_order{1,k,3}(cells_mask);
        matrix2plot = squeeze(median(all_activation(cells_mask, play_index_order, response_time_index), 3, 'omitmissing'));
        promedio_per_session = [];
        for sn=1:numel(session_list)
            sub_matrix = matrix2plot(ismember(y_labels, session_list{sn}),:);
            
            if size(sub_matrix,1)==1
                promedio = sub_matrix;
            elseif size(sub_matrix,1)>1
                 promedio = median(matrix2plot(ismember(y_labels, session_list{sn}),:));
            else
                promedio = nan(1,size(matrix2plot,2));
            end

        
            if numel(promedio)<size(matrix2plot,2)
                disp(numel(promedio))
                promedio = nan(1,size(matrix2plot,2));
            end
             promedio_per_session = [promedio_per_session;promedio'];

             if ci==1
                 stacked_labels  = [stacked_labels;labels(play_index_order)'];
                 staked_session_labels = [staked_session_labels; repmat(session_list(sn),numel(play_index_order),1)];
             end
        end 
        promedio_per_type = [promedio_per_type,promedio_per_session];
        true_median = median(matrix2plot, 'omitmissing');
        bar_responses(ci,k,:) = true_median;

        % --- permutation test (circular shift per neuron) ---
        [n_neurons,~] = size(matrix2plot);
        null_medians = nan(n_perm, n_behaviors);

        for p = 1:n_perm
            shifted = matrix2plot;
            for n = 1:n_neurons
                % shift_amount = randi(n_behaviors);
                % shifted(n,:) = circshift(shifted(n,:), shift_amount, 2);
                 shifted(n,:) = matrix2plot(n, randperm(n_behaviors));
            end
            null_medians(p,:) = median(shifted, 'omitmissing');
        end

        % --- compute p-values (two-tailed) ---
        pvals = mean(abs(null_medians) >= abs(true_median), 1);
        p_val_behavior_shuffle(ci,k,:) = pvals;

        % --- compute co-modulation p-values (one-tailed) ---
        n_significant = nan(n_perm,1);

        for p =  1:n_perm

            n_significant(p) = sum(mean(abs(null_medians) >= abs(null_medians(p,:)), 1)<0.05);
        end

        p_val_comodulation(k,ci) =  mean(sum(pvals<0.05)<n_significant);



        plot_matrix(matrix2plot, y_labels, c_lim, all_borders)
        title(cell_groups{ci,2})
    end



    all_lables = all_activation_order{1,k,1};

    shufl_median = nan(n_perm,n_behaviors,3);
    for p = 1:n_perm
        shuf_lables = randsample(all_lables, numel(all_lables));
        for ci = 1:3
            cells_mask = all_lables==ci;
            matrix2plot = squeeze(median(all_activation(cells_mask, play_index_order, response_time_index), 3, 'omitmissing'));
            shufl_median(p,:,ci) =   median(matrix2plot, 'omitmissing');
        end
    end


    for ci = 1:3
        for nb = 1:numel(area_list)
            null_dist = squeeze(shufl_median(:,:,ci));         % [n_perm × n_behaviors]

            % --- p-values ---
            p1 = mean(null_dist >= true_median(nb)', 1, 'omitnan');
            p2 = mean(abs(null_dist) >= abs(true_median(nb)'), 1, 'omitnan');

            % --- clip to avoid log(0) ---
            eps_clip = 1/(size(null_dist,1) + 1);
            p1 = min(max(p1, eps_clip), 1 - eps_clip);
            p2 = min(max(p2, eps_clip), 1 - eps_clip);

            % --- store ---
            p_val_labelshuffle_one_sided(ci,nb,:) = p1;
            p_val_labelshuffle_two_sided(ci,nb,:) = p2;
            logit_labelshuffle_one_sided(ci,nb,:) = log(p1 ./ (1 - p1));
            logit_labelshuffle_two_sided(ci,nb,:) = log(p2 ./ (1 - p2));
        end
    end




    promedio_per_area = [promedio_per_area;[num2cell(promedio_per_type),stacked_labels,staked_session_labels, repmat(area_list(k),numel(staked_session_labels),1)]];
    sgtitle(area_list{k})
     print(gcf,'-vector','-dsvg',[figure_folder_peakandtrough,'\', area_list{k}, 'LPAG peak and trough activation for each behavior.svg'])
     pause(.1)

end
promedio_per_area = cell2table(promedio_per_area);
promedio_per_area.Properties.VariableNames = {'Peak','Trough','NonEntrained','Behavior','Session', 'Area'};



%% A2) plotign all bars
% Bar plot of bar_responses for every area (rows) and cell type (columns).
% Red lines mark play / self-play / non-play borders; green shading is locomotive play.
figure('units','normalized','outerposition',[0 0 1 1]);
y_lim = [-2 3];
locomotive_border = [1.5  6.5; 8.5 13.5];
type_labels = {'Peak','Trough','NonModulated'};
for k = 1:numel(area_list)
    for type =1:3

        subplot(numel(area_list), 3, (k-1)*3 + type)
        bar(squeeze(bar_responses(type, k, :)), 'k')
       
        hold on
         for nb=1:numel(all_borders)
            plot(all_borders(nb)*[1 1]+.5, y_lim,  'r', 'LineWidth', 2);
         end

         for nb=1:size(locomotive_border,1)
             fill([locomotive_border(nb,[1 2 2 1])], y_lim([1 1 2 2]),'g', 'FaceAlpha', .1, 'EdgeColor', 'none' )
         end
         
         axis tight
         if type==1
             ylabel({area_list{k},'Median modualtion (zs)'})
         end
          xticks(1:numel(play_index_order))
         if k==1
             title(type_labels{type})
         elseif k== numel(area_list)
            
             xticklabels(labels(play_index_order))
         else
             xticklabels([])
         end
    end
end
%% A3) and plot only LPAG
% Same bars as A2, but only LPAG trough cells (area k = 3, type = 2).
figure
k=3;
type =2;
 bar(squeeze(bar_responses(type, k, :)), 'k')
       
        hold on
         for nb=1:numel(all_borders)
            plot(all_borders(nb)*[1 1]+.5, y_lim,  'r', 'LineWidth', 2);
         end

         for nb=1:size(locomotive_border,1)
             fill([locomotive_border(nb,[1 2 2 1])], y_lim([1 1 2 2]),'g', 'FaceAlpha', .1, 'EdgeColor', 'none' )
         end
         
         axis tight
         
             ylabel({area_list{k},'Median modualtion (zs)'})
     
          xticks(1:numel(play_index_order))
        
             title(type_labels{type})
        
            
             xticklabels(labels(play_index_order))
         
%%  A4) plot  p values
% Same layout as A2, but bars are filled when p < alpha_level and empty otherwise.
% p_val_2use chooses which permutation p-value colors the bars (and is reused in A5).
% Options from A1:
%   p_val_behavior_shuffle      — shuffle behaviors within each neuron (default)
%   p_val_labelshuffle_two_sided / _one_sided — shuffle cell-type labels
%   p_val_comodulation          — unusual number of significant behaviors
% bar_height is usually bar_responses; logit_labelshuffle_* is an alternative y-scale.

figure('units','normalized','outerposition',[0 0 1 1]);
y_lim = [-3 3];
locomotive_border = [1.5 6.5; 8.5 13.5];
% border_colors = {'','','',''}
type_labels = {'Peak','Trough','NonModulated'};
alpha_level = 0.05; % significance threshold
p_val_2use = p_val_behavior_shuffle; % p_val_behavior_shuffle | p_val_comodulation | p_val_labelshuffle_two_sided
bar_height = bar_responses; % logit_labelshuffle_one_sided; bar_responses;

for k = 1:numel(area_list)
    for type = 1:3
        subplot(numel(area_list), 3, (k-1)*3 + type)
        hold on

        % extract data and p-values
        vals = squeeze(bar_height(type, k, :));
        pvals = squeeze(p_val_2use(type, k, :));

        % determine which are significant
        sig_mask = pvals < alpha_level;

        % plot each bar individually
        for b = 1:numel(vals)
            if sig_mask(b)
                bar(b, vals(b), 'k', 'EdgeColor', 'none');  % filled black
            else
                bar(b, vals(b), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2); % empty with black edge
            end
        end

        % add borders
        for nb = 1:numel(all_borders)
            plot(all_borders(nb)*[1 1]+.5, y_lim, 'r', 'LineWidth', 2);
        end

        % highlight locomotive behaviors
        for nb = 1:size(locomotive_border,1)
            fill(locomotive_border(nb,[1 2 2 1]), y_lim([1 1 2 2]), ...
                'g', 'FaceAlpha', .1, 'EdgeColor', 'none');
        end

        % formatting
        axis tight
        ylim(y_lim)
        if type==1
            ylabel({area_list{k}, 'Median modulation (zs)'})
        end
        xticks(1:numel(play_index_order))
        if k==1
            title(type_labels{type})
        elseif k==numel(area_list)
            xticklabels(labels(play_index_order))
        else
            xticklabels([])
        end
    end
end

%% A5) now plot number of significnat behaviors
% For play vs non-play behaviors, stacked % of behaviors that are
% non-significant (black), significant negative (blue), or significant
% positive (red). Uses p_val_2use and bar_responses from A4 / A1.

 play_indexes = [1 4 14 18 20 24 6 8 5 15 19 21 25] ; no_play_indexes = [7 9 10 12 16 22 26 11 13 17 23 27];

counts_play = squeeze(sum(p_val_2use(:,:,ismember(play_index_order,play_indexes))<Inf,3));
play_counts = squeeze(sum(p_val_2use(:,:,ismember(play_index_order,play_indexes))<0.05,3));
play_counts_positive = squeeze(sum(p_val_2use(:,:,ismember(play_index_order,play_indexes))<0.05 & bar_responses(:,:,ismember(play_index_order,play_indexes))>0,3));
counts_noplay = squeeze(sum(p_val_2use(:,:,ismember(play_index_order,no_play_indexes))<Inf,3));
nonplay_counts = squeeze(sum(p_val_2use(:,:,ismember(play_index_order,no_play_indexes))<0.05,3));
nonplay_counts_positive = squeeze(sum(p_val_2use(:,:,ismember(play_index_order,no_play_indexes))<0.05 & bar_responses(:,:,ismember(play_index_order,no_play_indexes))>0 ,3));

figure
for  k = 1:numel(area_list)
    subplot(numel(area_list), 2, 2*(k-1) + 1)

    this_area_bar = 100*fliplr([counts_play(:,k)-play_counts(:,k) play_counts(:,k)- play_counts_positive(:,k) play_counts_positive(:,k)]./counts_play(:,k));
    this_area_bar = flipud(this_area_bar);
    b =barh(this_area_bar, 'stacked');
    b(1).FaceColor = [1 0 0];
    b(2).FaceColor = [0 0 1];
    b(3).FaceColor = [0 0 0];
    axis tight
    yticks(1:3)
    yticklabels({'NonEntrained','Trough','Peak'})
    if k==numel(area_list)
        xlabel('% of modualted behaviors')
    else
        xticklabels([])
    end
    if k==1
        title(['Play Behaviors (', num2str(numel(play_indexes)), ')'])
    end

    subplot(numel(area_list), 2, 2*(k-1) + 2)

    this_area_bar = 100*fliplr([counts_noplay(:,k)-nonplay_counts(:,k) nonplay_counts(:,k)- nonplay_counts_positive(:,k) nonplay_counts_positive(:,k)]./counts_noplay(:,k));
    this_area_bar = flipud(this_area_bar);
    b =barh(this_area_bar, 'stacked');
    b(1).FaceColor = [1 0 0];
    b(2).FaceColor = [0 0 1];
    b(3).FaceColor = [0 0 0];
    axis tight
    yticks(1:3)
    yticklabels([])
    if k==numel(area_list)
        xlabel('% of modualted behaviors')
    else
        xticklabels([])
    end
    if k==1
        title(['Non play Behaviors (', num2str(numel(no_play_indexes)), ')'])
    end
end
%% A6 )now pl;oting grouped behaviors (lik eplay no play etc..)
% Recompute median responses and permutation p-values for grouped behavior
% sets (play, no-play, self/other). Overwrites bar_responses and
% p_val_behavior_shuffle. group_means / group_pval are what A7 plots.
% n_perm is set to 1000 here (A1 used the checkpoint n_perm).

% ---- parameters ----
response_time = [0 5];
response_time_index = time2use >= response_time(1) & time2use <= response_time(2);

selected_behaviors = [1 4 5 6 7 8 9 14 15 18 19 20 21 24 25];
play_index_order   = [1 4 14 18 20 24 6 8 5 15 19 21 25 7 9 10 12 16 22 26 11 13 17 23 27];
play_border        = numel(selected_behaviors);
self_play_order    = 8;
sel_other_no_play  = 20;
all_borders        = [1 self_play_order play_border sel_other_no_play];
c_lim              = [-5 5];
matrix_type        = 4; % 2 or 4
n_perm             = 1000;

% ---- NEW: define behavior groups for grouped mean analysis ----
% Example: each cell array entry is a group of behavior indices
group_sets = {2:15, 16:25,2:8, 9:15, 16:20, 21:25}

bar_responses = nan(3, numel(area_list), numel(play_index_order));
p_val_behavior_shuffle         = nan(3, numel(area_list), numel(play_index_order));

% ---- NEW OUTPUTS for group analysis ----
group_means = nan(3, numel(area_list), numel(group_sets));
group_pval  = nan(3, numel(area_list), numel(group_sets));

% ---- helper function for plotting ----

% ---- main loop ----
for k = 1:numel(area_list)
    all_neurons = size(all_activation_order{1,k,2}, 1);
    all_activation = nan(all_neurons, size(all_activation_order,1), ...
                         size(all_activation_order{1,k,matrix_type}, 2));
    for bn = 1:size(all_activation_order,1)
        all_activation(:,bn,:) = all_activation_order{bn,k,matrix_type};
    end

    play_index = play_index_order;
    % play_index(selected_behaviors) = true;

    cell_groups = {
        all_activation_order{1,k,1}==1, 'Peak cells';
        all_activation_order{1,k,1}==2, 'Trough cells';
        all_activation_order{1,k,1}==3, 'Non-active cells';
    };

    for ci = 1:3
        cells_mask = cell_groups{ci,1};
        matrix2plot = squeeze(median(all_activation(cells_mask, play_index_order, response_time_index), 3, 'omitmissing'));
        true_median = median(matrix2plot, 'omitmissing');
        bar_responses(ci,k,:) = true_median;

        % ------------------------
        % 1️⃣ per-behavior permutation (same as before)
        % ------------------------
        [n_neurons, n_behaviors] = size(matrix2plot);
        null_medians = nan(n_perm, n_behaviors);
        for p = 1:n_perm
            shifted = matrix2plot;
            for n = 1:n_neurons
                shift_amount = randi(n_behaviors);
                shifted(n,:) = circshift(shifted(n,:), shift_amount, 2);
            end
            null_medians(p,:) = median(shifted, 'omitmissing');
        end
        p_val_behavior_shuffle(ci,k,:) = mean(abs(null_medians) >= abs(true_median), 1);

        % ------------------------
        % 2️⃣ NEW: grouped-behavior permutation
        % ------------------------
        for g = 1:numel(group_sets)
            idx = group_sets{g};
            n_group = numel(idx);

            % actual mean response (across behaviors)
            true_group_mean = mean(true_median(idx), 'omitmissing');
            group_means(ci,k,g) = true_group_mean;

            % permuted means: shuffle behavior labels (not circular)
            null_group_means = nan(n_perm, 1);
            for p = 1:n_perm
                rand_idx = randperm(n_behaviors, n_group);
                null_group_means(p) = mean(true_median(rand_idx), 'omitmissing');
            end

            % two-tailed p-value
            group_pval(ci,k,g) = mean(abs(null_group_means) >= abs(true_group_mean));
        end
    end
end
%% A7) and plot previous sectio
% Bars of group_means from A6; filled if group_pval < alpha_level.
figure
y_lim = [-2 3];
type_labels = {'Peak','Trough','NonModulated'};
alpha_level = 0.05; % significance threshold

% --- Define labels for your grouped behaviors ---
% (Must match number of entries in `group_sets`)
group_labels = ["play", "No play", "Self Play","Other Play", "Self No-play","Other No-play"];  

n_groups = numel(group_sets);

for k = 1:numel(area_list)
    for type = 1:3
        subplot(numel(area_list), 3, (k-1)*3 + type)
        hold on

        % extract data for this area/type
        vals  = squeeze(group_means(type, k, :));
        pvals = squeeze(group_pval(type, k, :));

        % plot bars for each group
        for g = 1:n_groups
            if pvals(g) < alpha_level
                bar(g, vals(g), 'k', 'EdgeColor', 'none');  % significant → filled
            else
                bar(g, vals(g), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.2); % non-sig → outline
            end
        end

        % add reference lines / formatting
        yline(0, 'k:')  % baseline
        ylim(y_lim)
        axis tight

        % axis labeling
        xticks(1:n_groups)
        xticklabels(group_labels)
        xtickangle(30)
        if type == 1
            ylabel({area_list{k}, 'Mean modulation (z)'})
        end
        if k == 1
            title(type_labels{type})
        end

        % optional: significance markers
        for g = 1:n_groups
            if pvals(g) < alpha_level
                text(g, vals(g) + 0.1, '*', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
            end
        end
    end
end

sgtitle('Grouped Behavior Responses Across Areas')


function plot_matrix(matrix2plot, y_labels, c_lim, all_borders)
    colormap(jet)
    imagesc(matrix2plot);
    hold on
    clim(c_lim)

    % Find group boundaries
    change_idx = find(~strcmp(y_labels(1:end-1), y_labels(2:end))) + 0.5;
    group_edges = [0.5; change_idx(:); numel(y_labels)+0.5];
    group_centers = movmean(group_edges, 2, 'Endpoints', 'discard');
    [unique_labels, ~] = unique(y_labels, 'stable');

    % Plot separators
    for i = 1:numel(change_idx)
        plot([1 size(matrix2plot,2)], [change_idx(i) change_idx(i)], ':w', 'LineWidth', 1);
    end
    for nb = 1:numel(all_borders)
        plot(all_borders(nb)*[1 1] + .5, [1 size(matrix2plot,1)], 'w', 'LineWidth', 2);
    end

    yticks(group_centers);
    yticklabels(unique_labels);
end
