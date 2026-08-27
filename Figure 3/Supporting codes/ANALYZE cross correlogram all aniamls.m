%% ANALYZE cross-correlograms (all animals)
% Plots spike-train cross-correlograms (CCGs) pooled across sessions,
% then splits pairs by area and by delta-phase cell type (peak / trough /
% unlocked). Play vs non-play; shuffle CCG is the within-pair control.

%% 1 Load cross-correlograms
% Per-session structs from Estimate cross correlogram all aniamls
% (play, non-play, entire session, and percentile/shuffle CCGs).
% Structs are under Data\Analysis results\Cross_correlogram.

data_root = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder = [data_root, '\Analysis results\Cross_correlogram'];
% mkdir(saving_folder)

bin_size_cc = 0.01;
hist_range_cc = [-3 3];

psth_edges_cc = hist_range_cc(1):bin_size_cc:hist_range_cc(2);

load([saving_folder,'\cross_corr_struct_structure_updated_with_non_playbouts_3sec.mat'],'cross_corr_struct');
load([saving_folder,'\cross_corr_struct_animal_names_updated_with_non_playbouts_3sec.mat'],'animal_names');

%% 2 Merge pairs across sessions
% Concatenate CCGs, shuffle CCGs, area labels, session, and cluster IDs.
entire_session_cc = [];
play_cc = [];
non_play_cc = [];
all_area_comb_cc = [];
session_id_cc = [];
clusters_id_cc = [];
all_cross_corr_play_pctl = [];
all_cross_corr_non_play_pctl = [];
for fn = 1:numel(cross_corr_struct)

    entire_session_cc            = [entire_session_cc;cross_corr_struct(fn).all_cross_corr_entire_session];
    all_cross_corr_play_pctl     = [all_cross_corr_play_pctl;cross_corr_struct(fn).all_cross_corr_play_pctl];
    all_cross_corr_non_play_pctl = [all_cross_corr_non_play_pctl;cross_corr_struct(fn).all_cross_corr_non_play_pctl];
    play_cc                      = cat(1,play_cc,cross_corr_struct(fn).all_cross_corr_play);
    non_play_cc                  = cat(1,non_play_cc,cross_corr_struct(fn).all_cross_corr_non_play);
    all_area_comb_cc             = [all_area_comb_cc;cross_corr_struct(fn).these_neurons_areas(cross_corr_struct(fn).cross_corr_comb)];
    session_id_cc                = [session_id_cc;repmat(animal_names(fn,1),size(cross_corr_struct(fn).cross_corr_comb,1),1)];
    clusters_id_cc               = [clusters_id_cc;cross_corr_struct(fn).cluster_info.cluster_id(cross_corr_struct(fn).cross_corr_comb)];

end

%% 3 Area–area combinations
% Unique directed pairs among SupCol, DLPAG, LPAG, VLPAG, DR.

area_combinations = cell2table(all_area_comb_cc);
area_combinations.Properties.VariableNames = {'Neuron1','Neuron2'};
unique_area_combinations = unique(area_combinations, 'rows');

relevant_areas = {'SupCol','DLPAG','LPAG','VLPAG','DR'};
unique_area_combinations = unique_area_combinations(ismember(unique_area_combinations.Neuron1,relevant_areas) & ismember(unique_area_combinations.Neuron2,relevant_areas),:);

%% 4 Pairwise CCG examples (play vs non-play)
% For each area pair: heatmap of shuffle-normalized play CCGs (sorted by
% lag-0 excess), mean ± CI (play red, non-play blue), and before vs after
% lag-0 scatter. Excess = (CCG − shuffle) / max(shuffle).
x_lim = [-1 1];
time_centers = .5*(psth_edges_cc(1:end-1)+psth_edges_cc(2:end));

alpha_level = 0.01/sum(time_centers>=x_lim(1) & time_centers<=x_lim(2));
order_range = [-1/3.5 1/3.5];

smoth_wind =10;
y_lim = [-.5 .75];

max_activation = 2;
square_axis  = [-.2 .2];
n_sp_col = 5;
nsp =1;
plot_bool = [true true];
plot_sign = false;
figure('units','normalized','outerposition',[0 0 1 1]);
colormap(jet(256))
for j=1:size(unique_area_combinations,1)
    if nsp>n_sp_col
        figure('units','normalized','outerposition',[0 0 1 1]);
        colormap(jet(256))

        nsp=1;
    end
    subplot(4,n_sp_col,nsp)
    area_combination_indexes = ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j));

    % % plot(time_centers,smooth(mean(all_cross_corr_norm(indexes,:))))
    % plot(time_centers,smooth(mean(all_cross_corr_norm(indexes,:), 'omitmissing'), 5), 'k')

    if plot_bool(1)
        % yyaxis right
        hold on
        this_cc = squeeze(play_cc(area_combination_indexes,:));
        for k = 1:size(this_cc,1)
            this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
        end

        this_cc_control = all_cross_corr_play_pctl(area_combination_indexes,:);
        for k = 1:size(this_cc_control,1)
            this_cc_control(k,:) = movmean( this_cc_control(k,:),smoth_wind);
        end
        direct_cc= this_cc;
        this_cc_diff = (this_cc-this_cc_control);
        % vector_norn = vecnorm(this_cc_control, 2, 2);
        vector_norn = max(this_cc_control, [], 2);
        vector_norn(vector_norn==0) = Inf;

        this_cc = diag(1./vector_norn)*this_cc_diff;
        this_cc(max(abs(this_cc),[],2)>max_activation,:) = NaN;
        this_play_cc = this_cc;
        delta_act = median(this_cc(:,time_centers>order_range(1) & time_centers<order_range(2)),2, 'omitmissing');
        before_activation_play = median(this_cc(:,time_centers>order_range(1) & time_centers<0),2, 'omitmissing');
        after_activation_play  = median(this_cc(:,time_centers>0 & time_centers<order_range(2)),2, 'omitmissing');
        [~, order] = sort(delta_act);
        imagesc(time_centers, 1:size(this_cc,1), this_cc(order,:))
        clim(2*y_lim)
        axis tight
        xlim(x_lim)
        title(num2str([sum(delta_act>0) sum(delta_act<0)]/numel(delta_act)))

        subplot(4,n_sp_col,nsp+n_sp_col)

        %
        % this_cc = ((this_cc-this_cc_control)-repmat(mean(this_cc_control,2, 'omitmissing'),1,size(this_cc_control,2)))./...
        %    repmat(std(this_cc_control,[],2, 'omitmissing'),1,size(this_cc_control,2));
        % this_cc(isinf(this_cc)) = NaN;
        if plot_sign
            this_cc_pos = this_cc(delta_act>0,:);
            [h, ~, ci] = ttest(this_cc_pos,0, 'Alpha', alpha_level);
            no_nan = ~any(isnan(ci));
            if sum(no_nan)>2
                fill([time_centers(no_nan) fliplr(time_centers(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'r', 'FaceAlpha',.25, 'EdgeColor','none')
            end
            hold on
            plot(time_centers,mean(this_cc_pos, 'omitmissing'), 'r')

            this_cc_neg = this_cc(delta_act<0,:);
            [~, ~, ci] = ttest(this_cc_neg,0, 'Alpha', alpha_level);
            no_nan = ~any(isnan(ci));
            if sum(no_nan)>2
                fill([time_centers(no_nan) fliplr(time_centers(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'b', 'FaceAlpha',.25, 'EdgeColor','none')
            end
            plot(time_centers,mean(this_cc_neg, 'omitmissing'), 'b')

        else
            hold on
            [h, ~, ci] = ttest(this_cc,0, 'Alpha', alpha_level);
            no_nan = ~any(isnan(ci));
            fill([time_centers(no_nan) fliplr(time_centers(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'r', 'FaceAlpha',.25, 'EdgeColor','none')
            plot(time_centers,mean(this_cc, 'omitmissing'), 'r')
        end

        % Find start and end indices of contiguous h==1 segments
        starts = find(diff([0 h]) == 1);
        ends   = find(diff([h 0]) == -1)-1;
        ends(ends==0)=1;
        % y-coordinates for rectangle height
        y_bottom = y_lim(2)*0.9;
        y_top    = y_lim(2);

        for k = 1:length(starts)
            xs = [time_centers(starts(k))  time_centers(ends(k)) ...
                time_centers(ends(k))   time_centers(starts(k))];
            ys = [y_bottom y_bottom y_top y_top];

            fill(xs, ys, 'r', 'FaceAlpha', .25, 'EdgeColor', 'none');
        end

    end

    if plot_bool(2)

        this_cc = squeeze(non_play_cc(area_combination_indexes,:));

        for k = 1:size(this_cc,1)
            this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
        end
        this_cc_control = all_cross_corr_non_play_pctl(area_combination_indexes,:);
        for k = 1:size(this_cc_control,1)
            this_cc_control(k,:) = movmean( this_cc_control(k,:),smoth_wind);
        end

        %   this_cc = ((this_cc-this_cc_control)-repmat(mean(this_cc_control,2, 'omitmissing'),1,size(this_cc_control,2)))./...
        %    repmat(std(this_cc_control,[],2, 'omitmissing'),1,size(this_cc_control,2));
        % this_cc(isinf(this_cc)) = NaN;
        this_cc_diff = (this_cc-this_cc_control);
        % vector_norn = vecnorm(this_cc_control, 2, 2);
        vector_norn = max(this_cc_control, [], 2);
        vector_norn(vector_norn==0) = Inf;

        this_cc = diag(1./vector_norn)*this_cc_diff;
        this_cc(abs(this_cc)>.2) = NaN;
        before_activation_non_play = median(this_cc(:,time_centers>order_range(1) & time_centers<0),2, 'omitmissing');
        after_activation_non_play  = median(this_cc(:,time_centers>0 & time_centers<order_range(2)),2, 'omitmissing');
        this_non_play_cc = this_cc;

        [h, ~, ci] = ttest(this_cc,0, 'Alpha', alpha_level);
        no_nan = ~any(isnan(ci));
        fill([time_centers(no_nan) fliplr(time_centers(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'b', 'FaceAlpha',.25, 'EdgeColor','none')
        plot(time_centers,mean(this_cc, 'omitmissing'), 'b')

        % Find start and end indices of contiguous h==1 segments
        starts = find(diff([0 h]) == 1);
        ends   = find(diff([h 0]) == -1)-1;
        ends(ends==0) = 1;
        % y-coordinates for rectangle height
        y_bottom = y_lim(2)*0.8;
        y_top    = y_lim(2)*.9;

        for k = 1:length(starts)
            xs = [time_centers(starts(k))  time_centers(ends(k)) ...
                time_centers(ends(k))   time_centers(starts(k))];
            ys = [y_bottom y_bottom y_top y_top];

            fill(xs, ys, 'b', 'FaceAlpha', .25, 'EdgeColor', 'none');
        end

        [h, ~, ci] = ttest(this_non_play_cc,this_play_cc, 'Alpha', alpha_level);

        % Find start and end indices of contiguous h==1 segments
        starts = find(diff([0 h]) == 1);
        ends   = find(diff([h 0]) == -1)-1;
        ends(ends==0) = 1;
        % y-coordinates for rectangle height
        y_bottom = y_lim(1)*0.9;
        y_top    = y_lim(1)*.8;

        for k = 1:length(starts)
            xs = [time_centers(starts(k))  time_centers(ends(k)) ...
                time_centers(ends(k))   time_centers(starts(k))];
            ys = [y_bottom y_bottom y_top y_top];

            fill(xs, ys, 'b', 'FaceAlpha', .25, 'EdgeColor', 'none');
        end
    end

    hold on
    plot([0 0],y_lim, 'g')
    plot([1 1]/3.5,y_lim, 'k')
    plot(-[1 1]/3.5,y_lim, 'k')
    plot(x_lim, [0 0], ':k', 'LineWidth',2)
    ylim(y_lim)
    xlim(x_lim)

    title(unique_area_combinations{j,1})
    ylabel([unique_area_combinations{j,2}{1}, ' excess activation'])

    subplot(4,n_sp_col,nsp+ 2*n_sp_col)
    hold on
    %      histogram(delta_act, linspace(-0.02, 0.02, 200), 'EdgeColor','none')
    %           histogram(delta_act, 200, 'FaceColor','k','EdgeColor','none')
    % p = signrank(delta_act);
    % y_lim_hist = ylim;
    % hold on
    % plot([1 1]*median(delta_act, 'omitmissing'),y_lim_hist, 'r')
    % title(p)
    plot(before_activation_play,after_activation_play, 'k.')
    axis tight
    these_axis = axis;
    % plot([min(these_axis) max(these_axis)],[min(these_axis) max(these_axis)], 'r')
    plot(square_axis, square_axis, 'r')
    p = signrank(before_activation_play,after_activation_play);
    title(['Play p=', num2str(p)])
    xlabel('Excess before')
    ylabel('Excess after')
    subplot(4,n_sp_col,nsp+3*n_sp_col)
    hold on

    plot(before_activation_non_play,after_activation_non_play, 'k.')
    axis tight
    these_axis = axis;
    % plot([min(these_axis) max(these_axis)],[min(these_axis) max(these_axis)], 'r')
    plot(square_axis, square_axis, 'r')
    p = signrank(before_activation_play,after_activation_play);
    title(['Non Play p=', num2str(p)])
    xlabel('Excess before')
    ylabel('Excess after')

    nsp=nsp+1;
end

%% 5 Mean CCG by area combination
% Probability-normalized mean CCG ± CI for play (red) and non-play (blue).
time_centers = .5*(psth_edges_cc(1:end-1)+psth_edges_cc(2:end));
figure
for j=1:size(unique_area_combinations,1)
    subplot(5,5,j)
    area_combination_indexes = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j)));
    % plot(time_centers,smooth(mean(all_cross_corr_norm(indexes,:))))

    this_cc = squeeze(play_cc(area_combination_indexes,:));
    for k = 1:size(this_cc,1)
        this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
    end
    this_cc = diag(1./sum(this_cc,2))*this_cc;
    [~, ~, ci] = ttest(this_cc);
    no_nan = any(isnan(ci));
    fill([time_centers(no_nan) fliplr(time_centers(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'k', 'FaceAlpha',.25, 'EdgeColor','none')
    plot(time_centers,mean(this_cc, 'omitmissing'), 'k')
    hold on
    y_lim = ylim;
    xlim([-.1 .1])

    % yyaxis right
    hold on
    this_cc = squeeze(play_cc(area_combination_indexes,:));
    for k = 1:size(this_cc,1)
        this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
    end
    this_cc = diag(1./sum(this_cc,2))*this_cc;
    [~, ~, ci] = ttest(this_cc);
    no_nan = ~any(isnan(ci));
    fill([time_centers(no_nan) fliplr(time_centers(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'r', 'FaceAlpha',.25, 'EdgeColor','none')
    plot(time_centers,mean(this_cc, 'omitmissing'), 'r')
    this_cc = squeeze(non_play_cc(area_combination_indexes ,:));

    for k = 1:size(this_cc,1)
        this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
    end
    this_cc = diag(1./sum(this_cc,2))*this_cc;
    [~, ~, ci] = ttest(this_cc);
    no_nan = ~any(isnan(ci));
    fill([time_centers(no_nan) fliplr(time_centers(no_nan))],[ci(1,no_nan) fliplr(ci(2,no_nan))],'b', 'FaceAlpha',.25, 'EdgeColor','none')
    plot(time_centers,mean(this_cc, 'omitmissing'), 'b')
    xlim(x_lim)
    hold on
    plot([0 0],y_lim, 'g')
    plot([1 1]/3.5,y_lim, 'k')
    plot(-[1 1]/3.5,y_lim, 'k')
    ylim(y_lim)

    title(unique_area_combinations{j,1})
    ylabel([unique_area_combinations{j,2}, ' activation'])

end

%% 6 Load peak / trough / unlocked labels
% Map each CCG pair onto all_neurons_TD, then label neurons by delta-phase
% locking (peak, trough, unlocked). Edit entrained_group* below to choose
% which cell-type combination section 7 plots.
saving_folder = [data_root, '\Analysis results\phase locking data'];

load([saving_folder,'\theta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
all_neurons_TD = all_neurons;
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner1'}))         = {'ThetaPartner1'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner2'}))         = {'ThetaPartner2'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Play'}))             = {'ThetaPlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'PrePlay'}))          = {'ThetaPrePlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'EntireSession'}))    = {'ThetaEntireSession'};

load([saving_folder,'\delta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'}))  =     {'isRt'  };
all_neurons_TD.DeltaPartner1                            = all_neurons.Partner1;
all_neurons_TD.DeltaPartner2                            = all_neurons.Partner2;
all_neurons_TD.DeltaEntireSession                       = all_neurons.EntireSession;
all_neurons_TD.DeltaPlay                                = all_neurons.Play;
all_neurons_TD.DeltaPrePlay                             = all_neurons.PrePlay;
all_neurons_TD.Exited                                   = nan(size(all_neurons_TD,1),1);
all_neurons_TD.Inhibited                                = nan(size(all_neurons_TD,1),1);

nRows = size(clusters_id_cc,1);
idx_pairs_cc = zeros(nRows,2);

for i = 1:nRows
    s = session_id_cc{i};
    id_pair = clusters_id_cc(i,:);
    for j = 1:2
        idx_pairs_cc(i,j) = find( ...
            strcmp(all_neurons_TD.session, s) & ...
            all_neurons_TD.cluster_id == id_pair(j), ...
            1);  % only one match expected
    end
end

%% 7 CCG by cell type and area pair
% Shuffle-normalized play CCG for the selected cell-type combination
% (red) vs the comparison combination (black), one subplot per area pair.

trough          = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & (all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
peak            = all_neurons_TD.DeltaEntireSession.PPCPval<0.01 & ~(all_neurons_TD.DeltaEntireSession.PreferedAngle>pi/2 | all_neurons_TD.DeltaEntireSession.PreferedAngle<-pi/2);
non_entrained   = all_neurons_TD.DeltaEntireSession.PPCPval>.1;
alpha_level = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% here you select with type of combination to plot %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entrained_group1 = trough;
entrained_group2 = trough;
comparison_group1 = non_entrained;
comparison_group2 = non_entrained;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_centers = .5*(psth_edges_cc(1:end-1)+psth_edges_cc(2:end));
z_scored_limit = Inf;
figure
smoth_wind = 1;
smoth_wind2 = 10;
y_lim = [-1 4];
parametric_cond = true;
x_lim = [-2 2];
plot_ci = true;
plot_logit = ~plot_ci;
for j=1:size(unique_area_combinations,1)

    subplot(5,5,j)
    hold on

    area_combination_indexes       = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
        (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)));
    cell_type_comb_trough          = (entrained_group1(idx_pairs_cc(:,1)) &  entrained_group2(idx_pairs_cc(:,2)) &  area_combination_indexes );
    cell_type_comb_non_entrained   = ( comparison_group1(idx_pairs_cc(:,1)) &  comparison_group2(idx_pairs_cc(:,2)) &  area_combination_indexes);

    this_cc = squeeze(play_cc(cell_type_comb_non_entrained,:));
    for k = 1:size(this_cc,1)
        this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
    end

    this_cc_control = all_cross_corr_play_pctl(cell_type_comb_non_entrained,:);
    for k = 1:size(this_cc_control,1)
        this_cc_control(k,:) = movmean( this_cc_control(k,:),smoth_wind);
    end

    this_cc_zs          = this_cc;
    this_cc_control_zs  = this_cc_control;
    for k = 1:size(this_cc,1)
        this_cc_zs(k,:) =(this_cc(k,:)-mean(this_cc_control(k,:)))/std(this_cc_control(k,:));
        this_cc_control_zs(k,:) =(this_cc_control_zs(k,:)-mean(this_cc_control(k,:)))/std(this_cc_control(k,:));
    end

    this_cc = this_cc_zs-this_cc_control_zs;
    this_cc(any(abs(this_cc)>z_scored_limit,2),:) = NaN;

    if parametric_cond
        p_non_entrained = nan(1,size(this_cc,2));

        for t_index=1:numel(p_non_entrained)
            if sum(~isnan(this_cc(:,t_index)))>5
                p_non_entrained(t_index) = signrank(this_cc(:,t_index));
            end
        end

    else
        [~, p_non_entrained, ci] = ttest(this_cc,0, 'Alpha', alpha_level);
    end

    if plot_ci
        plot(time_centers,smooth(median(this_cc, 'omitmissing'), smoth_wind2), 'k')
        ylim(y_lim)
    else
        logit_val = -log10(p_non_entrained).*sign(movmean(mean(this_cc, 'omitmissing'),20));
        plot(time_centers,logit_val, 'k')
        plot(x_lim, [2 2], ':k')
        plot(x_lim, -[2 2], ':k')
    end

    hold on
    this_cc = squeeze(play_cc(cell_type_comb_trough,:));
    for k = 1:size(this_cc,1)
        this_cc(k,:) = movmean( this_cc(k,:),smoth_wind);
    end
    this_cc_control = all_cross_corr_play_pctl(cell_type_comb_trough,:);
    for k = 1:size(this_cc_control,1)
        this_cc_control(k,:) = movmean( this_cc_control(k,:),smoth_wind);
    end

    this_cc_zs = this_cc;
    this_cc_control_zs = this_cc_control;
    for k = 1:size(this_cc,1)
        this_cc_zs(k,:) =(this_cc(k,:)-mean(this_cc_control(k,:)))/std(this_cc_control(k,:));
        this_cc_control_zs(k,:) =(this_cc_control_zs(k,:)-mean(this_cc_control(k,:)))/std(this_cc_control(k,:));
    end

    this_cc = this_cc_zs-this_cc_control_zs;
    this_cc(any(abs(this_cc)>z_scored_limit,2),:) = NaN;

    if parametric_cond
        p_entrained = nan(1,size(this_cc,2));

        for t_index=1:numel(p_entrained)
            if sum(~isnan(this_cc(:,t_index)))>5
                p_entrained(t_index) = signrank(this_cc(:,t_index));
            end
        end
        p_entrained = min(p_entrained,1);
    else
        [~, p_entrained, ci] = ttest(this_cc,0, 'Alpha', alpha_level);
    end

    if plot_ci
        plot(time_centers,smooth(median(this_cc, 'omitmissing'),smoth_wind2), 'r')
        ylim(y_lim)
    else

        logit_val = -log10(p_entrained).*sign(movmean(mean(this_cc, 'omitmissing'),20));
        plot(time_centers,logit_val, 'b')
        plot(x_lim, [2 2], ':k')
        plot(x_lim, -[2 2], ':k')
    end

    title([unique_area_combinations{j,1}{1}, ' #',num2str(numel(cell_type_comb_non_entrained)), '/',num2str(numel(cell_type_comb_trough)) ] )
    ylabel(unique_area_combinations{j,2}{1})
    xlim(x_lim)
    pause(.1)

end

%% 8 Prepare CCG for one area–area pair
% z-score play or non-play CCGs to the shuffle for trough, peak, and
% unlocked pairs. j indexes unique_area_combinations (e.g. 14 = LPAG–SupCol).

cc2analysie = play_cc;
cc2analyse_control = all_cross_corr_play_pctl;
% cc2analysie = non_play_cc;
% cc2analyse_control = all_cross_corr_non_play_pctl;

j =20;
smoth_wind = 1;
area_combination_indexes = ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j));
cell_type_comb_trough          = find(trough(idx_pairs_cc(:,1)) &  trough(idx_pairs_cc(:,2)) & area_combination_indexes);
cell_type_comb_non_entrained   = find(non_entrained(idx_pairs_cc(:,1)) & non_entrained(idx_pairs_cc(:,2))  & area_combination_indexes) ;
cell_type_comb_peak            = find(peak(idx_pairs_cc(:,1)) &  peak(idx_pairs_cc(:,2))  & area_combination_indexes);

%Trough cells
trough_this_cc_control     = cc2analyse_control( cell_type_comb_trough,:);
trough_this_cc             = squeeze(cc2analysie(cell_type_comb_trough,:));
for k = 1:size(trough_this_cc,1)
    trough_this_cc(k,:) = movmean( trough_this_cc(k,:),smoth_wind);
end
for k = 1:size(trough_this_cc_control,1)
    trough_this_cc_control(k,:) = movmean( trough_this_cc_control(k,:),smoth_wind);
end

trough_this_cc_zs          =  trough_this_cc;
trough_this_cc_control_zs  =  trough_this_cc_control;
for k = 1:size( trough_this_cc,1)
    trough_this_cc_zs(k,:) = ( trough_this_cc(k,:) - mean( trough_this_cc_control(k,:)))/std( trough_this_cc_control(k,:));
    trough_this_cc_control_zs(k,:) = ( trough_this_cc_control(k,:) - mean( trough_this_cc_control(k,:)))/std( trough_this_cc_control(k,:));
end

%Peak cells
peak_this_cc_control     = cc2analyse_control(cell_type_comb_peak,:);
peak_this_cc             = squeeze(cc2analysie(cell_type_comb_peak,:));
for k = 1:size(peak_this_cc,1)
    peak_this_cc(k,:) = movmean( peak_this_cc(k,:),smoth_wind);
end
for k = 1:size(peak_this_cc_control,1)
    peak_this_cc_control(k,:) = movmean( peak_this_cc_control(k,:),smoth_wind);
end

peak_this_cc_zs          =  peak_this_cc;
peak_this_cc_control_zs  =  peak_this_cc_control;
for k = 1:size( peak_this_cc,1)
    peak_this_cc_zs(k,:) = ( peak_this_cc(k,:) - mean( peak_this_cc_control(k,:)))/std( peak_this_cc_control(k,:));
    peak_this_cc_control_zs(k,:) = ( peak_this_cc_control(k,:) - mean( peak_this_cc_control(k,:)))/std( peak_this_cc_control(k,:));
end

%Non entrained cells
ne_this_cc_control         = cc2analyse_control(cell_type_comb_non_entrained,:);
ne_this_cc                 = squeeze(cc2analysie(cell_type_comb_non_entrained,:));
for k = 1:size(ne_this_cc,1)
    ne_this_cc(k,:) = movmean( ne_this_cc(k,:),smoth_wind);
end
for k = 1:size(ne_this_cc_control,1)
    ne_this_cc_control(k,:) = movmean( ne_this_cc_control(k,:),smoth_wind);
end

ne_this_cc_zs              =  ne_this_cc;
ne_this_cc_control_zs      =  ne_this_cc_control;
for k = 1:size( ne_this_cc,1)
    ne_this_cc_zs(k,:) = ( ne_this_cc(k,:) - mean( ne_this_cc_control(k,:)))/std( ne_this_cc_control(k,:));
    ne_this_cc_control_zs(k,:) = ( ne_this_cc_control(k,:) - mean( ne_this_cc_control(k,:)))/std( ne_this_cc_control(k,:));
end

%% 9 Single-pair CCG examples (trough, peak, unlocked)
% Ten example pairs per cell type: real CCG vs shuffle (black).

j =14;
raw_n = 50;
% raw_n = 10;
x_lim = [-.5 .5];
entrained_group1 = trough;
entrained_group2 = trough;
bin_size_cc = mean(diff(time_centers));
area_combination_indexes = ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j));
cell_type_comb_trough          = find(trough(idx_pairs_cc(:,1)) &  trough(idx_pairs_cc(:,2)) & area_combination_indexes);
figure('units','normalized','outerposition',[0 0 .3 1]);
for kk=1:10
    subplot(5,2,kk)
    index1 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(raw_n+kk-1),1) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(raw_n+kk-1)}));
    index2 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(raw_n+kk-1),2) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(raw_n+kk-1)}));
    all_neurons_TD([index1,index2],:)
    if ~isempty(index1) &  ~isempty(index2)
        plot(time_centers, movmean( trough_this_cc_zs(raw_n+kk-1,:),10), 'b')
        hold on
        plot(time_centers, movmean( trough_this_cc_control_zs(raw_n+kk-1,:),10), 'k')
        xlim(x_lim)
        title([all_neurons_TD.session{index1}, ' ID', num2str(all_neurons_TD.cluster_id(index1))])
        ylabel(num2str(all_neurons_TD.cluster_id(index2)))
    end
end
sgtitle('Torugh cells')
raw_n = 10;

figure('units','normalized','outerposition',[.33 0 .3 1]);
for kk=1:10
    subplot(5,2,kk)
    index1 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_peak(raw_n+kk-1),1) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_peak(raw_n+kk-1)}));
    index2 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_peak(raw_n+kk-1),2) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_peak(raw_n+kk-1)}));
    if ~isempty(index1) &  ~isempty(index2)
        plot(time_centers, movmean( peak_this_cc_zs(raw_n+kk-1,:),10), 'r')
        hold on
        plot(time_centers, movmean( peak_this_cc_control_zs(raw_n+kk-1,:),10), 'k')
        xlim(x_lim)
        title([all_neurons_TD.session{index1}, ' ID', num2str(all_neurons_TD.cluster_id(index1))])
        ylabel(num2str(all_neurons_TD.cluster_id(index2)))
    end
end
sgtitle('Peak cells')

raw_n = 10;

figure('units','normalized','outerposition',[.66 0 .3 1]);
for kk=1:10
    subplot(5,2,kk)
    index1 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_non_entrained(raw_n+kk-1),1) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_non_entrained(raw_n+kk-1)}));
    index2 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_non_entrained(raw_n+kk-1),2) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_non_entrained(raw_n+kk-1)}));
    if ~isempty(index1) &  ~isempty(index2)
        plot(time_centers, movmean( ne_this_cc_zs(raw_n+kk-1,:),10), 'g')
        hold on
        plot(time_centers, movmean( ne_this_cc_control_zs(raw_n+kk-1,:),10), 'k')
        xlim(x_lim)
        title([all_neurons_TD.session{index1}, ' ID', num2str(all_neurons_TD.cluster_id(index1))])
        ylabel(num2str(all_neurons_TD.cluster_id(index2)))
    end
end
sgtitle('non-entrained cells')

%% 10 Trough-cell CCG vs shuffle (selected pair)
% Heatmaps of trough–trough CCGs, shuffle, and difference; median traces
% and signed-rank p over lag. j = 14 is LPAG vs SupCol.

alpha4rank_test = 0.00000001;
p_width = 0.01;

x_lim = [-.5 .5];

y_lim = [0 1];
y_lim2 = [0 .2];

figure

subplot(5,1,1)
colormap(1-gray)

[~, order] = sort(max(trough_this_cc,[],2));
imagesc(time_centers, 1:size(trough_this_cc,1),trough_this_cc(order,:))
clim([-2 2])
xlim(x_lim)
title([unique_area_combinations{j,1}{1}, ' #',num2str(numel(cell_type_comb_non_entrained)), '/',num2str(numel(cell_type_comb_trough)) ] )
ylabel(unique_area_combinations{j,2}{1})

subplot(5,1,2)
colormap(1-gray)
imagesc(time_centers, 1:size(trough_this_cc,1),trough_this_cc_control(order,:))
clim([-2 2])
xlim(x_lim)

subplot(5,1,3)
colormap(1-gray)
imagesc(time_centers, 1:size(trough_this_cc,1),trough_this_cc(order,:)-trough_this_cc_control(order,:))
clim([-2 2]/2)
xlim(x_lim)

subplot(5,1,4)

plot(time_centers, median(trough_this_cc, 'omitmissing'), 'b')
hold on
plot(time_centers, median(trough_this_cc_control, 'omitmissing'), 'k')
ylim(y_lim)
xlim(x_lim)

subplot(5,1,5)
rank_p_value = nan(1,size(trough_this_cc,2));
z_over_time = nan(1,size(trough_this_cc,2));

for time_idx = 1:numel(rank_p_value)
    [ rank_p_value(time_idx), ~, stats] = signrank(trough_this_cc(:,time_idx)-trough_this_cc_control(:,time_idx));
    n = sum((trough_this_cc(:,time_idx) - trough_this_cc_control(:,time_idx)) ~= 0);
    z_over_time(time_idx) = stats.zval/sqrt(n);
end
plot(time_centers, median(trough_this_cc-trough_this_cc_control, 'omitmissing'), 'b')
hold on

start_end = [find(diff([0 rank_p_value<alpha4rank_test 0])==1)' (find(diff([0 rank_p_value<alpha4rank_test 0])==-1)-1)'];

for per_n = 1:size(start_end,1)

    fill([time_centers(start_end(per_n,:)) fliplr(time_centers(start_end(per_n,:)))], y_lim2(2)*[1 1 1 1] + [0 0 -p_width -p_width], 'r')
end
xlim(x_lim)
ylim(y_lim2)

figure
plot(time_centers,z_over_time)
title([unique_area_combinations{j,1}{1}, ' to ',unique_area_combinations{j,2}{1} ] )
legend({'Play', 'non-Play'})

%% 11 Unlocked-cell CCG vs shuffle (selected pair)
% Same layout as section 10, for unlocked–unlocked pairs.

alpha4rank_test = 0.01;
p_width = 0.02;
j =14;

x_lim = [-.5 .5];

y_lim = [-.01 .01];
y_lim2 = [-.1 .1];

figure

subplot(5,1,1)
colormap(1-gray)

[~, order] = sort(max(ne_this_cc,[],2));
imagesc(time_centers, 1:size(ne_this_cc,1),ne_this_cc(order,:))
clim([-2 2])
xlim(x_lim)
title([unique_area_combinations{j,1}{1}, ' #',num2str(numel(cell_type_comb_non_entrained)), '/',num2str(numel(cell_type_comb_trough)) ] )
ylabel(unique_area_combinations{j,2}{1})

subplot(5,1,2)
colormap(1-gray)
imagesc(time_centers, 1:size(ne_this_cc,1),ne_this_cc_control(order,:))
clim([-2 2])
xlim(x_lim)

subplot(5,1,3)
colormap(1-gray)
imagesc(time_centers, 1:size(ne_this_cc,1),ne_this_cc(order,:)-ne_this_cc_control(order,:))
clim([-2 2]/2)
xlim(x_lim)

subplot(5,1,4)

plot(time_centers, mean(ne_this_cc, 'omitmissing'), 'b')
hold on
plot(time_centers, mean(ne_this_cc_control, 'omitmissing'), 'k')

xlim(x_lim)

subplot(5,1,5)
rank_p_value = nan(1,size(ne_this_cc,2));

for time_idx = 1:numel(rank_p_value)
    rank_p_value(time_idx) = signrank(ne_this_cc(:,time_idx)-ne_this_cc_control(:,time_idx));
end
plot(time_centers, mean(ne_this_cc-ne_this_cc_control, 'omitmissing'), 'b')
hold on

start_end = [find(diff([0 rank_p_value<alpha4rank_test 0])==1)' (find(diff([0 rank_p_value<alpha4rank_test 0])==-1)-1)'];

for per_n = 1:size(start_end,1)

    fill([time_centers(start_end(per_n,:)) fliplr(time_centers(start_end(per_n,:)))], y_lim2(2)*[1 1 1 1] + [0 0 -p_width -p_width], 'r')
end
xlim(x_lim)
ylim(y_lim2)

%% 12 All area pairs: signed-rank z over lag
% For trough, peak, and unlocked: z of (CCG − shuffle) vs lag, play (magenta)
% vs non-play (black). Stores combinations{cell_type, play/nonplay, pair}.

n=5;
[I,J] = find(triu(true(n)));   % upper triangle including diagonal
idx = sub2ind([n n], I, J);

combinations = cell(3,2,size(unique_area_combinations,1));

CELLTYPE_NAMES = {'TROUGH','PEAK','NON ENTRAINED'};
for cell_type = 1:3
    figure
    smoth_wind = 1;
    for cc_type = 1:2

        if cc_type==1
            cc2analysie = play_cc;
            cc2analyse_control = all_cross_corr_play_pctl;
        else
            cc2analysie = non_play_cc;
            cc2analyse_control = all_cross_corr_non_play_pctl;
        end

        for j =1:size(unique_area_combinations,1)

            area_combination_indexes = (ismember(area_combinations.Neuron1,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron2,unique_area_combinations.Neuron2(j))) | ...
                (ismember(area_combinations.Neuron2,unique_area_combinations.Neuron1(j)) & ismember(area_combinations.Neuron1,unique_area_combinations.Neuron2(j))) ;
            cell_type_comb_trough          = find(trough(idx_pairs_cc(:,1)) &  trough(idx_pairs_cc(:,2)) & area_combination_indexes);
            cell_type_comb_non_entrained   = find(non_entrained(idx_pairs_cc(:,1)) & non_entrained(idx_pairs_cc(:,2))  & area_combination_indexes) ;
            cell_type_comb_peak            = find(peak(idx_pairs_cc(:,1)) &  peak(idx_pairs_cc(:,2))  & area_combination_indexes);
            if cell_type==1
                %Trough cells
                trough_this_cc_control     = cc2analyse_control( cell_type_comb_trough,:);
                trough_this_cc             = squeeze(cc2analysie(cell_type_comb_trough,:));
                for k = 1:size(trough_this_cc,1)
                    trough_this_cc(k,:) = movmean( trough_this_cc(k,:),smoth_wind);
                end
                for k = 1:size(trough_this_cc_control,1)
                    trough_this_cc_control(k,:) = movmean( trough_this_cc_control(k,:),smoth_wind);
                end

                this_cc4ANAL            = trough_this_cc;
                this_cc4ANAL_control    = trough_this_cc_control;

            elseif cell_type==2
                %Peak cells
                peak_this_cc_control     = cc2analyse_control(cell_type_comb_peak,:);
                peak_this_cc             = squeeze(cc2analysie(cell_type_comb_peak,:));
                for k = 1:size(peak_this_cc,1)
                    peak_this_cc(k,:) = movmean( peak_this_cc(k,:),smoth_wind);
                end
                for k = 1:size(peak_this_cc_control,1)
                    peak_this_cc_control(k,:) = movmean( peak_this_cc_control(k,:),smoth_wind);
                end

                this_cc4ANAL            = peak_this_cc;
                this_cc4ANAL_control    = peak_this_cc_control;
            elseif cell_type==3
                %Non entrained cells
                ne_this_cc_control         = cc2analyse_control(cell_type_comb_non_entrained,:);
                ne_this_cc                 = squeeze(cc2analysie(cell_type_comb_non_entrained,:));
                for k = 1:size(ne_this_cc,1)
                    ne_this_cc(k,:) = movmean( ne_this_cc(k,:),smoth_wind);
                end
                for k = 1:size(ne_this_cc_control,1)
                    ne_this_cc_control(k,:) = movmean( ne_this_cc_control(k,:),smoth_wind);
                end

                this_cc4ANAL            = ne_this_cc;
                this_cc4ANAL_control    = ne_this_cc_control;

            end
            if min(size(this_cc4ANAL))>1

                rank_p_value = nan(1,size(this_cc4ANAL,2));
                z_over_time = nan(1,size(this_cc4ANAL_control,2));

                for time_idx = 1:numel(rank_p_value)
                    [ rank_p_value(time_idx), ~, stats] = signrank(this_cc4ANAL(:,time_idx)-this_cc4ANAL_control(:,time_idx));
                    n = sum((this_cc4ANAL(:,time_idx) - this_cc4ANAL_control(:,time_idx)) ~= 0);

                    if ismember('zval', fields(stats))
                        z_over_time(time_idx) = stats.zval/sqrt(n);
                    end
                end

                z_over_time = movmean(z_over_time, 20);

                subplot(5,5,j)
                hold on
                if cc_type==1
                    current_color = 'm';
                else
                    current_color = 'k';
                end
                plot(time_centers,z_over_time,current_color)
                combinations{cell_type, cc_type, j} = z_over_time;
            end
            title([unique_area_combinations{j,1}{1}, ' to ',unique_area_combinations{j,2}{1} ] )

        end
    end
    sgtitle(CELLTYPE_NAMES{cell_type})
    pause(.1)
end

%% 13 Play minus non-play z (unique area pairs)
% Difference of the section-12 z traces (play − non-play), upper triangle
% of the 5×5 area matrix. Rows = trough / peak / unlocked.

figure
mean_change = nan(3, numel(idx));
x_labels = cell(numel(idx), 1);
for ctype = 1:3
    for j = 1:numel(idx)
        subplot(3, numel(idx), j + (ctype-1)*numel(idx))
        plot(time_centers, movmean(combinations{ctype,1,idx(j)},20) - movmean(combinations{ctype,2,idx(j)},20))
        mean_change(ctype, j) = mean(movmean(combinations{ctype,1,idx(j)},20) - movmean(combinations{ctype,2,idx(j)},20));
        hold on
        plot([-3 3], [0 0], ':k')
        title([unique_area_combinations{idx(j),1}{1}, ' ', unique_area_combinations{idx(j),2}{1}])
        ylim([-.8 .8])
        x_labels{j} = [unique_area_combinations{idx(j),1}{1}, ' ', unique_area_combinations{idx(j),2}{1}];
    end
end

%% 14 Mean play–non-play Δz per area pair
% Bar of the time-averaged difference from section 13, one row per cell type.

figure
for j = 1:3
    subplot(3, 1, j)
    bar(mean_change(j,:))
    ylim([-.2 .2])
    xticks(1:numel(idx))
    xticklabels(x_labels)
end

%% 15 Browse trough-cell pairs for an example
% Grid of trough–trough CCGs (real vs shuffle) with mean CCG > 1, to pick
% a pair for section 16.
coincident_level = 0.005;
raw_n = 150;
x_lim = [-.25 .25];
bin_size_cc = mean(diff(time_centers));

figure('units','normalized','outerposition',[0 0 .5 1]);
nsp = 1;
list_of_neuron = find(mean(trough_this_cc,2)>1)';
fig_n=1;

for kk_n=1:numel(list_of_neuron)
    kk = list_of_neuron(kk_n);
    if nsp>25
        nsp=1;
        sgtitle(num2str(fig_n))
        fig_n = fig_n+1;
        figure('units','normalized','outerposition',[0 0 .5 1]);

    end
    subplot(5,5,nsp)
    index1 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),1) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
    index2 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),2) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
    % all_neurons_TD([index1,index2],:)
    if ~isempty(index1) &  ~isempty(index2)
        plot(time_centers, trough_this_cc(kk,:), 'b')
        hold on
        plot(time_centers, trough_this_cc_control(kk,:), 'k')
        xlim(x_lim)
        title([all_neurons_TD.session{index1}, ' ID', num2str(all_neurons_TD.cluster_id(index1)), ' cc#', num2str(kk_n)])
        ylabel({[' ID ',num2str(all_neurons_TD.cluster_id(index2))],['Ploting index: ' num2str(kk_n)]})
        ylim tight
        y_lim = ylim;
        fill([-1 1 1 -1]*coincident_level,y_lim([1 1 2 2]), 'g', 'FaceAlpha',.5, 'EdgeColor','none')

    end
    nsp = nsp+1;
end

%% 16 Plot a chosen pair (GENERATE_CROSS_CORR_EXAMPLE)
% Uses list_of_neuron(kk_n) from section 15: session + cluster IDs → example
% CCG at higher time precision.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% wrote ploting index as selected among all subplots in previous section %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% some examples that work %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kk=list_of_neuron(60);
% kk=list_of_neuron(17);
kk=list_of_neuron(70);
time_precision = 0.005;

npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';

bin_size_cc = 0.01;
hist_range_cc = [-3 3];

psth_edges_cc = hist_range_cc(1):bin_size_cc:hist_range_cc(2);

index1 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),1) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
index2 = find(all_neurons_TD.cluster_id==clusters_id_cc(cell_type_comb_trough(kk),2) & ismember(all_neurons_TD.session,session_id_cc{cell_type_comb_trough(kk)}));
animal_id = all_neurons_TD.session{index1};
id_1    = all_neurons_TD.cluster_id(index1);
id_2     = all_neurons_TD.cluster_id(index2);

animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];

animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];

animal2analize = animal_id;
animal_list = animal_list(ismember(animal_file_names,animal2analize));
GENERATE_CROSS_CORR_EXAMPLE([npx_Raw_Data, '\', animal_list.name],bin_size_cc, hist_range_cc,[id_1 id_2],time_precision);

pause(.1)

%% 17 Browse unlocked-cell pairs
% Same grid as section 15, for unlocked–unlocked pairs (mean CCG > 0.1).

raw_n = 150;
x_lim = [-.5 .5];
bin_size_cc = mean(diff(time_centers));

figure('units','normalized','outerposition',[0 0 .3 1]);
nsp = 1;
list_of_neuron = find(mean(ne_this_cc, 2) > .1)';
fig_n = 1;

for kk_n = 1:numel(list_of_neuron)
    kk = list_of_neuron(kk_n);
    if nsp > 25
        nsp = 1;
        sgtitle(num2str(fig_n))
        fig_n = fig_n + 1;
        figure
    end
    subplot(5, 5, nsp)
    index1 = find(all_neurons_TD.cluster_id == clusters_id_cc(cell_type_comb_trough(kk),1) & ismember(all_neurons_TD.session, session_id_cc{cell_type_comb_trough(kk)}));
    index2 = find(all_neurons_TD.cluster_id == clusters_id_cc(cell_type_comb_trough(kk),2) & ismember(all_neurons_TD.session, session_id_cc{cell_type_comb_trough(kk)}));
    if ~isempty(index1) && ~isempty(index2)
        plot(time_centers, ne_this_cc(kk,:), 'b')
        hold on
        plot(time_centers, ne_this_cc_control(kk,:), 'k')
        xlim(x_lim)
        title([all_neurons_TD.session{index1}, ' ID', num2str(all_neurons_TD.cluster_id(index1)), ' cc#', num2str(kk_n)])
        ylabel(num2str(all_neurons_TD.cluster_id(index2)))
    end
    nsp = nsp + 1;
end
