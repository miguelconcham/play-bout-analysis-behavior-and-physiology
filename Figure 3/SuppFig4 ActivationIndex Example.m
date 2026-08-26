%% (Supp fig 4 a-c) EXAMPLE PLOT all neurons per area sorting images by activity in response_time
% Run independently: execute this section only (loads the section-8 file).
checkpoint_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\activation index';
example_checkpoint = [checkpoint_folder, '\activation_index_example_plot.mat'];
disp('Loading example-plot checkpoint')
load(example_checkpoint);
disp('Checkpoint loaded')

color_Codes = {[1 0 0],[0 0 1],[0 .25 0]};
border_width = 5;
n_perm = 10000;
y_lim_cor = [-1 1];
response_time = [0 5];
histogram_edge = -2:0.25:8;
alpha_level = 0.01;
non_entrained_lvl =0.1;
x_lim = [-5 10];
face_alpha = .2;
y_lim = [-2 4];
c_lim = [-2 10];
sw_lim = [-4 4];
color_by_group = {'r','b','k'};
group_labels = {'Peak','Trough','Non Entrained'};
smooth_window = 5;
warped_time     = (((1:60)/20)*5) - 5 ;
non_wraped_time = (1:50)/5 - 5;
baseline  = [-Inf 0];

onset_time      = [0 0];
offset_time     = [5 5];
all_neurons_TD.Modulated =all_neurons_TD.Exited==1 |  all_neurons_TD.Inhibited==1;
alpha_ttest     = 0.1;
response_type   = 'All';  %options: Modulated  Exited  Inhibited All
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Here you select the behvior you want to plot %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psth2use        = 'PsthWarpedCHOther'; %options: PsthOnset  PsthOffset PsthWarped PsthOnlyPB PsthWarpedOther PsthWarpedSelf
                                    % PsthWarpedCHSelf PsthWarpedCHOther
                                    % PsthWarpedPOSelf PsthWarpedPOOther
                                    % PsthWarpedPWIOther PsthWarpedPWISelf
                                    % etc (full list it psth_map(:,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq2use        = 'DeltaEntireSession'; %options:   ThetaEntireSession DeltaEntireSession
area_list       = {'SupCol' 'DLPAG'	'LPAG'	'VLPAG' 'DR' };
if strcmp(psth2use, 'PsthOnlyPB')
    time2use = mig_edges_centers;
elseif strcmp(psth2use, 'PsthOnset') || strcmp(psth2use, 'PsthOffset')
    time2use= non_wraped_time;
else
    time2use = warped_time;
end
colored_ranges = [-5 -4; 9 10];
color_range_index =any(time2use>=colored_ranges(:,1) & time2use<=(colored_ranges(:,2)),1);
baseline_index      = time2use>=baseline(1) & time2use<=baseline(2);
response_time_index = time2use>=response_time(1) & time2use<=response_time(2);




psth_figure=figure('units','normalized','outerposition',[0 0 .5 1]);
corr_figure=figure('units','normalized','outerposition',[.5 .75 .5 .25]);
session_list = unique(all_neurons_TD.session(~isnan(all_neurons_TD.Exited)));
session_index =ismember(all_neurons_TD.session,session_list );

for an=1:numel(area_list)
    


    if strcmp(response_type, 'All')
        response_index = true(size(all_neurons_TD, 1),1);
    else
        response_index  = all_neurons_TD.(response_type)==1;
    end

    entreined       = all_neurons_TD.(freq2use).PPCPval<=alpha_level ;
    not_entrained   = all_neurons_TD.(freq2use).PPCPval>non_entrained_lvl ;

    area_index      = ismember(all_neurons_TD.area, area_list{an}) & session_index;
    angletype_peak  = ~(all_neurons_TD.(freq2use).PreferedAngle>=pi/2 | all_neurons_TD.(freq2use).PreferedAngle<-pi/2);
    no_nan          = ~isnan(all_neurons_TD.Inhibited);


    peak_index          = area_index & no_nan & response_index & entreined & angletype_peak  ;
    trough_index        = area_index & no_nan & response_index & entreined & ~angletype_peak ;
    nonentrained_index  = area_index & no_nan & response_index & not_entrained ;
    all_indexes         = zeros(size(nonentrained_index));
    all_indexes(peak_index)=1;
    all_indexes(trough_index)=2;
    all_indexes(nonentrained_index)=3;


    
    peak_psth      = all_neurons_TD.(psth2use)(peak_index, :);
    % Baseline correction
    for j=1:size(peak_psth,1)
        peak_psth(j,:) = smooth( peak_psth(j,:),smooth_window);
        if  std( peak_psth(j,baseline_index), 'omitmissing')>0.01
            peak_psth(j,:) = ( peak_psth(j,:) - mean( peak_psth(j,baseline_index), 'omitmissing'))/ std( peak_psth(j,baseline_index), 'omitmissing');
        else
            peak_psth(j,:) = ( peak_psth(j,:) - mean( peak_psth(j,:), 'omitmissing'))/ std( peak_psth(j,:), 'omitmissing');


        end

    end

    trough_psth  =  all_neurons_TD.(psth2use)(trough_index, :);
    % Baseline correction
    for j=1:size(trough_psth,1)
        trough_psth(j,:) = smooth( trough_psth(j,:),smooth_window);
        if  std( trough_psth(j,baseline_index), 'omitmissing')>0.01
            trough_psth(j,:) = ( trough_psth(j,:) - mean( trough_psth(j,baseline_index), 'omitmissing'))/ std( trough_psth(j,baseline_index), 'omitmissing');
        else
            trough_psth(j,:) = ( trough_psth(j,:) - mean( trough_psth(j,:), 'omitmissing'))/ std( trough_psth(j,:), 'omitmissing');


        end
    end

    nonentrained_psth  =  all_neurons_TD.(psth2use)(nonentrained_index, :);
    % Baseline correction
    for j=1:size(nonentrained_psth,1)
        nonentrained_psth(j,:) = smooth( nonentrained_psth(j,:),smooth_window);
        if  std( nonentrained_psth(j,baseline_index), 'omitmissing')>0.01
            nonentrained_psth(j,:) = ( nonentrained_psth(j,:) - mean( nonentrained_psth(j,baseline_index), 'omitmissing'))/ std( nonentrained_psth(j,baseline_index), 'omitmissing');
        else
            nonentrained_psth(j,:) = ( nonentrained_psth(j,:) - mean( nonentrained_psth(j,:), 'omitmissing'))/ std( nonentrained_psth(j,:), 'omitmissing');


        end
    end

    all_psth = all_neurons_TD.(psth2use)(area_index & no_nan & response_index, :);
    % Baseline correction
    for j=1:size(all_psth,1)
        all_psth(j,:) = smooth( all_psth(j,:),smooth_window);
        if  std( all_psth(j,baseline_index), 'omitmissing')>0.01
            all_psth(j,:) = ( all_psth(j,:) - mean( all_psth(j,baseline_index), 'omitmissing'))/ std( all_psth(j,baseline_index), 'omitmissing');
        else
            all_psth(j,:) = ( all_psth(j,:) - mean( all_psth(j,:), 'omitmissing'))/ std( all_psth(j,:), 'omitmissing');

        end
    end

    all_psth_indexes = all_indexes(area_index & no_nan & response_index);


    

    null_psth_peak = nan(n_perm, size(peak_psth,2));
    for pn = 1:n_perm

        sub_selection = all_psth(randperm(size(all_psth,1),size(peak_psth,1)),:);

        null_psth_peak(pn,:) = median(sub_selection, 'omitmissing');
    end
    % peak_pctl_activation = 100*mean(null_psth>mean(peak_psth, 'omitmissing'));
    peak_pctl_activation = mean(null_psth_peak<median(peak_psth, 'omitmissing'));
    peak_pctl_activation(peak_pctl_activation==0) = 1/n_perm;
    peak_pctl_activation(peak_pctl_activation==1) = (n_perm-1)/n_perm;
    peak_pctl_activation = log10(peak_pctl_activation./(1-peak_pctl_activation));
    peak_pctl_inhibition = mean(null_psth_peak>mean(peak_psth, 'omitmissing'));

    null_psth_trough = nan(n_perm, size(trough_psth,2));
    for pn = 1:n_perm

        sub_selection = all_psth(randperm(size(all_psth,1),size(trough_psth,1)),:);

        null_psth_trough(pn,:) = median(sub_selection, 'omitmissing');
    end
    % trough_pctl_activation = 100*mean(null_psth>mean(trough_psth, 'omitmissing'));
    trough_pctl_activation = mean(null_psth_trough<median(trough_psth, 'omitmissing'));

    trough_pctl_activation(trough_pctl_activation==0) = 1/n_perm;
    trough_pctl_activation(trough_pctl_activation==1) = (n_perm-1)/n_perm;
    trough_pctl_activation = log10(trough_pctl_activation./(1-trough_pctl_activation));
    trough_pctl_inhibition = mean(null_psth_trough>median(trough_psth, 'omitmissing'));

    null_psth_nonentrained = nan(n_perm, size(nonentrained_psth,2));
    for pn = 1:n_perm

        sub_selection = all_psth(randperm(size(all_psth,1),size(nonentrained_psth,1)),:);

        null_psth_nonentrained(pn,:) = median(sub_selection, 'omitmissing');
    end
    % nonentrained_pctl_activation = 100*median(null_psth>median(nonentrained_psth, 'omitmissing'));
    nonentrained_pctl_activation = mean(null_psth_nonentrained<median(nonentrained_psth, 'omitmissing'));
    nonentrained_pctl_activation(nonentrained_pctl_activation==0) = 1/n_perm;
    nonentrained_pctl_activation(nonentrained_pctl_activation==1) = (n_perm-1)/n_perm;
    nonentrained_pctl_activation = log10(nonentrained_pctl_activation./(1-nonentrained_pctl_activation));
    nonentrained_pctl_inhibition = mean(null_psth_nonentrained>median(nonentrained_psth, 'omitmissing'));


    group1_levl      = all_neurons_TD.(freq2use).PPC(peak_index);
    group2_levl      = all_neurons_TD.(freq2use).PPC(trough_index);
    group3_levl      = all_neurons_TD.(freq2use).PPC(nonentrained_index);
    all_level        = mean(all_psth(:, response_time_index),2);



    [~, order] = sort(all_level);
    figure(psth_figure)
    %peak neurons imagesc
    subplot(6, numel(area_list),an + [0 1]*numel(area_list) )   

    % psth_rgb = repmat(all_psth*0 ,1,1,3);
    % for j=1:3
    % psth_rgb(all_psth_indexes==j,:,color_Codes{j}==1)=repmat(all_psth(all_psth_indexes==j,:),1,1,sum(color_Codes{j}));
    % psth_rgb(all_psth_indexes==j,color_range_index,color_Codes{j}==1)=c_lim(2);    
    % end
    % 
    % psth_rgb(psth_rgb>c_lim(2))=c_lim(2);
    % psth_rgb(psth_rgb<c_lim(1))=c_lim(1);
    % psth_rgb = (psth_rgb-c_lim(1) )/c_lim(2);
    % all_matrix2plot = psth_rgb(order,:,:) ;
    % imagesc(time2use, 1:size(psth_rgb,1), all_matrix2plot)

    psth_rgb = zeros(size(all_psth,1), size(all_psth,2), 3);

    for j = 1:3
        % Fill image values
        psth_rgb(all_psth_indexes==j,:,ceil(color_Codes{j})==1) = ...
            repmat(all_psth(all_psth_indexes==j,:), 1, 1, sum(ceil(color_Codes{j})));
    end

    % Clamp and normalize
    psth_rgb(psth_rgb > c_lim(2)) = c_lim(2);
    psth_rgb(psth_rgb < c_lim(1)) = c_lim(1);
    psth_rgb = (psth_rgb - c_lim(1)) / (c_lim(2) - c_lim(1));

    % Create left border (category color bar)
    nRows = size(psth_rgb,1);
              % number of columns for the border
    border_rgb = zeros(nRows, 1, 3);  % single-column border

    for j = 1:3
        border_rgb(all_psth_indexes==j, 1, ceil(color_Codes{j})==1) = 1;
    end

    % Expand border to desired width
    border_rgb = repmat(border_rgb, 1, border_width, 1);

    % Concatenate border + image
    all_matrix2plot = cat(2, border_rgb, psth_rgb);
    all_matrix2plot = all_matrix2plot(order,:,:);

    % Plot
    border_x = time2use(1) - (border_width:-1:1) * .25;
xvals = [border_x, time2use];
    imagesc(xvals, ...
        1:size(all_matrix2plot,1), ...
        all_matrix2plot ...
        )
    axis tight
    hold on
    axis xy
    
    plot(onset_time, [.5 size(psth_rgb,1)+.5], 'w')
    plot(offset_time, [.5 size(psth_rgb,1)+.5], 'w')
    % clim(c_lim)
    xlim([xvals(2) x_lim(2)])
    xticklabels([])
    yticks([])
    % ylabel(group_labels{1})
    title([area_list{an}, '  ' response_type])



    %plot mean population responses
    subplot(6, numel(area_list),an + 2*numel(area_list))
    pctl2plot =  prctile(null_psth_peak,[5 95]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    hold on
    pctl2plot =  prctile(null_psth_peak,[10 90]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    pctl2plot =  prctile(null_psth_peak,[1 99]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    hold on
    plot(time2use, median(null_psth_peak, 'omitmissing'), 'Color', 'k')
    plot(time2use, median(peak_psth, 'omitmissing'), 'Color', color_Codes{1}, 'LineWidth',2)
    plot(onset_time, y_lim, ':k')
    plot(offset_time, y_lim, ':k')
    ylim(y_lim)
    xlim([xvals(2) x_lim(2)])


    %
    subplot(6, numel(area_list),an + 3*numel(area_list))
    pctl2plot =  prctile(null_psth_trough,[5 95]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    hold on
    pctl2plot =  prctile(null_psth_trough,[10 90]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    pctl2plot =  prctile(null_psth_trough,[1 99]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    hold on
    plot(time2use, median(null_psth_trough, 'omitmissing'), 'Color', 'k')
    plot(time2use, median(trough_psth, 'omitmissing'), 'Color', color_Codes{2}, 'LineWidth',2)
    plot(onset_time, y_lim, ':k')
    plot(offset_time, y_lim, ':k')
    ylim(y_lim)
    xlim([xvals(2) x_lim(2)])



    subplot(6, numel(area_list),an + 4*numel(area_list))
    pctl2plot =  prctile(null_psth_nonentrained,[5 95]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    hold on
    pctl2plot =  prctile(null_psth_nonentrained,[10 90]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    pctl2plot =  prctile(null_psth_nonentrained,[1 99]);
    fill([time2use fliplr(time2use)], [pctl2plot(1,:) fliplr(pctl2plot(2,:))],'k', 'FaceAlpha',face_alpha, 'EdgeColor','none')
    hold on
    plot(time2use, median(null_psth_nonentrained), 'Color', 'k')
    plot(time2use, median(nonentrained_psth, 'omitmissing'), 'Color', color_Codes{3}, 'LineWidth',2)
    plot(onset_time, y_lim, ':k')
    plot(offset_time, y_lim, ':k')
   
    ylim(y_lim)
    xlim([xvals(2) x_lim(2)])


    peak_pop = median(peak_pctl_activation(response_time_index));
    trough_pop = median(trough_pctl_activation(response_time_index));
    nonentrained_pop = median(nonentrained_pctl_activation(response_time_index));

    subplot(6, numel(area_list),an + 5*numel(area_list))
    plot(time2use, nonentrained_pctl_activation, 'Color',color_Codes{3})
    hold on
    plot(time2use, trough_pctl_activation, 'Color', color_Codes{2})
    plot(time2use, peak_pctl_activation, 'Color', color_Codes{1})
    plot(onset_time, sw_lim, ':k')
    plot(offset_time,sw_lim, ':k')
     plot([xvals(2) x_lim(2)], [2 2], ':k')
    plot([xvals(2) x_lim(2)], -[2 2], ':k')
    title(num2str(round([peak_pop trough_pop nonentrained_pop]-50)/50))
    ylim(sw_lim)
    xlim([xvals(2) x_lim(2)])
    xticklabels([])


    sgtitle([psth2use,  ' // ' ,freq2use])




end
