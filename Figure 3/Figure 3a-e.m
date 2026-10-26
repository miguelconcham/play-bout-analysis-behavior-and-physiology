%% Figure 3 — phase locking (peak / trough / unlocked)
% Preferred phase, locking strength, and rate for cells grouped as
% peak-locked, trough-locked, or unlocked. Default load is delta; switch
% the block in section 3 for theta / beta / gamma. Analysis .mat files are
% under Data\.

%% 1 Folders
data_root = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data';
saving_folder = [data_root, '\Analysis results\phase locking data'];
figures_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure 3 Updated';

%% 2 FIg 3a Peak and trough schematic
% Sine wave split at 0: red = peak half-cycle, blue = trough half-cycle.
figure
x = (1:2000)/1000;
y_red = sin(2*pi*((1:2000)/1000)-pi/2 );
y_red(y_red<0) = NaN;
y_blue = sin(2*pi*((1:2000)/1000)-pi/2 );
y_blue(y_blue>=0) = NaN;


plot(360*x - 180,y_red, 'r' )
hold on
plot(360*x - 180,y_blue, 'b' )
xticks([- 180 -90 0 90 180 270 360])


%% 3 Load phase-coupling data
% Concatenate per-session phase_struct into all_neurons. Uncomment one
% structure+names pair below to plot that frequency band.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% here you select the type of locking you want to plot %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load delta or theta
load([saving_folder,'\delta_phase_couplig_structure_updated_with_non_playbouts.mat'],'phase_struct');
load([saving_folder,'\delta_phase_couplig_animal_names_updated_with_non_playbouts.mat'],'animal_names');

% load([saving_folder,'\theta_phase_couplig_structure_updated_with_non_playbouts.mat'],'phase_struct');
% load([saving_folder,'\theta_phase_couplig_animal_names_updated_with_non_playbouts.mat'],'animal_names');
 % 
 % load([saving_folder,'\low_gamma_phase_couplig_structure_updated_with_non_playbouts_peak_peaks.mat'],'phase_struct');
 % load([saving_folder,'\low_gamma_phase_couplig_animal_names_updated_with_non_playbouts_peak_peaks.mat'],'animal_names');
 % 
 % load([saving_folder,'\betha_phase_couplig_structure_updated_with_non_playbouts_peak_peaks.mat'],'phase_struct');
 % load([saving_folder,'\betha_phase_couplig_animal_names_updated_with_non_playbouts_peak_peaks.mat'],'animal_names');
% 
%  load([saving_folder,'\gamma_couplig_structure_updated_with_non_playbouts_peak_peaks.mat'],'phase_struct');
% load([saving_folder,'\gamma_couplig_animal_names_updated_with_non_playbouts_peak_peaks.mat'],'animal_names');
% 




% delta_bin_size = 0.01;
phase_prop_names = {'PreferedAngle','MVL','MVLPval','PPC','PPCPval','MeanRate', 'Id'};


all_session_phase_stats = [];
all_session_psth        = [];
entire_recording_psth = [];

all_neurons = [];
for ns = 1:numel(phase_struct)

    all_session_phase_stats = cat(2,all_session_phase_stats, phase_struct(ns).session_phase_stats(1:2,:,:));
    all_session_psth = cat(2,all_session_psth, phase_struct(ns).session_psth(1:2,:,:));
    entire_recording_psth = cat(2,entire_recording_psth, phase_struct(ns).entire_recording_psth);
    entre_session_stats = phase_struct(ns).entire_recording_phase_stats;
    this_session_cluster_info = phase_struct(ns).cluster_info;


    sub_Table =  array2table(squeeze(phase_struct(ns).session_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.Partner1 =sub_Table;

    sub_Table =  array2table(squeeze(phase_struct(ns).session_phase_stats(2,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.Partner2 =sub_Table;

    sub_Table =  array2table(squeeze(phase_struct(ns).play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.Play =sub_Table;

    sub_Table =  array2table(squeeze(phase_struct(ns).pre_play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.PrePlay =sub_Table;

    sub_Table =  array2table(squeeze(phase_struct(ns).non_play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.NonPlay =sub_Table;

    sub_Table =  array2table(squeeze(phase_struct(ns).pre_non_play_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.PreNonPlay =sub_Table;

    sub_Table =  array2table(squeeze(phase_struct(ns).entire_recording_phase_stats(1,:,:)));
    sub_Table.Properties.VariableNames = phase_prop_names;
    this_session_cluster_info.EntireSession =sub_Table;

    this_session_cluster_info.session = repmat(animal_names(ns,1),size(this_session_cluster_info,1),1);

    all_neurons = [all_neurons; this_session_cluster_info];
end
%%
% save([saving_folder,'\gamma_all_neurons_v2.mat'],'all_neurons');
%% (Not shown) 4 Entrainment per area and PSTHs
% Heatmaps / mean PSTHs of locked vs unlocked cells in each structure.
psth_edges = phase_struct(1).edges_freq;
psth_centers = .5*(psth_edges(1:end-1) + psth_edges(2:end));

% x_lim = [-.025 .025];
x_lim= [psth_centers(1) psth_centers(end)];


center_index = psth_centers>=-.25 & psth_centers<=.25;
selected_bin_size=0.001;
freq2ilter = 40;
c_lim =[-10 10] ;
y_lim = [-5 5]
alpha = 0.01;
figure
area_list = {'SupCol' 'DLPAG'	'LPAG'	'VLPAG' 'DR' 'isRt',{'SupCol' 'DLPAG'	'LPAG'	'VLPAG' 'DR'}};

modulation_range_entrained = cell(numel(area_list) ,4);
modulation_range_non_entrained = cell(numel(area_list) ,4);

modulation_value = cell(numel(area_list) ,1);
modulation_value_stacked = [];
pctl_list = [25 50 75];
pctls = nan(numel(area_list),numel(pctl_list));
median_modulaton = nan(numel(area_list),1);
for an=1:numel(area_list)
    area_index   = ismember(all_neurons.area, area_list{an});
    entrained_index = all_neurons.EntireSession.PPCPval<alpha & ~isnan(all_neurons.EntireSession.PPC);

    psth_entrained = squeeze(entire_recording_psth(1,entrained_index & area_index,:));
    if size(psth_entrained,2)==1
        psth_entrained = psth_entrained';
    end
    entrained_angles = all_neurons.EntireSession.PreferedAngle(entrained_index & area_index);
    entrained_mvl = all_neurons.EntireSession.MVL(entrained_index & area_index);



    for j=1:size(psth_entrained,1)
        psth_entrained(j,:) = smooth(psth_entrained(j,:),round((1/freq2ilter)/selected_bin_size));
        psth_entrained(j,:) = 100*(psth_entrained(j,:) - mean(psth_entrained(j,:),'omitmissing'))/mean(psth_entrained(j,:),'omitmissing');
    end

    modulation_range_entrained{an,1} = range(psth_entrained(:,center_index),2);
    modulation_range_entrained{an,2} = entrained_angles;
    modulation_range_entrained{an,3} = entrained_mvl;
    modulation_range_entrained{an,4} = psth_entrained;

    modulation_value{an} = range(psth_entrained,2);
    pctls(an,:) = prctile(range(psth_entrained,2), pctl_list);
    median_modulaton(an) = median(range(psth_entrained,2));
    modulation_value_stacked = [modulation_value_stacked;[range(psth_entrained,2)*0+an range(psth_entrained,2)]];


    psth_non_entrained = squeeze(entire_recording_psth(1,~entrained_index & area_index,:));
    non_entrained_angles = all_neurons.EntireSession.PreferedAngle(~entrained_index & area_index);
    non_entrained_mvl = all_neurons.EntireSession.MVL(~entrained_index & area_index);
    for j=1:size(psth_non_entrained,1)
        psth_non_entrained(j,:) = smooth(psth_non_entrained(j,:),round((1/freq2ilter)/selected_bin_size));
        psth_non_entrained(j,:) = 100*(psth_non_entrained(j,:) - mean(psth_non_entrained(j,:),'omitmissing'))/mean(psth_non_entrained(j,:),'omitmissing');
    end

    modulation_range_non_entrained{an,1} = range(psth_non_entrained(:,center_index),2);
    modulation_range_non_entrained{an,2} = non_entrained_angles;
    modulation_range_non_entrained{an,3} = non_entrained_mvl;
    modulation_range_non_entrained{an,4} = psth_non_entrained;

    entrained_angles(isnan(entrained_angles)) = 0;
    [sorted_angles, order] = sort(entrained_angles);
    matrix2plot = psth_entrained(order,:);
    subplot(4,numel(area_list) ,an)
    hold on
    y_ticks = 180*sorted_angles/pi;
    imagesc(psth_centers,y_ticks ,matrix2plot)
    hold on
    axis xy
    [pos_90, loc_90]    = min(abs((180*sorted_angles/pi )-90));   
    [pos_270, loc_270]  = min(abs((180*sorted_angles/pi)+90));
    [pos_0, loc_0]       = min(abs((180*sorted_angles/pi)));
    clim(c_lim)
    xlim(x_lim)
    ylim([-180 180])
    plot([psth_centers(1) psth_centers(end)],y_ticks([loc_0 loc_0]), 'w')
    plot([psth_centers(1) psth_centers(end)],y_ticks([loc_270 loc_270]), 'w')
    plot([psth_centers(1) psth_centers(end)],y_ticks([loc_90 loc_90]), 'w')
    % yticks([-180 y_ticks([loc_270 loc_0 loc_90])' 180])
    yticklabels({'-180','-90','0','90','180'})
    title(area_list{an})

    subplot(4,numel(area_list) ,numel(area_list) + an)
    selection = max(abs((matrix2plot-repmat(mean(matrix2plot,2),1,size(matrix2plot,2)))./repmat(std(matrix2plot,[],2),1,size(matrix2plot,2))),[],2)<8;

    index = sorted_angles>pi/2 | sorted_angles<-pi/2;
    [~, ~, ci] = ttest(matrix2plot(index & selection,:));
    if length(ci)>2
    fill([psth_centers fliplr(psth_centers)], [ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(psth_centers, mean(matrix2plot(index & selection,:)), 'b', 'LineWidth',2)
    end

    [~, ~, ci] = ttest(matrix2plot(~index & selection,:));
    if length(ci)>2
    fill([psth_centers fliplr(psth_centers)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(psth_centers, mean(matrix2plot(~index & selection,:)), 'r', 'LineWidth',2)
    end
    ylim(y_lim)
    xlim(x_lim)
      title({num2str(100*[sum(index) sum(~index)]/numel(index)),num2str([sum(index) sum(~index)])})

    subplot(4,numel(area_list) ,2*numel(area_list) +an)
    hold on
    non_entrained_angles(isnan(non_entrained_angles)) =0;
    [sorted_angles, order] = sort(non_entrained_angles);
    matrix2plot = psth_non_entrained(order,:);

    imagesc(psth_centers, 180*sorted_angles/pi,matrix2plot)

    ylim([-180 180])
    plot([psth_centers(1) psth_centers(end)],[0 0], 'w')
    plot([psth_centers(1) psth_centers(end)],-[90 90], 'w')
    plot([psth_centers(1) psth_centers(end)],[90 90], 'w')
    axis xy
    clim(c_lim)
    yticks([-180 -90 0 90 180])
    xlim(x_lim)

    subplot(4,numel(area_list) ,3*numel(area_list) +an)
    selection = max(abs((matrix2plot-repmat(mean(matrix2plot,2),1,size(matrix2plot,2)))./repmat(std(matrix2plot,[],2),1,size(matrix2plot,2))),[],2)<Inf;
    index = sorted_angles>-pi & sorted_angles<0;
    [~, ~, ci] = ttest(matrix2plot(index & selection,:));
    fill([psth_centers fliplr(psth_centers)], [ci(1,:) fliplr(ci(2,:))], 'b', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(psth_centers, mean(matrix2plot(index & selection,:)), 'b', 'LineWidth',2)
  


    [~, ~, ci] = ttest(matrix2plot(~index & selection,:));
    fill([psth_centers fliplr(psth_centers)], [ci(1,:) fliplr(ci(2,:))], 'r', 'FaceAlpha',.25, 'EdgeColor','none')
    hold on
    plot(psth_centers, mean(matrix2plot(~index & selection,:)), 'r', 'LineWidth',2)
    ylim(y_lim)
    xlim(x_lim)



    pause(.1)
end

%% 5 Oscillation frequency of locked neurons
% Peak spacing and a sine fit on each locked cell’s PSTH (set min_freq
% higher for gamma so slow rate changes do not hide the oscillation).
% min_freq = 35; % change to 0 for lower frequencies, like delta and theta,
% for gamma the slow modualtions in firing rate can conceal the high 
% frequency modualtion,  and then a higher tresh (like 35) is needed;
min_freq = 0; 
time_rrange = [-.2 .2]; % also change range for the relevant time scale
locked_neurons = psth_entrained;
bin_size = mean(diff(psth_centers));
smooth_window = round(.005/bin_size);

time_index = psth_centers>=time_rrange(1) & psth_centers<=time_rrange(2);
selected_time = psth_centers(time_index);
frequency_range = [1 5];
estiamted_frequencies = nan(size(locked_neurons,1),1);
measured_frequency = nan(size(locked_neurons,1),1);
modulation_amp = nan(size(locked_neurons,1),1);
peak_distance = nan(size(locked_neurons,1),1);
for j=1:size(locked_neurons,1)
    y = movmean(locked_neurons(j,time_index),smooth_window);
    x = selected_time;
    [~,peaks_mean_positive] = findpeaks(y, 'MinPeakDistance',min((1/frequency_range(2))/bin_size,78), 'MinPeakProminence', std(y));
    [~,peaks_mean_negative] = findpeaks(-y, 'MinPeakDistance',min((1/frequency_range(2))/bin_size,78), 'MinPeakProminence', std(y));
    consecutive_peaks = sort([peaks_mean_positive,peaks_mean_negative]);
    per = 2*mean(diff(consecutive_peaks))*bin_size;
    measured_frequency(j) = 1/per;
    if 1/per<min_freq
        per = mean(frequency_range);
    end

    if numel(consecutive_peaks) <2 || per==0
        per = mean(frequency_range);
    end
    peak_distance(j) = 1/per;
    mdl = fittype('a*sin(b*x + c) + d*x + e','indep','x');
    [fittedmdl,gof] = fit(x',y',mdl,'start',[rand(),1/(per/(2*pi)),rand(), rand(), mean(y)]) ;
    % fited_y = fittedmdl(x);
    estiamted_frequencies(j) = fittedmdl.b/(2*pi);
    modulation_amp(j) = fittedmdl.a/fittedmdl.b;
end
%% (Reproted in text) 6 Frequency and modulation histograms
% Distributions of fitted frequency, modulation depth, and seed frequency.
title('estimated frequencies')
xlabel('Hz')
figure
subplot(1,3,1)

histogram(estiamted_frequencies, 0:.1:20,'Normalization','percentage')
title('Estimated frequency')
xlim([0 10])
subplot(1,3,2)
histogram(abs(modulation_amp), 0:.1:100,'Normalization','percentage')
title('proportion of modualtion from mean rate')
xlim([0 5])

subplot(1,3,3)
histogram(peak_distance, 'Normalization','percentage')
title('Seed Frequency used for estimating freq')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot saving line example %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figures_folder2 = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure call supplementary'
% print(gcf,'-vector','-dsvg',[figures_folder, '\phaselocking delt aall areas together theta.svg'])

%% 7 (Fig 3c) Bimodal preferred-angle distribution
% Histogram of preferred phase for significantly locked cells in one area.
% Mixture of two von Mises; bootstrap LR is saved (do not rerun unless needed).

area_list = unique(all_neurons.area);

freq2use        = 'EntireSession'; %options:   ThetaEntireSession DeltaEntireSession
alpha_tresh = 0.01;
an=2;
area_index = ismember(all_neurons.area,area_list{an} );
entrainment_index = all_neurons.(freq2use).PPCPval<alpha_tresh;
this_Area_Angles = all_neurons.(freq2use).PreferedAngle(~isnan(all_neurons.(freq2use).PreferedAngle) & entrainment_index);
this_area_names = all_neurons.area(~isnan(all_neurons.(freq2use).PreferedAngle) & entrainment_index);
figure
polarhistogram(this_Area_Angles, -pi:(pi/32):pi)
hold on

[mu, kappa, w, LL] = vm2_mixture_EM(this_Area_Angles+pi, 10000);
LL2 = LL(end);
r_lim = rlim;
polarplot(mu([1 1])-pi,r_lim, 'r')
hold on
polarplot(mu([2 2])-pi,r_lim, 'b')


[theta, kapa_unimodal] = circ_vmpar(this_Area_Angles);
LL1 = sum(log(circ_vmpdf(this_Area_Angles, theta, kapa_unimodal)));

LR = 2*(LL2 - LL1);
legend('Phase distribution','estiamted mean phase peak','estimated mean phase trough' )

%% Run to estiamte estiatsitc only, its already saved dont run again unless needed
% figure_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure 7 Inputs';
% 
% N = numel(this_Area_Angles);
% nBoot = 10000;
% LR_boot = zeros(nBoot,1);
% 
% for b = 1:nBoot
%     boot_sample = circ_vmrnd(theta, kapa_unimodal, N);
% 
%     % Fit 2-component mixture
%     [~, ~, ~, LL_boot] = vm2_mixture_EM(boot_sample, 100, 1e-6, false); % fewer iter for speed
%     LL2_boot = LL_boot(end);
% 
%     % Log-likelihood under unimodal null
%     LL1_boot = sum(log(circ_vmpdf(boot_sample, theta, kapa_unimodal)));
% 
%     LR_boot(b) = 2*(LL2_boot - LL1_boot);
% end
% 
% p_value = mean(LR_boot >= LR);
% 
% save([figure_folder, '\LR_boot_theta.mat'], 'LR_boot')
%% 8 (Reported in text) Mixture fit vs bootstrap LR
% Overlay two-component circular mixture on the angle histogram; compare
% observed LR to the saved bootstrap (LR_boot.mat = delta, LR_boot_theta = theta).
figure_folder = [data_root, '\Figure codes\Figure 7 Inputs'];
% figures_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure 5 Inputs';
load([figure_folder, '\LR_boot.mat'], 'LR_boot') % for delta
% load([figure_folder, '\LR_boot_theta.mat'], 'LR_boot') % for theta
anged_edges = -pi:(pi/32):pi;
angle_counts = histcounts(this_Area_Angles, anged_edges);
angle_counts = angle_counts/sum(angle_counts);

angle_centers = .5*(anged_edges(1:end-1)+anged_edges(2:end));

mixtured_model  = w .* circ_vmpdf(angle_centers, mu(1)-pi, kappa(1)) + ...
      (1-w) .* circ_vmpdf(angle_centers, mu(2)-pi, kappa(2));
mixtured_model = mixtured_model'/sum(mixtured_model);
figure
plot([angle_centers ,angle_centers+2*pi],[angle_counts angle_counts], 'k')
hold on
plot([angle_centers ,angle_centers+2*pi],[mixtured_model mixtured_model], 'b')


figure
subplot(2,1,1)
histogram(LR_boot, 250, 'FaceColor', 'k', 'EdgeColor', 'none')
y_lim = ylim;
hold on
plot([LR LR],y_lim,'r')


xscale log
xlim([.1 LR+10])

mixtured_model  = w .* circ_vmpdf(anged_edges, mu(1)-pi, kappa(1)) + ...
      (1-w) .* circ_vmpdf(anged_edges, mu(2)-pi, kappa(2));
mixtured_model = mixtured_model'/sum(mixtured_model);

subplot(2,1,2)
polarplot(anged_edges,[angle_counts,angle_counts(1)], 'k' )
hold on
polarplot(anged_edges,mixtured_model, 'r' )
r_lim = rlim;

polarplot([mu(1) mu(1)]-pi,r_lim, 'b' )

polarplot([mu(2) mu(2)]-pi,r_lim, 'r' )
%%

print(gcf,'-vector','-dsvg',[figures_folder, '\two distribution plot theta.svg'])


%% 9 Mean preferred angle per structure
% Watson–Williams test on trough-preferring cells across PAG / SC.
% Convert labels to grouping variable

angles2test = (this_Area_Angles>pi/2 | this_Area_Angles<-pi/2) & ismember(this_area_names,{'DLPAG','LPAG','VLPAG','SupCol'});
angles = this_Area_Angles(angles2test);
labels = this_area_names(angles2test);

[group_names, ~, group_idx] = unique(labels);
nGroups = numel(group_names);

[p, table] = circ_wwtest(angles, group_idx);
disp(p)

%% (not reported) 10 Pairwise circular tests, no significnat difference in mean angle per structure
% Pairwise circ_wwtest between structures.

pairs = nchoosek(1:nGroups, 2);
pvals = nan(size(pairs,1),1);

for i = 1:size(pairs,1)
    g1 = pairs(i,1);
    g2 = pairs(i,2);
    
    idx = group_idx==g1 | group_idx==g2;
    
    pvals(i) = circ_wwtest(angles(idx), group_idx(idx));
end

% Bonferroni correction
pvals_corr = min(pvals * size(pairs,1), 1);
disp([pairs, pvals_corr])

%% 11 Plot mean preferred angle per structure
circular_violin_plot(angles, group_idx, true); %the function create a figure 
x_tick_orde = str2double(xticklabels);

xticklabels(group_names(x_tick_orde))


 angular_range = circular_violin_plot_percentiles(angles, group_idx, [20 80], true)
 x_tick_orde = str2double(xticklabels);

xticklabels(group_names(x_tick_orde))


%% saving if needed
figures_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure coincident phases';

print(gcf,'-vector','-dsvg',[figures_folder, '\phaselocking delt  angle violin plots.svg'])



%%
figures_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure coincident phases';

print(gcf,'-vector','-dsvg',[figures_folder, '\phaselocking delta  phases percentiles and time equivalents.svg'])

%% 12 (Not used any more) Percent locked / trough / peak (older counts)
area_list = {'SupCol','DLPAG','LPAG','VLPAG', 'DR'}
number_per_area = nan(numel(area_list),4);
non_entrained_tresh = 0.1;
for an=1:numel(area_list)
    area_index = ismember(all_neurons.area,area_list{an} );
    entrainment_index = all_neurons.(freq2use).PPCPval<alpha_tresh;
    this_Area_Angles = all_neurons.(freq2use).PreferedAngle(~isnan(all_neurons.(freq2use).PreferedAngle) & entrainment_index & area_index);

    number_per_area(an,1) = sum(this_Area_Angles>pi/2 | this_Area_Angles<-pi/2);
     number_per_area(an,2) = sum(~(this_Area_Angles>pi/2 | this_Area_Angles<-pi/2));
     number_per_area(an,3) = sum(~entrainment_index & area_index & ~isnan(all_neurons.(freq2use).PreferedAngle));
     number_per_area(an,4) = sum( area_index & ~isnan(all_neurons.(freq2use).PreferedAngle));

end


figure
subplot(4,1,1)
bar(diag(1./number_per_area(:,4))*[number_per_area(:,4)-number_per_area(:,3) number_per_area(:,3)], 'stacked')
legend({'Entrained','Non entrained'})
ylabel('%')
subplot(4,1,2:4)
bar(diag(1./sum(number_per_area(:,1:2),2))*number_per_area(:,1:2), 'stacked')
xticklabels(area_list)
ylabel('%')
legend({'Trough','Peak'})




figure
subplot(4,1,1)
bar([number_per_area(:,4)-number_per_area(:,3) number_per_area(:,3)])
subplot(4,1,2:4)
bar(number_per_area(:,1:2))
xticklabels(area_list)

%% 13 (Fig 3d) Percent per area (updated 2026-02-06)
% Stacked bars: trough vs peak vs unlocked, by structure.
freq2use        = 'EntireSession'; %options:   ThetaEntireSession DeltaEntireSession
alpha_tresh = 0.01;
area_list = {'SupCol','DLPAG','LPAG','VLPAG', 'DR'}
number_per_area = nan(numel(area_list),5);
non_entrained_tresh = 0.1;
for an=1:numel(area_list)
    area_index = ismember(all_neurons.area,area_list{an} );
    entrainment_index = all_neurons.(freq2use).PPCPval<alpha_tresh;
    unlocked_index =  all_neurons.(freq2use).PPCPval>=non_entrained_tresh;
    this_Area_Angles = all_neurons.(freq2use).PreferedAngle(~isnan(all_neurons.(freq2use).PreferedAngle) & entrainment_index & area_index);

    number_per_area(an,1) = sum(this_Area_Angles>pi/2 | this_Area_Angles<-pi/2);
     number_per_area(an,2) = sum(~(this_Area_Angles>pi/2 | this_Area_Angles<-pi/2));
     number_per_area(an,3) = sum(entrainment_index & area_index & ~isnan(all_neurons.(freq2use).PreferedAngle));
     number_per_area(an,4) = sum( area_index & ~isnan(all_neurons.(freq2use).PreferedAngle));
      number_per_area(an,5) = sum(unlocked_index & area_index & ~isnan(all_neurons.(freq2use).PreferedAngle));

end

counts = [number_per_area(:,4)-number_per_area(:,1)-number_per_area(:,2)-number_per_area(:,5) number_per_area(:,1) number_per_area(:,2) number_per_area(:,5)];
figure
percentages = diag(1./number_per_area(:,4))*[number_per_area(:,4)-number_per_area(:,1)-number_per_area(:,2)-number_per_area(:,5) number_per_area(:,1) number_per_area(:,2) number_per_area(:,5)];
percentages = percentages(:,[2 3 4 1]);
counts = counts(:,[2 3 4 1]);
b =bar((percentages), 'stacked');
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
b(3).FaceColor = 'g';
b(4).FaceColor = 'k';
legend({'Trough-locked','Peak-locked','Unlocked','Unclassified'})
ylabel('%')




colNames ={'Trough-locked','Peak-locked','Unlocked','Unclassified'}
rowNames = {'SC','DLPAG','LPAG', 'VLPAG', 'DR', 'Total%'};

f = figure;

uitable(f, ...
    'Data', [counts;sum(counts)/sum(sum(counts))], ...
    'RowName', rowNames, ...
    'ColumnName', colNames, ...
    'Position', [20 20 600 150]);

%%
print(gcf,'-vector','-dsvg',[figures_folder, '\phaselocking theta angle and modulation.svg'])

%% 14 (Fig 3e) Partner-session locking and rate
% Compare two conditions (Partner1 vs Partner2 by default; Play vs PrePlay
% or non-play alternatives are commented). Square plots: paired cells.
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
% area_list = unique(all_neurons.area)';

% area_list =  {'4N'	'DLPAG'	'DMPAG'	'DR'	'InfCol'	'LPAG'	'LSD'	'LSI'	'SupCol'	'VLPAG'	'isRt'	'mlf'};
% area_list =  {	'DLPAG'	'DMPAG'	'DR'	'LPAG'		'SupCol'	'VLPAG'	'isRt'	'mlf'};
% area_list = {'SupCol'  'DLPAG'	'LPAG' 'VLPAG'	'DR', 'isRt'	};
y_lim = [10^-6 5]
area_list =  {	'DLPAG'	'DMPAG'	'DR'	'LPAG'		'SupCol'	'VLPAG'	};
col = strcmp(phase_prop_names, 'PPC');

concatenated_phase_differences = [];
mean_phase_differences = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% here you select the conditions you want to plot %%%%%
%%%% either as on main or supp matierial             %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% condition_1 = 'Partner1';
% condition_2 = 'Partner2';

condition_1 = 'NonPlay';
condition_2 = 'PreNonPlay';
% condition_1 = 'Play';
% condition_2 = 'PrePlay';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for an=1:numel(area_list)
    index = strcmp(all_neurons.area, area_list{an}) & all_neurons.EntireSession.PPCPval<alpha &  ~isnan(all_neurons.(condition_1).PPC + all_neurons.(condition_2).PPC);
    x_jit = (rand(sum(index),2) -.5)*.25;
    matrix2plot = [all_neurons.(condition_1).PPC(index) all_neurons.(condition_2).PPC(index)  ];

    data = [diff(matrix2plot')' ones(size(matrix2plot,1),1)*an ] ;
    if size(data,1)>1
        mean_phase_differences = [mean_phase_differences;mean(data, 'omitmissing')];
    else
        mean_phase_differences = [mean_phase_differences;data];
    end

    concatenated_phase_differences = [concatenated_phase_differences;[diff(matrix2plot')' ones(size(matrix2plot,1),1)*an all_neurons.EntireSession.PreferedAngle(index) matrix2plot]];

end

%% 15 (Fig 3e) Paired comparisons (square plots)
% Per-cell Partner1 vs Partner2 (or Play vs PrePlay): PPC and rate, trough vs peak.

  alpha=0.05;
  figure('units','normalized','outerposition',[0 0 1 1]);
  y_lim_ppc = [-0.005 0.07];
  y_lim_Rate = [-100 100];
  line_width = .25;
  names = {'Trough','Peak'};
  bar_per_area = [];
  for an=1:numel(area_list)
      subplot(2,numel(area_list)+1,an)
      index           = concatenated_phase_differences(:,2)==an;
      angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;


      p =ranksum(concatenated_phase_differences(angle_index & index,1),concatenated_phase_differences(~angle_index & index,1));
      [p1,h,stats] = signrank(concatenated_phase_differences(angle_index & index,1));
      % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))




      y2plot = atanh(concatenated_phase_differences(angle_index & index, end-1:end));
      plot(y2plot(:,1),y2plot(:,2), 'k.')
      hold on
      this_medians = median(y2plot);
      plot(y_lim_ppc, y_lim_ppc, 'r')

      axis([y_lim_ppc y_lim_ppc])

      plot(this_medians(1), this_medians(2), 'xr')
      if  ismember('zval',fields(stats))
          z_val = stats.zval;
      else
          z_val = NaN;
      end
      title([area_list{an},' ' num2str(round(p1,4)), ' ', num2str(z_val), ' ',num2str(stats.signedrank)])


      axis square
      if an==1
          ylabel({'Trough', condition_2})
      end
      xlabel(condition_1)




      subplot(2,numel(area_list)+1,numel(area_list)+1 +an)

      y2plot = atanh(concatenated_phase_differences(~angle_index & index, end-1:end));
      plot(y2plot(:,1),y2plot(:,2), 'k.')
      hold on
      this_medians = median(y2plot);
      plot(this_medians(1), this_medians(2), 'xr')
      plot(y_lim_ppc, y_lim_ppc, 'r')

      axis([y_lim_ppc y_lim_ppc])

      [p2,h,stats] = signrank(concatenated_phase_differences(~angle_index & index, end-1),concatenated_phase_differences(~angle_index & index, end));
      % [h,p2]= ttest(concatenated_phase_differences(~angle_index & index, end-1)+.2,concatenated_phase_differences(~angle_index & index, end)+.2);

      this_medians = median(concatenated_phase_differences(~angle_index & index, end-1:end));
      if  ismember('zval',fields(stats))
          z_val = stats.zval;
      else
          z_val = NaN;
      end
      title([area_list{an},' ' num2str(round(p2,4)), ' ', num2str(z_val), ' ',num2str(stats.signedrank)])

      axis square
      if an==1
          ylabel({'Peak', condition_2})
      end
      xlabel(condition_1)





      % ylim(y_lim_ppc)
      % xticks(1:4)
      % xticklabels({[names{1}, ' ', condition_1],[names{1}, ' ', condition_2],[names{2} , ' ', condition_1],[names{2}, ' ', condition_2]})
      %




  end


  subplot(2,numel(area_list)+1,numel(area_list)+1)
  index           = true(size(concatenated_phase_differences(:,2)));
  angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;



  [p1,h,stats] = signrank(concatenated_phase_differences(angle_index & index, end),concatenated_phase_differences(angle_index & index, end-1));
  % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))



  y2plot = atanh(concatenated_phase_differences(angle_index & index, end-1:end));
  plot(y2plot(:,1),y2plot(:,2), 'k.')
  hold on
  this_medians = median(y2plot);
  plot(y_lim_ppc, y_lim_ppc, 'r')

  axis([y_lim_ppc y_lim_ppc])
  plot(this_medians(1), this_medians(2), 'xr')
  if  ismember('zval',fields(stats))
      z_val = stats.zval;
  else
      z_val = NaN;
  end
  title(['ALL AREAS ', num2str(round(p1,6)), ' ', num2str(z_val), ' ',num2str(stats.signedrank), ' ', num2str(sum(angle_index & index))])
  axis square
  subplot(2,numel(area_list)+1,2*(numel(area_list)+1))
  index           = true(size(concatenated_phase_differences(:,2)));
  angle_index     = concatenated_phase_differences(:,3)<-pi/2 | concatenated_phase_differences(:,3)>pi/2;



  [p1,h,stats] = signrank(concatenated_phase_differences(~angle_index & index, end),concatenated_phase_differences(~angle_index & index, end-1));
  % [h,p1]= ttest(concatenated_phase_differences(angle_index & index,end-1),concatenated_phase_differences(angle_index & index,end))



  y2plot = atanh(concatenated_phase_differences(~angle_index & index, end-1:end));
  plot(y2plot(:,1),y2plot(:,2), 'k.')
  hold on
  this_medians = median(y2plot);
  plot(y_lim_ppc, y_lim_ppc, 'r')

  axis([y_lim_ppc y_lim_ppc])
  plot(this_medians(1), this_medians(2), 'xr')
  if ismember('zval',fields(stats))
      z_val = stats.zval;
  else
      z_val = NaN;
  end
  title(['ALL AREAS ' num2str(round(p1,6)), ' ', num2str(z_val), ' ',num2str(stats.signedrank), ' ', num2str(sum(~angle_index & index))])

  axis square


%% 16 (Fig 3a-b) Example neuron (LFP + spikes)
% Highest-MVL cell in a chosen structure/session: raw LFP, delta/theta
% peaks, and spike phase. NPX raw data are under Data\NPX data\NPX raw data.


npx_Raw_Data    = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
saving_folder   = [data_root, '\Analysis results\phase locking data'];
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];


animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];
% animal2exclude = {'B4D4 0826 Dual'};
animal2exclude = {''};
animal_list(ismember(animal_file_names,animal2exclude)) = [];
animal_names ={};
n_strctut = 1;

freq_range_1    = [1 5];
freq_range_2    = [6 12];
sr              = 2500;
filter_order    = 2000;


phase_struct = [];
% Parameters for delta
Hd_freq = designfilt('bandpassfir', ...
    'FilterOrder', filter_order, ...
    'CutoffFrequency1', freq_range_1(1), ...
    'CutoffFrequency2', freq_range_1(2), ...
    'SampleRate', sr, ...
    'DesignMethod', 'window', ...
    'Window', 'hamming');
bin_size_freq = 0.01;


% Parameters for theta
% Hd_freq = designfilt('bandpassfir', ...
% 'FilterOrder', filter_order, ...
% 'CutoffFrequency1', freq_range_2(1), ...
% 'CutoffFrequency2', freq_range_2(2), ...
% 'SampleRate', sr, ...
% 'DesignMethod', 'window', ...
% 'Window', 'hamming');
% bin_size_freq = 0.001;
%% possible sessions to chose (here you plot the neurons with highest locking value per session, within the chosen structre (area2select))
alpha = 0.01;
area2select = 'SupCol';
area_index   = strcmp(all_neurons.area, area2select);
entrained_index = all_neurons.EntireSession.PPCPval<alpha;
session_list = unique(all_neurons.session(area_index & entrained_index));
max_mvl_table = cell(numel(session_list),5);
id_loked = {};

figure
for ss=1:numel(session_list)
    subplot(4,3,ss)

    session_index= find(strcmp(all_neurons.session, session_list(ss)));
    [ max_mvl, loc] = max(all_neurons.EntireSession.PPC(session_index));

    max_mvl_table(ss,:) = {max_mvl,all_neurons.cluster_id(session_index(loc)),all_neurons.EntireSession.MeanRate(session_index(loc)),session_list{ss},session_index(loc)};
    plot(psth_centers,squeeze(entire_recording_psth(1,session_index(loc),:))/bin_size_freq, 'k')
end
max_mvl_table= cell2table(max_mvl_table);
max_mvl_table.Properties.VariableNames = {'PPC','ClusterID','MeanRate','Session', 'PsthRawIndex'};
disp(max_mvl_table)



%% estimate neurons (select neurons from the displayed table, write below the raw from max_mvl_table you want to select)

row2select = 1; % for delta 1 for theta 10 (to replicate paper figure)
plot_bool = true;
this_neuron_phase_struct = GENERATE_PHASE_COUPLING_NEURON_ID([npx_Raw_Data, '\', max_mvl_table.Session{row2select}],Hd_freq,bin_size_freq , max_mvl_table.ClusterID(row2select), plot_bool );
%%
ylim tight
print(gcf,'-vector','-dsvg',[figures_folder, '\neuron example psth vectorized theta.svg'])


%%
y_lim = [-2000 2000];
edges_tn = this_neuron_phase_struct.edges_freq;
peak_delta_psth = this_neuron_phase_struct.all_psth{1,1};
psth_centers_tn = (edges_tn(2:end) +edges_tn(1:end-1))/2;
x_lim = [-.2 .2];
figure
colormap(1-gray)
subplot(4,1,1:2)
peakswithspijkes = unique(this_neuron_phase_struct.spikes4raster{1,2}(:,2));
peak_number = nan(max(peakswithspijkes),1);
peak_number(peakswithspijkes) = 1:numel(peakswithspijkes);
plot(this_neuron_phase_struct.spikes4raster{1,2}(:,1),peak_number(this_neuron_phase_struct.spikes4raster{1,2}(:,2)), '.k')

axis tight
xlim(x_lim)
yticks([])
pbaspect([9.799 13.884 9.799]/9.799)
xticklabels([])


subplot(4,1,3:4)
plot(psth_centers_tn,mean(peak_delta_psth)/bin_size_freq, 'k')
xlim(x_lim)

pbaspect([9.799 13.884 9.799]/9.799)

%%
print(gcf,'-vector','-dsvg',[figures_folder, '\neuron example psth theta.svg'])

%% now plot neurons

figure
x_lim = [217 221]; % for theta
% x_lim = [286 292]; % for delta
subplot(2,5,1:4)
hold off
real_lfp = this_neuron_phase_struct.all_lfp{1,1}-mean(this_neuron_phase_struct.all_lfp{1,1});
real_lfp = smooth(real_lfp, (1/50)*2500);
real_time = this_neuron_phase_struct.all_lfp{1,3};
lfp_phases = this_neuron_phase_struct.all_lfp{1,4};
this_neuron_spikes = this_neuron_phase_struct.all_spikes{1}';
this_neuron_phases = this_neuron_phase_struct.all_phases{1,1};

index2plot_lfp  = real_time>=x_lim(1) & real_time<=x_lim(2);
index2plot_spikes = this_neuron_spikes>=x_lim(1) & this_neuron_spikes<=x_lim(2);
filtered_lfp = this_neuron_phase_struct.all_lfp{1,2}-mean(this_neuron_phase_struct.all_lfp{1,2});
plot(real_time(index2plot_lfp),real_lfp(index2plot_lfp) , 'k')
hold on

plot(real_time(index2plot_lfp), filtered_lfp(index2plot_lfp), 'r', 'LineWidth',2)
plot([this_neuron_spikes(index2plot_spikes);this_neuron_spikes(index2plot_spikes)],y_lim, 'b')
xlim(x_lim)
dela_peak_times     = this_neuron_phase_struct.all_loc_times{1,1};
delta_peak_indexes  = ismember(real_time,dela_peak_times);
delta_peak_phases = lfp_phases(delta_peak_indexes);

subplot(2,5,5)
polarhistogram(lfp_phases, -pi:(pi/32):pi, 'FaceColor', 'k','EdgeColor','none', 'Normalization','percentage')
hold on
polarhistogram(this_neuron_phases, -pi:(pi/32):pi, 'FaceColor', 'r', 'EdgeColor','none', 'Normalization','percentage')
hold on
r_lim = rlim;
polarplot([1 1]*circ_mean(delta_peak_phases'), r_lim,  'b')



sub_time = real_time(index2plot_lfp);
dela_peak_times_this_range = dela_peak_times(dela_peak_times>=x_lim(1) & dela_peak_times<=x_lim(2));
positive_phases =filtered_lfp(index2plot_lfp);
negative_indexes = lfp_phases(index2plot_lfp)<-pi/2 | lfp_phases(index2plot_lfp)>pi/2;
positive_phases( negative_indexes) = NaN;

negative_phases = filtered_lfp(index2plot_lfp);
negative_phases(~negative_indexes) = NaN;




subplot(2,5,5 + (1:4))

plot(sub_time,negative_phases, 'b')
hold on
plot(sub_time,positive_phases, 'r')
y_lim = [-1000 1000];
plot([dela_peak_times_this_range;dela_peak_times_this_range], y_lim, 'r')

subplot(2,5,10)
polarhistogram(delta_peak_phases, -pi:(pi/32):pi,'FaceColor', 'r', 'EdgeColor','none')


%%
print(gcf,'-vector','-dsvg',[figures_folder, '\neuron example lfp and spikes theta.svg'])


%% 17 (Former cross-freq-coupling figure) Theta and delta locking together
% Reload theta, gamma, and delta all_neurons tables and plot both bands.

saving_folder = [data_root, '\Analysis results\phase locking data'];

load([saving_folder,'\theta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };
all_neurons_TD = all_neurons;
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner1'}))         = {'ThetaPartner1'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Partner2'}))         = {'ThetaPartner2'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'Play'}))             = {'ThetaPlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'PrePlay'}))          = {'ThetaPrePlay'};
all_neurons_TD.Properties.VariableNames(ismember(all_neurons_TD.Properties.VariableNames,{'EntireSession'}))    = {'ThetaEntireSession'};




load([saving_folder,'\gamma_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'})) =     {'isRt'  };

all_neurons_TD.GammaPartner1                            = all_neurons.Partner1;
all_neurons_TD.GammaPartner2                            = all_neurons.Partner2;
all_neurons_TD.GammaEntireSession                       = all_neurons.EntireSession;
all_neurons_TD.GammaPlay                                = all_neurons.Play;
all_neurons_TD.GammaPrePlay                             = all_neurons.PrePlay;



load([saving_folder,'\delta_all_neurons_v2.mat'],'all_neurons');
all_neurons.area(ismember(all_neurons.area, {'isRT'}))  =     {'isRt'  };
all_neurons_TD.DeltaPartner1                            = all_neurons.Partner1;
all_neurons_TD.DeltaPartner2                            = all_neurons.Partner2;
all_neurons_TD.DeltaEntireSession                       = all_neurons.EntireSession;
all_neurons_TD.DeltaPlay                                = all_neurons.Play;
all_neurons_TD.DeltaPrePlay                             = all_neurons.PrePlay;
all_neurons_TD.Exited                                   = nan(size(all_neurons_TD,1),1);
all_neurons_TD.Inhibited                                = nan(size(all_neurons_TD,1),1);





areas = {'SupCol','DLPAG','LPAG','VLPAG', 'DR'};


proportions_all_areaas = nan(numel(areas),4);

for an=1:numel(areas)

no_nan = ~isnan(all_neurons_TD.DeltaEntireSession.PPCPval) & ~isnan(all_neurons_TD.ThetaEntireSession.PPCPval);
this_area_index = ismember(all_neurons_TD.area    , areas{an}) & no_nan;

deltaandtheta   = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)<=0.01 & all_neurons_TD.ThetaEntireSession.PPCPval(this_area_index)<=0.01;
onlydelta       = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)<=0.01 & all_neurons_TD.ThetaEntireSession.PPCPval(this_area_index)>0.01;
onlytheta        = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)>0.01 & all_neurons_TD.ThetaEntireSession.PPCPval(this_area_index)<=0.01;
unlocked       = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)>0.01 & all_neurons_TD.ThetaEntireSession.PPCPval(this_area_index)>0.01;

proportions_all_areaas(an,:) = [sum(deltaandtheta) sum(onlydelta) sum(onlytheta) sum(unlocked)]/sum(this_area_index);



end


figure

bar(proportions_all_areaas, 'stacked')
legend({'dleta and theta-locked','delta-locked theta-unlocked','theta-locked delta-unlocked','theta-unlocked delta-unlocked'}, 'Location','northoutside')





areas = {'SupCol','DLPAG','LPAG','VLPAG', 'DR'};


proportions_all_areaas = nan(numel(areas),4);

for an=1:numel(areas)

no_nan = ~isnan(all_neurons_TD.DeltaEntireSession.PPCPval) & ~isnan(all_neurons_TD.ThetaEntireSession.PPCPval);
this_area_index = ismember(all_neurons_TD.area    , areas{an}) & no_nan;

deltaandgamma   = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)<=0.01 & all_neurons_TD.GammaEntireSession.PPCPval(this_area_index)<=0.01;
onlydelta       = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)<=0.01 & all_neurons_TD.GammaEntireSession.PPCPval(this_area_index)>0.01;
onlygamma        = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)>0.01 & all_neurons_TD.GammaEntireSession.PPCPval(this_area_index)<=0.01;
unlocked       = all_neurons_TD.DeltaEntireSession.PPCPval(this_area_index)>0.01 & all_neurons_TD.GammaEntireSession.PPCPval(this_area_index)>0.01;

proportions_all_areaas(an,:) = [sum(deltaandgamma) sum(onlydelta) sum(onlygamma) sum(unlocked)]/sum(this_area_index);



end


figure

bar(proportions_all_areaas, 'stacked')
% legend({'dleta and gamma-locked','delta-locked gamma-unlocked','gamma-locked delta-unlocked','gamma-unlocked delta-unlocked'}, 'Location','northoutside')
legend({'dleta and gamma-locked','delta-locked gamma-unlocked','gamma-locked delta-unlocked','gamma-unlocked delta-unlocked'}, 'Location','northoutside')


%% 18 LPAG-centered preferred phase
% Re-references the preferred angles loaded in section 3 to the mid-LPAG
% channel. Sections 1-17 are unchanged. This block only reads phase_struct,
% animal_names, and data_root.
%
% Coherence estimates are the full pairwise maps from Estimate_coherence_maps.m
% (GENERATE_COHERNECE_MAPS_STRUCT), saved as:
%   coherence_pli_struct_2026.mat
%   coherence_pli_structure_animal_names2026.mat
% under Analysis results\psth power by frequency and behavior.
% The phase lag of channel X is the column of that matrix against the
% mid-LPAG site on the same probe (phase(X) - phase(mid-LPAG)).
% Delta uses mean_phaseLag_freq1. Theta uses mean_phaseLag_freq2.
% The centered angle is
%
%   wrap(preferred angle - channel phase lag)
%
% The matrix row is the 1-based LFP row, matched to cluster_info.ch.
% Peak versus trough uses the section 13 split: trough outside
% [-pi/2, pi/2], peak inside it. Locked neurons with no mid-LPAG lag are
% left out of the new figures and of the retention counts.
%
% Run sections 1 and 3 first. lpag_phase_lag_band has to match that load.

lpag_phase_lag_band = 'delta'; % 'delta' or 'theta'

if ~exist('phase_struct', 'var') || ~exist('animal_names', 'var') || ~exist('data_root', 'var')
    error('Figure3:MissingPhaseData', ...
        'Run sections 1 and 3 before this section so phase_struct, animal_names, and data_root exist.');
end
if numel(phase_struct) ~= size(animal_names, 1)
    error('Figure3:SessionCount', ...
        'phase_struct has %d sessions and animal_names has %d rows.', ...
        numel(phase_struct), size(animal_names, 1));
end

% Pairwise coherence files written by Estimate_coherence_maps.m
lpag_coherence_folder = [data_root, '\Analysis results\psth power by frequency and behavior'];
lpag_structure_file = [lpag_coherence_folder, '\coherence_pli_struct_2026.mat'];
lpag_names_file = [lpag_coherence_folder, '\coherence_pli_structure_animal_names2026.mat'];
disp('Coherence estimates from Estimate_coherence_maps.m')
disp(['  ', lpag_structure_file])
disp(['  ', lpag_names_file])

if exist(lpag_structure_file, 'file') ~= 2 || exist(lpag_names_file, 'file') ~= 2
    error('Figure3:MissingCoherence', ...
        ['Pairwise coherence files were not found. Run Estimate_coherence_maps.m. ', ...
        'Expected %s and %s.'], lpag_structure_file, lpag_names_file);
end

switch lower(lpag_phase_lag_band)
    case 'delta'
        lpag_lag_field = 'mean_phaseLag_freq1';
    case 'theta'
        lpag_lag_field = 'mean_phaseLag_freq2';
    otherwise
        error('Figure3:PhaseLagBand', ...
            'lpag_phase_lag_band must be ''delta'' or ''theta''.');
end

lpag_repo_root = fileparts(data_root);
addpath(fullfile(lpag_repo_root, 'Custom functions'));
addpath(fullfile(lpag_repo_root, 'Figure 2', 'Supporting codes v2'));
addpath(fullfile(lpag_repo_root, 'Figure 3', 'Supporting codes'));

disp('Loading pairwise coherence maps')
lpag_coherence = load(lpag_structure_file, 'pli_struct');
lpag_pli = lpag_coherence.pli_struct;
lpag_names_loaded = load(lpag_names_file, 'animal_names');
lpag_coherence_names = lpag_names_loaded.animal_names;
if numel(lpag_pli) ~= size(lpag_coherence_names, 1)
    warning('Figure3:CoherenceNames', ...
        'coherence_pli_struct_2026.mat has %d sessions and the names file has %d rows.', ...
        numel(lpag_pli), size(lpag_coherence_names, 1));
end

% One column per probe: phase of each same-probe channel versus mid-LPAG.
lpag_traces = extract_delta_phase_vs_lpag_distance(lpag_pli, lpag_coherence_names, lpag_lag_field);
clear lpag_pli lpag_coherence
if isempty(lpag_traces)
    error('Figure3:NoPhaseLags', ...
        'No mid-LPAG phase-lag column in %s.', lpag_structure_file);
end

lpag_lag_by_channel = containers.Map('KeyType', 'char', 'ValueType', 'double');
lpag_coherence_sessions = cell(numel(lpag_traces), 1);
lpag_ref_lags = [];
for lpag_k = 1:numel(lpag_traces)
    lpag_session_name = strtrim(char(lpag_traces(lpag_k).session_name));
    lpag_coherence_sessions{lpag_k} = lpag_session_name;
    lpag_rows = round(double(lpag_traces(lpag_k).channel_index(:)));
    lpag_lags = double(lpag_traces(lpag_k).phase_rad(:));
    lpag_n = min(numel(lpag_rows), numel(lpag_lags));
    for lpag_i = 1:lpag_n
        if ~isfinite(lpag_rows(lpag_i)) || ~isfinite(lpag_lags(lpag_i))
            continue
        end
        lpag_key = sprintf('%s|%d', lpag_session_name, lpag_rows(lpag_i));
        if isKey(lpag_lag_by_channel, lpag_key)
            lpag_previous = lpag_lag_by_channel(lpag_key);
            if abs(angle(exp(1i * (lpag_previous - lpag_lags(lpag_i))))) > 1e-3
                warning('Figure3:DuplicateChannelLag', ...
                    'Session %s channel %d has two different phase lags. Keeping the first.', ...
                    lpag_session_name, lpag_rows(lpag_i));
            end
        else
            lpag_lag_by_channel(lpag_key) = lpag_lags(lpag_i);
        end
    end
    lpag_ref = round(double(lpag_traces(lpag_k).lpag_center_index));
    lpag_ref_pos = find(lpag_rows(1:lpag_n) == lpag_ref, 1);
    if ~isempty(lpag_ref_pos) && isfinite(lpag_lags(lpag_ref_pos))
        lpag_ref_lags(end + 1, 1) = lpag_lags(lpag_ref_pos); %#ok<SAGROW>
    end
end
clear lpag_traces
if lpag_lag_by_channel.Count == 0
    error('Figure3:NoPhaseLags', ...
        'No mid-LPAG phase lags were read from %s.', lpag_structure_file);
end
fprintf('Read %d channel phase lags from %s\n', lpag_lag_by_channel.Count, lpag_structure_file);
if ~isempty(lpag_ref_lags)
    fprintf('Mid-LPAG reference |phase lag| max = %.3g rad (expected near 0).\n', ...
        max(abs(lpag_ref_lags)));
end

lpag_n_sessions = numel(phase_struct);
lpag_counts = zeros(lpag_n_sessions, 1);
for lpag_ns = 1:lpag_n_sessions
    if ~isfield(phase_struct(lpag_ns), 'cluster_info') || ~isfield(phase_struct(lpag_ns), 'entire_recording_phase_stats')
        error('Figure3:PhaseStructFields', ...
            'Session %d is missing cluster_info or entire_recording_phase_stats.', lpag_ns);
    end
    lpag_info_names = phase_struct(lpag_ns).cluster_info.Properties.VariableNames;
    if ~ismember('ch', lpag_info_names) || ~ismember('area', lpag_info_names)
        error('Figure3:NoChannel', ...
            'cluster_info in session %d has no ch or area column.', lpag_ns);
    end
    lpag_counts(lpag_ns) = size(phase_struct(lpag_ns).cluster_info, 1);
end

lpag_n_neurons = sum(lpag_counts);
lpag_session = cell(lpag_n_neurons, 1);
lpag_area = cell(lpag_n_neurons, 1);
lpag_channel = nan(lpag_n_neurons, 1);
lpag_preferred_angle = nan(lpag_n_neurons, 1);
lpag_ppc_p = nan(lpag_n_neurons, 1);
lpag_cursor = 1;
for lpag_ns = 1:lpag_n_sessions
    lpag_n_this = lpag_counts(lpag_ns);
    lpag_idx = lpag_cursor:(lpag_cursor + lpag_n_this - 1);
    lpag_info = phase_struct(lpag_ns).cluster_info;
    lpag_stats = phase_struct(lpag_ns).entire_recording_phase_stats;
    if ndims(lpag_stats) == 3
        lpag_preferred_part = reshape(lpag_stats(1, :, 1), [], 1);
        lpag_ppc_part = reshape(lpag_stats(1, :, 5), [], 1);
    elseif ismatrix(lpag_stats) && size(lpag_stats, 2) >= 5 && size(lpag_stats, 1) == lpag_n_this
        lpag_preferred_part = lpag_stats(:, 1);
        lpag_ppc_part = lpag_stats(:, 5);
    else
        error('Figure3:PhaseStatsShape', ...
            'entire_recording_phase_stats in session %d has an unexpected size.', lpag_ns);
    end
    if numel(lpag_preferred_part) ~= lpag_n_this || numel(lpag_ppc_part) ~= lpag_n_this
        error('Figure3:PhaseStatsCount', ...
            'Session %d phase stats do not match cluster_info (%d neurons).', ...
            lpag_ns, lpag_n_this);
    end

    if iscell(animal_names)
        lpag_session_name = animal_names{lpag_ns, 1};
    elseif isstring(animal_names)
        lpag_session_name = animal_names(lpag_ns, 1);
    else
        lpag_session_name = animal_names(lpag_ns, :);
    end
    if iscell(lpag_session_name)
        lpag_session_name = lpag_session_name{1};
    end
    lpag_session_name = strtrim(char(lpag_session_name));

    lpag_area_column = lpag_info.area;
    lpag_area_part = cell(lpag_n_this, 1);
    if iscell(lpag_area_column)
        for lpag_a = 1:lpag_n_this
            lpag_value = lpag_area_column{lpag_a};
            while iscell(lpag_value) && numel(lpag_value) == 1
                lpag_value = lpag_value{1};
            end
            if isempty(lpag_value)
                lpag_area_part{lpag_a} = '';
            else
                lpag_area_part{lpag_a} = strtrim(char(string(lpag_value)));
            end
        end
    else
        lpag_area_part = cellstr(strtrim(string(lpag_area_column)));
    end

    lpag_session(lpag_idx) = repmat({lpag_session_name}, lpag_n_this, 1);
    lpag_area(lpag_idx) = lpag_area_part;
    lpag_channel(lpag_idx) = double(lpag_info.ch(:));
    lpag_preferred_angle(lpag_idx) = lpag_preferred_part;
    lpag_ppc_p(lpag_idx) = lpag_ppc_part;
    lpag_cursor = lpag_cursor + lpag_n_this;
end

lpag_phase_lag = nan(lpag_n_neurons, 1);
for lpag_n = 1:lpag_n_neurons
    if ~isfinite(lpag_channel(lpag_n))
        continue
    end
    lpag_key = sprintf('%s|%d', lpag_session{lpag_n}, round(lpag_channel(lpag_n)));
    if isKey(lpag_lag_by_channel, lpag_key)
        lpag_phase_lag(lpag_n) = lpag_lag_by_channel(lpag_key);
    end
end
lpag_centered_angle = mod(lpag_preferred_angle - lpag_phase_lag + pi, 2 * pi) - pi;

lpag_has_angle = isfinite(lpag_preferred_angle);
lpag_has_lag = isfinite(lpag_phase_lag) & lpag_has_angle;
lpag_locked = lpag_has_angle & lpag_ppc_p < 0.01;
lpag_unlocked = lpag_has_angle & lpag_ppc_p >= 0.1;
lpag_original_trough = lpag_locked & lpag_has_lag & ...
    (lpag_preferred_angle > pi / 2 | lpag_preferred_angle < -pi / 2);
lpag_original_peak = lpag_locked & lpag_has_lag & ~lpag_original_trough;
lpag_new_trough = lpag_has_lag & isfinite(lpag_centered_angle) & ...
    (lpag_centered_angle > pi / 2 | lpag_centered_angle < -pi / 2);
lpag_new_peak = lpag_has_lag & isfinite(lpag_centered_angle) & ~lpag_new_trough;

fprintf(['LPAG-centered phase (%s): %d neurons, %d locked, %d locked with a ', ...
    'mid-LPAG lag, %d locked without a lag.\n'], ...
    lpag_phase_lag_band, lpag_n_neurons, sum(lpag_locked), ...
    sum(lpag_locked & lpag_has_lag), sum(lpag_locked & ~lpag_has_lag));
if sum(lpag_locked & lpag_has_lag) < 0.5 * sum(lpag_locked)
    warning('Figure3:PartialMatch', ...
        ['Only %d of %d locked neurons matched an LFP row in coherence_pli_struct_2026.mat. ', ...
        'Preferred angles use cluster_info.ch. Phase lags use that same LFP row in the pairwise matrix.'], ...
        sum(lpag_locked & lpag_has_lag), sum(lpag_locked));
end
if ~any(lpag_locked & lpag_has_lag)
    lpag_phase_sessions = unique(lpag_session, 'stable');
    lpag_lag_sessions = unique(lpag_coherence_sessions, 'stable');
    error('Figure3:NoPhaseLagMatch', ...
        ['No locked neuron matched a mid-LPAG phase lag. ', ...
        'Phase-coupling sessions include: %s. Coherence sessions include: %s. ', ...
        'Session folder names and LFP channel rows have to agree.'], ...
        strjoin(lpag_phase_sessions(1:min(5, numel(lpag_phase_sessions))), ' | '), ...
        strjoin(lpag_lag_sessions(1:min(5, numel(lpag_lag_sessions))), ' | '));
end

lpag_centered_neurons = table(lpag_session, lpag_area, lpag_channel, ...
    lpag_preferred_angle, lpag_phase_lag, lpag_centered_angle, lpag_ppc_p, ...
    'VariableNames', {'session', 'area', 'channel', 'preferred_angle', ...
    'phase_lag', 'lpag_centered_angle', 'ppc_p'});

%% 18a Phase distribution, channel-centered and mid-LPAG-centered
lpag_plot_locked = lpag_locked & lpag_has_lag;
lpag_angle_edges = -pi:(pi / 32):pi;
lpag_distribution_fig = figure('Name', 'LPAG-centered preferred phase');
subplot(1, 2, 1)
polarhistogram(lpag_preferred_angle(lpag_plot_locked), lpag_angle_edges)
title('Channel-centered')
subplot(1, 2, 2)
polarhistogram(lpag_centered_angle(lpag_plot_locked), lpag_angle_edges)
title('Mid-LPAG-centered')
lpag_rmax = 0;
for lpag_panel = 1:2
    subplot(1, 2, lpag_panel)
    lpag_this_r = rlim;
    lpag_rmax = max(lpag_rmax, lpag_this_r(2));
end
for lpag_panel = 1:2
    subplot(1, 2, lpag_panel)
    rlim([0 lpag_rmax])
end

%% 18b Two-component mixture on the mid-LPAG-centered angles
lpag_locked_angles = lpag_centered_angle(lpag_plot_locked);
[lpag_mu, lpag_kappa, lpag_w, lpag_LL] = vm2_mixture_EM(lpag_locked_angles + pi, 10000);
lpag_displayed_mu = lpag_mu - pi;
[~, lpag_peak_component] = min(abs(lpag_displayed_mu));
lpag_component_order = [lpag_peak_component, 3 - lpag_peak_component];
[lpag_theta_uni, lpag_kappa_uni] = circ_vmpar(lpag_locked_angles);
lpag_LL1 = sum(log(circ_vmpdf(lpag_locked_angles, lpag_theta_uni, lpag_kappa_uni)));
lpag_LR = 2 * (lpag_LL(end) - lpag_LL1);
fprintf('LPAG-centered mixture LR (two components vs unimodal) = %.4f\n', lpag_LR);

lpag_angle_counts = histcounts(lpag_locked_angles, lpag_angle_edges);
lpag_angle_counts = lpag_angle_counts(:)' / sum(lpag_angle_counts);
lpag_angle_centers = 0.5 * (lpag_angle_edges(1:end - 1) + lpag_angle_edges(2:end));
lpag_angle_centers = lpag_angle_centers(:)';
lpag_mix_centers = lpag_w .* circ_vmpdf(lpag_angle_centers(:), lpag_mu(1) - pi, lpag_kappa(1)) + ...
    (1 - lpag_w) .* circ_vmpdf(lpag_angle_centers(:), lpag_mu(2) - pi, lpag_kappa(2));
lpag_mix_centers = lpag_mix_centers(:)' / sum(lpag_mix_centers);

lpag_mixture_fig = figure('Name', 'LPAG-centered mixture');
subplot(2, 1, 1)
plot([lpag_angle_centers, lpag_angle_centers + 2 * pi], [lpag_angle_counts, lpag_angle_counts], 'k')
hold on
plot([lpag_angle_centers, lpag_angle_centers + 2 * pi], [lpag_mix_centers, lpag_mix_centers], 'b')
xlabel('Phase (rad)')
ylabel('Probability')
legend({'Mid-LPAG-centered phases', 'Two-component von Mises'}, 'Location', 'northeast')
xlim([-pi, 3 * pi])

lpag_mix_edges = lpag_w .* circ_vmpdf(lpag_angle_edges(:), lpag_mu(1) - pi, lpag_kappa(1)) + ...
    (1 - lpag_w) .* circ_vmpdf(lpag_angle_edges(:), lpag_mu(2) - pi, lpag_kappa(2));
lpag_mix_edges = lpag_mix_edges(:)' / sum(lpag_mix_edges);
subplot(2, 1, 2)
polarplot(lpag_angle_edges, [lpag_angle_counts, lpag_angle_counts(1)], 'k')
hold on
polarplot(lpag_angle_edges, lpag_mix_edges, 'r')
lpag_r_lim = rlim;
polarplot([1 1] * lpag_displayed_mu(lpag_component_order(1)), lpag_r_lim, 'r')
polarplot([1 1] * lpag_displayed_mu(lpag_component_order(2)), lpag_r_lim, 'b')
legend({'Phase distribution', 'Mixture', 'Mean phase peak', 'Mean phase trough'})

%% 18c Mean mid-LPAG-centered angle of trough neurons, per structure
lpag_trough_regions = {'DLPAG', 'LPAG', 'VLPAG', 'SupCol'};
lpag_trough_mask = lpag_new_trough & lpag_locked & ismember(lpag_area, lpag_trough_regions);
lpag_trough_angles = lpag_centered_angle(lpag_trough_mask);
lpag_trough_labels = lpag_area(lpag_trough_mask);
lpag_violin_fig = [];
lpag_percentile_fig = [];
lpag_group_names = {};
lpag_group_idx = [];
if numel(lpag_trough_angles) > 0
    [lpag_group_names, ~, lpag_group_idx] = unique(lpag_trough_labels);
end
if numel(lpag_trough_angles) == 0 || numel(lpag_group_names) < 2
    warning('Figure3:TroughViolin', ...
        'Not enough mid-LPAG-centered trough neurons to compare structures.');
else
    circular_violin_plot(lpag_trough_angles, lpag_group_idx, true);
    lpag_tick_order = str2double(xticklabels);
    xticklabels(lpag_group_names(lpag_tick_order))
    title('Mid-LPAG-centered trough neurons')
    lpag_violin_fig = gcf;

    lpag_angular_range = circular_violin_plot_percentiles(lpag_trough_angles, lpag_group_idx, [20 80], true);
    lpag_tick_order = str2double(xticklabels);
    xticklabels(lpag_group_names(lpag_tick_order))
    title('Mid-LPAG-centered trough neurons')
    lpag_percentile_fig = gcf;
    disp(lpag_angular_range)
end

%% 18d Trough / peak / unlocked percentages from the mid-LPAG-centered angles
lpag_area_list = {'SupCol', 'DLPAG', 'LPAG', 'VLPAG', 'DR'};
lpag_area_labels = {'SC', 'DLPAG', 'LPAG', 'VLPAG', 'DR'};
lpag_number_per_area = nan(numel(lpag_area_list), 5);
for lpag_an = 1:numel(lpag_area_list)
    lpag_in_area = strcmp(lpag_area, lpag_area_list{lpag_an}) & lpag_has_lag;
    lpag_area_angles = lpag_centered_angle(lpag_in_area & lpag_locked);
    lpag_area_finite = isfinite(lpag_area_angles);
    lpag_number_per_area(lpag_an, 1) = sum(lpag_area_angles > pi / 2 | lpag_area_angles < -pi / 2);
    lpag_number_per_area(lpag_an, 2) = sum(lpag_area_finite & ~(lpag_area_angles > pi / 2 | lpag_area_angles < -pi / 2));
    lpag_number_per_area(lpag_an, 3) = sum(lpag_in_area & lpag_locked);
    lpag_number_per_area(lpag_an, 4) = sum(lpag_in_area);
    lpag_number_per_area(lpag_an, 5) = sum(lpag_in_area & lpag_unlocked);
end

lpag_class_counts = [lpag_number_per_area(:, 1), lpag_number_per_area(:, 2), ...
    lpag_number_per_area(:, 5), ...
    lpag_number_per_area(:, 4) - lpag_number_per_area(:, 1) - lpag_number_per_area(:, 2) - lpag_number_per_area(:, 5)];
lpag_class_percent = nan(size(lpag_class_counts));
for lpag_an = 1:numel(lpag_area_list)
    if lpag_number_per_area(lpag_an, 4) > 0
        lpag_class_percent(lpag_an, :) = 100 * lpag_class_counts(lpag_an, :) / lpag_number_per_area(lpag_an, 4);
    end
end

lpag_percent_fig = figure('Name', 'LPAG-centered trough and peak percentages');
lpag_bars = bar(lpag_class_percent, 'stacked');
lpag_bars(1).FaceColor = 'b';
lpag_bars(2).FaceColor = 'r';
lpag_bars(3).FaceColor = 'g';
lpag_bars(4).FaceColor = 'k';
xticklabels(lpag_area_labels)
ylim([0 100])
legend({'Trough-locked', 'Peak-locked', 'Unlocked', 'Unclassified'})
ylabel('Percent')
title('Mid-LPAG-centered classification')

lpag_total_counts = sum(lpag_class_counts, 1);
if sum(lpag_total_counts) > 0
    lpag_total_row = lpag_total_counts / sum(lpag_total_counts);
else
    lpag_total_row = nan(1, size(lpag_class_counts, 2));
end
lpag_class_names = {'Trough-locked', 'Peak-locked', 'Unlocked', 'Unclassified'};
lpag_percent_table = figure('Name', 'LPAG-centered class counts');
uitable(lpag_percent_table, ...
    'Data', [lpag_class_counts; lpag_total_row], ...
    'RowName', [lpag_area_labels, {'Total%'}], ...
    'ColumnName', lpag_class_names, ...
    'Position', [20 20 640 170]);

%% 18e Neurons that stayed trough or stayed peak after the re-reference
lpag_n_regions = numel(lpag_area_list);
lpag_ret_trough_n = zeros(lpag_n_regions + 1, 1);
lpag_ret_trough_stay = zeros(lpag_n_regions + 1, 1);
lpag_ret_peak_n = zeros(lpag_n_regions + 1, 1);
lpag_ret_peak_stay = zeros(lpag_n_regions + 1, 1);
for lpag_an = 1:lpag_n_regions
    lpag_in_area = strcmp(lpag_area, lpag_area_list{lpag_an});
    lpag_ret_trough_n(lpag_an) = sum(lpag_original_trough & lpag_in_area);
    lpag_ret_trough_stay(lpag_an) = sum(lpag_original_trough & lpag_new_trough & lpag_in_area);
    lpag_ret_peak_n(lpag_an) = sum(lpag_original_peak & lpag_in_area);
    lpag_ret_peak_stay(lpag_an) = sum(lpag_original_peak & lpag_new_peak & lpag_in_area);
end
lpag_ret_trough_n(end) = sum(lpag_ret_trough_n(1:lpag_n_regions));
lpag_ret_trough_stay(end) = sum(lpag_ret_trough_stay(1:lpag_n_regions));
lpag_ret_peak_n(end) = sum(lpag_ret_peak_n(1:lpag_n_regions));
lpag_ret_peak_stay(end) = sum(lpag_ret_peak_stay(1:lpag_n_regions));

lpag_pct_trough_stay = 100 * lpag_ret_trough_stay ./ lpag_ret_trough_n;
lpag_pct_peak_stay = 100 * lpag_ret_peak_stay ./ lpag_ret_peak_n;
lpag_pct_trough_stay(lpag_ret_trough_n == 0) = NaN;
lpag_pct_peak_stay(lpag_ret_peak_n == 0) = NaN;

lpag_retention = table([lpag_area_labels, {'Total'}]', lpag_ret_trough_n, ...
    lpag_ret_trough_stay, lpag_pct_trough_stay, lpag_ret_peak_n, ...
    lpag_ret_peak_stay, lpag_pct_peak_stay, ...
    'VariableNames', {'region', 'n_trough', 'n_trough_still_trough', ...
    'pct_trough_still_trough', 'n_peak', 'n_peak_still_peak', 'pct_peak_still_peak'});
disp(lpag_retention)
fprintf('%-8s %8s %14s %8s %14s\n', 'Region', 'N trough', '% still trough', 'N peak', '% still peak');
for lpag_an = 1:(lpag_n_regions + 1)
    fprintf('%-8s %8d %14.1f %8d %14.1f\n', lpag_retention.region{lpag_an}, ...
        lpag_retention.n_trough(lpag_an), lpag_retention.pct_trough_still_trough(lpag_an), ...
        lpag_retention.n_peak(lpag_an), lpag_retention.pct_peak_still_peak(lpag_an));
end

lpag_retention_fig = figure('Name', 'Trough and peak labels kept after LPAG centering');
bar([lpag_pct_trough_stay, lpag_pct_peak_stay]);
xticklabels([lpag_area_labels, {'Total'}])
ylim([0 100])
ylabel('Percent that kept the label')
legend({'Trough neurons still trough', 'Peak neurons still peak'}, 'Location', 'northeast')
title('Label retention after mid-LPAG re-referencing')

lpag_retention_table = figure('Name', 'Trough and peak retention counts');
uitable(lpag_retention_table, ...
    'Data', [lpag_ret_trough_n, lpag_ret_trough_stay, lpag_pct_trough_stay, ...
    lpag_ret_peak_n, lpag_ret_peak_stay, lpag_pct_peak_stay], ...
    'RowName', [lpag_area_labels, {'Total'}], ...
    'ColumnName', {'N trough', 'Still trough', '% still trough', ...
    'N peak', 'Still peak', '% still peak'}, ...
    'Position', [20 20 860 180]);

%% 18f Save the mid-LPAG-centered figures and the neuron table
lpag_figures_folder = fullfile(data_root, 'Figure codes', 'Figure 3 Updated');
if exist(lpag_figures_folder, 'dir') ~= 7
    mkdir(lpag_figures_folder);
end
print(lpag_distribution_fig, '-vector', '-dsvg', ...
    fullfile(lpag_figures_folder, 'lpag_centered_phase_distribution.svg'));
print(lpag_mixture_fig, '-vector', '-dsvg', ...
    fullfile(lpag_figures_folder, 'lpag_centered_phase_mixture.svg'));
if ~isempty(lpag_violin_fig)
    print(lpag_violin_fig, '-vector', '-dsvg', ...
        fullfile(lpag_figures_folder, 'lpag_centered_trough_angle_violin.svg'));
    print(lpag_percentile_fig, '-vector', '-dsvg', ...
        fullfile(lpag_figures_folder, 'lpag_centered_trough_angle_percentiles.svg'));
end
print(lpag_percent_fig, '-vector', '-dsvg', ...
    fullfile(lpag_figures_folder, 'lpag_centered_trough_peak_percent.svg'));
print(lpag_retention_fig, '-vector', '-dsvg', ...
    fullfile(lpag_figures_folder, 'lpag_centered_label_retention.svg'));

lpag_saving_folder = fullfile(data_root, 'Analysis results', 'phase locking data');
if exist(lpag_saving_folder, 'dir') ~= 7
    mkdir(lpag_saving_folder);
end
save(fullfile(lpag_saving_folder, 'lpag_centered_preferred_angles.mat'), ...
    'lpag_centered_neurons', 'lpag_retention', 'lpag_phase_lag_band', '-v7.3');

if numel(lpag_trough_angles) > 0 && numel(lpag_group_names) >= 2
    [lpag_ww_p, lpag_ww_table] = circ_wwtest(lpag_trough_angles, lpag_group_idx);
    disp(lpag_ww_p)
    disp(lpag_ww_table)
end



