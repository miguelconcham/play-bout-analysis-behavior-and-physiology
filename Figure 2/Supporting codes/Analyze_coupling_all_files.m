




npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\psth power by frequency and behavior';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];


freq_range_1    = [.1 5];
% freq_range_2    = [6 12];
freq_range_2    = [35 90];
sr              = 2500;
filter_order    = 2000;


Hd_freq1 = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_1(1), ...
'CutoffFrequency2', freq_range_1(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');

Hd_freq2 = designfilt('bandpassfir', ...
'FilterOrder', filter_order, ...
'CutoffFrequency1', freq_range_2(1), ...
'CutoffFrequency2', freq_range_2(2), ...
'SampleRate', sr, ...
'DesignMethod', 'window', ...
'Window', 'hamming');




%% load if needed
disp('loading')
% load([saving_folder,'\couplig_structure_V3.mat'],'psth_structure');
% load([saving_folder,'\animal_names_coupling_V3.mat'],'animal_names');


load([saving_folder,'\couplig_structure_delta_gamma.mat'],'psth_structure');
load([saving_folder,'\animal_names_coupling_delta_gamma.mat'],'animal_names');
disp('ready')



%% merging

all_mean_AMP    = [];
all_mean_norm   = [];
all_real_MI     = [];
all_shuffled_MI = [];
all_real_r     = [];
all_shuffled_r = [];

 nBins = 72; 
 edges = linspace(-pi, pi, nBins+1);
edges = (edges(1:end-1)+edges(2:end))/2;
mean_angle = nan(numel(psth_structure),1);
figure

for ns =1:numel(psth_structure)
    subplot(5,3,ns)
    all_mean_AMP    = [all_mean_AMP;psth_structure(ns).meanAmp_real];
    % all_mean_norm   = [all_mean_norm;(psth_structure(ns).meanAmp_real - psth_structure(ns).mean_ampl)/psth_structure(ns).std_ampl];


    [mode_ang, ~, ~, ~] = circ_mode_chat(psth_structure(ns).centered_distr, 72, 'SigmaBins', 2, 'Plot', false);

    pause(.1)
    mean_angle(ns) = mode_ang;

    a = mode_ang;  % shift in radians
 
    % Wrap to [-pi, pi]

    x_wrap      = wrapToPi(edges);
    x_shifted   = wrapToPi(x_wrap - a);
    y           = zscore(psth_structure(ns).meanAmp_real);
    polarhistogram(psth_structure(ns).centered_distr,72)
    r_lim = rlim;
    hold on
    polarplot([a a], r_lim, 'r')
    title(animal_names{ns,1})
    % Interpolate y onto original x grid (circularly extended)
    % Extend domain to handle wrap-around
    xx = [x_wrap-2*pi, x_wrap, x_wrap+2*pi];
    yy = [y, y, y];

    y_shifted = interp1(xx, yy, x_shifted, 'linear');
    % y_shifted = movmean(y_shifted,2);
    all_mean_norm   = [all_mean_norm;y_shifted];


    this_session_shuffled_MI    = psth_structure(ns).MI_perm;
    this_session_real_MI        = psth_structure(ns).MI_real;
    this_session_real_MI        = (this_session_real_MI - mean(this_session_shuffled_MI))/std(this_session_shuffled_MI);
    this_session_shuffled_MI    = zscore(this_session_shuffled_MI);
    all_shuffled_MI             = [all_shuffled_MI;this_session_shuffled_MI];
    all_real_MI                 = [all_real_MI;this_session_real_MI];




    this_session_shuffled_r     = psth_structure(ns).r_perm;
    this_session_real_r         = psth_structure(ns).r_real;
    this_session_real_r         = (this_session_real_r - mean(this_session_shuffled_r))/std(this_session_shuffled_r);
    this_session_shuffled_r     = zscore(this_session_shuffled_r);
    all_real_r                  = [all_real_r;this_session_real_r];
    all_shuffled_r              = [all_shuffled_r;this_session_shuffled_r];


end

%%


figure
y_lim = [-2 2]
 nBins = 72; 
 % session2remove = {'B4D4 0826 Dual',2};
  session2remove = {'B1D1 1007 Dual',1}; % the detection of delta peak didnt work out (check distribution of delta peas in first subplot)

electrode_index = ~((ismember([animal_names{:,2}],session2remove{2})' &  ismember(animal_names(:,1),session2remove{1})) | any(isnan(all_mean_norm),2));
 edges = linspace(-pi, pi, nBins);
 zero_deg_mod = mean(all_mean_norm(:,edges>= -.5 & edges<=.5 ),2);
 % [~, mod_order] = sort(mean_angle);
 
  [~, mod_order] = sort(zero_deg_mod);
 electrode_index = electrode_index(mod_order);
figure_figure =  figure('units','normalized','outerposition',[0 .5 1 .5]);
subplot(1,5,1)
imagesc(180*edges/pi, 1:sum(electrode_index),all_mean_norm(mod_order(electrode_index),:))
axis xy
xticks([-180 -90 0 90 180])
yticks(1:sum(electrode_index))
yticklabels(animal_names(mod_order(electrode_index),1))
colorbar('location', 	'northoutside')
clim(y_lim)
xlabel('Delta phase (deg)')

subplot(1,5,2)
[~, ~, ci] = ttest(all_mean_norm(mod_order(electrode_index),:));
no_nan = ~any(isnan(ci));
fill(180*[edges, fliplr(edges)]/pi, [ci(1,no_nan) fliplr(ci(2, no_nan))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
hold on
plot(180*edges/pi, all_mean_norm(mod_order(electrode_index),:), 'k:')
plot(180*edges/pi, mean(all_mean_norm(mod_order(electrode_index),:), 'omitmissing'), 'k', 'LineWidth',2)
ylim(y_lim)
ylabel('Theta Power')
xlabel('Delta Phase (deg)')

subplot(1,5,3)
hold off
histogram(all_shuffled_MI, 100, 'EdgeColor','none','FaceColor','k', 'Normalization','percentage')
xscale log
ylim([0 2.5])
ylabel('%')
yyaxis right
hold on
histogram(all_real_MI, numel(all_real_MI), 'EdgeColor','none','FaceColor','r')
xticks([1 10 100 1000 10000])
ylim([0 2.5])
xlabel({'Modualtion Index', '(z-scored)'})
ylabel('Count')




%% load single session for correlation example example

fn =4;
current_dir = [npx_Raw_Data, '\', animal_list(fn).name];
area_limit_table    = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\Area_limits_GoodLooking.xlsx';
% npx_raw_data = 
animal_code         = strsplit(current_dir, '\');
animal_code         = animal_code{end};
animal_code_params  = strsplit(animal_code, ' ');
animal_batch        = animal_code_params{1};
repeated_animal     = animal_code_params{3};

%%  load lfp from current dir
disp('LOADING LFP')
if exist([current_dir,'\','LFP_PAG.mat'], 'file')==2

    NPX_Type        = 2;
    load([current_dir,'\','LFP_PAG.mat'], 'LFP')
elseif exist([current_dir,'\','LFP_PAG.dat'], 'file')==2
    NPX_Type        = 1;
    file_pointer    = fopen([current_dir,'\','LFP_PAG.dat'], 'r');
    LFP             = fread(file_pointer,'int16');
    LFP             = reshape(LFP, 384, numel(LFP)/384);
end


disp('LFP LOADED')
%% select_mid_pag_channel
disp('Loading Channel Map')
hard_coded_x_coords = [8 40;258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
this_animal = ['Batch', animal_batch(2), repeated_animal,animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName,this_animal),:);

if NPX_Type == 1

    PAG_channels = area_limit{ismember(area_limit.area, {'LPAG'}), {'ch_start', 'ch_end'}};
    PAG_channels = str2double(PAG_channels);
    channel_Range = [min(PAG_channels(:)) max(PAG_channels(:))];
    mid_PAG_channel = round(mean(channel_Range));
else
    load([current_dir,'\ChannelMap.mat'], 'xcoords', 'ycoords','chanMap')
    Y_Range = area_limit{ismember(area_limit.area, {'LPAG'}), {'ProbeNum','depth_start', 'depth_end'}};
    mid_PAG_channel = nan(size(Y_Range,1),1);
    figure
    plot(xcoords,ycoords, 'k.')
    hold on


    for j=1:size(Y_Range,1)
        this_indexes = ycoords>=Y_Range(j,2) & ycoords<=Y_Range(j,3) & ismember(xcoords,hard_coded_x_coords(Y_Range(j,1),:));

        all_locs = [xcoords(this_indexes) ycoords(this_indexes)];
        plot(all_locs(:,1),all_locs(:,2), 'r.')

        mean_loc = mean(all_locs);
        [~, closest_channel]= min(sum(abs([xcoords ycoords]-repmat(mean_loc,numel(ycoords),1)),2));
        plot(xcoords(closest_channel), ycoords(closest_channel), 'xb')
        mid_PAG_channel(j) = chanMap(closest_channel);

    end

end
%% obtain LFP AND FILTER DATA



PAG_LFP         = double(LFP(mid_PAG_channel,:));
clear LFP
ch_n=1;
disp('Filtering Signal')

filtered_signal_freq1 = filtfilt(Hd_freq1.Coefficients, 1, PAG_LFP(ch_n,:));
hiblert_data1 = hilbert(filtered_signal_freq1);

filtered_signal_freq2 = filtfilt(Hd_freq2.Coefficients, 1, PAG_LFP(ch_n,:));
hiblert_data2 = hilbert(filtered_signal_freq2);


phase_data1         = angle(hiblert_data1);
amplitud_data1      = abs(hiblert_data1);
amplitud_data2      = abs(hiblert_data2);

%%
figure(figure_figure)
max_n_points = 200;
max_delta   = 1600;
points2plot = 1000;
% index = randsample(numel(amplitud_data1), points2plot);
n_bins = 100;
[counts, ~, bn_index] = histcounts(amplitud_data1,linspace(min(amplitud_data1), max_delta,n_bins+1));

valsperbin = min(round(.5*min(counts)),max_n_points);;

subsampled_data = nan(valsperbin*n_bins,1);

for j=1:n_bins
    subsampled_data((1 + valsperbin*(j-1)):(valsperbin*(j))) = randsample(find(bn_index==j),valsperbin);
end

data2corr = [amplitud_data1(subsampled_data)',amplitud_data2(subsampled_data)'];
plot(data2corr(:,1),data2corr(:,2), '.' )
hold on
fit_lm = fitlm(amplitud_data1(subsampled_data)',amplitud_data2(subsampled_data)');







central_values = abs(zscore(data2corr(:,1)))<2 & abs(zscore(data2corr(:,2)))<2;
subplot(1,5,4)
plot(data2corr(central_values,1),data2corr(central_values,2), 'k.', 'MarkerSize',.1 )
[c,p] = corr(data2corr(central_values,:));
ylabel('Theta Power')
xlabel('Delta Power')
title(p(1,2))


subplot(1,5,5)
hold off
histogram(all_shuffled_r, 100, 'EdgeColor','none','FaceColor','k', 'Normalization','percentage')
xscale log
hold on
ylim([0 4])
ylabel('%')
yyaxis right
histogram(all_real_r, numel(all_real_r), 'EdgeColor','none','FaceColor','r')
xticks([0.1 1 10 100 1000 10000])
xlim tight
ylabel('count')
xlabel({'R2', '(z-scored)'})
ylim([0 4])


%% save resulint figure

figure_dir='\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Codes\Figure codes\Data cros s freuqnecy coupling';
print(gcf,'-vector','-dsvg',[figure_dir,'\gamma delta entreinment.svg'])