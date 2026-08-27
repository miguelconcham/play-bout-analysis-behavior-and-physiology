%% Analyze_freq_coupling
% Cross-frequency coupling analysis (delta phase vs gamma amplitude).

this_file = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\CrossFreq coupling.m';
repo_root = fileparts(fileparts(this_file));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');
npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\NPX data\NPX raw data';
area_limit_table = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Area_limits_GoodLooking.xlsx';

freq_range_1 = [.1 5];
% freq_range_2 = [35 90]; % for gamma
freq_range_2 = [6 12]; % for theta
sr = 2500;
filter_order = 2000;

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

%% Load precomputed coupling structures
disp('loading')

load([saving_folder,'\couplig_structure_delta_theta.mat'],'psth_structure');
load([saving_folder,'\animal_names_coupling_delta_theta.mat'],'animal_names');
% load([saving_folder,'\couplig_structure_delta_gamma.mat'],'psth_structure');
% load([saving_folder,'\animal_names_coupling_delta_gamma.mat'],'animal_names');
disp('ready')

%% Merge sessions — phase-align amplitude distributions
nBins = 72;
edges = linspace(-pi, pi, nBins + 1);
edges = (edges(1:end-1) + edges(2:end)) / 2;
mean_angle = nan(numel(psth_structure), 1);

all_mean_norm = [];
all_real_MI = [];
all_shuffled_MI = [];
all_real_r = [];
all_shuffled_r = [];

figure
for ns = 1:numel(psth_structure)
    subplot(5, 3, ns)
    [mode_ang, ~, ~, ~] = circ_mode_chat(psth_structure(ns).centered_distr, 72, 'SigmaBins', 2, 'Plot', false);
    mean_angle(ns) = mode_ang;

    x_wrap = wrapToPi(edges);
    x_shifted = wrapToPi(x_wrap - mode_ang);
    y = zscore(psth_structure(ns).meanAmp_real);
    polarhistogram(psth_structure(ns).centered_distr, 72)
    r_lim = rlim;
    hold on
    polarplot([mode_ang mode_ang], r_lim, 'r')
    title(animal_names{ns, 1})

    xx = [x_wrap - 2*pi, x_wrap, x_wrap + 2*pi];
    yy = [y, y, y];
    y_shifted = interp1(xx, yy, x_shifted, 'linear');
    all_mean_norm = [all_mean_norm; y_shifted];

    this_session_shuffled_MI = psth_structure(ns).MI_perm;
    this_session_real_MI = psth_structure(ns).MI_real;
    this_session_real_MI = (this_session_real_MI - mean(this_session_shuffled_MI)) / std(this_session_shuffled_MI);
    this_session_shuffled_MI = zscore(this_session_shuffled_MI);
    all_shuffled_MI = [all_shuffled_MI; this_session_shuffled_MI];
    all_real_MI = [all_real_MI; this_session_real_MI];

    this_session_shuffled_r = psth_structure(ns).r_perm;
    this_session_real_r = psth_structure(ns).r_real;
    this_session_real_r = (this_session_real_r - mean(this_session_shuffled_r)) / std(this_session_shuffled_r);
    this_session_shuffled_r = zscore(this_session_shuffled_r);
    all_real_r = [all_real_r; this_session_real_r];
    all_shuffled_r = [all_shuffled_r; this_session_shuffled_r];
end

%% Summary figure — phase-amplitude map, MI and R2 null distributions
y_lim = [-2 2];
session2remove = {'B1D1 1007 Dual', 1};
electrode_index = ~((ismember([animal_names{:, 2}], session2remove{2})' & ismember(animal_names(:, 1), session2remove{1})) | any(isnan(all_mean_norm), 2));

edges_plot = linspace(-pi, pi, nBins);
zero_deg_mod = mean(all_mean_norm(:, edges_plot >= -.5 & edges_plot <= .5), 2);
[~, mod_order] = sort(zero_deg_mod);
electrode_index = electrode_index(mod_order);

figure_figure = figure('units', 'normalized', 'outerposition', [0 .5 1 .5]);

subplot(1, 5, 1)
imagesc(180 * edges_plot / pi, 1:sum(electrode_index), all_mean_norm(mod_order(electrode_index), :))
axis xy
xticks([-180 -90 0 90 180])
yticks(1:sum(electrode_index))
yticklabels(animal_names(mod_order(electrode_index), 1))
colorbar('location', 'northoutside')
clim(y_lim)
xlabel('Delta phase (deg)')

subplot(1, 5, 2)
[~, ~, ci] = ttest(all_mean_norm(mod_order(electrode_index), :));
no_nan = ~any(isnan(ci));
fill(180 * [edges_plot, fliplr(edges_plot)] / pi, [ci(1, no_nan) fliplr(ci(2, no_nan))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
hold on
plot(180 * edges_plot / pi, all_mean_norm(mod_order(electrode_index), :), 'k:')
plot(180 * edges_plot / pi, mean(all_mean_norm(mod_order(electrode_index), :), 'omitmissing'), 'k', 'LineWidth', 2)
ylim(y_lim)
ylabel('Theta Power')
xlabel('Delta Phase (deg)')

subplot(1, 5, 3)
histogram(all_shuffled_MI, 100, 'EdgeColor', 'none', 'FaceColor', 'k', 'Normalization', 'percentage')
xscale log
ylim([0 2.5])
ylabel('%')
yyaxis right
hold on
histogram(all_real_MI, numel(all_real_MI), 'EdgeColor', 'none', 'FaceColor', 'r')
xticks([1 10 100 1000 10000])
ylim([0 2.5])
xlabel({'Modulation Index', '(z-scored)'})
ylabel('Count')

%% Single-session correlation example (LPAG channel)
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];
fn = 4;
current_dir = fullfile(npx_Raw_Data, animal_list(fn).name);

animal_code = strsplit(current_dir, '\');
animal_code = animal_code{end};
animal_code_params = strsplit(animal_code, ' ');
animal_batch = animal_code_params{1};
repeated_animal = animal_code_params{3};

disp('LOADING LFP')
if exist(fullfile(current_dir, 'LFP_PAG.mat'), 'file') == 2
    NPX_Type = 2;
    load(fullfile(current_dir, 'LFP_PAG.mat'), 'LFP')
elseif exist(fullfile(current_dir, 'LFP_PAG.dat'), 'file') == 2
    NPX_Type = 1;
    file_pointer = fopen(fullfile(current_dir, 'LFP_PAG.dat'), 'r');
    LFP = fread(file_pointer, 'int16');
    LFP = reshape(LFP, 384, numel(LFP) / 384);
end
disp('LFP LOADED')

hard_coded_x_coords = [8 40; 258 290; 508 540; 758 790];
area_limit = readtable(area_limit_table);
if strcmp(repeated_animal, 'Single2')
    this_animal = ['Batch', animal_batch(2), repeated_animal];
else
    this_animal = ['Batch', animal_batch(2), repeated_animal, animal_batch(4)];
end
area_limit = area_limit(ismember(area_limit.AnimalName, this_animal), :);

if NPX_Type == 1
    PAG_channels = area_limit{ismember(area_limit.area, {'LPAG'}), {'ch_start', 'ch_end'}};
    PAG_channels = str2double(PAG_channels);
    channel_Range = [min(PAG_channels(:)) max(PAG_channels(:))];
    mid_PAG_channel = round(mean(channel_Range));
else
    load(fullfile(current_dir, 'chann_map_PAG.mat'), 'xcoords', 'ycoords', 'chanMap')
    Y_Range = area_limit{ismember(area_limit.area, {'LPAG'}), {'ProbeNum', 'depth_start', 'depth_end'}};
    mid_PAG_channel = nan(size(Y_Range, 1), 1);
    for j = 1:size(Y_Range, 1)
        this_indexes = ycoords >= Y_Range(j, 2) & ycoords <= Y_Range(j, 3) & ismember(xcoords, hard_coded_x_coords(Y_Range(j, 1), :));
        all_locs = [xcoords(this_indexes) ycoords(this_indexes)];
        mean_loc = mean(all_locs);
        [~, closest_channel] = min(sum(abs([xcoords ycoords] - repmat(mean_loc, numel(ycoords), 1)), 2));
        mid_PAG_channel(j) = chanMap(closest_channel);
    end
end

PAG_LFP = double(LFP(mid_PAG_channel, :));
clear LFP
filtered_signal_freq1 = filtfilt(Hd_freq1.Coefficients, 1, PAG_LFP(1, :));
hiblert_data1 = hilbert(filtered_signal_freq1);
filtered_signal_freq2 = filtfilt(Hd_freq2.Coefficients, 1, PAG_LFP(1, :));
hiblert_data2 = hilbert(filtered_signal_freq2);
amplitud_data1 = abs(hiblert_data1);
amplitud_data2 = abs(hiblert_data2);

figure(figure_figure)
max_n_points = 200;
max_delta = 1600;
n_bins = 100;
[counts, ~, bn_index] = histcounts(amplitud_data1, linspace(min(amplitud_data1), max_delta, n_bins + 1));
valsperbin = min(round(.5 * min(counts)), max_n_points);

subsampled_data = nan(valsperbin * n_bins, 1);
for j = 1:n_bins
    subsampled_data((1 + valsperbin * (j - 1)):(valsperbin * j)) = randsample(find(bn_index == j), valsperbin);
end

data2corr = [amplitud_data1(subsampled_data)', amplitud_data2(subsampled_data)'];
central_values = abs(zscore(data2corr(:, 1))) < 2 & abs(zscore(data2corr(:, 2))) < 2;

subplot(1, 5, 4)
plot(data2corr(central_values, 1), data2corr(central_values, 2), 'k.', 'MarkerSize', .1)
[c, p] = corr(data2corr(central_values, :));
ylabel('Theta Power')
xlabel('Delta Power')
title(p(1, 2))

subplot(1, 5, 5)
histogram(all_shuffled_r, 100, 'EdgeColor', 'none', 'FaceColor', 'k', 'Normalization', 'percentage')
xscale log
hold on
ylim([0 4])
ylabel('%')
yyaxis right
histogram(all_real_r, numel(all_real_r), 'EdgeColor', 'none', 'FaceColor', 'r')
xticks([0.1 1 10 100 1000 10000])
xlim tight
ylabel('count')
xlabel({'R2', '(z-scored)'})
ylim([0 4])
