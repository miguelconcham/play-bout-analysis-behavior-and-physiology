


npx_Raw_Data = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\NPX data\NPX raw data';
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';
animal_list = dir(npx_Raw_Data);
animal_list(1:2) = [];


animal_file_names =  cellfun(@(x) ['B', x],strsplit([animal_list.name], 'B'), 'UniformOutput',false)';
animal_file_names(1) = [];

freq_range_1    = [.1 6];
freq_range_2    = [6 12];
sr              = 2500;
filter_order    = 2000;


psth_structure = [];
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

%%
figure
center = [];
center_zs = [];
allResults = struct([]);

PROBES = {'NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX2','NPX2'};
for ns =1:numel(psth_structure)
subplot(4,3,ns)
y_pos = psth_structure(ns).channel_map(:,2);
y_pos(1) = [];
[~, order] = sort(y_pos);
mean_amp=psth_structure(ns).meanAmp_real;

center = [center;movmean((mean(mean_amp(:,28:44),2)),10)'];
eges = psth_structure(ns).edges;
centers = .5*(edges(1:end-1) + edges(2:end));

for ch=1:384
    mean_amp(ch,:) = (mean_amp(ch,:) - mean(mean_amp(ch,:), 'omitmissing'))/std(mean_amp(ch,:), 'omitmissing');
end
centers_this_probe = movmean((mean(mean_amp(:,28:44),2, 'omitmissing')),10)';
ARRAY1 = psth_structure(ns).channel_map(2:end,:);
amp_difference = psth_structure(ns).power_hig_freq;
MI_difference = zscore(psth_structure(ns).MI_real(:,1));
phase_difference = psth_structure(ns).original_distr;
[~,loc] = min(abs(phase_difference-pi));
CELLARRAY2 = psth_structure(ns).areas_by_channel(2:end);
ARRAY3 = zscore(centers_this_probe);
% ARRAY3 = zscore(amp_difference);
% ARRAY3 = MI_difference;
% ARRAY3 = phase_difference+pi;
nbins = 20;
probeType = PROBES{ns};
[wVals, wAreas, wY, borders] = ...
    wrap_probe_values(ARRAY1, CELLARRAY2, ARRAY3, nbins, probeType);
allResults(ns).wrappedVals = wVals;
    allResults(ns).wrappedAreas = wAreas;
    allResults(ns).wrappedY = wY;
    allResults(ns).borders = borders;
center_zs = [center_zs;centers_this_probe];
imagesc(centers, 1:384, mean_amp(order,:))
axis xy
clim([-2 2])
title(animal_names{ns,1})
end

%%
[alignedVals, alignedAreas, alignedY] = align_by_area(allResults, {'LPAG','VLPAG', 'DLPAG'});
[valsResampled, yResampled] = resample_segments_simple(alignedVals, alignedY, nbins);
figure; hold on;
stacked_val = []
for i = 1:numel(alignedVals)
    stacked_val = [stacked_val; valsResampled{i}];
    plot(yResampled{i}, valsResampled{i},':k', 'DisplayName', sprintf('Exp %d', i));
end
xline(0,'k--','Target Area'); % target area centered
xlabel('Relative Wrapped Position');
ylabel('Value');
legend show

%%
figure
% plot(stacked_val', ':k')
hold on
[~,~,ci]=ttest(stacked_val);
no_nan = ~any(isnan(ci));
fill([yResampled{1}(no_nan) fliplr(yResampled{1}(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
plot(yResampled{1},mean(stacked_val, 'omitmissing'), 'k')

%%
figure
% plot(stacked_val', ':k')
hold on
% [~,~,ci]=ttest(stacked_val);
% no_nan = ~any(isnan(ci));
% fill([yResampled{1}(no_nan) fliplr(yResampled{1}(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha',.25, 'EdgeColor','none')
plot(yResampled{1},mod(circ_mean(stacked_val)/pi,1), 'k')