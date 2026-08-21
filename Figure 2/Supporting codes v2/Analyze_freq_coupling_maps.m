%% Analyze_freq_coupling_maps
% Spatial maps of cross-frequency coupling aligned by PAG sub-area.

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';

%% Load precomputed coupling map structures
disp('loading')
load([saving_folder,'\couplig_structure_V2.mat'],'psth_structure');
load([saving_folder,'\animal_names_coupling_V2.mat'],'animal_names');
disp('ready')

%% Per-session phase-amplitude maps on probe
PROBES = {'NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX1','NPX2','NPX2'};
allResults = struct([]);
nbins = 20;

figure
for ns = 1:numel(psth_structure)
    subplot(4, 3, ns)
    y_pos = psth_structure(ns).channel_map(:, 2);
    y_pos(1) = [];
    [~, order] = sort(y_pos);
    mean_amp = psth_structure(ns).meanAmp_real;

    edges = psth_structure(ns).edges;
    centers = .5 * (edges(1:end-1) + edges(2:end));

    for ch = 1:384
        mean_amp(ch, :) = (mean_amp(ch, :) - mean(mean_amp(ch, :), 'omitmissing')) / std(mean_amp(ch, :), 'omitmissing');
    end

    centers_this_probe = movmean((mean(mean_amp(:, 28:44), 2, 'omitmissing')), 10)';
    ARRAY1 = psth_structure(ns).channel_map(2:end, :);
    CELLARRAY2 = psth_structure(ns).areas_by_channel(2:end);
    ARRAY3 = zscore(centers_this_probe);
    probeType = PROBES{ns};

    [wVals, wAreas, wY, borders] = wrap_probe_values(ARRAY1, CELLARRAY2, ARRAY3, nbins, probeType);
    allResults(ns).wrappedVals = wVals;
    allResults(ns).wrappedAreas = wAreas;
    allResults(ns).wrappedY = wY;
    allResults(ns).borders = borders;

    imagesc(centers, 1:384, mean_amp(order, :))
    axis xy
    clim([-2 2])
    title(animal_names{ns, 1})
end

%% Align and resample by LPAG / VLPAG / DLPAG
[alignedVals, alignedAreas, alignedY] = align_by_area(allResults, {'LPAG', 'VLPAG', 'DLPAG'});
[valsResampled, yResampled] = resample_segments_simple(alignedVals, alignedY, nbins);

figure
hold on
stacked_val = [];
for i = 1:numel(alignedVals)
    stacked_val = [stacked_val; valsResampled{i}];
    plot(yResampled{i}, valsResampled{i}, ':k', 'DisplayName', sprintf('Exp %d', i));
end
xline(0, 'k--', 'Target Area');
xlabel('Relative Wrapped Position');
ylabel('Value');
legend show

%% Group mean with confidence band
figure
hold on
[~, ~, ci] = ttest(stacked_val);
no_nan = ~any(isnan(ci));
fill([yResampled{1}(no_nan) fliplr(yResampled{1}(no_nan))], [ci(1, no_nan) fliplr(ci(2, no_nan))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
plot(yResampled{1}, mean(stacked_val, 'omitmissing'), 'k')

%% Circular mean across sessions
figure
hold on
plot(yResampled{1}, mod(circ_mean(stacked_val) / pi, 1), 'k')
