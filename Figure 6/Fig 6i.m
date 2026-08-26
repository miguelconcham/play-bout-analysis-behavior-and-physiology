
%% LOAD DATA
figure_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Figure codes\Figure Mutual informatin inputs';

saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\Theta psth';

disp('loading')
%%%%%%%%%%%%%%%% The one below contains the frequency ranges as onf Fig 6i
load([saving_folder,'\mi_structure_20bins_delta_by_frequency_extended.mat'],'mi_structure');
%%%%%%%%%%%%%%%% The one blow contains frequency ranges of 1 hz
% load([saving_folder,'\mi_structure_20bins_delta_by_frequency_full_extended.mat'],'mi_structure');

%% Label each frequency band manualy (or automatically using last commented lines)

nAnimals = numel(mi_structure);

% Extract frequency labels
% Plot MI z-scored values across frequency bands
% Uses actual frequency-band centers on x-axis

nAnimals = numel(mi_structure);

% Extract frequency ranges

freq_ranges = mi_structure(1).freq_pow_range;

freq_band_names = { ...
    'Delta (1-4 Hz)', ...
    'Low Theta (5-6 Hz)', ...
    'Theta (7-12 Hz)', ...
    'Beta (12-20 Hz)', ...
    'Low Gamma (20-40 Hz)'};


freq_band_names = { ...
    'Delta (1-4 Hz)', ...
    '(4-6 Hz)', ...
    '(5-6 Hz)', ...
    '(6 - 12 Hz)', ...
    '(7 - 12 Hz)',...
    '(12 - 20 Hz)', ...
    '(20 - 40 Hz)'};

% freq_band_names = arrayfun(@(i) sprintf('%d %d', freq_ranges(i,1), freq_ranges(i,2)), ...
%                        (1:size(freq_ranges,1))', ...
%                        'UniformOutput', false);
% % center frequency of each band;
freq_centers = mean(mi_structure(1).freq_pow_range,2);


%% Build matrix for ploting mi

% rows = animals/sessions
% cols = frequency bands
MI_zscore = nan(nAnimals, numel(freq_centers));
p_values = nan(nAnimals, numel(freq_centers));

for i = 1:nAnimals

    MI_zscore(i,:) = (mi_structure(i).mi_global_per_freq(:,2));
    p_values(i,:) = mi_structure(i).mi_global_per_freq(:,3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Population stats %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_mi = mean(MI_zscore,1,'omitnan');

 var_meas = std(MI_zscore,[],1,'omitnan') ./ ...
          sqrt(sum(~isnan(MI_zscore),1));

tcrit = tinv(0.975, size(MI_zscore,1)-1);
CI = tcrit * var_meas;
var_meas = CI;


% var_meas= std(MI_zscore,[],1,'omitnan') ;

%% Fig 6I  (p values when the mi is significantly above zero)
freq_centers = mean(freq_ranges,2);

figure
hold on
% freq_centers(3) = freq_centers(3)+.25;
[~,sorted_order] = sort(freq_centers);

% individual animals
for i = 1:nAnimals

    plot(freq_centers(sorted_order), ...
         MI_zscore(i,:), ...
         '-o', ...
         'Color', [.7 .7 .7], ...
         'LineWidth',1)

end

% population mean
errorbar(freq_centers(sorted_order), ...
         mean_mi, ...
         var_meas, ...
         '-ko', ...
         'LineWidth',3, ...
         'MarkerFaceColor','k', ...
         'MarkerSize',8)

% Formatting

xlabel('Frequency (Hz)')
ylabel('MI z-score')

title('Population synchrony spectrum')

xlim([min(freq_centers)-1 max(freq_centers)+1])
xticks(freq_centers(sorted_order))
xticklabels(freq_band_names(sorted_order))

box off
set(gca,'FontSize',14)

y_lim = ylim;

for col = 1:numel(freq_centers)
    if signrank(MI_zscore(:,col))<0.05
        text(freq_centers(sorted_order(col)), y_lim(2)+1, '*')
    end
end
ylim([y_lim(1) y_lim(2)+2])
plot([0 15], [0 0], ':k')

plot([0 15], [2 2], ':k')


%% Fig 6I: Paired comparison against first frequency band (need to be added manually to the figure)
% compares every band to band #1 across animals/sessions

close all

nAnimals = numel(mi_structure);

% Frequency info

freq_ranges  = mi_structure(1).freq_pow_range;
freq_centers = mean(freq_ranges,2);

%%%%%%%%%%%%%%%%%%%%%
%%%% Build matrix %%%
%%%%%%%%%%%%%%%%%%%%%

MI_zscore = nan(nAnimals, numel(freq_centers));

for i = 1:nAnimals

    MI_zscore(i,:) = mi_structure(i).mi_global_per_freq(:,2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Reference group = first frequency band %%%%%%
%%%%% Select other column if needed          %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


reference_band = MI_zscore(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compuite paired statistics %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p_values = nan(numel(freq_centers)-1,1);
t_stats  = nan(numel(freq_centers)-1,1);

disp(' ')
disp('Paired t-tests vs first frequency band')
disp('--------------------------------------')

for f = 2:numel(freq_centers)

    comparison_band = MI_zscore(:,f);

    [~,p,~,stats] = ttest(comparison_band- reference_band);

    p_values(f-1) = p;
    t_stats(f-1)  = stats.tstat;

    fprintf('%s vs %s --> t(%d)=%.3f, p=%.5f\n', ...
        freq_band_names{f}, ...
        freq_band_names{1}, ...
        stats.df, ...
        stats.tstat, ...
        p);

end
