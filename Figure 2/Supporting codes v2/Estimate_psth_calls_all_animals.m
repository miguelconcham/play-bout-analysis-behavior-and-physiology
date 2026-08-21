%% Estimate_psth_calls_all_animals
% Driver script: compute call-triggered PSTH across all animals, then analyze/plot.
% Calls GENERATE_PSTH_CALLS for estimation; subsequent sections load saved data
% for merging, plotting, and mixed-model statistics.

%% Animal list and paths
list_of_animals = {'B1D1 1013 Dual','B1S3 1008 Single','B1S3 1009 Single', ...
    'B2S2 1110 Single2','B2S2 1111 Single2','B3D2 1130 Dual','B4S2 0825 Single'};
saving_folder = '\\experimentfs.bccn-berlin.pri\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\DataSets\Analysis results\Theta psth';

%% Spectrogram parameters (gamma — adjust for theta/delta as needed)
f           = 35:1:100;
freq_range  = [35 90];
wind_params = [.250 .245];

%% Estimate call-triggered PSTH
psth_structure = [];
animal_names = [];

for fn = 1:numel(list_of_animals)
    new_struct = GENERATE_PSTH_CALLS(list_of_animals{fn}, f, freq_range, wind_params);

    if isempty(psth_structure)
        psth_structure = new_struct;
    else
        start_idx = numel(psth_structure) + 1;
        for sub_j = 1:numel(new_struct)
            psth_structure(start_idx + sub_j - 1) = new_struct(sub_j);
        end
    end

    animal_names = [animal_names; ...
        [repmat({list_of_animals{fn}}, numel(new_struct), 1), num2cell(1:numel(new_struct))']];
end

%% Save results (uncomment when needed)
% disp('saving')
% save(fullfile(saving_folder, 'psth_structure_call_gamma.mat'), 'psth_structure');
% save(fullfile(saving_folder, 'animal_names_call_gamma.mat'), 'animal_names');

%% Load saved data for analysis (delta band example)
f           = .1:.05:6;
freq_range  = [.1 5];
wind_params = [.250 .245];

load(fullfile(saving_folder, 'psth_structure_call_delta.mat'), 'psth_structure');
load(fullfile(saving_folder, 'animal_names_call_delta.mat'), 'animal_names');

%% Merge data across animals
smooth_wind          = 20;
baseline_range       = [-2 0];
animal_label         = {'B1D1','B1S3','B2S2','B3D2', 'B4S2'};
bin_size             = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
psth_ranges          = psth_structure(1).hist_range;
time                 = psth_ranges(1):bin_size:psth_ranges(2)+bin_size;
baseline_index       = time<baseline_range(2) & time>baseline_range(1);

all_psth_onset       = [];
all_psth_offset      = [];
all_onset_regressors = [];
all_offset_regressors = [];
all_Calls            = [];
electrode_index      = [];
animal_index         = [];

for j = 1:numel(psth_structure)
    if contains(animal_names{j}, animal_label)
        animal_num    = find(cell2mat(cellfun(@(x) contains(animal_names{j}, x), animal_label, 'UniformOutput', false)));
        electrode_num = animal_names{j,2};

        this_psth_onset = psth_structure(j).call_onset;
        animal_index    = [animal_index; repmat(animal_num, size(this_psth_onset,1), 1)];
        electrode_index = [electrode_index; ones(size(this_psth_onset,1), 1)*electrode_num];

        for trial = 1:size(this_psth_onset,1)
            this_psth_onset(trial,:) = (this_psth_onset(trial,:) - mean(this_psth_onset(trial,baseline_index))) / std(this_psth_onset(trial,baseline_index));
            this_psth_onset(trial,:) = movmean(this_psth_onset(trial,:), smooth_wind);
        end
        all_psth_onset = [all_psth_onset; this_psth_onset];

        this_psth_onset  = psth_structure(j).call_onset;
        this_psth_offset = psth_structure(j).call_offset;
        for trial = 1:size(this_psth_offset,1)
            this_psth_offset(trial,:) = (this_psth_offset(trial,:) - mean(this_psth_onset(trial,baseline_index))) / std(this_psth_onset(trial,baseline_index));
            this_psth_offset(trial,:) = movmean(this_psth_offset(trial,:), smooth_wind);
        end
        all_psth_offset = [all_psth_offset; this_psth_offset];

        all_onset_regressors  = [all_onset_regressors; psth_structure(j).call_onset_regressor];
        all_offset_regressors = [all_offset_regressors; psth_structure(j).call_onset_regressor];
        all_Calls             = [all_Calls; psth_structure(j).CallStats];
    end
end

call_lengths = all_Calls.CallLengths;
[sorted_call_lengths, order] = sort(call_lengths);

%% Plot — per-animal call PSTH
X_lim = [-.5 1];
figure
min_length = .0;
respones_array = [];

for an = 1:numel(animal_label)
    animal_bool    = animal_index == an;
    length_bool    = call_lengths > min_length;
    electrode_bool = electrode_index == 1;
    [sorted_call_lengths, order] = sort(call_lengths(animal_bool & length_bool & electrode_bool,:));

    subplot(5, numel(animal_label), (1:numel(animal_label):2*numel(animal_label)) + an-1)
    array = all_psth_onset(animal_bool & length_bool & electrode_bool,:);
    imagesc(time, 1:numel(sorted_call_lengths), array(order,:))
    xlim(X_lim)
    clim([-2 2])
    axis xy
    hold on
    plot([0 0], [1 numel(sorted_call_lengths)], 'w')
    plot(sorted_call_lengths, 1:numel(sorted_call_lengths), 'w')
    title(animal_label{an})

    subplot(5, numel(animal_label), ((2*numel(animal_label) + 1):numel(animal_label):5*numel(animal_label)) + an-1)
    [~, ~, ci] = ttest(array);
    no_nan = ~any(isnan(ci));
    fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
    hold on
    plot(time, trimmean(array, 5), 'k')
    yyaxis right
    plot(time, mean(all_onset_regressors(animal_bool & length_bool,:)), 'r')
    respones_array = [respones_array; trimmean(array, 5)];
    xlim(X_lim)
end

%% Plot — all sessions combined (no mixed model)
figure
X_lim = [-2 2];
min_length = 0.05;
length_bool    = call_lengths > min_length;
electrode_bool = electrode_index == 1;
[sorted_call_lengths, order] = sort(call_lengths(length_bool & electrode_bool,:));

subplot(5,1,1:3)
array = all_psth_onset(length_bool & electrode_bool,:);
imagesc(time, 1:numel(sorted_call_lengths), array(order,:))
xlim(X_lim)
clim([-2 2])
axis xy
hold on
plot([0 0], [1 numel(sorted_call_lengths)], 'w')
plot(sorted_call_lengths, 1:numel(sorted_call_lengths), 'w')

subplot(5,1,4:5)
[~, ~, ci] = ttest(array);
no_nan = ~any(isnan(ci));
fill([time(no_nan) fliplr(time(no_nan))], [ci(1,no_nan) fliplr(ci(2,no_nan))], 'k', 'FaceAlpha', .25, 'EdgeColor', 'none')
hold on
plot(time, mean(array, 'omitmissing'), 'k')
yyaxis right
plot(time, mean(all_onset_regressors(length_bool,:)), 'r')
xlim(X_lim)

%% Mixed-model — per time bin
length_limit  = .0;
indextinclude = call_lengths >= length_limit;
power         = all_psth_onset(indextinclude,:);
subject_idx   = animal_index(indextinclude);
time_range    = [-1 2];
limted_time   = time;
valid_bins    = find(limted_time >= time_range(1) & limted_time <= time_range(2));
power         = power(:, valid_bins);
limted_time   = limted_time(valid_bins);

[nTrials, nTime] = size(power);
est   = nan(nTime, 1);
se    = nan(nTime, 1);
pvals = nan(nTime, 1);
d     = nan(nTime, 1);
ci    = nan(nTime, 2);

[trial_idx, time_idx] = ndgrid(1:nTrials, 1:nTime);
tbl = table;
tbl.Power = power(:);
data2include = abs(tbl.Power) < 2.5;
tbl.Subject = categorical(subject_idx(trial_idx(:)));
time_matrix = repmat(limted_time, nTrials, 1);
tbl.Time = time_matrix(:);
tbl = tbl(data2include,:);

for i = 1:nTime
    tbl_t = tbl(tbl.Time == limted_time(i), :);
    if sum(~isnan(tbl_t.Power)) > 10
        lme = fitlme(tbl_t, 'Power ~ 1 + (1|Subject)');
        est(i)   = lme.Coefficients.Estimate(1);
        se(i)    = lme.Coefficients.SE(1);
        pvals(i) = lme.Coefficients.pValue(1);
        ci_t = coefCI(lme);
        ci(i,:) = ci_t(1,:);
        sd_within = std(tbl_t.Power, 'omitmissing');
        d(i) = est(i) / sd_within;
    end
end

pvals_fdr = mafdr(pvals, 'BHFDR', true);
results.time      = limted_time(:);
results.est       = est;
results.se        = se;
results.ci        = ci;
results.pvals     = pvals;
results.pvals_fdr = pvals_fdr;
results.d         = d;

save(fullfile(saving_folder, 'results_call_updated_gamma.mat'), 'results', 'respones_array', 'time');

%% Plot — mixed-model results
load(fullfile(saving_folder, 'results_call_updated_delta.mat'), 'results');

alpha = 0.05;
forced_x_lim = [-1 2];
limited_time = results.time;
est = results.est;
ci = results.ci;
pvals_fdr = results.pvals_fdr;
d = results.d;

figure;
subplot(2,1,1); hold on;
plot(time, respones_array, 'k:')
no_nan = ~any(isnan(ci'));
fill([limited_time(no_nan); flipud(limited_time(no_nan))], [ci(no_nan,1); flipud(ci(no_nan,2))], ...
    [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(limited_time, est, 'b', 'LineWidth', 2);
sig_idx = pvals_fdr < alpha;
plot(limited_time(sig_idx), est(sig_idx), 'r*', 'MarkerSize', 6);
ylabel('Mean Power');
title('Mixed-Effects Power (CI + FDR-corrected sig)');
grid on;
xlim(forced_x_lim)

subplot(2,1,2);
plot(limited_time, d, 'k', 'LineWidth', 2);
ylabel('Effect size (Cohen''s d)');
xlabel('Time (s)');
grid on;
xlim(forced_x_lim)
