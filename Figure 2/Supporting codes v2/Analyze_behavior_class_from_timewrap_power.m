%% Analyze_behavior_class_from_timewrap_power
% Can band power tell two Type2 behavior groups apart?
% Uses the time-wrapped PSTHs already saved by
% Estimate_psth_all_behaviors_theta_gamma.m. It does not reload LFP.
%
% Each row of classification_table is one behavior event from the repeated
% animal, electrode 1. Column 1, mean_power, is the mean log10 band power
% over the time-wrapped behavior. Column 2, behavior_class, is 1 for
% behavior_class_1 and 2 for behavior_class_2. behavior_name and
% session_name are kept so each row can be traced back to a Type2 event.
%
% classifier_name:
%   'lda'          fitcdiscr, the default
%   'svm'          fitcsvm, linear kernel
%   'naive_bayes'  fitcnb
%   'tree'         fitctree
%   'logit'        fitglm, binomial
%
% behavior_class_2 is the positive class. A cross-validated call is then
% true positive, true negative, false positive, or false negative.
% The last figure compares those four counts for delta, theta, and gamma.

%% Choose the band, the two behavior groups, and the classifier
band_name = 'gamma';
% 'delta', 'theta', or 'gamma'

behavior_class_1 = {'Pounce'};
behavior_class_2 = {'PounceI'};
% Type2 names. Several names in one cell share that class.
% Valid names: Pounce, CC, Boxing, Evasion, Pin, Escape, CB, CD,
% Grooming, Rearing, Sniffing, Scratching, PounceI, Bite.

classifier_name = 'lda';
% 'lda', 'svm', 'naive_bayes', 'tree', 'logit'

% Time-wrapped axis. (0, 5] is the behavior itself. Pre-onset is (-5, 0].
power_time_window = [0 5];
n_folds_requested = 10;
n_permutations = 200;
rng_seed = 1;

%% Paths
this_file = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Figure 2\Supporting codes v2\Analyze_behavior_class_from_timewrap_power.m';
repo_root = fileparts(fileparts(fileparts(this_file)));
data_root = fullfile(repo_root, 'Data');
saving_folder = fullfile(data_root, 'Analysis results', 'psth power by frequency and behavior');

all_behaviors = {'Pounce', 'CC', 'Boxing', 'Evasion', 'Pin', 'Escape', 'CB', 'CD', 'Grooming', 'Rearing', 'Sniffing', 'Scratching', 'PounceI', 'Bite'};
animal_label = {'B1D1', 'B1S3', 'B2S2', 'B3D2', 'B4S2', 'B4D4'};
allowed_classifiers = {'lda', 'svm', 'naive_bayes', 'tree', 'logit'};

behavior_class_1 = cellstr(behavior_class_1);
behavior_class_2 = cellstr(behavior_class_2);
classifier_name = lower(classifier_name);
if ~any(strcmp(classifier_name, allowed_classifiers))
    error('classifier_name must be one of: %s.', strjoin(allowed_classifiers, ', '));
end
if any(ismember(behavior_class_1, behavior_class_2))
    error('A behavior is listed in both behavior_class_1 and behavior_class_2.');
end
unknown_behaviors = setdiff([behavior_class_1, behavior_class_2], all_behaviors);
if ~isempty(unknown_behaviors)
    error('Unknown Type2 name: %s. Valid names: %s.', ...
        strjoin(unknown_behaviors, ', '), strjoin(all_behaviors, ', '));
end

%% Mean power of each behavior event
classification_table = behavior_event_table(saving_folder, band_name, all_behaviors, ...
    animal_label, behavior_class_1, behavior_class_2, power_time_window);

n_class_1 = nnz(classification_table.behavior_class == 1);
n_class_2 = nnz(classification_table.behavior_class == 2);
if n_class_1 < 2 || n_class_2 < 2
    error('Need at least 2 electrode-1 events in each class. Class 1 has %d, class 2 has %d.', ...
        n_class_1, n_class_2);
end
n_folds = min([n_folds_requested, n_class_1, n_class_2]);

%% Cross-validated classification
rng(rng_seed);
[cv_prediction, cv_scored] = cross_validated_predictions( ...
    classification_table.mean_power, classification_table.behavior_class, ...
    classifier_name, n_folds);
if ~any(cv_scored)
    error('No cross-validation fold was fit.');
end

scored_class = classification_table.behavior_class(cv_scored);
scored_prediction = cv_prediction(cv_scored);
cv_accuracy = mean(scored_prediction == scored_class);
cv_balanced_accuracy = balanced_accuracy(scored_class, scored_prediction);
chance_accuracy = max(mean(scored_class == 1), mean(scored_class == 2));
confusion_counts = confusionmat(scored_class, scored_prediction, 'Order', [1 2]);
% Class 2 is the positive class. Order is TN, FP, FN, TP.
outcome_count = confusion_outcome_counts(scored_class, scored_prediction);
true_negative = outcome_count(1);
false_positive = outcome_count(2);
false_negative = outcome_count(3);
true_positive = outcome_count(4);
class_1_label = strjoin(behavior_class_1, ', ');
class_2_label = strjoin(behavior_class_2, ', ');
outcome_name = {'true_negative'; 'false_positive'; 'false_negative'; 'true_positive'};
outcome_meaning = {
    sprintf('%s predicted as %s', class_1_label, class_1_label)
    sprintf('%s predicted as %s', class_1_label, class_2_label)
    sprintf('%s predicted as %s', class_2_label, class_1_label)
    sprintf('%s predicted as %s', class_2_label, class_2_label)
    };
confusion_outcomes = table(outcome_name, outcome_count(:), outcome_meaning, ...
    'VariableNames', {'outcome', 'count', 'meaning'});

permutation_balanced_accuracy = nan(n_permutations, 1);
for permutation_n = 1:n_permutations
    shuffled_class = classification_table.behavior_class(randperm(height(classification_table)));
    [shuffled_prediction, shuffled_scored] = cross_validated_predictions( ...
        classification_table.mean_power, shuffled_class, classifier_name, n_folds);
    permutation_balanced_accuracy(permutation_n) = balanced_accuracy( ...
        shuffled_class(shuffled_scored), shuffled_prediction(shuffled_scored));
end
if n_permutations > 0
    permutation_p = (nnz(permutation_balanced_accuracy >= cv_balanced_accuracy) + 1) / ...
        (n_permutations + 1);
else
    permutation_p = NaN;
end

[session_prediction, session_scored] = leave_one_session_out_predictions( ...
    classification_table.mean_power, classification_table.behavior_class, ...
    classification_table.session_name, classifier_name);
if any(session_scored)
    session_balanced_accuracy = balanced_accuracy( ...
        classification_table.behavior_class(session_scored), ...
        session_prediction(session_scored));
else
    session_balanced_accuracy = NaN;
end

classification_result = struct();
classification_result.band_name = band_name;
classification_result.behavior_class_1 = behavior_class_1;
classification_result.behavior_class_2 = behavior_class_2;
classification_result.classifier_name = classifier_name;
classification_result.power_time_window = power_time_window;
classification_result.n_class_1 = n_class_1;
classification_result.n_class_2 = n_class_2;
classification_result.n_folds = n_folds;
classification_result.cv_accuracy = cv_accuracy;
classification_result.cv_balanced_accuracy = cv_balanced_accuracy;
classification_result.chance_accuracy = chance_accuracy;
classification_result.confusion_counts = confusion_counts;
classification_result.confusion_outcomes = confusion_outcomes;
classification_result.true_negative = true_negative;
classification_result.false_positive = false_positive;
classification_result.false_negative = false_negative;
classification_result.true_positive = true_positive;
classification_result.permutation_balanced_accuracy = permutation_balanced_accuracy;
classification_result.permutation_p = permutation_p;
classification_result.session_balanced_accuracy = session_balanced_accuracy;
classification_result.cv_prediction = cv_prediction;

fprintf('\n%s power, %s\n', band_name, classifier_name)
fprintf('Class 1 (%s): %d events\n', strjoin(behavior_class_1, ', '), n_class_1)
fprintf('Class 2 (%s): %d events\n', strjoin(behavior_class_2, ', '), n_class_2)
fprintf('Stratified %d-fold accuracy: %.3f\n', n_folds, cv_accuracy)
fprintf('Stratified %d-fold balanced accuracy: %.3f\n', n_folds, cv_balanced_accuracy)
fprintf('Majority-class chance: %.3f\n', chance_accuracy)
fprintf('Permutation p of balanced accuracy: %.4f (%d shuffles)\n', permutation_p, n_permutations)
fprintf('Leave-one-session-out balanced accuracy: %.3f\n', session_balanced_accuracy)
fprintf('Positive class is %s.\n', class_2_label)
disp(confusion_outcomes)
disp(classification_table(1:min(8, height(classification_table)), :))

outcome_figure_name = classification_figure_name(band_name, behavior_class_1, behavior_class_2, ...
    [classifier_name '_outcomes']);
figure
set(gcf, 'WindowState', 'maximized', 'NumberTitle', 'off', 'Name', outcome_figure_name)
outcome_bars = bar(outcome_count);
outcome_bars.FaceColor = 'flat';
outcome_bars.CData = [0.20 0.45 0.70; 0.85 0.35 0.20; 0.85 0.35 0.20; 0.20 0.45 0.70];
xticklabels({'True negative', 'False positive', 'False negative', 'True positive'})
ylabel('Cross-validated events')
title(sprintf('%s, %s. Positive class: %s', band_name, classifier_name, class_2_label))
ylim([0, max(outcome_count) * 1.15 + 1])
text(1:4, outcome_count, compose('%d', outcome_count), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
set(gca, 'TickDir', 'out')

%% Plot power and the cross-validated confusion matrix
true_label = repmat({class_1_label}, height(classification_table), 1);
true_label(classification_table.behavior_class == 2) = {class_2_label};
predicted_label = repmat({''}, height(classification_table), 1);
predicted_label(cv_scored & cv_prediction == 1) = {class_1_label};
predicted_label(cv_scored & cv_prediction == 2) = {class_2_label};

figure_name = classification_figure_name(band_name, behavior_class_1, behavior_class_2, classifier_name);
figure
set(gcf, 'WindowState', 'maximized', 'NumberTitle', 'off', 'Name', figure_name)

subplot(1, 2, 1)
boxchart(classification_table.behavior_class, classification_table.mean_power);
hold on
scatter(classification_table.behavior_class + (rand(height(classification_table), 1) - 0.5) * 0.25, ...
    classification_table.mean_power, 12, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({class_1_label, class_2_label})
ylabel('Mean log10 power')
title(sprintf('%s, %s', band_name, classifier_name))
set(gca, 'TickDir', 'out')

subplot(1, 2, 2)
confusion_chart = confusionchart( ...
    categorical(true_label(cv_scored), {class_1_label, class_2_label}), ...
    categorical(predicted_label(cv_scored), {class_1_label, class_2_label}), ...
    'RowSummary', 'row-normalized');
confusion_chart.Title = sprintf('Balanced accuracy %.3f, permutation p %.4f', ...
    cv_balanced_accuracy, permutation_p);

%% Compare confusion outcomes across delta, theta, and gamma
bands_to_compare = {'delta', 'theta', 'gamma'};
band_outcome_counts = nan(numel(bands_to_compare), 4);
for band_n = 1:numel(bands_to_compare)
    this_band = bands_to_compare{band_n};
    if strcmp(this_band, band_name)
        band_outcome_counts(band_n, :) = outcome_count(:).';
        continue
    end
    try
        this_table = behavior_event_table(saving_folder, this_band, all_behaviors, ...
            animal_label, behavior_class_1, behavior_class_2, power_time_window);
    catch band_error
        warning('Could not load %s: %s', this_band, band_error.message);
        continue
    end
    this_n_class_1 = nnz(this_table.behavior_class == 1);
    this_n_class_2 = nnz(this_table.behavior_class == 2);
    if this_n_class_1 < 2 || this_n_class_2 < 2
        warning('%s has too few events for classification.', this_band);
        continue
    end
    this_n_folds = min([n_folds_requested, this_n_class_1, this_n_class_2]);
    rng(rng_seed);
    [this_prediction, this_scored] = cross_validated_predictions( ...
        this_table.mean_power, this_table.behavior_class, classifier_name, this_n_folds);
    if ~any(this_scored)
        warning('No cross-validation fold was fit for %s.', this_band);
        continue
    end
    band_outcome_counts(band_n, :) = confusion_outcome_counts( ...
        this_table.behavior_class(this_scored), this_prediction(this_scored));
end

band_confusion_outcomes = array2table(band_outcome_counts, ...
    'VariableNames', outcome_name', 'RowNames', bands_to_compare);
disp(band_confusion_outcomes)

comparison_figure_name = classification_figure_name('delta_theta_gamma', ...
    behavior_class_1, behavior_class_2, [classifier_name '_outcomes']);
figure
set(gcf, 'WindowState', 'maximized', 'NumberTitle', 'off', 'Name', comparison_figure_name)
comparison_bars = bar(band_outcome_counts.');
comparison_bars(1).FaceColor = [0.35 0.35 0.35];
comparison_bars(2).FaceColor = [0.20 0.55 0.35];
comparison_bars(3).FaceColor = [0.75 0.35 0.55];
xticklabels({'True negative', 'False positive', 'False negative', 'True positive'})
ylabel('Cross-validated events')
legend(bands_to_compare, 'Location', 'best')
title(sprintf('%s. Positive class: %s', classifier_name, class_2_label))
set(gca, 'TickDir', 'out')

function [prediction, scored] = cross_validated_predictions(power, class_label, classifier_name, n_folds)
prediction = nan(size(class_label));
partition = cvpartition(class_label, 'KFold', n_folds, 'Stratify', true);
for fold_n = 1:n_folds
    train_index = training(partition, fold_n);
    test_index = test(partition, fold_n);
    try
        model = train_behavior_classifier(power(train_index), class_label(train_index), classifier_name);
        prediction(test_index) = predict_behavior_class( ...
            model, power(test_index), classifier_name);
    catch fold_error
        warning('Fold %d was not fit: %s', fold_n, fold_error.message);
    end
end
scored = isfinite(prediction);
end

function [prediction, scored] = leave_one_session_out_predictions(power, class_label, session_label, classifier_name)
prediction = nan(size(class_label));
sessions = unique(session_label);
for session_n = 1:numel(sessions)
    test_index = strcmp(session_label, sessions{session_n});
    train_index = ~test_index;
    if numel(unique(class_label(train_index))) < 2
        continue
    end
    try
        model = train_behavior_classifier(power(train_index), class_label(train_index), classifier_name);
        prediction(test_index) = predict_behavior_class( ...
            model, power(test_index), classifier_name);
    catch session_error
        warning('Session %s was not fit: %s', sessions{session_n}, session_error.message);
    end
end
scored = isfinite(prediction);
end

function model = train_behavior_classifier(power, class_label, classifier_name)
power = power(:);
class_label = class_label(:);
switch classifier_name
    case 'lda'
        model = fitcdiscr(power, class_label, 'DiscrimType', 'linear');
    case 'svm'
        model = fitcsvm(power, class_label, 'KernelFunction', 'linear', 'Standardize', true);
    case 'naive_bayes'
        model = fitcnb(power, class_label);
    case 'tree'
        model = fitctree(power, class_label);
    case 'logit'
        model = fitglm(power, class_label == 2, 'Distribution', 'binomial');
end
end

function prediction = predict_behavior_class(model, power, classifier_name)
power = power(:);
if strcmp(classifier_name, 'logit')
    probability_class_2 = predict(model, power);
    prediction = ones(size(power));
    prediction(probability_class_2 >= 0.5) = 2;
else
    prediction = predict(model, power);
end
prediction = prediction(:);
end

function value = balanced_accuracy(true_class, predicted_class)
recall_class_1 = mean(predicted_class(true_class == 1) == 1);
recall_class_2 = mean(predicted_class(true_class == 2) == 2);
value = mean([recall_class_1, recall_class_2]);
end

function figure_name = classification_figure_name(band_name, class_1, class_2, classifier_name)
class_1_name = strjoin(class_1, '_');
class_2_name = strjoin(class_2, '_');
figure_name = ['classify_' band_name '_' class_1_name '_vs_' class_2_name '_' classifier_name];
figure_name = regexprep(figure_name, '[^A-Za-z0-9_]', '');
end

function counts = confusion_outcome_counts(true_class, predicted_class)
% Class 2 is the positive class.
% counts is [true_negative, false_positive, false_negative, true_positive].
confusion_counts = confusionmat(true_class, predicted_class, 'Order', [1 2]);
counts = [confusion_counts(1, 1), confusion_counts(1, 2), ...
    confusion_counts(2, 1), confusion_counts(2, 2)];
end

function event_table = behavior_event_table(saving_folder, band_name, all_behaviors, ...
        animal_label, behavior_class_1, behavior_class_2, power_time_window)
structure_file = ['psth_structure_' band_name '_all_behaviors_self_other_partner.mat'];
names_file = ['animal_names_' band_name '_all_behaviors_self_other_partner.mat'];
disp(['loading ' band_name])
load(fullfile(saving_folder, structure_file), 'psth_structure');
load(fullfile(saving_folder, names_file), 'animal_names');
disp('ready')

bin_size = psth_structure(1).wind_length - psth_structure(1).wind_overlap;
wrap_range = psth_structure(1).range_time_wrap;
n_bins_time_wrap = psth_structure(1).n_bins_time_wrap;
pre_index = 1:round(abs(wrap_range(1)) / bin_size);
post_index = 1:round(abs(wrap_range(2)) / bin_size);
time_wrap_time = [(pre_index * bin_size) + wrap_range(1), ...
    linspace(bin_size, wrap_range(2), n_bins_time_wrap), ...
    wrap_range(2) + post_index * bin_size];

requested_behaviors = unique([behavior_class_1, behavior_class_2], 'stable');
requested_index = zeros(1, numel(requested_behaviors));
for name_n = 1:numel(requested_behaviors)
    requested_index(name_n) = find(strcmp(all_behaviors, requested_behaviors{name_n}));
end

psth_timewrap_self = cell(1, numel(all_behaviors));
session_index_self = cell(1, numel(all_behaviors));
electrode_index_self = cell(1, numel(all_behaviors));
session_name_self = cell(1, numel(all_behaviors));
for behavior_n = 1:numel(all_behaviors)
    session_name_self{behavior_n} = cell(0, 1);
end

for j = 1:numel(psth_structure)
    if ~contains(animal_names{j, 1}, animal_label)
        continue
    end
    animal_num = find(cell2mat(cellfun(@(name) contains(animal_names{j, 1}, name), ...
        animal_label, 'UniformOutput', false)));
    electrode_num = animal_names{j, 2};
    this_session_name = animal_label{animal_num(1)};

    for b = requested_index
        this_psth = psth_structure(j).all_behavior_psth_self{3, b};
        n_events = size(this_psth, 1);
        if n_events == 0
            continue
        end
        psth_timewrap_self{b} = [psth_timewrap_self{b}; this_psth];
        session_index_self{b} = [session_index_self{b}; animal_num(1) * ones(n_events, 1)];
        electrode_index_self{b} = [electrode_index_self{b}; electrode_num * ones(n_events, 1)];
        session_name_self{b} = [session_name_self{b}; repmat({this_session_name}, n_events, 1)];
    end
end

mean_power = [];
behavior_class = [];
behavior_name = cell(0, 1);
session_name = cell(0, 1);
for class_id = 1:2
    if class_id == 1
        class_behaviors = behavior_class_1;
    else
        class_behaviors = behavior_class_2;
    end
    for name_n = 1:numel(class_behaviors)
        behavior_index = find(strcmp(all_behaviors, class_behaviors{name_n}));
        psth = psth_timewrap_self{behavior_index};
        electrode = electrode_index_self{behavior_index};
        event_session_name = session_name_self{behavior_index};
        if isempty(psth)
            continue
        end
        n_use = min([size(psth, 1), numel(electrode), numel(event_session_name)]);
        psth = psth(1:n_use, :);
        electrode = electrode(1:n_use);
        event_session_name = event_session_name(1:n_use);
        keep = electrode == 1;
        psth = psth(keep, :);
        event_session_name = event_session_name(keep);
        if isempty(psth)
            continue
        end
        n_time = min(size(psth, 2), numel(time_wrap_time));
        in_window = false(1, size(psth, 2));
        in_window(1:n_time) = time_wrap_time(1:n_time) > power_time_window(1) & ...
            time_wrap_time(1:n_time) <= power_time_window(2);
        if ~any(in_window)
            error('power_time_window [%g %g] contains no time-wrapped samples.', ...
                power_time_window(1), power_time_window(2));
        end
        this_mean = mean(psth(:, in_window), 2, 'omitnan');
        valid = isfinite(this_mean);
        mean_power = [mean_power; this_mean(valid)]; %#ok<AGROW>
        behavior_class = [behavior_class; class_id * ones(nnz(valid), 1)]; %#ok<AGROW>
        behavior_name = [behavior_name; repmat(class_behaviors(name_n), nnz(valid), 1)]; %#ok<AGROW>
        session_name = [session_name; event_session_name(valid)]; %#ok<AGROW>
    end
end

event_table = table(mean_power, behavior_class, behavior_name, session_name);
end
