%% Analyze_behavior_class_from_timewrap_power
% Can band power tell two Type2 behavior groups apart?
% Uses the time-wrapped PSTHs already saved by
% Estimate_psth_all_behaviors_theta_gamma.m. It does not reload LFP.
%
% First approach: one model, trained on events pooled from every kept
% animal and every session. Events are not split into a separate model
% per animal or per session.
%
% Each row of classification_table is one behavior event from the repeated
% animal, electrode 1. Column 1, mean_power, is the mean baseline z-score
% over the time-wrapped behavior. Each event is z-scored to its own
% pre-onset baseline before that mean is taken. Column 2, behavior_class,
% is 1 for behavior_class_1 and 2 for behavior_class_2. behavior_name,
% animal_name, and session_name trace each row back to a Type2 event.
%
% classifier_name:
%   'lda'          fitcdiscr, the default
%   'svm'          fitcsvm, linear kernel
%   'naive_bayes'  fitcnb
%   'tree'         fitctree
%   'logit'        fitglm, binomial
%
% Cross-validation holds out single events from that pooled table. Each
% training fold still mixes the remaining sessions and animals, then keeps
% the same number of events from each class. Held-out events keep their
% natural counts. behavior_class_2 is the positive class. Confusion values
% are the percent of those single events, and the four percentages of one
% band sum to 100. A later figure compares delta, theta, and gamma on
% that same pooled scale. One further figure splits the same calls by the
% true behavior: behavior A classified as A, A as B, B as B, and B as A.
% Each behavior's two rates sum to 100. Each of those rates is tested
% against a label-shuffle null. A shuffle keeps the class counts, retrains
% the same cross-validated classifier, and records the rate. p is the
% upper tail of that null: the share of shuffles whose rate is at least
% as large as the observed rate. Leave-one-session-out is printed
% afterward as a later check. It is not the first model.

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

% Time-wrapped axis. (0, 5] is the behavior itself.
% Each event is z-scored to its own samples inside the baseline window.
power_time_window = [0 5];
baseline_time_window = [-5 0];
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
    animal_label, behavior_class_1, behavior_class_2, power_time_window, baseline_time_window);

n_class_1 = nnz(classification_table.behavior_class == 1);
n_class_2 = nnz(classification_table.behavior_class == 2);
if n_class_1 < 2 || n_class_2 < 2
    error('Need at least 2 electrode-1 events in each class. Class 1 has %d, class 2 has %d.', ...
        n_class_1, n_class_2);
end
n_folds = min([n_folds_requested, n_class_1, n_class_2]);
pooled_animals = unique(classification_table.animal_name, 'stable');
pooled_sessions = unique(classification_table.session_name, 'stable');
n_animals_pooled = numel(pooled_animals);
n_sessions_pooled = numel(pooled_sessions);

%% Cross-validated classification
% First approach. classification_table holds every electrode-1 event from
% every kept animal and every session. Each training fold uses that pool
% minus the held-out events, then subsamples to equal class counts.
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
n_single_events = numel(scored_class);
outcome_percent = 100 * outcome_count(:) / n_single_events;
confusion_outcomes = table(outcome_name, outcome_count(:), outcome_percent, outcome_meaning, ...
    'VariableNames', {'outcome', 'count', 'percent_of_single_events', 'meaning'});
% Percent of each true behavior. Order is A as A, A as B, B as A, B as B.
within_behavior_percent = within_behavior_percent_from_counts(outcome_count);

permutation_balanced_accuracy = nan(n_permutations, 1);
% Columns match within_behavior_percent: A as A, A as B, B as A, B as B.
permutation_within_behavior_percent = nan(n_permutations, 4);
for permutation_n = 1:n_permutations
    shuffled_class = classification_table.behavior_class(randperm(height(classification_table)));
    [shuffled_prediction, shuffled_scored] = cross_validated_predictions( ...
        classification_table.mean_power, shuffled_class, classifier_name, n_folds);
    permutation_balanced_accuracy(permutation_n) = balanced_accuracy( ...
        shuffled_class(shuffled_scored), shuffled_prediction(shuffled_scored));
    if ~any(shuffled_scored)
        continue
    end
    permutation_within_behavior_percent(permutation_n, :) = within_behavior_percent_from_counts( ...
        confusion_outcome_counts(shuffled_class(shuffled_scored), shuffled_prediction(shuffled_scored)));
end
if n_permutations > 0
    permutation_p = (nnz(permutation_balanced_accuracy >= cv_balanced_accuracy) + 1) / ...
        (n_permutations + 1);
else
    permutation_p = NaN;
end
within_behavior_null_mean = mean(permutation_within_behavior_percent, 1, 'omitnan');
within_behavior_permutation_p = permutation_p_above( ...
    permutation_within_behavior_percent, within_behavior_percent);

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
classification_result.baseline_time_window = baseline_time_window;
classification_result.pooling = 'all_sessions_all_animals';
classification_result.n_animals = n_animals_pooled;
classification_result.n_sessions = n_sessions_pooled;
classification_result.pooled_animals = pooled_animals;
classification_result.pooled_sessions = pooled_sessions;
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
classification_result.n_single_events = n_single_events;
classification_result.outcome_percent = outcome_percent;
classification_result.within_behavior_percent = within_behavior_percent;
classification_result.permutation_balanced_accuracy = permutation_balanced_accuracy;
classification_result.permutation_p = permutation_p;
classification_result.permutation_within_behavior_percent = permutation_within_behavior_percent;
classification_result.within_behavior_null_mean = within_behavior_null_mean;
classification_result.within_behavior_permutation_p = within_behavior_permutation_p;
classification_result.session_balanced_accuracy = session_balanced_accuracy;
classification_result.cv_prediction = cv_prediction;

fprintf('\nFirst approach: train on events pooled from all sessions and all animals.\n')
fprintf('%s power, %s\n', band_name, classifier_name)
fprintf('Pooled animals (%d): %s\n', n_animals_pooled, strjoin(pooled_animals, ', '))
fprintf('Pooled sessions (%d): %s\n', n_sessions_pooled, strjoin(pooled_sessions, ', '))
fprintf('Class 1 (%s): %d single events\n', strjoin(behavior_class_1, ', '), n_class_1)
fprintf('Class 2 (%s): %d single events\n', strjoin(behavior_class_2, ', '), n_class_2)
fprintf('Held-out single events: %d\n', n_single_events)
fprintf('Stratified %d-fold accuracy on the pooled events: %.1f%% of single events\n', n_folds, 100 * cv_accuracy)
fprintf('Stratified %d-fold balanced accuracy: %.1f%%\n', n_folds, 100 * cv_balanced_accuracy)
fprintf('Majority-class chance: %.1f%%\n', 100 * chance_accuracy)
fprintf('Permutation p of balanced accuracy: %.4f (%d shuffles)\n', permutation_p, n_permutations)
fprintf('Positive class is %s. Each training fold uses the pooled events, minus the held-out events, with equal class counts.\n', class_2_label)
fprintf('Confusion percentages sum to 100 over held-out single events.\n')
fprintf('Within each behavior, the two rates below sum to 100.\n')
fprintf('Each p is the upper tail of the label-shuffle null for that rate (%d shuffles).\n', n_permutations)
fprintf('  %s classified as %s: %.1f%% (null mean %.1f%%, p=%.4f)\n', ...
    class_1_label, class_1_label, within_behavior_percent(1), ...
    within_behavior_null_mean(1), within_behavior_permutation_p(1))
fprintf('  %s classified as %s: %.1f%% (null mean %.1f%%, p=%.4f)\n', ...
    class_1_label, class_2_label, within_behavior_percent(2), ...
    within_behavior_null_mean(2), within_behavior_permutation_p(2))
fprintf('  %s classified as %s: %.1f%% (null mean %.1f%%, p=%.4f)\n', ...
    class_2_label, class_2_label, within_behavior_percent(4), ...
    within_behavior_null_mean(4), within_behavior_permutation_p(4))
fprintf('  %s classified as %s: %.1f%% (null mean %.1f%%, p=%.4f)\n', ...
    class_2_label, class_1_label, within_behavior_percent(3), ...
    within_behavior_null_mean(3), within_behavior_permutation_p(3))
fprintf('Later check, not the first model: leave-one-session-out balanced accuracy %.1f%%\n', ...
    100 * session_balanced_accuracy)
fprintf('%s baseline z-score: class 1 %.3f (n=%d), class 2 %.3f (n=%d)\n', ...
    band_name, mean(classification_table.mean_power(classification_table.behavior_class == 1)), n_class_1, ...
    mean(classification_table.mean_power(classification_table.behavior_class == 2)), n_class_2)
disp(confusion_outcomes)
disp(classification_table(1:min(8, height(classification_table)), :))

outcome_figure_name = classification_figure_name(band_name, behavior_class_1, behavior_class_2, ...
    [classifier_name '_outcomes']);
figure
set(gcf, 'WindowState', 'maximized', 'NumberTitle', 'off', 'Name', outcome_figure_name)
outcome_bars = bar(outcome_percent);
outcome_bars.FaceColor = 'flat';
outcome_bars.CData = [0.20 0.45 0.70; 0.85 0.35 0.20; 0.85 0.35 0.20; 0.20 0.45 0.70];
xticklabels({'True negative', 'False positive', 'False negative', 'True positive'})
ylabel('Percent of single events')
title(sprintf('%s, %s. Pooled, %d animals, %d sessions, %d events. Positive class: %s', ...
    band_name, classifier_name, n_animals_pooled, n_sessions_pooled, n_single_events, class_2_label))
ylim([0 100])
text(1:4, outcome_percent, compose('%.1f%%', outcome_percent), ...
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
ylabel('Baseline z-score of one event')
title(sprintf('%s, %s, one event per point', band_name, classifier_name))
set(gca, 'TickDir', 'out')

subplot(1, 2, 2)
confusion_chart = confusionchart( ...
    categorical(true_label(cv_scored), {class_1_label, class_2_label}), ...
    categorical(predicted_label(cv_scored), {class_1_label, class_2_label}), ...
    'Normalization', 'total-normalized');
confusion_chart.Title = sprintf('Pooled events. Percent of %d single events. Balanced accuracy %.1f%%', ...
    n_single_events, 100 * cv_balanced_accuracy);

%% Compare confusion outcomes across delta, theta, and gamma
bands_to_compare = {'delta', 'theta', 'gamma'};
band_outcome_counts = nan(numel(bands_to_compare), 4);
band_outcome_percent = nan(numel(bands_to_compare), 4);
band_within_behavior_percent = nan(numel(bands_to_compare), 4);
band_within_behavior_null = nan(numel(bands_to_compare), n_permutations, 4);
band_within_behavior_null_mean = nan(numel(bands_to_compare), 4);
band_within_behavior_p = nan(numel(bands_to_compare), 4);
for band_n = 1:numel(bands_to_compare)
    this_band = bands_to_compare{band_n};
    if strcmp(this_band, band_name)
        this_table = classification_table;
        this_counts = outcome_count(:).';
    else
        try
            this_table = behavior_event_table(saving_folder, this_band, all_behaviors, ...
                animal_label, behavior_class_1, behavior_class_2, power_time_window, baseline_time_window);
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
        this_counts = confusion_outcome_counts( ...
            this_table.behavior_class(this_scored), this_prediction(this_scored));
    end
    class_1_power = this_table.mean_power(this_table.behavior_class == 1);
    class_2_power = this_table.mean_power(this_table.behavior_class == 2);
    fprintf('%s pooled across %d animals and %d sessions. Baseline z-score: class 1 %.3f (n=%d), class 2 %.3f (n=%d)\n', ...
        this_band, numel(unique(this_table.animal_name)), numel(unique(this_table.session_name)), ...
        mean(class_1_power), numel(class_1_power), mean(class_2_power), numel(class_2_power))
    band_outcome_counts(band_n, :) = this_counts;
    band_outcome_percent(band_n, :) = 100 * this_counts / sum(this_counts);
    band_within_behavior_percent(band_n, :) = within_behavior_percent_from_counts(this_counts);
    if strcmp(this_band, band_name)
        this_null = permutation_within_behavior_percent;
    else
        fprintf('Label-shuffle null for %s (%d permutations)\n', this_band, n_permutations)
        this_null = permutation_within_behavior_rates( ...
            this_table.mean_power, this_table.behavior_class, ...
            classifier_name, this_n_folds, n_permutations);
    end
    band_within_behavior_null(band_n, :, :) = reshape( ...
        this_null, [1, size(this_null, 1), size(this_null, 2)]);
    band_within_behavior_null_mean(band_n, :) = mean(this_null, 1, 'omitnan');
    band_within_behavior_p(band_n, :) = permutation_p_above( ...
        this_null, band_within_behavior_percent(band_n, :));
end

band_confusion_outcomes = array2table(band_outcome_percent, ...
    'VariableNames', outcome_name', 'RowNames', bands_to_compare);
disp('Percent of single events')
disp(band_confusion_outcomes)

comparison_figure_name = classification_figure_name('delta_theta_gamma', ...
    behavior_class_1, behavior_class_2, [classifier_name '_outcomes']);
figure
set(gcf, 'WindowState', 'maximized', 'NumberTitle', 'off', 'Name', comparison_figure_name)
comparison_bars = bar(band_outcome_percent.');
comparison_bars(1).FaceColor = [0.35 0.35 0.35];
comparison_bars(2).FaceColor = [0.20 0.55 0.35];
comparison_bars(3).FaceColor = [0.75 0.35 0.55];
xticklabels({'True negative', 'False positive', 'False negative', 'True positive'})
ylabel('Percent of single events')
ylim([0 100])
legend(bands_to_compare, 'Location', 'best')
title(sprintf('%s. Pooled across sessions and animals. Equal training counts. Positive class: %s', ...
    classifier_name, class_2_label))
set(gca, 'TickDir', 'out')

%% Confusion rates separated by the true behavior
% Columns are A as A, A as B, B as A, B as B. Each behavior's two rates
% sum to 100.
classification_result.band_within_behavior_percent = band_within_behavior_percent;
classification_result.band_within_behavior_null = band_within_behavior_null;
classification_result.band_within_behavior_null_mean = band_within_behavior_null_mean;
classification_result.band_within_behavior_p_above_null = band_within_behavior_p;
classification_result.bands_compared = bands_to_compare;
rate_variable_names = {'class_1_as_class_1', 'class_1_as_class_2', ...
    'class_2_as_class_1', 'class_2_as_class_2'};
fprintf('Percent of each true behavior. Class 1 is %s. Class 2 is %s.\n', ...
    class_1_label, class_2_label)
fprintf('Columns: class 1 as class 1, class 1 as class 2, class 2 as class 1, class 2 as class 2.\n')
within_behavior_rates = array2table(band_within_behavior_percent, ...
    'VariableNames', rate_variable_names, 'RowNames', bands_to_compare);
disp(within_behavior_rates)
fprintf('Null mean percent from label shuffles. Same columns.\n')
disp(array2table(band_within_behavior_null_mean, ...
    'VariableNames', rate_variable_names, 'RowNames', bands_to_compare))
fprintf('One-sided p that the observed percent is above the label-shuffle null. Same columns.\n')
fprintf('Class 1 as class 2 is 100 minus class 1 as class 1. The same pairing holds for class 2, so each pair is two tails of one null.\n')
disp(array2table(band_within_behavior_p, ...
    'VariableNames', rate_variable_names, 'RowNames', bands_to_compare))

by_behavior_figure_name = classification_figure_name('delta_theta_gamma', ...
    behavior_class_1, behavior_class_2, [classifier_name '_by_behavior']);
figure
set(gcf, 'WindowState', 'maximized', 'NumberTitle', 'off', 'Name', by_behavior_figure_name)
band_colors = [0.35 0.35 0.35; 0.20 0.55 0.35; 0.75 0.35 0.55];

subplot(1, 2, 1)
class_1_bars = bar(band_within_behavior_percent(:, 1:2).');
color_classification_bars(class_1_bars, band_colors);
mark_within_behavior_null(class_1_bars, ...
    band_within_behavior_null_mean(:, 1:2), band_within_behavior_p(:, 1:2));
xticklabels({
    sprintf('%s as %s', class_1_label, class_1_label)
    sprintf('%s as %s', class_1_label, class_2_label)
    })
xtickangle(20)
ylabel('Percent of this behavior''s events')
ylim([0 115])
title(sprintf('True behavior: %s', class_1_label))
legend(bands_to_compare, 'Location', 'best')
set(gca, 'TickDir', 'out')

subplot(1, 2, 2)
class_2_bars = bar(band_within_behavior_percent(:, [4 3]).');
color_classification_bars(class_2_bars, band_colors);
mark_within_behavior_null(class_2_bars, ...
    band_within_behavior_null_mean(:, [4 3]), band_within_behavior_p(:, [4 3]));
xticklabels({
    sprintf('%s as %s', class_2_label, class_2_label)
    sprintf('%s as %s', class_2_label, class_1_label)
    })
xtickangle(20)
ylabel('Percent of this behavior''s events')
ylim([0 115])
title(sprintf('True behavior: %s', class_2_label))
legend(bands_to_compare, 'Location', 'best')
set(gca, 'TickDir', 'out')
sgtitle(sprintf(['%s. Pooled across sessions and animals. Black dots are label-shuffle means. ' ...
    'Stars mark a rate above that null (one star p < 0.05, two p < 0.01, three p < 0.001).'], ...
    classifier_name))

function [prediction, scored] = cross_validated_predictions(power, class_label, classifier_name, n_folds)
prediction = nan(size(class_label));
partition = cvpartition(class_label, 'KFold', n_folds, 'Stratify', true);
for fold_n = 1:n_folds
    train_index = balanced_training_index(training(partition, fold_n), class_label);
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
    train_index = balanced_training_index(~test_index, class_label);
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
        model = fitcdiscr(power, class_label, 'DiscrimType', 'linear', 'Prior', 'uniform');
    case 'svm'
        model = fitcsvm(power, class_label, 'KernelFunction', 'linear', 'Standardize', true);
    case 'naive_bayes'
        model = fitcnb(power, class_label, 'Prior', 'uniform');
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

function balanced_index = balanced_training_index(train_index, class_label)
% Keep every event from the smaller class, and the same number from the larger class.
positions = find(train_index);
labels = class_label(positions);
class_1_positions = positions(labels == 1);
class_2_positions = positions(labels == 2);
n_keep = min(numel(class_1_positions), numel(class_2_positions));
balanced_index = false(size(train_index));
if n_keep < 1
    return
end
class_1_positions = class_1_positions(randperm(numel(class_1_positions), n_keep));
class_2_positions = class_2_positions(randperm(numel(class_2_positions), n_keep));
balanced_index(class_1_positions) = true;
balanced_index(class_2_positions) = true;
end

function value = balanced_accuracy(true_class, predicted_class)
recall_class_1 = mean(predicted_class(true_class == 1) == 1);
recall_class_2 = mean(predicted_class(true_class == 2) == 2);
value = mean([recall_class_1, recall_class_2]);
end

function figure_name = classification_figure_name(band_name, class_1, class_2, classifier_name)
class_1_name = strjoin(class_1, '_');
class_2_name = strjoin(class_2, '_');
figure_name = ['classify_' band_name '_' class_1_name '_vs_' class_2_name '_' classifier_name '_pooled'];
figure_name = regexprep(figure_name, '[^A-Za-z0-9_]', '');
end

function counts = confusion_outcome_counts(true_class, predicted_class)
% Class 2 is the positive class.
% counts is [true_negative, false_positive, false_negative, true_positive].
confusion_counts = confusionmat(true_class, predicted_class, 'Order', [1 2]);
counts = [confusion_counts(1, 1), confusion_counts(1, 2), ...
    confusion_counts(2, 1), confusion_counts(2, 2)];
end

function percent = within_behavior_percent_from_counts(counts)
% counts is [A as A, A as B, B as A, B as B].
% percent is those four calls as a percent of the true behavior.
counts = counts(:).';
n_class_1 = counts(1) + counts(2);
n_class_2 = counts(3) + counts(4);
percent = nan(1, 4);
if n_class_1 > 0
    percent(1:2) = 100 * counts(1:2) / n_class_1;
end
if n_class_2 > 0
    percent(3:4) = 100 * counts(3:4) / n_class_2;
end
end

function color_classification_bars(bar_handles, band_colors)
for band_n = 1:min(numel(bar_handles), size(band_colors, 1))
    bar_handles(band_n).FaceColor = band_colors(band_n, :);
end
end

function null_rates = permutation_within_behavior_rates(power, class_label, classifier_name, n_folds, n_permutations)
% Each row is one label shuffle: A as A, A as B, B as A, B as B.
null_rates = nan(n_permutations, 4);
for permutation_n = 1:n_permutations
    shuffled_class = class_label(randperm(numel(class_label)));
    [shuffled_prediction, shuffled_scored] = cross_validated_predictions( ...
        power, shuffled_class, classifier_name, n_folds);
    if ~any(shuffled_scored)
        continue
    end
    null_rates(permutation_n, :) = within_behavior_percent_from_counts( ...
        confusion_outcome_counts(shuffled_class(shuffled_scored), shuffled_prediction(shuffled_scored)));
end
end

function p_value = permutation_p_above(null_values, observed)
% One-sided upper tail. null_values is n_permutations by n_rates.
% Shuffles that did not produce a rate are left out of that tail.
n_permutations = size(null_values, 1);
observed = observed(:).';
p_value = nan(1, numel(observed));
if n_permutations == 0
    return
end
for rate_n = 1:numel(observed)
    if ~isfinite(observed(rate_n))
        continue
    end
    null_rate = null_values(:, rate_n);
    null_rate = null_rate(isfinite(null_rate));
    n_valid = numel(null_rate);
    if n_valid == 0
        continue
    end
    p_value(rate_n) = (nnz(null_rate >= observed(rate_n)) + 1) / (n_valid + 1);
end
end

function mark_within_behavior_null(bar_handles, null_mean, p_values)
% null_mean and p_values are n_bands by n_groups, one row per bar series.
hold on
for band_n = 1:min(numel(bar_handles), size(null_mean, 1))
    if ~isprop(bar_handles(band_n), 'XEndPoints')
        continue
    end
    x_end = bar_handles(band_n).XEndPoints;
    y_end = bar_handles(band_n).YEndPoints;
    for group_n = 1:numel(x_end)
        if group_n > size(null_mean, 2)
            break
        end
        if isfinite(null_mean(band_n, group_n))
            plot(x_end(group_n), null_mean(band_n, group_n), 'k.', ...
                'MarkerSize', 16, 'HandleVisibility', 'off');
        end
        if ~isfinite(p_values(band_n, group_n)) || ~isfinite(y_end(group_n))
            continue
        end
        if p_values(band_n, group_n) < 0.001
            star_label = '***';
        elseif p_values(band_n, group_n) < 0.01
            star_label = '**';
        elseif p_values(band_n, group_n) < 0.05
            star_label = '*';
        else
            star_label = '';
        end
        if ~isempty(star_label)
            text(x_end(group_n), y_end(group_n), star_label, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 12, 'HandleVisibility', 'off');
        end
    end
end
end

function event_table = behavior_event_table(saving_folder, band_name, all_behaviors, ...
        animal_label, behavior_class_1, behavior_class_2, power_time_window, baseline_time_window)
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
animal_name_self = cell(1, numel(all_behaviors));
session_name_self = cell(1, numel(all_behaviors));
for behavior_n = 1:numel(all_behaviors)
    animal_name_self{behavior_n} = cell(0, 1);
    session_name_self{behavior_n} = cell(0, 1);
end

for j = 1:numel(psth_structure)
    if ~contains(animal_names{j, 1}, animal_label)
        continue
    end
    animal_num = find(cell2mat(cellfun(@(name) contains(animal_names{j, 1}, name), ...
        animal_label, 'UniformOutput', false)));
    electrode_num = animal_names{j, 2};
    this_animal_name = animal_label{animal_num(1)};
    this_session_name = char(animal_names{j, 1});

    for b = requested_index
        this_psth = psth_structure(j).all_behavior_psth_self{3, b};
        n_events = size(this_psth, 1);
        if n_events == 0
            continue
        end
        psth_timewrap_self{b} = [psth_timewrap_self{b}; this_psth];
        session_index_self{b} = [session_index_self{b}; animal_num(1) * ones(n_events, 1)];
        electrode_index_self{b} = [electrode_index_self{b}; electrode_num * ones(n_events, 1)];
        animal_name_self{b} = [animal_name_self{b}; repmat({this_animal_name}, n_events, 1)];
        session_name_self{b} = [session_name_self{b}; repmat({this_session_name}, n_events, 1)];
    end
end

mean_power = [];
behavior_class = [];
behavior_name = cell(0, 1);
animal_name = cell(0, 1);
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
        event_animal_name = animal_name_self{behavior_index};
        event_session_name = session_name_self{behavior_index};
        if isempty(psth)
            continue
        end
        n_use = min([size(psth, 1), numel(electrode), numel(event_animal_name), numel(event_session_name)]);
        psth = psth(1:n_use, :);
        electrode = electrode(1:n_use);
        event_animal_name = event_animal_name(1:n_use);
        event_session_name = event_session_name(1:n_use);
        keep = electrode == 1;
        psth = psth(keep, :);
        event_animal_name = event_animal_name(keep);
        event_session_name = event_session_name(keep);
        if isempty(psth)
            continue
        end
        n_time = min(size(psth, 2), numel(time_wrap_time));
        in_window = false(1, size(psth, 2));
        in_baseline = false(1, size(psth, 2));
        in_window(1:n_time) = time_wrap_time(1:n_time) > power_time_window(1) & ...
            time_wrap_time(1:n_time) <= power_time_window(2);
        in_baseline(1:n_time) = time_wrap_time(1:n_time) >= baseline_time_window(1) & ...
            time_wrap_time(1:n_time) <= baseline_time_window(2);
        if ~any(in_window)
            error('power_time_window [%g %g] contains no time-wrapped samples.', ...
                power_time_window(1), power_time_window(2));
        end
        if ~any(in_baseline)
            error('baseline_time_window [%g %g] contains no time-wrapped samples.', ...
                baseline_time_window(1), baseline_time_window(2));
        end
        baseline = psth(:, in_baseline);
        baseline_mean = mean(baseline, 2, 'omitnan');
        baseline_sd = std(baseline, 0, 2, 'omitnan');
        baseline_sd(baseline_sd == 0) = NaN;
        z_behavior = (psth(:, in_window) - baseline_mean) ./ baseline_sd;
        this_mean = mean(z_behavior, 2, 'omitnan');
        valid = isfinite(this_mean);
        mean_power = [mean_power; this_mean(valid)]; %#ok<AGROW>
        behavior_class = [behavior_class; class_id * ones(nnz(valid), 1)]; %#ok<AGROW>
        behavior_name = [behavior_name; repmat(class_behaviors(name_n), nnz(valid), 1)]; %#ok<AGROW>
        animal_name = [animal_name; event_animal_name(valid)]; %#ok<AGROW>
        session_name = [session_name; event_session_name(valid)]; %#ok<AGROW>
    end
end

event_table = table(mean_power, behavior_class, behavior_name, animal_name, session_name);
end
