%% Analyze_speed_peak_frequency
% Delta and theta FOOOF peak frequency against speed.
%
% Loads the structure written by Estimate_speed_peak_frequency.m. Each
% session keeps its own correlation. The regressions are fit twice per
% channel and session, then once on every session together:
%
%   PeakHz ~ Speed
%   log(PeakHz) ~ Speed
%   PeakHz ~ log(Speed)
%
% The pooled model is outcome ~ predictor + (1|Animal) + (1|Session),
% with maximum likelihood. Session is the animal, partner, and session
% index, so channels of one behavioral session share that intercept.
% Per-session slopes use raw pixels per second. Pooled PeakHz and
% log(PeakHz) use a within-animal speed z-score, because a pixel scale
% multiplies a linear slope and the animal intercept only shifts the
% baseline. Pooled PeakHz ~ log(Speed) keeps pixels per second, because
% that scale is an animal intercept on the log-speed axis.

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Analyze_speed_peak_frequency.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

parameters = speed_peak_options();
saving_folder = fullfile(parameters.data_root, 'Analysis results', ...
    'speed peak frequency');
input_file = fullfile(saving_folder, 'speed_peak_frequency.mat');
if ~isfile(input_file)
    error('Analyze_speed_peak_frequency:missing', ...
        'Run Estimate_speed_peak_frequency first. Missing %s.', input_file);
end

loaded = load(input_file, 'recordings', 'analysis_parameters');
recordings = loaded.recordings;
analysis_parameters = loaded.analysis_parameters;
[peak_table, ~] = speed_peak_session_table(recordings);
if isempty(peak_table)
    error('Analyze_speed_peak_frequency:empty', ...
        'The saved recordings have no speed bins.');
end

fprintf('Parameters stored with this estimate:\n');
disp(analysis_parameters);
fprintf(['Speed unit: %s\nDelta %.1f-%.1f Hz, theta %.1f-%.1f Hz\n' ...
    'Bin fractions %.0f-%.0f%% of each session\n%s\n'], ...
    analysis_parameters.speed_units, ...
    analysis_parameters.delta_band(1), analysis_parameters.delta_band(2), ...
    analysis_parameters.theta_band(1), analysis_parameters.theta_band(2), ...
    100 * analysis_parameters.speed_bin_fraction_limits(1), ...
    100 * analysis_parameters.speed_bin_fraction_limits(2), ...
    analysis_parameters.spectrum_speed_definition);

methods = {'Welch', 'Fourier'};
delta_columns = {'WelchDeltaCF', 'FourierDeltaCF'};
theta_columns = {'WelchThetaCF', 'FourierThetaCF'};
session_keys = unique(peak_table(:, {'Animal', 'Partner', 'SessionIndex', ...
    'Probe', 'LPAGChannel', 'SessionSource'}), 'rows', 'stable');

correlation_rows = cell(height(session_keys) * numel(methods), 18);
row_n = 0;
for session_n = 1:height(session_keys)
    this_session = peak_table(session_row_mask(peak_table, session_keys, session_n), :);
    for method_n = 1:numel(methods)
        stats = speed_peak_correlation( ...
            this_session.SpectrumSpeed, ...
            this_session.(delta_columns{method_n}), ...
            this_session.(theta_columns{method_n}), ...
            analysis_parameters.min_bins_for_correlation);
        row_n = row_n + 1;
        correlation_rows(row_n, :) = { ...
            session_keys.Animal{session_n}, ...
            session_keys.Partner(session_n), ...
            session_keys.SessionIndex(session_n), ...
            session_keys.SessionSource{session_n}, ...
            session_keys.Probe(session_n), ...
            session_keys.LPAGChannel(session_n), ...
            methods{method_n}, ...
            stats.n_delta, stats.pearson_delta, stats.pearson_delta_p, ...
            stats.spearman_delta, stats.n_theta, stats.pearson_theta, ...
            stats.pearson_theta_p, stats.spearman_theta, stats.n_delta_theta, ...
            stats.pearson_delta_theta, stats.pearson_delta_theta_p};
        fprintf(['%s P%d session %d probe %g channel %g %s: ' ...
            'delta r = %+.3f (n = %d), theta r = %+.3f (n = %d), ' ...
            'delta-theta r = %+.3f (n = %d)\n'], ...
            session_keys.Animal{session_n}, session_keys.Partner(session_n), ...
            session_keys.SessionIndex(session_n), session_keys.Probe(session_n), ...
            session_keys.LPAGChannel(session_n), methods{method_n}, ...
            stats.pearson_delta, stats.n_delta, ...
            stats.pearson_theta, stats.n_theta, ...
            stats.pearson_delta_theta, stats.n_delta_theta);
    end
end

correlation_table = cell2table(correlation_rows, 'VariableNames', { ...
    'Animal', 'Partner', 'SessionIndex', 'SessionSource', 'Probe', ...
    'LPAGChannel', 'Method', 'NDelta', 'PearsonDelta', 'PearsonDeltaP', ...
    'SpearmanDelta', 'NTheta', 'PearsonTheta', 'PearsonThetaP', ...
    'SpearmanTheta', 'NDeltaTheta', 'PearsonDeltaTheta', 'PearsonDeltaThetaP'});

for method_n = 1:numel(methods)
    method_rows = strcmp(correlation_table.Method, methods{method_n});
    report_session_lme(correlation_table(method_rows, :), 'PearsonDelta', ...
        [methods{method_n} ' delta peak vs speed']);
    report_session_lme(correlation_table(method_rows, :), 'PearsonTheta', ...
        [methods{method_n} ' theta peak vs speed']);
    report_session_lme(correlation_table(method_rows, :), 'PearsonDeltaTheta', ...
        [methods{method_n} ' delta peak vs theta peak']);
end

model_specification = speed_peak_model_specification();
fprintf('\nPeak-frequency regressions\n%s\n', model_specification.why_zscore_linear_speed);
fprintf('Comparable AIC: %s\n', model_specification.aic_hz);
session_model_table = fit_session_models(peak_table, session_keys, ...
    methods, delta_columns, theta_columns, analysis_parameters.min_bins_for_correlation);
pooled_model_table = fit_pooled_models(peak_table, methods, ...
    delta_columns, theta_columns, analysis_parameters.min_bins_for_correlation);

plot_peak_vs_speed(peak_table, session_keys, delta_columns, theta_columns, methods);
plot_session_correlations(correlation_table);
plot_peak_model_fits(peak_table, session_keys, methods, delta_columns, ...
    theta_columns, analysis_parameters.min_bins_for_correlation);

if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
analysis_file = fullfile(saving_folder, 'speed_peak_frequency_analysis.mat');
save(analysis_file, 'peak_table', 'correlation_table', 'session_model_table', ...
    'pooled_model_table', 'model_specification', 'analysis_parameters', '-v7.3');
writetable(peak_table, fullfile(saving_folder, 'speed_peak_by_session.csv'));
writetable(correlation_table, fullfile(saving_folder, 'speed_peak_correlations_recomputed.csv'));
writetable(session_model_table, fullfile(saving_folder, 'speed_peak_session_models.csv'));
writetable(pooled_model_table, fullfile(saving_folder, 'speed_peak_pooled_models.csv'));
fprintf('Saved %s\n', analysis_file);

function report_session_lme(correlation_table, value_name, label)
values = correlation_table.(value_name);
animals = correlation_table.Animal;
keep = isfinite(values);
values = values(keep);
animals = animals(keep);
fprintf('%s: %d sessions with a correlation\n', label, numel(values));
if numel(values) < 2
    return
end
[~, p_value, ci, stats] = ttest(values);
fprintf('  one-sample t: mean %+.3f, 95%% CI [%+.3f, %+.3f], t(%d) = %.2f, p = %.4f\n', ...
    mean(values), ci(1), ci(2), stats.df, stats.tstat, p_value);
if numel(unique(animals)) < 2
    fprintf('  LME skipped: only one animal\n');
    return
end
lme = fitlme(table(values, categorical(animals), ...
    'VariableNames', {'PearsonR', 'Animal'}), 'PearsonR ~ 1 + (1|Animal)');
coef = fixedEffects(lme);
ci_lme = coefCI(lme);
fprintf('  LME PearsonR ~ 1 + (1|Animal): mean %+.3f, 95%% CI [%+.3f, %+.3f], p = %.4f\n', ...
    coef, ci_lme(1), ci_lme(2), lme.Coefficients.pValue(1));
end

function plot_peak_vs_speed(peak_table, session_keys, delta_columns, theta_columns, methods)
n_sessions = height(session_keys);
colors = lines(max(n_sessions, 1));
labels = cell(n_sessions, 1);
for session_n = 1:n_sessions
    labels{session_n} = sprintf('%s P%d s%d ch%g', ...
        session_keys.Animal{session_n}, session_keys.Partner(session_n), ...
        session_keys.SessionIndex(session_n), session_keys.LPAGChannel(session_n));
end

for method_n = 1:numel(methods)
    figure('Name', [methods{method_n} ' peak frequency vs speed'])
    subplot(1, 2, 1)
    plot_one_peak(peak_table, session_keys, delta_columns{method_n}, colors, ...
        [methods{method_n} ' delta peak'])
    subplot(1, 2, 2)
    plot_one_peak(peak_table, session_keys, theta_columns{method_n}, colors, ...
        [methods{method_n} ' theta peak'])
    if n_sessions <= 12
        legend(labels, 'Location', 'eastoutside')
    end
end
end

function mask = session_row_mask(peak_table, session_keys, session_n)
animal = session_keys.Animal{session_n};
mask = strcmp(peak_table.Animal, animal) & ...
    peak_table.Partner == session_keys.Partner(session_n) & ...
    peak_table.SessionIndex == session_keys.SessionIndex(session_n) & ...
    peak_table.Probe == session_keys.Probe(session_n) & ...
    peak_table.LPAGChannel == session_keys.LPAGChannel(session_n);
end

function plot_one_peak(peak_table, session_keys, value_name, colors, plot_title)
hold on
for session_n = 1:height(session_keys)
    this_session = peak_table(session_row_mask(peak_table, session_keys, session_n), :);
    keep = isfinite(this_session.SpectrumSpeed) & isfinite(this_session.(value_name));
    if ~any(keep)
        continue
    end
    [speed, order] = sort(this_session.SpectrumSpeed(keep));
    frequency = this_session.(value_name);
    frequency = frequency(keep);
    plot(speed, frequency(order), '-o', 'Color', colors(session_n, :), ...
        'MarkerSize', 4, 'LineWidth', 1)
end
xlabel('Spectrum speed (pixels/s)')
ylabel('Peak frequency (Hz)')
title(plot_title)
box off
end

function plot_session_correlations(correlation_table)
methods = {'Welch', 'Fourier'};
value_names = {'PearsonDelta', 'PearsonTheta', 'PearsonDeltaTheta'};
panel_titles = {'Delta vs speed', 'Theta vs speed', 'Delta vs theta'};
for method_n = 1:numel(methods)
    method_rows = correlation_table(strcmp(correlation_table.Method, methods{method_n}), :);
    labels = cell(height(method_rows), 1);
    for row_n = 1:height(method_rows)
        labels{row_n} = sprintf('%s P%d s%d ch%g', ...
            method_rows.Animal{row_n}, method_rows.Partner(row_n), ...
            method_rows.SessionIndex(row_n), method_rows.LPAGChannel(row_n));
    end
    figure('Name', [methods{method_n} ' within-session correlations'])
    for panel_n = 1:numel(value_names)
        subplot(1, 3, panel_n)
        barh(method_rows.(value_names{panel_n}))
        set(gca, 'YTick', 1:height(method_rows), 'YTickLabel', labels, 'YDir', 'reverse')
        xline(0, 'k:')
        xlim([-1 1])
        xlabel('Pearson r')
        title(panel_titles{panel_n})
    end
end
end

function session_model_table = fit_session_models(peak_table, session_keys, ...
        methods, delta_columns, theta_columns, min_n)
bands = {'Delta', 'Theta'};
names = {'Animal', 'Partner', 'SessionIndex', 'SessionSource', 'Probe', ...
    'LPAGChannel', 'Method', 'Band', 'Model', 'Formula', 'N', ...
    'Intercept', 'Slope', 'InterceptSE', 'SlopeSE', 'InterceptP', ...
    'SlopeP', 'RSquared', 'Sigma', 'AIC', 'AICOnHz'};
rows = cell(height(session_keys) * numel(methods) * numel(bands) * 3, numel(names));
row_n = 0;
for session_n = 1:height(session_keys)
    this_session = peak_table(session_row_mask(peak_table, session_keys, session_n), :);
    for method_n = 1:numel(methods)
        peak_columns = {delta_columns{method_n}, theta_columns{method_n}};
        for band_n = 1:numel(bands)
            models = speed_peak_models(this_session.SpectrumSpeed, ...
                this_session.(peak_columns{band_n}), min_n);
            for model_n = 1:numel(models)
                model = models(model_n);
                row_n = row_n + 1;
                rows(row_n, :) = { ...
                    session_keys.Animal{session_n}, ...
                    session_keys.Partner(session_n), ...
                    session_keys.SessionIndex(session_n), ...
                    session_keys.SessionSource{session_n}, ...
                    session_keys.Probe(session_n), ...
                    session_keys.LPAGChannel(session_n), ...
                    methods{method_n}, bands{band_n}, model.name, model.formula, ...
                    model.n, model.intercept, model.slope, model.intercept_se, ...
                    model.slope_se, model.intercept_p, model.slope_p, ...
                    model.r_squared, model.sigma, model.aic, model.aic_hz};
                fprintf(['%s P%d session %d channel %g %s %s %s: ' ...
                    'slope %+.4g, p = %.4g, R2 = %.3f, AIC_Hz = %.1f, n = %d\n'], ...
                    session_keys.Animal{session_n}, session_keys.Partner(session_n), ...
                    session_keys.SessionIndex(session_n), ...
                    session_keys.LPAGChannel(session_n), methods{method_n}, ...
                    bands{band_n}, model.formula, model.slope, model.slope_p, ...
                    model.r_squared, model.aic_hz, model.n);
            end
        end
    end
end
session_model_table = cell2table(rows(1:row_n, :), 'VariableNames', names);
end

function pooled_model_table = fit_pooled_models(peak_table, methods, ...
        delta_columns, theta_columns, min_n)
bands = {'Delta', 'Theta'};
model_defs = { ...
    'linear', 'PeakHz ~ SpeedZ_within_animal', 'peak', 'speed_z', false, ...
        'within-animal z-score of pixels/s'; ...
    'log_frequency', 'log(PeakHz) ~ SpeedZ_within_animal', 'log_peak', 'speed_z', true, ...
        'within-animal z-score of pixels/s'; ...
    'log_speed', 'PeakHz ~ log(Speed)', 'peak', 'log_speed', false, ...
        'natural log of pixels/s'};
names = {'Method', 'Band', 'Model', 'Formula', 'SpeedNormalization', ...
    'FitMethod', 'N', 'NAnimals', 'NSessions', 'Intercept', 'Slope', ...
    'InterceptSE', 'SlopeSE', 'InterceptP', 'SlopeP', 'InterceptCILow', ...
    'InterceptCIHigh', 'SlopeCILow', 'SlopeCIHigh', 'AIC', 'AICOnHz', ...
    'LogLikelihood', 'AnimalSD', 'SessionSD', 'ResidualSD', ...
    'PeakRatioPerSpeedSD', 'Converged', 'Message'};
rows = cell(numel(methods) * numel(bands) * size(model_defs, 1), numel(names));
row_n = 0;
[animal, session_key, bin_id, speed] = pooled_speed_keys(peak_table);
for method_n = 1:numel(methods)
    peak_columns = {delta_columns{method_n}, theta_columns{method_n}};
    for band_n = 1:numel(bands)
        obs = speed_peak_pooled_observations(animal, session_key, bin_id, ...
            speed, peak_table.(peak_columns{band_n}));
        choice = speed_peak_pooled_formula(obs.animal, obs.session_key);
        for model_n = 1:size(model_defs, 1)
            outcome = obs.(model_defs{model_n, 3});
            predictor = obs.(model_defs{model_n, 4});
            fit = fit_pooled_one(obs.animal, obs.session_key, outcome, ...
                predictor, choice, obs.log_peak, model_defs{model_n, 5}, min_n);
            formula = model_defs{model_n, 2};
            if ~isempty(fit.random_terms)
                formula = [formula ' + ' fit.random_terms]; %#ok<AGROW>
            end
            row_n = row_n + 1;
            rows(row_n, :) = { ...
                methods{method_n}, bands{band_n}, model_defs{model_n, 1}, ...
                formula, model_defs{model_n, 6}, fit.fit_method, fit.n, ...
                fit.n_animals, fit.n_sessions, fit.intercept, fit.slope, ...
                fit.intercept_se, fit.slope_se, fit.intercept_p, fit.slope_p, ...
                fit.intercept_ci(1), fit.intercept_ci(2), fit.slope_ci(1), ...
                fit.slope_ci(2), fit.aic, fit.aic_hz, fit.log_likelihood, ...
                fit.animal_sd, fit.session_sd, fit.residual_sd, ...
                fit.peak_ratio, double(fit.converged), fit.message};
            fprintf(['Pooled %s %s %s: slope %+.4g (SE %.4g), p = %.4g, ' ...
                'AIC_Hz = %.1f, n = %d, animals = %d, sessions = %d\n  %s\n'], ...
                methods{method_n}, bands{band_n}, formula, fit.slope, ...
                fit.slope_se, fit.slope_p, fit.aic_hz, fit.n, fit.n_animals, ...
                fit.n_sessions, fit.message);
        end
    end
end
pooled_model_table = cell2table(rows(1:row_n, :), 'VariableNames', names);
end

function [animal, session_key, bin_id, speed] = pooled_speed_keys(peak_table)
n = height(peak_table);
animal = peak_table.Animal;
session_key = cell(n, 1);
for row_n = 1:n
    session_key{row_n} = sprintf('%s|P%d|S%d', animal{row_n}, ...
        peak_table.Partner(row_n), peak_table.SessionIndex(row_n));
end
bin_id = peak_table.Bin;
speed = peak_table.SpectrumSpeed;
end

function fit = fit_pooled_one(animal, session_key, outcome, predictor, ...
        choice, log_peak, outcome_is_log, min_n)
fit = empty_pooled_fit();
fit.n = numel(outcome);
fit.n_animals = choice.n_animals;
fit.n_sessions = choice.n_sessions;
fit.message = choice.note;
if fit.n < min_n
    fit.message = sprintf(['Fewer than %d rows with positive speed, a ' ...
        'positive peak, and a within-animal speed scale. %s'], min_n, choice.note);
    return
end
if isempty(choice.random_terms)
    fit = ols_pooled_fit(fit, predictor, outcome, log_peak, outcome_is_log, min_n);
    return
end
candidates = pooled_formula_candidates(choice.random_terms);
last_error = '';
for candidate_n = 1:numel(candidates)
    try
        fit = lme_pooled_fit(fit, animal, session_key, outcome, predictor, ...
            candidates{candidate_n}, log_peak, outcome_is_log);
        if candidate_n > 1
            fit.message = sprintf('%s Fell back to %s after: %s', ...
                choice.note, candidates{candidate_n}, last_error);
        end
        return
    catch err
        last_error = err.message;
    end
end
fit.message = sprintf('%s Mixed model failed: %s', choice.note, last_error);
fit = ols_pooled_fit(fit, predictor, outcome, log_peak, outcome_is_log, min_n);
fit.message = sprintf('%s Ordinary least squares was used instead.', fit.message);
end

function candidates = pooled_formula_candidates(random_terms)
full = ['Peak ~ Predictor + ' random_terms];
candidates = {full};
if contains(random_terms, '(1|Animal)') && contains(random_terms, '(1|Session)')
    candidates = [candidates, { ...
        'Peak ~ Predictor + (1|Animal)', ...
        'Peak ~ Predictor + (1|Session)'}]; %#ok<AGROW>
end
end

function fit = lme_pooled_fit(fit, animal, session_key, outcome, predictor, ...
        formula, log_peak, outcome_is_log)
model_table = table(outcome, predictor, categorical(animal), categorical(session_key), ...
    'VariableNames', {'Peak', 'Predictor', 'Animal', 'Session'});
lme = fitlme(model_table, formula, 'FitMethod', 'ML');
[intercept, intercept_se, intercept_p, intercept_ci] = lme_coefficient(lme, '(Intercept)');
[slope, slope_se, slope_p, slope_ci] = lme_coefficient(lme, 'Predictor');
fit.converged = isfinite(slope);
fit.fit_method = 'ML';
fit.intercept = intercept;
fit.slope = slope;
fit.intercept_se = intercept_se;
fit.slope_se = slope_se;
fit.intercept_p = intercept_p;
fit.slope_p = slope_p;
fit.intercept_ci = intercept_ci;
fit.slope_ci = slope_ci;
fit.aic = double(lme.ModelCriterion.AIC);
fit.log_likelihood = double(lme.LogLikelihood);
fit.aic_hz = aic_on_hz(fit.aic, log_peak, outcome_is_log);
fit.peak_ratio = peak_ratio_of(slope, outcome_is_log);
fit.random_terms = random_terms_from_formula(formula);
[fit.animal_sd, fit.session_sd, fit.residual_sd] = random_effect_sds(lme, formula);
end

function fit = ols_pooled_fit(fit, predictor, outcome, log_peak, outcome_is_log, min_n)
ols = speed_peak_ols(predictor, outcome, min_n);
fit.fit_method = 'OLS';
fit.converged = ols.ok;
fit.n = ols.n;
fit.intercept = ols.intercept;
fit.slope = ols.slope;
fit.intercept_se = ols.intercept_se;
fit.slope_se = ols.slope_se;
fit.intercept_p = ols.intercept_p;
fit.slope_p = ols.slope_p;
if ols.ok
    t_intercept = tinv(0.975, ols.residual_df);
    t_slope = t_intercept;
    fit.intercept_ci = ols.intercept + [-1 1] * t_intercept * ols.intercept_se;
    fit.slope_ci = ols.slope + [-1 1] * t_slope * ols.slope_se;
end
fit.aic = ols.aic;
fit.log_likelihood = ols.log_likelihood;
fit.aic_hz = aic_on_hz(ols.aic, log_peak, outcome_is_log);
fit.peak_ratio = peak_ratio_of(ols.slope, outcome_is_log);
fit.residual_sd = ols.sigma;
fit.random_terms = '';
end

function [estimate, se, p_value, ci] = lme_coefficient(lme, name)
coef = lme.Coefficients;
names = cellstr(coef.Name);
match = strcmp(names, name);
if ~any(match)
    error('Analyze_speed_peak_frequency:coefficient', ...
        'The mixed model has no %s coefficient.', name);
end
estimate = coef.Estimate(match);
se = coef.SE(match);
p_value = coef.pValue(match);
ci = [NaN NaN];
intervals = coefCI(lme);
ci(1) = intervals(match, 1);
ci(2) = intervals(match, 2);
end

function random_terms = random_terms_from_formula(formula)
prefix = 'Peak ~ Predictor + ';
if strncmp(formula, prefix, numel(prefix))
    random_terms = formula(numel(prefix) + 1:end);
else
    random_terms = '';
end
end

function [animal_sd, session_sd, residual_sd] = random_effect_sds(lme, formula)
animal_sd = NaN;
session_sd = NaN;
[psi, mse] = covarianceParameters(lme);
groups = {};
if contains(formula, '(1|Animal)')
    groups{end + 1} = 'Animal'; %#ok<AGROW>
end
if contains(formula, '(1|Session)')
    groups{end + 1} = 'Session'; %#ok<AGROW>
end
for group_n = 1:numel(groups)
    if iscell(psi)
        variance = psi{group_n};
    else
        variance = psi;
    end
    sd = sqrt(max(double(variance(1)), 0));
    if strcmp(groups{group_n}, 'Animal')
        animal_sd = sd;
    else
        session_sd = sd;
    end
end
residual_sd = sqrt(max(double(mse(1)), 0));
end

function aic_hz = aic_on_hz(aic, log_peak, outcome_is_log)
aic_hz = aic;
if outcome_is_log && isfinite(aic)
    aic_hz = aic + 2 * sum(log_peak(:));
end
end

function ratio = peak_ratio_of(slope, outcome_is_log)
ratio = NaN;
if outcome_is_log && isfinite(slope)
    ratio = exp(slope);
end
end

function fit = empty_pooled_fit()
fit = struct( ...
    'converged', false, ...
    'fit_method', '', ...
    'n', 0, ...
    'n_animals', 0, ...
    'n_sessions', 0, ...
    'intercept', NaN, ...
    'slope', NaN, ...
    'intercept_se', NaN, ...
    'slope_se', NaN, ...
    'intercept_p', NaN, ...
    'slope_p', NaN, ...
    'intercept_ci', [NaN NaN], ...
    'slope_ci', [NaN NaN], ...
    'aic', NaN, ...
    'aic_hz', NaN, ...
    'log_likelihood', NaN, ...
    'animal_sd', NaN, ...
    'session_sd', NaN, ...
    'residual_sd', NaN, ...
    'peak_ratio', NaN, ...
    'random_terms', '', ...
    'message', '');
end

function plot_peak_model_fits(peak_table, session_keys, methods, ...
        delta_columns, theta_columns, min_n)
n_sessions = height(session_keys);
colors = lines(max(n_sessions, 1));
labels = cell(n_sessions, 1);
for session_n = 1:n_sessions
    labels{session_n} = sprintf('%s P%d s%d ch%g', ...
        session_keys.Animal{session_n}, session_keys.Partner(session_n), ...
        session_keys.SessionIndex(session_n), session_keys.LPAGChannel(session_n));
end
bands = {'Delta', 'Theta'};
for method_n = 1:numel(methods)
    peak_columns = {delta_columns{method_n}, theta_columns{method_n}};
    figure('Name', [methods{method_n} ' peak-frequency models'])
    for band_n = 1:numel(bands)
        subplot(2, 2, band_n)
        plot_session_fits(peak_table, session_keys, peak_columns{band_n}, ...
            colors, labels, min_n, 'speed')
        title(sprintf('%s: solid PeakHz ~ Speed, dashed log(PeakHz) ~ Speed', ...
            bands{band_n}))
        subplot(2, 2, band_n + 2)
        plot_session_fits(peak_table, session_keys, peak_columns{band_n}, ...
            colors, labels, min_n, 'log_speed')
        title(sprintf('%s: PeakHz ~ log(Speed)', bands{band_n}))
    end
    if n_sessions <= 12 && ~isempty(findobj(gca, 'Type', 'line'))
        legend('show', 'Location', 'eastoutside')
    end
end
end

function plot_session_fits(peak_table, session_keys, value_name, colors, ...
        labels, min_n, x_scale)
hold on
for session_n = 1:height(session_keys)
    this_session = peak_table(session_row_mask(peak_table, session_keys, session_n), :);
    speed = this_session.SpectrumSpeed;
    peak = this_session.(value_name);
    keep = isfinite(speed) & isfinite(peak) & speed > 0 & peak > 0;
    if ~any(keep)
        continue
    end
    speed = speed(keep);
    peak = peak(keep);
    models = speed_peak_models(speed, peak, min_n);
    if strcmp(x_scale, 'speed')
        x = speed;
        plot(x, peak, 'o', 'Color', colors(session_n, :), ...
            'MarkerSize', 4, 'DisplayName', labels{session_n})
        draw_fit_line(models, 'linear', speed, colors(session_n, :), '-')
        draw_log_frequency_curve(models, speed, colors(session_n, :))
        xlabel('Spectrum speed (pixels/s)')
    else
        plot(log(speed), peak, 'o', 'Color', colors(session_n, :), ...
            'MarkerSize', 4, 'DisplayName', labels{session_n})
        draw_fit_line(models, 'log_speed', log(speed), colors(session_n, :), '-')
        xlabel('log spectrum speed (log pixels/s)')
    end
end
ylabel('Peak frequency (Hz)')
box off
end

function draw_fit_line(models, model_name, x, color, style)
model = models(strcmp({models.name}, model_name));
if ~model.ok || max(x) <= min(x)
    return
end
x_line = linspace(min(x), max(x), 40);
y_line = model.intercept + model.slope * x_line;
plot(x_line, y_line, style, 'Color', color, 'LineWidth', 1, ...
    'HandleVisibility', 'off')
end

function draw_log_frequency_curve(models, speed, color)
model = models(strcmp({models.name}, 'log_frequency'));
if ~model.ok || max(speed) <= min(speed)
    return
end
speed_line = linspace(min(speed), max(speed), 40);
peak_line = exp(model.intercept + model.slope * speed_line);
plot(speed_line, peak_line, '--', 'Color', color, 'LineWidth', 1, ...
    'HandleVisibility', 'off')
end
