%% Analyze_speed_peak_frequency
% Delta and theta FOOOF peak frequency against speed, with every session
% kept separate.
%
% Loads the structure written by Estimate_speed_peak_frequency.m. The
% within-session correlations are recomputed from the saved peak
% frequencies and spectrum speeds. A mixed model then asks whether those
% session-level correlations share a non-zero mean. It does not pool
% speed bins across sessions before the correlation.

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

plot_peak_vs_speed(peak_table, session_keys, delta_columns, theta_columns, methods);
plot_session_correlations(correlation_table);

if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
analysis_file = fullfile(saving_folder, 'speed_peak_frequency_analysis.mat');
save(analysis_file, 'peak_table', 'correlation_table', 'analysis_parameters', '-v7.3');
writetable(peak_table, fullfile(saving_folder, 'speed_peak_by_session.csv'));
writetable(correlation_table, fullfile(saving_folder, 'speed_peak_correlations_recomputed.csv'));
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
