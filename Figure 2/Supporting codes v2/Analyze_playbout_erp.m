%% Analyze_playbout_erp
% Grand-average play-bout ERP after saturated-trial rejection.
%
% Load the struct written by Estimate_playbout_erp.m. For the first mid-LPAG
% site in each recording:
%   1. Discard trials that reach the int16 rails (abs(sample) >= 32767).
%   2. Subtract each remaining trial's mean on [-2, 0) from that trial.
%   3. Average the kept trials into one ERP per recording.
% At every time point, fit
%
%   Voltage ~ 1 + (1|Animal)
%
% The intercept is the grand mean. The figure draws that mean and its
% confidence interval. Recordings from the same animal share the random
% intercept.
%
% Run add_repo_paths.m once per MATLAB session before this script.

%% Paths and rejection settings
options = playbout_erp_defaults();
options.saturation_limit = 32767;
options.baseline_window = [-2 0];
options.electrode_index = 1;
options.lme_step_samples = 1;
options.alpha = 0.05;

saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\playbout erp';
erp_file = fullfile(saving_folder, 'playbout_erp.mat');

%% Clean session ERPs
disp('loading')
loaded = load(erp_file, 'playbout_erp');
playbout_erp = loaded.playbout_erp;
assembled = session_erps_for_lme(playbout_erp, options);

if assembled.lfp_abs_max < options.saturation_limit
    warning('Analyze_playbout_erp:saturation', ...
        ['Largest |LFP| sample is %g and the saturation limit is %g. ' ...
        'No trial can meet the int16 rail. If LFP_PAG was rescaled, ' ...
        'set options.saturation_limit in the stored sample units.'], ...
        assembled.lfp_abs_max, options.saturation_limit);
end
fprintf(['Electrode %g. Complete trials %d, saturated trials removed %d, ' ...
    'trials kept %d, recordings in the model %d.\n'], ...
    assembled.electrode_index, assembled.n_trials_complete, ...
    assembled.n_trials_saturated, assembled.n_trials_kept, ...
    assembled.n_recordings_used);

step = options.lme_step_samples;
if ~isscalar(step) || step < 1 || step ~= round(step)
    error('Analyze_playbout_erp:step', 'lme_step_samples must be a positive integer.');
end
fit_time = assembled.time_s(1:step:end);
fit_erp = assembled.session_erp(:, 1:step:end);

%% Grand average
results = grand_average_erp_lme(fit_erp, assembled.animal_id, fit_time, options.alpha);
results.session_id = assembled.session_id;
results.n_trials_complete = assembled.n_trials_complete;
results.n_trials_saturated = assembled.n_trials_saturated;
results.n_trials_kept = assembled.n_trials_kept;
results.saturation_limit = assembled.saturation_limit;
results.baseline_window = assembled.baseline_window;
results.electrode_index = assembled.electrode_index;
results.lme_step_samples = step;
results.window_s = playbout_erp(1).window_s;
results.lfp_units = playbout_erp(1).lfp_units;

%% Mean and confidence interval
figure('Color', 'w', 'Name', 'Play-bout ERP');
hold on
t = results.time_s(:).';
session_lines = plot(t, results.session_erp.');
set(session_lines, 'Color', [0.82 0.82 0.82], 'LineWidth', 0.6);
ci_patch = fill_ci_runs(t, results.ci(:, 1), results.ci(:, 2), results.mean);
mean_line = plot(t, results.mean, 'k', 'LineWidth', 1.5);
xline(0, ':k');
yline(0, ':k');
xlim([t(1), t(end)]);
xlabel('Time from play-bout onset (s)');
ylabel('LFP (baseline-corrected)');
title(sprintf('Play-bout ERP grand average (%d recordings, %d trials)', ...
    results.n_sessions, results.n_trials_kept));
ci_label = sprintf('Pointwise %g%% CI', (1 - options.alpha) * 100);
if isempty(session_lines) || isempty(ci_patch)
    legend(mean_line, 'Mixed-model mean', 'Location', 'best');
else
    legend([session_lines(1), ci_patch(1), mean_line], ...
        {'Recording ERP', ci_label, 'Mixed-model mean'}, 'Location', 'best');
end
box off

if ~isfolder(saving_folder)
    mkdir(saving_folder);
end
results_file = fullfile(saving_folder, 'playbout_erp_grand_average.mat');
save(results_file, 'results', '-v7.3');
fprintf('Saved %s\n', results_file);
try
    savefig(gcf, fullfile(saving_folder, 'playbout_erp_grand_average.fig'));
catch save_error
    warning('Analyze_playbout_erp:figure', ...
        'The figure was not written: %s', save_error.message);
end

function ci_patch = fill_ci_runs(time_s, ci_low, ci_high, mean_line)
finite = isfinite(time_s(:)) & isfinite(ci_low(:)) & ...
    isfinite(ci_high(:)) & isfinite(mean_line(:));
edges = diff([false; finite; false]);
run_start = find(edges == 1);
run_end = find(edges == -1) - 1;
ci_patch = gobjects(numel(run_start), 1);
for run_n = 1:numel(run_start)
    idx = run_start(run_n):run_end(run_n);
    this_time = time_s(idx);
    this_low = ci_low(idx);
    this_high = ci_high(idx);
    this_time = this_time(:).';
    this_low = this_low(:).';
    this_high = this_high(:).';
    ci_patch(run_n) = fill( ...
        [this_time, fliplr(this_time)], ...
        [this_low, fliplr(this_high)], ...
        [0.55 0.55 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
end
end
