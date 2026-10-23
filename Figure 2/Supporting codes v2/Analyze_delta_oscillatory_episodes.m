%% Analyze_delta_oscillatory_episodes
% Plots and the session summary for candidate continuous delta episodes.
%
% Run Estimate_delta_oscillatory_episodes.m first. This script reads
% Data/Analysis results/delta oscillatory episodes/delta_oscillatory_episodes.mat
% and writes:
%   session_summary.csv
%   figures/<recording>_probe<P>_channel<C>.png
%       minimally processed LFP, 1-5 Hz trace, Hilbert amplitude,
%       detected episode, peaks, troughs (cycle boundaries), and the
%       trough-to-trough period histogram
%   figures/cycle_period_histogram_pooled.png
%   figures/session_summary.png
%       episode duration, trough-cycle count, period variability,
%       and rise-decay / peak-trough symmetry
%
% The histogram axis is 50-1000 ms. A 1-5 Hz cycle is 200-1000 ms
% (1/5 Hz to 1/1 Hz); that band is shaded. Periods near 50 ms are faster
% than the delta band and are shown so they are not mistaken for 5 Hz.
% When trough cycles are saved, the histogram uses those periods.
% Older files that only have peak-to-peak intervals still plot those.
%
% Episodes are candidates. Filtering can ring around a transient.

%% Load
show_figures = false;

this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Analyze_delta_oscillatory_episodes.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

options = delta_episode_options();
saving_folder = fullfile(options.data_root, 'Analysis results', ...
    'delta oscillatory episodes');
input_file = fullfile(saving_folder, 'delta_oscillatory_episodes.mat');
figure_folder = fullfile(saving_folder, 'figures');

fprintf('Loading %s\n', input_file);
load(input_file, 'delta_episodes');
if isempty(delta_episodes)
    error('The saved delta_episodes structure is empty.');
end
fprintf('%s\n', delta_episodes(1).parameters.interpretation);

if ~isfolder(figure_folder)
    mkdir(figure_folder);
end

session_summary = delta_episode_summary_table(delta_episodes);
writetable(session_summary, fullfile(saving_folder, 'session_summary.csv'));
disp(session_summary);
print_session_lines(session_summary);

%% One figure per mid-LPAG channel
for recording_n = 1:numel(delta_episodes)
    recording = delta_episodes(recording_n);
    figure_handle = new_figure(show_figures, [1200 900]);
    plot_recording(recording);
    figure_name = sprintf('%s_probe%g_channel%g.png', ...
        file_token(recording.animal_name), recording.lpag_probe_number, ...
        recording.lpag_channel);
    print(figure_handle, fullfile(figure_folder, figure_name), '-dpng', '-r140');
    if ~show_figures
        close(figure_handle);
    end
end

%% Pooled cycle-period histogram
[periods_s, period_kind] = all_cycle_periods(delta_episodes);
figure_handle = new_figure(show_figures, [900 500]);
parameters = delta_episodes(1).parameters;
plot_cycle_histogram(periods_s, parameters.delta_freq, ...
    parameters.histogram_ms, parameters.histogram_bin_ms);
title(cycle_histogram_title(periods_s, parameters.delta_freq, ...
    parameters.histogram_ms));
xlabel(period_axis_label(period_kind));
ylabel('Cycles');
print(figure_handle, fullfile(figure_folder, 'cycle_period_histogram_pooled.png'), ...
    '-dpng', '-r140');
if ~show_figures
    close(figure_handle);
end

%% Session summary
figure_handle = new_figure(show_figures, [1400 980]);
plot_session_summary(delta_episodes, session_summary);
print(figure_handle, fullfile(figure_folder, 'session_summary.png'), ...
    '-dpng', '-r140');
if ~show_figures
    close(figure_handle);
end
fprintf('Figures written to %s\n', figure_folder);

%% -------------------------------------------------------------------------
function plot_recording(recording)
parameters = recording.parameters;
example = recording.example;
[periods_s, period_kind] = recording_cycle_periods(recording);

subplot(4, 1, 1);
plot_trace_panel(example, 'lfp_minimal', [0.25 0.25 0.25]);
hold on
plot_cycle_marks(example, 'lfp_minimal');
ylabel('Minimal LFP');
title(sprintf('%s probe %g channel %g', recording.animal_name, ...
    recording.lpag_probe_number, recording.lpag_channel), 'Interpreter', 'none');

subplot(4, 1, 2);
plot_trace_panel(example, 'lfp_delta', [0.75 0.15 0.15]);
hold on
if example.has_episode && ~isempty(example.peak_index)
    peak_index = example.peak_index(example.peak_index >= 1 & ...
        example.peak_index <= numel(example.lfp_delta));
    plot(example.time_s(peak_index), double(example.lfp_delta(peak_index)), ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
end
plot_trough_marks(example, 'lfp_delta');
ylabel('Delta 1-5 Hz');
title('Circles: 1-5 Hz peaks. Triangles: troughs (cycle boundaries).');

subplot(4, 1, 3);
plot_trace_panel(example, 'amplitude', [0.1 0.25 0.7]);
hold on
if ~isempty(example.time_s) && isfinite(example.threshold)
    plot(example.time_s([1 end]), [example.threshold example.threshold], 'k--');
end
ylabel('Hilbert amplitude');
xlabel('Recording time (s)');
if example.has_episode
    text(0.01, 0.9, episode_annotation(example), ...
        'Units', 'normalized', 'FontSize', 9);
else
    text(0.01, 0.9, 'No candidate episode in this recording', ...
        'Units', 'normalized', 'FontSize', 9);
end

subplot(4, 1, 4);
plot_cycle_histogram(periods_s, parameters.delta_freq, ...
    parameters.histogram_ms, parameters.histogram_bin_ms);
title(cycle_histogram_title(periods_s, parameters.delta_freq, ...
    parameters.histogram_ms));
xlabel(period_axis_label(period_kind));
ylabel('Cycles');
end

function note = episode_annotation(example)
note = sprintf(['full episode %.2f s, %d peak-to-peak cycles, ' ...
    'median %.0f ms, IQR %.0f ms'], ...
    example.episode_duration_s, example.n_cycles, ...
    1000 * example.median_cycle_period_s, ...
    1000 * example.cycle_period_iqr_s);
if isfield(example, 'n_trough_cycles')
    sustained = 'no';
    if example.has_sustained_sequence
        sustained = 'yes';
    end
    note = sprintf(['%s\n%d trough cycles, median %.0f ms, CV %.2f, ' ...
        'sustained sequence %s'], note, example.n_trough_cycles, ...
        1000 * example.median_trough_period_s, example.trough_period_cv, ...
        sustained);
end
end

function plot_cycle_marks(example, field_name)
plot_trough_marks(example, field_name);
if ~example.has_episode || ~isfield(example, 'cycle_peak_index') || ...
        isempty(example.cycle_peak_index) || isempty(example.time_s)
    return
end
peak_index = example.cycle_peak_index(example.cycle_peak_index >= 1 & ...
    example.cycle_peak_index <= numel(example.(field_name)));
if isempty(peak_index)
    return
end
plot(example.time_s(peak_index), double(example.(field_name)(peak_index)), ...
    '^', 'Color', [0.85 0.4 0.1], 'MarkerFaceColor', [0.85 0.4 0.1], ...
    'MarkerSize', 5);
end

function plot_trough_marks(example, field_name)
if ~example.has_episode || ~isfield(example, 'trough_index') || ...
        isempty(example.trough_index) || isempty(example.time_s)
    return
end
trough_index = example.trough_index(example.trough_index >= 1 & ...
    example.trough_index <= numel(example.(field_name)));
if isempty(trough_index)
    return
end
plot(example.time_s(trough_index), double(example.(field_name)(trough_index)), ...
    'v', 'Color', [0.1 0.35 0.75], 'MarkerFaceColor', [0.1 0.35 0.75], ...
    'MarkerSize', 5);
end

function label = period_axis_label(period_kind)
if strcmp(period_kind, 'trough')
    label = 'Trough-to-trough cycle period (ms)';
else
    label = 'Peak-to-peak cycle duration (ms)';
end
end

function plot_trace_panel(example, field_name, color)
if isempty(example.time_s)
    text(0.5, 0.5, 'No excerpt saved', 'HorizontalAlignment', 'center');
    axis off
    return
end
plot(example.time_s, double(example.(field_name)), 'Color', color);
hold on
shade_episode(example);
xlim(example.time_s([1 end]));
box off
end

function shade_episode(example)
if ~example.has_episode || isempty(example.episode_sample_range) || ...
        isempty(example.time_s)
    return
end
sample_range = example.episode_sample_range;
if size(sample_range, 2) ~= 2 || sample_range(1) > numel(example.time_s)
    return
end
sample_range(1) = max(1, sample_range(1));
sample_range(2) = min(numel(example.time_s), sample_range(2));
if sample_range(2) < sample_range(1)
    return
end
y_limits = ylim;
patch_x = example.time_s([sample_range(1), sample_range(2), ...
    sample_range(2), sample_range(1)]);
patch(patch_x, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
    [0.85 0.2 0.2], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
uistack(findobj(gca, 'Type', 'patch'), 'bottom');
ylim(y_limits);
end

function plot_cycle_histogram(periods_s, delta_freq, histogram_ms, bin_ms)
periods_ms = periods_s(:) * 1000;
edges = histogram_ms(1):bin_ms:histogram_ms(2);
if numel(edges) < 2
    edges = histogram_ms(:).';
end
histogram(periods_ms, edges, 'FaceColor', [0.2 0.35 0.55], 'EdgeColor', 'none');
hold on
y_limits = ylim;
if y_limits(2) <= 0
    y_limits(2) = 1;
end
band_ms = [1000 / delta_freq(2), 1000 / delta_freq(1)];
patch([band_ms(1) band_ms(2) band_ms(2) band_ms(1)], ...
    [0 0 y_limits(2) y_limits(2)], [0.55 0.75 0.95], ...
    'FaceAlpha', 0.35, 'EdgeColor', 'none');
uistack(findobj(gca, 'Type', 'patch'), 'bottom');
ylim(y_limits);
xlim(histogram_ms);
box off
end

function title_text = cycle_histogram_title(periods_s, delta_freq, histogram_ms)
periods_ms = periods_s(:) * 1000;
band_ms = [1000 / delta_freq(2), 1000 / delta_freq(1)];
n_in_band = nnz(periods_ms >= band_ms(1) & periods_ms <= band_ms(2));
n_outside = nnz(periods_ms < histogram_ms(1) | periods_ms > histogram_ms(2));
title_text = sprintf(['%d cycles; %d in 1-5 Hz (%.0f-%.0f ms); ' ...
    '%d outside the %.0f-%.0f ms axis'], ...
    numel(periods_ms), n_in_band, band_ms(1), band_ms(2), ...
    n_outside, histogram_ms(1), histogram_ms(2));
end

function plot_session_summary(recordings, session_summary)
n_recordings = height(session_summary);
labels = recording_labels(session_summary);
subplot(3, 3, 1);
if ismember('n_episodes_with_sustained_sequence', session_summary.Properties.VariableNames) && ...
        any(isfinite(session_summary.n_episodes_with_sustained_sequence))
    sustained = session_summary.n_episodes_with_sustained_sequence;
    sustained(~isfinite(sustained)) = 0;
    counts = [session_summary.n_episodes(:), sustained(:)];
    if size(counts, 1) == 1
        bar(1, counts);
    else
        bar(counts);
    end
    legend({'Hilbert episodes', 'With a sustained sequence'}, ...
        'Location', 'best', 'FontSize', 8);
else
    bar(session_summary.n_episodes);
end
label_recordings(labels);
ylabel('Episodes');
title('Candidate episodes');

metric_fields = {'duration_s', 'n_trough_cycles', 'median_trough_period_s', ...
    'trough_period_cv'};
metric_fallbacks = {'duration_s', 'n_cycles', 'median_cycle_period_s', ...
    'cycle_period_iqr_s'};
metric_labels = {'Episode duration (s)', 'Trough cycles per episode', ...
    'Median trough period (s)', 'Trough-period CV'};
for metric_n = 1:numel(metric_fields)
    subplot(3, 3, metric_n + 1);
    hold on
    for recording_n = 1:n_recordings
        episodes = recordings(recording_n).episodes;
        if isempty(episodes)
            continue
        end
        values = episode_metric(episodes, metric_fields{metric_n}, ...
            metric_fallbacks{metric_n});
        plot(recording_n * ones(numel(values), 1), values, '.', ...
            'Color', [0.2 0.4 0.7], 'MarkerSize', 10);
        finite_values = values(isfinite(values));
        if isempty(finite_values)
            median_value = NaN;
        else
            median_value = median(finite_values);
        end
        plot(recording_n, median_value, 'ks', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 5);
    end
    label_recordings(labels);
    ylabel(metric_labels{metric_n});
    if metric_n == 1
        title('Squares mark the median across episodes');
    end
    if metric_n == 3
        plot(xlim, [0.2 0.2], 'Color', [0.4 0.4 0.4]);
        plot(xlim, [1 1], 'Color', [0.4 0.4 0.4]);
    end
    if metric_n == 4
        title('No CV cutoff is applied');
    end
    box off
end

subplot(3, 3, 6);
bar(session_summary.episode_time_fraction);
label_recordings(labels);
ylabel('Fraction of analyzed time');
title('Time inside candidate episodes');

[symmetry_rd, symmetry_pt] = all_symmetry(recordings);
subplot(3, 3, 7);
plot_symmetry_histogram(symmetry_rd, [0.75 0.3 0.2]);
xlabel('Rise-decay symmetry');
title('Rise / cycle period');
subplot(3, 3, 8);
plot_symmetry_histogram(symmetry_pt, [0.2 0.45 0.7]);
xlabel('Peak-trough symmetry');
title('Peak / (peak + previous trough)');
subplot(3, 3, 9);
plot_criteria_note(recordings(1).parameters);
end

function values = episode_metric(episodes, primary, fallback)
if isfield(episodes, primary)
    values = [episodes.(primary)]';
else
    values = [episodes.(fallback)]';
end
end

function [rdsym, ptsym] = all_symmetry(recordings)
rdsym = zeros(0, 1);
ptsym = zeros(0, 1);
if ~isfield(recordings, 'cycles')
    return
end
for recording_n = 1:numel(recordings)
    cycles = recordings(recording_n).cycles;
    if isempty(cycles)
        continue
    end
    rdsym = [rdsym; [cycles.rdsym]']; %#ok<AGROW>
    ptsym = [ptsym; [cycles.ptsym]']; %#ok<AGROW>
end
end

function plot_symmetry_histogram(values, color)
values = values(isfinite(values));
edges = 0:0.05:1;
if isempty(values)
    text(0.5, 0.5, 'No cycles', 'HorizontalAlignment', 'center');
    axis off
    return
end
histogram(values, edges, 'FaceColor', color, 'EdgeColor', 'none');
hold on
y_limits = ylim;
plot([0.5 0.5], y_limits, 'k--');
xlim([0 1]);
ylabel('Cycles');
box off
end

function plot_criteria_note(parameters)
axis off
if ~isfield(parameters, 'amp_consistency_threshold')
    text(0, 0.5, 'Cycle-consistency thresholds were not saved with this file.');
    return
end
if isfinite(parameters.cycle_amplitude_percentile)
    amplitude_line = sprintf('amplitude >= %gth percentile of the interval', ...
        parameters.cycle_amplitude_percentile);
else
    amplitude_line = 'no amplitude percentile';
end
note = strjoin({ ...
    sprintf('Sustained sequence (%s)', parameters.cycle_threshold_preset), ...
    sprintf('amp consistency >= %g', parameters.amp_consistency_threshold), ...
    sprintf('period consistency >= %g', parameters.period_consistency_threshold), ...
    sprintf('monotonicity >= %g', parameters.monotonicity_threshold), ...
    amplitude_line, ...
    sprintf('>= %g consecutive cycles', parameters.min_consecutive_cycles), ...
    'CV and symmetry are not gated'}, newline);
text(0, 1, note, 'VerticalAlignment', 'top', 'FontSize', 9, ...
    'Interpreter', 'none');
end

function labels = recording_labels(session_summary)
labels = cell(height(session_summary), 1);
for row_n = 1:height(session_summary)
    labels{row_n} = sprintf('%s p%g', session_summary.animal_name{row_n}, ...
        session_summary.probe_number(row_n));
end
end

function label_recordings(labels)
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
    'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');
xlim([0.5, numel(labels) + 0.5]);
box off
end

function [periods_s, kind] = recording_cycle_periods(recording)
if isfield(recording, 'cycles')
    kind = 'trough';
    if isempty(recording.cycles)
        periods_s = zeros(0, 1);
    else
        periods_s = [recording.cycles.period_s]';
    end
    return
end
kind = 'peak';
periods_s = zeros(0, 1);
for episode_n = 1:numel(recording.episodes)
    periods_s = [periods_s; recording.episodes(episode_n).cycle_periods_s(:)]; %#ok<AGROW>
end
end

function [periods_s, kind] = all_cycle_periods(recordings)
periods_s = zeros(0, 1);
kind = 'peak';
if isfield(recordings, 'cycles')
    kind = 'trough';
end
for recording_n = 1:numel(recordings)
    [recording_periods, ~] = recording_cycle_periods(recordings(recording_n));
    periods_s = [periods_s; recording_periods(:)]; %#ok<AGROW>
end
end

function print_session_lines(session_summary)
fprintf(['\nSession summary. Medians are taken across episodes. ' ...
    'Trough period and CV are the cycle-by-cycle measures.\n']);
for row_n = 1:height(session_summary)
    fprintf(['%-24s probe %4g ch %4g  episodes %4d  sustained %4g  ' ...
        'dur %6.2f s  trough cycles %5.1f  period %6.0f ms  CV %5.2f  ' ...
        'time %5.1f%%\n'], ...
        session_summary.animal_name{row_n}, ...
        session_summary.probe_number(row_n), ...
        session_summary.lpag_channel(row_n), ...
        session_summary.n_episodes(row_n), ...
        session_summary.n_episodes_with_sustained_sequence(row_n), ...
        session_summary.median_duration_s(row_n), ...
        session_summary.median_n_trough_cycles(row_n), ...
        1000 * session_summary.median_trough_period_s(row_n), ...
        session_summary.median_trough_period_cv(row_n), ...
        100 * session_summary.episode_time_fraction(row_n));
end
end

function figure_handle = new_figure(show_figures, size_pixels)
if show_figures
    visibility = 'on';
else
    visibility = 'off';
end
figure_handle = figure('Color', 'w', 'Visible', visibility, ...
    'Position', [80 80 size_pixels]);
set(figure_handle, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');
end

function name = file_token(text)
name = regexprep(char(text), '[^\w\-]+', '_');
name = regexprep(name, '_+', '_');
name = regexprep(name, '^_|_$', '');
if isempty(name)
    name = 'recording';
end
end
