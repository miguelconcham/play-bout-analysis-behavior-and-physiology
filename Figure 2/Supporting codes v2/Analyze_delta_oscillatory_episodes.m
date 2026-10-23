%% Analyze_delta_oscillatory_episodes
% Plots and the session summary for candidate continuous delta episodes.
%
% Run Estimate_delta_oscillatory_episodes.m first. This script reads
% Data/Analysis results/delta oscillatory episodes/delta_oscillatory_episodes.mat
% and writes:
%   session_summary.csv
%   figures/<recording>_probe<P>_channel<C>.png
%       minimally processed LFP, 1-5 Hz trace, Hilbert amplitude,
%       detected episode, positive peaks, and the cycle-period histogram
%   figures/cycle_period_histogram_pooled.png
%   figures/session_summary.png
%
% The histogram axis is 50-1000 ms. A 1-5 Hz cycle is 200-1000 ms
% (1/5 Hz to 1/1 Hz); that band is shaded. Periods near 50 ms are faster
% than the delta band and are shown so they are not mistaken for 5 Hz.
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
periods_s = all_cycle_periods(delta_episodes);
figure_handle = new_figure(show_figures, [900 500]);
parameters = delta_episodes(1).parameters;
plot_cycle_histogram(periods_s, parameters.delta_freq, ...
    parameters.histogram_ms, parameters.histogram_bin_ms);
title(cycle_histogram_title(periods_s, parameters.delta_freq, ...
    parameters.histogram_ms));
xlabel('Peak-to-peak cycle duration (ms)');
ylabel('Cycles');
print(figure_handle, fullfile(figure_folder, 'cycle_period_histogram_pooled.png'), ...
    '-dpng', '-r140');
if ~show_figures
    close(figure_handle);
end

%% Session summary
figure_handle = new_figure(show_figures, [1400 800]);
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
periods_s = recording_cycle_periods(recording);

subplot(4, 1, 1);
plot_trace_panel(example, 'lfp_minimal', [0.25 0.25 0.25]);
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
ylabel('Delta 1-5 Hz');

subplot(4, 1, 3);
plot_trace_panel(example, 'amplitude', [0.1 0.25 0.7]);
hold on
if ~isempty(example.time_s) && isfinite(example.threshold)
    plot(example.time_s([1 end]), [example.threshold example.threshold], 'k--');
end
ylabel('Hilbert amplitude');
xlabel('Recording time (s)');
if example.has_episode
    text(0.01, 0.9, sprintf(['full episode %.2f s, %d cycles, ' ...
        'median period %.0f ms, IQR %.0f ms'], ...
        example.episode_duration_s, example.n_cycles, ...
        1000 * example.median_cycle_period_s, ...
        1000 * example.cycle_period_iqr_s), ...
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
xlabel('Peak-to-peak cycle duration (ms)');
ylabel('Cycles');
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
subplot(2, 3, 1);
bar(session_summary.n_episodes);
label_recordings(labels);
ylabel('Episodes');
title('Candidate episodes');

metric_fields = {'duration_s', 'n_cycles', 'median_cycle_period_s', ...
    'cycle_period_iqr_s'};
metric_labels = {'Episode duration (s)', 'Cycles per episode', ...
    'Median cycle period (s)', 'Cycle-period IQR (s)'};
for metric_n = 1:numel(metric_fields)
    subplot(2, 3, metric_n + 1);
    hold on
    for recording_n = 1:n_recordings
        episodes = recordings(recording_n).episodes;
        if isempty(episodes)
            continue
        end
        values = [episodes.(metric_fields{metric_n})]';
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
    box off
end

subplot(2, 3, 6);
bar(session_summary.episode_time_fraction);
label_recordings(labels);
ylabel('Fraction of analyzed time');
title('Time inside candidate episodes');
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

function periods_s = recording_cycle_periods(recording)
periods_s = zeros(0, 1);
for episode_n = 1:numel(recording.episodes)
    periods_s = [periods_s; recording.episodes(episode_n).cycle_periods_s(:)]; %#ok<AGROW>
end
end

function periods_s = all_cycle_periods(recordings)
periods_s = zeros(0, 1);
for recording_n = 1:numel(recordings)
    periods_s = [periods_s; recording_cycle_periods(recordings(recording_n))]; %#ok<AGROW>
end
end

function print_session_lines(session_summary)
fprintf(['\nSession summary. Median duration, cycle count, period, and ' ...
    'period IQR are taken across episodes.\n']);
for row_n = 1:height(session_summary)
    fprintf(['%-24s probe %4g ch %4g  episodes %4d  ' ...
        'dur %6.2f s  cycles %5.1f  period %6.0f ms  IQR %5.0f ms  ' ...
        'time %5.1f%%\n'], ...
        session_summary.animal_name{row_n}, ...
        session_summary.probe_number(row_n), ...
        session_summary.lpag_channel(row_n), ...
        session_summary.n_episodes(row_n), ...
        session_summary.median_duration_s(row_n), ...
        session_summary.median_n_cycles(row_n), ...
        1000 * session_summary.median_cycle_period_s(row_n), ...
        1000 * session_summary.median_cycle_period_iqr_s(row_n), ...
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
