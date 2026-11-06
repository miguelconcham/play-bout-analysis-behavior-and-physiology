%% Analyse_play_session_lpag_irasa
% Plot session-wise mid-LPAG IRASA spectra saved by
% Estimate_play_session_lpag_irasa.
%
% The saved file keeps every window. This script averages those windows
% in linear power, one trace per Partners session. The raw trace is the
% Welch power. The fractal and periodic traces are the IRASA components.
% Periodic power can be negative. The printed exponent is the power-law
% fit of that session-mean fractal spectrum, over the saved fit range.
%
% Figures are written next to the loaded .mat file and left open.
% Set close_figures_after_saving to true to close each figure after it
% is saved.

%% Figure options
close_figures_after_saving = false;

saving_folder = '\\experimentfs.home.bccn-berlin.de\experiment\PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\play bout analysis behavior and physiology\Data\Analysis results\psth power by frequency and behavior';
input_file = fullfile(saving_folder, 'play_session_lpag_irasa.mat');

disp(['Loading ', input_file])
load(input_file, 'play_session_lpag_irasa', 'animal_names');
if isempty(play_session_lpag_irasa)
    error('The saved play_session_lpag_irasa structure is empty.');
end

n_recordings = numel(play_session_lpag_irasa);
n_sessions = sum(arrayfun(@(x) numel(x.sessions), play_session_lpag_irasa));
frequency = play_session_lpag_irasa(1).frequency(:).';
n_frequency = numel(frequency);

fractal_log_power = nan(n_sessions, n_frequency);
oscillatory_power = nan(n_sessions, n_frequency);
raw_frequency = play_session_lpag_irasa(1).sessions(1).raw.frequency(:).';
raw_log_power = nan(n_sessions, numel(raw_frequency));
exponent = nan(n_sessions, 1);
log_intercept = nan(n_sessions, 1);
broadband_log_power = nan(n_sessions, 1);
n_windows = nan(n_sessions, 1);
session_name = cell(n_sessions, 1);
recording_index = nan(n_sessions, 1);
session_index = nan(n_sessions, 1);
lpag_channel = nan(n_sessions, 1);

output_session_n = 0;
for recording_n = 1:n_recordings
    this_recording = play_session_lpag_irasa(recording_n);
    if numel(this_recording.frequency) ~= n_frequency || ...
            max(abs(this_recording.frequency(:).' - frequency)) > 1e-6
        error('Frequency vectors differ between saved IRASA sessions.');
    end
    for recording_session_n = 1:numel(this_recording.sessions)
        output_session_n = output_session_n + 1;
        this_session = this_recording.sessions(recording_session_n);
        session_name{output_session_n} = sprintf('%s session %d', ...
            this_recording.animal_name, recording_session_n);
        recording_index(output_session_n) = recording_n;
        session_index(output_session_n) = recording_session_n;
        lpag_channel(output_session_n) = this_recording.lpag_channel;
        n_windows(output_session_n) = this_session.n_windows;
        if numel(this_session.raw.frequency) ~= numel(raw_frequency) || ...
                max(abs(this_session.raw.frequency(:).' - raw_frequency)) > 1e-6
            error('Raw frequency vectors differ between saved IRASA sessions.');
        end
        raw_mean = mean(this_session.raw.welch, 1, 'omitnan');
        fractal_mean = mean(this_session.fractal, 1, 'omitnan').';
        periodic_mean = mean(this_session.periodic, 1, 'omitnan');
        session_fit = irasa_power_law_fit(struct( ...
            'frequency', this_session.frequency, ...
            'fractal', fractal_mean), this_session.fit_range);
        raw_log_power(output_session_n, :) = ...
            log10(max(raw_mean, realmin('double')));
        fractal_log_power(output_session_n, :) = ...
            log10(max(fractal_mean(:).', realmin('double')));
        oscillatory_power(output_session_n, :) = periodic_mean(:).';
        exponent(output_session_n) = session_fit.exponent;
        log_intercept(output_session_n) = session_fit.log_intercept;
        broadband_log_power(output_session_n) = session_fit.broadband_log_power;
    end
end

irasa_session_table = table(session_name, recording_index, session_index, ...
    lpag_channel, n_windows, exponent, log_intercept, broadband_log_power, ...
    'VariableNames', {'Session', 'RecordingIndex', 'SessionIndex', ...
    'LPAGChannel', 'NWindows', 'Exponent', 'LogIntercept', ...
    'BroadbandLogPower'});
disp(irasa_session_table)

raw_figure = figure('Name', 'Play-session LPAG raw and IRASA spectra');
subplot(1, 2, 1)
session_lines = plot(raw_frequency, raw_log_power.', ...
    'Color', [0.8 0.8 0.8]);
hold on
mean_line = plot(raw_frequency, mean(raw_log_power, 1, 'omitnan'), ...
    'k', 'LineWidth', 2);
xlabel('Frequency (Hz)')
ylabel('log_{10} power')
title('Raw Welch power, mean of windows')
legend([session_lines(1), mean_line], ...
    {'Partners sessions', 'Mean'}, 'Location', 'best')
box off

subplot(1, 2, 2)
session_lines = plot(frequency, fractal_log_power.', ...
    'Color', [0.8 0.8 0.8]);
hold on
mean_line = plot(frequency, mean(fractal_log_power, 1, 'omitnan'), ...
    'k', 'LineWidth', 2);
xlabel('Frequency (Hz)')
ylabel('log_{10} power')
title(sprintf('IRASA fractal, mean exponent %.2f', ...
    mean(exponent, 'omitnan')))
legend([session_lines(1), mean_line], ...
    {'Partners sessions', 'Mean'}, 'Location', 'best')
box off
save_irasa_figure(raw_figure, saving_folder, ...
    'play_session_lpag_irasa_raw_and_fractal', close_figures_after_saving);

periodic_figure = figure('Name', 'Play-session LPAG IRASA periodic spectra');
hold on
session_lines = plot(frequency, oscillatory_power.', ...
    'Color', [0.8 0.8 0.8]);
mean_line = plot(frequency, mean(oscillatory_power, 1, 'omitnan'), ...
    'k', 'LineWidth', 2);
yline(0, ':k');
xline(1, '--', 'Color', [0.4 0.4 0.4]);
xline(5, '--', 'Color', [0.4 0.4 0.4]);
xline(10, '--', 'Color', [0.4 0.4 0.4]);
xlabel('Frequency (Hz)')
ylabel('Periodic power')
title('IRASA periodic spectrum, mean of windows')
legend([session_lines(1), mean_line], ...
    {'Partners sessions', 'Mean'}, 'Location', 'best')
box off
save_irasa_figure(periodic_figure, saving_folder, ...
    'play_session_lpag_irasa_periodic', close_figures_after_saving);

irasa_results = struct();
irasa_results.frequency = frequency(:);
irasa_results.raw_frequency = raw_frequency(:);
irasa_results.session_table = irasa_session_table;
irasa_results.raw_log_power = raw_log_power;
irasa_results.fractal_log_power = fractal_log_power;
irasa_results.periodic_power = oscillatory_power;
irasa_results.animal_names = animal_names;
disp('Play-session LPAG IRASA analysis complete.')

function save_irasa_figure(figure_handle, saving_folder, figure_name, close_after_saving)
figure_base = fullfile(saving_folder, figure_name);
savefig(figure_handle, [figure_base, '.fig']);
saveas(figure_handle, [figure_base, '.png']);
disp(['Saved ', figure_base, '.fig'])
if close_after_saving
    close(figure_handle);
end
end
