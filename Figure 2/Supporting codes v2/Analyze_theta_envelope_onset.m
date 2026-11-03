%% Analyze_theta_envelope_onset
% Animal and single-behavior test of the theta envelope before play onset.
%
% Run Estimate_theta_envelope_delta.m first. Then run this file from the top,
% or run the sections in order. Each plotting section opens its own figure.
%
% Gray line: one animal with at least 15 play onsets. Black line: mean of
% those animals. Dashed line: the same mean using only onsets whose preceding
% 2 s contain no other play behavior. Blue band and the red line mark
% -2 to 0 s and play onset.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes v2\Analyze_theta_envelope_onset.m'];
end
repo_root = fileparts(fileparts(fileparts(this_file)));
run(fullfile(repo_root, 'add_repo_paths.m'));

options = theta_envelope_delta_options();
result_file = fullfile(options.output_folder, options.output_filename);

%% Load
fprintf('Loading %s\n', result_file);
loaded = load(result_file, 'recordings', 'parameters');
recordings = loaded.recordings;
if isempty(recordings)
    error('No recordings in %s.', result_file);
end
if isfield(loaded, 'parameters') && isfield(loaded.parameters, 'min_events')
    min_events = loaded.parameters.min_events;
else
    min_events = 15;
end
grouped = theta_envelope_group_onsets(recordings, min_events);
included = find(grouped.included);

%% Print animal and group tests
fprintf(['Onset-locked theta envelope. %d animals, %d with at least %d onsets.\n' ...
    'Confirmatory locked amplitude p %.4f.\n' ...
    'Confirmatory phase consistency p %.4f.\n' ...
    'Leave-one-out correlation p %.4f.\n' ...
    'Single-behavior phase fraction p %.4f.\n'], ...
    numel(grouped.animals), sum(grouped.included), min_events, ...
    grouped.locked_p, grouped.ppc_p, grouped.loo_p, grouped.phase_fraction_p);
for animal_n = 1:numel(grouped.animals)
    fprintf(['  %s: %d onsets, included %d, locked p %.4f, ' ...
        'phase fraction %.2f\n'], ...
        grouped.animals{animal_n}, grouped.n_trials(animal_n), ...
        grouped.included(animal_n), grouped.animal_locked_p(animal_n), ...
        grouped.phase_fraction(animal_n));
end

%% Waveform aligned to play onset
y_values = grouped.mean_epoch(grouped.included, :);
y_values = y_values(isfinite(y_values));
if isempty(y_values)
    span = [-1 1];
else
    span = [min(y_values), max(y_values)];
    if span(1) == span(2)
        span = span + [-1 1];
    end
end

figure('Name', 'Theta envelope at play onset', 'Color', 'w');
hold on
patch([-2 0 0 -2], [span(1) span(1) span(2) span(2)], [0.85 0.9 1], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
animal_line = plot(nan, nan, 'Color', [0.75 0.75 0.75], 'DisplayName', 'Animal');
for animal_n = included(:).'
    plot(grouped.epoch_time, grouped.mean_epoch(animal_n, :), ...
        'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
end
mean_line = plot(grouped.epoch_time, grouped.grand_mean_epoch, 'k', ...
    'LineWidth', 1.5, 'DisplayName', 'Mean');
clean_line = plot(grouped.epoch_time, grouped.grand_clean_epoch, 'k--', ...
    'LineWidth', 1.2, 'DisplayName', 'Clean pre-window');
onset_line = plot([0 0], span, 'r', 'DisplayName', 'Play onset');
hold off
xlabel('Time from play onset (s)')
ylabel('Theta envelope (z)')
title('Animal means')
legend([animal_line, mean_line, clean_line, onset_line], 'Location', 'northeast')

%% Spectrum of the average waveform
fine_frequencies = recordings(1).frequencies_fine;
fine_amplitude = nan(numel(grouped.animals), numel(fine_frequencies));
pre = grouped.epoch_time >= -2 & grouped.epoch_time <= 0;
pre_time = grouped.epoch_time(pre);
for animal_n = 1:numel(grouped.animals)
    snippet = grouped.mean_epoch(animal_n, pre);
    if isempty(snippet) || any(~isfinite(snippet))
        continue
    end
    coefficient = theta_envelope_fourier(snippet, pre_time, fine_frequencies);
    fine_amplitude(animal_n, :) = abs(coefficient);
end
if any(grouped.included)
    fine_mean = nan(1, numel(fine_frequencies));
    selected = fine_amplitude(grouped.included, :);
    for column = 1:numel(fine_mean)
        finite = selected(isfinite(selected(:, column)), column);
        if ~isempty(finite)
            fine_mean(column) = mean(finite);
        end
    end
else
    fine_mean = nan(1, numel(fine_frequencies));
end

figure('Name', 'Spectrum of the pre-play average', 'Color', 'w');
hold on
animal_line = plot(nan, nan, 'Color', [0.75 0.75 0.75], 'DisplayName', 'Animal');
for animal_n = included(:).'
    plot(fine_frequencies, fine_amplitude(animal_n, :), ...
        'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
end
mean_line = plot(fine_frequencies, fine_mean, 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Mean');
hold off
xlabel('Frequency (Hz)')
ylabel('Amplitude of the average')
title(sprintf('0.2 Hz shape of -2 to 0 s, 2-3 Hz locked p = %.4f', grouped.locked_p))
legend([animal_line, mean_line], 'Location', 'northeast')

%% Phase consistency across behaviors
if any(grouped.included)
    ppc_mean = nan(1, numel(grouped.frequencies));
    selected = grouped.ppc(grouped.included, :);
    for column = 1:numel(ppc_mean)
        finite = selected(isfinite(selected(:, column)), column);
        if ~isempty(finite)
            ppc_mean(column) = mean(finite);
        end
    end
else
    ppc_mean = nan(1, numel(grouped.frequencies));
end

figure('Name', 'Phase consistency before play', 'Color', 'w');
hold on
animal_line = plot(nan, nan, 'Color', [0.75 0.75 0.75], 'DisplayName', 'Animal');
for animal_n = included(:).'
    plot(grouped.frequencies, grouped.ppc(animal_n, :), ...
        'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
end
mean_line = plot(grouped.frequencies, ppc_mean, 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Mean');
hold off
xlabel('Frequency (Hz)')
ylabel('Pairwise phase consistency')
title(sprintf('2-3 Hz consistency p = %.4f', grouped.ppc_p))
legend([animal_line, mean_line], 'Location', 'northeast')

%% Single-behavior correlations
loo = [];
for recording_n = 1:numel(recordings)
    loo = [loo; recordings(recording_n).onset.loo_r(:)]; %#ok<AGROW>
end
loo = loo(isfinite(loo));

figure('Name', 'Single behaviors before play', 'Color', 'w');
if isempty(loo)
    text(0.5, 0.5, 'No single-behavior correlations', 'HorizontalAlignment', 'center')
    axis off
else
    histogram(loo, 20, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none')
    xlabel('Leave-one-out correlation')
    ylabel('Behaviors')
end
title(sprintf('Single behaviors, p = %.4f', grouped.loo_p))

%% Save animal table
summary = table(grouped.animals, grouped.n_trials, grouped.included, ...
    grouped.confirmatory_locked, grouped.animal_locked_p, ...
    grouped.confirmatory_ppc, grouped.animal_ppc_p, ...
    grouped.loo_mean_r, grouped.animal_loo_p, ...
    grouped.phase_fraction, grouped.animal_phase_p, ...
    'VariableNames', {'animal', 'n_onsets', 'included', ...
    'locked_amplitude', 'locked_p', 'phase_consistency', 'phase_consistency_p', ...
    'leave_one_out_r', 'leave_one_out_p', 'phase_fraction', 'phase_fraction_p'});
[folder, ~, ~] = fileparts(result_file);
writetable(summary, fullfile(folder, 'theta_envelope_onset_animals.csv'));
disp(summary)
