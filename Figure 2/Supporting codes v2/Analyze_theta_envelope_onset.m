function Analyze_theta_envelope_onset(result_file)
%ANALYZE_THETA_ENVELOPE_ONSET  Animal and single-behavior test of pre-play modulation.
%
% The waveform plot is the animal mean of the theta envelope aligned to each
% play onset. The gray line is every included animal. The black line is their
% mean. The dashed line is the same average restricted to onsets whose
% preceding 2 s contain no other play behavior.
%
% The printed p-values compare the confirmatory 2-3 Hz band with the
% onset-shift null. Locked amplitude is the size of that band in the animal
% average. Pairwise phase consistency and the phase fraction ask whether
% single behaviors share that phase. Leave-one-out correlation asks whether
% a single pre-onset window resembles the other windows from that recording.

if nargin < 1 || isempty(result_file)
    options = theta_envelope_delta_options();
    result_file = fullfile(options.output_folder, options.output_filename);
end
loaded = load(result_file, 'recordings', 'parameters');
recordings = loaded.recordings;
if isempty(recordings)
    error('Analyze_theta_envelope_onset:empty', ...
        'No recordings in %s.', result_file);
end
if isfield(loaded, 'parameters') && isfield(loaded.parameters, 'min_events')
    min_events = loaded.parameters.min_events;
else
    min_events = 15;
end
grouped = theta_envelope_group_onsets(recordings, min_events);

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

figure('Name', 'Theta envelope at play onset', 'Color', 'w');
subplot(2, 2, 1)
hold on
included = find(grouped.included);
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
shade = patch([-2 0 0 -2], [span(1) span(1) span(2) span(2)], [0.85 0.9 1], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
for animal_n = included(:).'
    plot(grouped.epoch_time, grouped.mean_epoch(animal_n, :), 'Color', [0.75 0.75 0.75])
end
plot(grouped.epoch_time, grouped.grand_mean_epoch, 'k', 'LineWidth', 1.5)
plot(grouped.epoch_time, grouped.grand_clean_epoch, 'k--', 'LineWidth', 1.2)
plot([0 0], span, 'r')
uistack(shade, 'bottom')
hold off
xlabel('Time from play onset (s)')
ylabel('Theta envelope (z)')
title('Animal means')

subplot(2, 2, 2)
[fine_frequencies, fine_amplitude] = average_waveform_spectrum( ...
    grouped, recordings(1).frequencies_fine);
hold on
for animal_n = included(:).'
    plot(fine_frequencies, fine_amplitude(animal_n, :), 'Color', [0.75 0.75 0.75])
end
if any(grouped.included)
    plot(fine_frequencies, finite_mean(fine_amplitude(grouped.included, :)), ...
        'k', 'LineWidth', 1.5)
end
hold off
xlabel('Frequency (Hz)')
ylabel('Amplitude of the average')
title(sprintf('0.2 Hz shape, 2-3 Hz locked p = %.4f', grouped.locked_p))

subplot(2, 2, 3)
hold on
for animal_n = included(:).'
    plot(grouped.frequencies, grouped.ppc(animal_n, :), 'Color', [0.75 0.75 0.75])
end
if any(grouped.included)
    plot(grouped.frequencies, finite_mean(grouped.ppc(grouped.included, :)), ...
        'k', 'LineWidth', 1.5)
end
hold off
xlabel('Frequency (Hz)')
ylabel('Pairwise phase consistency')
title(sprintf('2-3 Hz consistency p = %.4f', grouped.ppc_p))

subplot(2, 2, 4)
loo = [];
for recording_n = 1:numel(recordings)
    loo = [loo; recordings(recording_n).onset.loo_r(:)]; %#ok<AGROW>
end
loo = loo(isfinite(loo));
if isempty(loo)
    text(0.5, 0.5, 'No single-behavior correlations', 'HorizontalAlignment', 'center')
else
    histogram(loo, 20, 'FaceColor', [0.2 0.2 0.2], 'EdgeColor', 'none')
    xlabel('Leave-one-out correlation')
    ylabel('Behaviors')
end
title(sprintf('Single behaviors, p = %.4f', grouped.loo_p))

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
end

function [frequencies, amplitude] = average_waveform_spectrum(grouped, frequencies)
n_animals = numel(grouped.animals);
amplitude = nan(n_animals, numel(frequencies));
time = grouped.epoch_time(:);
pre = time >= -2 & time <= 0;
if ~any(pre)
    return
end
pre_time = time(pre);
for animal_n = 1:n_animals
    snippet = grouped.mean_epoch(animal_n, pre.');
    if isempty(snippet) || any(~isfinite(snippet))
        continue
    end
    coefficient = theta_envelope_fourier(snippet(:).', pre_time, frequencies);
    amplitude(animal_n, :) = abs(coefficient);
end
end

function value = finite_mean(matrix)
value = nan(1, size(matrix, 2));
for column = 1:size(matrix, 2)
    finite = matrix(isfinite(matrix(:, column)), column);
    if ~isempty(finite)
        value(column) = mean(finite);
    end
end
end
