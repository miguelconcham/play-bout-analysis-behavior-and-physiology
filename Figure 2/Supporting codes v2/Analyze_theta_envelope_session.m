function Analyze_theta_envelope_session(result_file)
%ANALYZE_THETA_ENVELOPE_SESSION  Animal-level spectrum and coherence tests.
%
% The plotted curves are surrogate z-scores on the 0.2 Hz grid. The printed
% test is a sign-flip of the animal means inside the confirmatory band, plus
% a cluster sign-flip on the 0.5 Hz grid.

if nargin < 1 || isempty(result_file)
    options = theta_envelope_delta_options();
    result_file = fullfile(options.output_folder, options.output_filename);
end
loaded = load(result_file, 'recordings', 'parameters');
recordings = loaded.recordings;
if isempty(recordings)
    error('Analyze_theta_envelope_session:empty', ...
        'No recordings in %s.', result_file);
end
if isfield(loaded, 'parameters') && isfield(loaded.parameters, 'n_group_flip')
    n_flip = loaded.parameters.n_group_flip;
else
    n_flip = 5000;
end

grouped = theta_envelope_group_sessions(recordings);
power_scalar = row_mean(grouped.power_z_test, grouped.confirmatory_mask);
coherence_scalar = row_mean(grouped.coherence_z_test, grouped.confirmatory_mask);
power_flip = theta_envelope_signflip(power_scalar, n_flip);
coherence_flip = theta_envelope_signflip(coherence_scalar, n_flip);
power_cluster = theta_envelope_cluster_signflip(grouped.power_z_test, n_flip);
coherence_cluster = theta_envelope_cluster_signflip(grouped.coherence_z_test, n_flip);
power_overlap = overlapping_cluster(power_cluster, grouped.confirmatory_mask);
coherence_overlap = overlapping_cluster(coherence_cluster, grouped.confirmatory_mask);

fprintf(['Session theta-envelope test. %d animals.\n' ...
    'Confirmatory power z mean %.3f, sign-flip p %.4f.\n' ...
    'Confirmatory coherence z mean %.3f, sign-flip p %.4f.\n'], ...
    power_flip.n_animals, power_flip.mean, power_flip.p_greater, ...
    coherence_flip.mean, coherence_flip.p_greater);
print_cluster('Power', power_overlap, grouped.frequencies_test);
print_cluster('Coherence', coherence_overlap, grouped.frequencies_test);
fprintf('Animal envelope peaks (Hz): %s\n', ...
    sprintf('%.2f ', grouped.peak_power_hz));
fprintf('Animal coherence peaks (Hz): %s\n', ...
    sprintf('%.2f ', grouped.peak_coherence_hz));

figure('Name', 'Theta envelope session', 'Color', 'w');
subplot(2, 1, 1)
plot_animal_curves(grouped.frequencies_fine, grouped.band_power_fine, ...
    grouped.peak_power_hz);
xlabel('Frequency (Hz)')
ylabel('Envelope power')
title(sprintf('Envelope spectrum, confirmatory z p = %.4f', power_flip.p_greater))
subplot(2, 1, 2)
plot_animal_curves(grouped.frequencies_fine, grouped.band_coherence_fine, ...
    grouped.peak_coherence_hz);
xlabel('Frequency (Hz)')
ylabel('Coherence')
title(sprintf('Envelope-LFP coherence, confirmatory z p = %.4f', coherence_flip.p_greater))

summary = table(grouped.animals, grouped.n_recordings, grouped.peak_power_hz, ...
    grouped.peak_coherence_hz, power_scalar, coherence_scalar, ...
    'VariableNames', {'animal', 'n_recordings', 'peak_power_hz', ...
    'peak_coherence_hz', 'confirmatory_power_z', 'confirmatory_coherence_z'});
[folder, ~, ~] = fileparts(result_file);
writetable(summary, fullfile(folder, 'theta_envelope_session_animals.csv'));
disp(summary)
end

function value = row_mean(matrix, mask)
value = nan(size(matrix, 1), 1);
for row = 1:size(matrix, 1)
    selected = matrix(row, logical(mask));
    selected = selected(isfinite(selected));
    if ~isempty(selected)
        value(row) = mean(selected);
    end
end
end

function chosen = overlapping_cluster(result, mask)
chosen = [];
best = -inf;
if ~isfield(result, 'clusters')
    return
end
for c = 1:numel(result.clusters)
    idx = result.clusters(c).start:result.clusters(c).stop;
    if any(mask(idx)) && result.clusters(c).sum > best
        chosen = result.clusters(c);
        best = chosen.sum;
    end
end
end

function print_cluster(label, cluster, frequencies)
if isempty(cluster)
    fprintf('%s: no positive cluster overlaps 2-3 Hz.\n', label);
    return
end
fprintf('%s cluster %.2f-%.2f Hz, sum %.3f, cluster p %.4f.\n', ...
    label, frequencies(cluster.start), frequencies(cluster.stop), ...
    cluster.sum, cluster.p);
end

function plot_animal_curves(frequencies, animal_z, peaks)
y_values = animal_z(isfinite(animal_z));
if isempty(y_values)
    span = [-1 1];
else
    span = [min(y_values), max(y_values)];
    if span(1) == span(2)
        span = span + [-1 1];
    end
end
hold on
shade = patch([2 3 3 2], [span(1) span(1) span(2) span(2)], [0.85 0.9 1], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');
for row = 1:size(animal_z, 1)
    plot(frequencies, animal_z(row, :), 'Color', [0.7 0.7 0.7])
end
plot(frequencies, finite_mean(animal_z), 'k', 'LineWidth', 1.5)
uistack(shade, 'bottom')
hold off
xlim([frequencies(1), frequencies(end)])
if any(isfinite(peaks))
    text(0.02, 0.9, sprintf('peaks %.2f Hz', mean(peaks(isfinite(peaks)))), ...
        'Units', 'normalized')
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
