%% Analyze_theta_envelope_session
% Animal-level envelope spectrum and envelope-LFP coherence.
%
% Run Estimate_theta_envelope_delta.m first. Then run this file from the top,
% or run the sections in order. Each plotting section opens its own figure.
%
% Gray line: one animal. Black line: mean across animals. Blue band: 2-3 Hz.
% The curves are envelope power and coherence. The printed p-values are the
% animal sign-flip of the 2-3 Hz surrogate z-score.

%% Paths
this_file = mfilename('fullpath');
if isempty(this_file)
    this_file = ['\\experimentfs.home.bccn-berlin.de\experiment\' ...
        'PlayNeuralData\NPX-OPTO PLAY NMM\PlayBout Analysis\' ...
        'play bout analysis behavior and physiology\Figure 2\' ...
        'Supporting codes\Theta envelope\Analyze_theta_envelope_session.m'];
end
repo_root = fileparts(fileparts(fileparts(fileparts(this_file))));
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
if isfield(loaded, 'parameters') && isfield(loaded.parameters, 'n_group_flip')
    n_flip = loaded.parameters.n_group_flip;
else
    n_flip = 5000;
end
grouped = theta_envelope_group_sessions(recordings);

%% Confirmatory sign-flip and frequency clusters
n_animals = size(grouped.power_z_test, 1);
power_scalar = nan(n_animals, 1);
coherence_scalar = nan(n_animals, 1);
for animal_n = 1:n_animals
    power_values = grouped.power_z_test(animal_n, logical(grouped.confirmatory_mask));
    power_values = power_values(isfinite(power_values));
    coherence_values = grouped.coherence_z_test(animal_n, logical(grouped.confirmatory_mask));
    coherence_values = coherence_values(isfinite(coherence_values));
    if ~isempty(power_values)
        power_scalar(animal_n) = mean(power_values);
    end
    if ~isempty(coherence_values)
        coherence_scalar(animal_n) = mean(coherence_values);
    end
end
power_flip = theta_envelope_signflip(power_scalar, n_flip);
coherence_flip = theta_envelope_signflip(coherence_scalar, n_flip);
power_cluster = theta_envelope_cluster_signflip(grouped.power_z_test, n_flip);
coherence_cluster = theta_envelope_cluster_signflip(grouped.coherence_z_test, n_flip);

power_overlap = [];
best_sum = -inf;
if isfield(power_cluster, 'clusters')
    for cluster_n = 1:numel(power_cluster.clusters)
        cluster = power_cluster.clusters(cluster_n);
        in_band = any(grouped.confirmatory_mask(cluster.start:cluster.stop));
        if in_band && cluster.sum > best_sum
            power_overlap = cluster;
            best_sum = cluster.sum;
        end
    end
end
coherence_overlap = [];
best_sum = -inf;
if isfield(coherence_cluster, 'clusters')
    for cluster_n = 1:numel(coherence_cluster.clusters)
        cluster = coherence_cluster.clusters(cluster_n);
        in_band = any(grouped.confirmatory_mask(cluster.start:cluster.stop));
        if in_band && cluster.sum > best_sum
            coherence_overlap = cluster;
            best_sum = cluster.sum;
        end
    end
end

fprintf(['Session theta-envelope test. %d animals.\n' ...
    'Confirmatory power z mean %.3f, sign-flip p %.4f.\n' ...
    'Confirmatory coherence z mean %.3f, sign-flip p %.4f.\n'], ...
    power_flip.n_animals, power_flip.mean, power_flip.p_greater, ...
    coherence_flip.mean, coherence_flip.p_greater);
if isempty(power_overlap)
    fprintf('Power: no positive cluster overlaps 2-3 Hz.\n');
else
    fprintf('Power cluster %.2f-%.2f Hz, sum %.3f, cluster p %.4f.\n', ...
        grouped.frequencies_test(power_overlap.start), ...
        grouped.frequencies_test(power_overlap.stop), ...
        power_overlap.sum, power_overlap.p);
end
if isempty(coherence_overlap)
    fprintf('Coherence: no positive cluster overlaps 2-3 Hz.\n');
else
    fprintf('Coherence cluster %.2f-%.2f Hz, sum %.3f, cluster p %.4f.\n', ...
        grouped.frequencies_test(coherence_overlap.start), ...
        grouped.frequencies_test(coherence_overlap.stop), ...
        coherence_overlap.sum, coherence_overlap.p);
end
fprintf('Animal envelope peaks (Hz): %s\n', sprintf('%.2f ', grouped.peak_power_hz));
fprintf('Animal coherence peaks (Hz): %s\n', sprintf('%.2f ', grouped.peak_coherence_hz));

%% Envelope spectrum
y_values = grouped.band_power_fine(isfinite(grouped.band_power_fine));
if isempty(y_values)
    span = [-1 1];
else
    span = [min(y_values), max(y_values)];
    if span(1) == span(2)
        span = span + [-1 1];
    end
end
power_mean = nan(1, size(grouped.band_power_fine, 2));
for column = 1:numel(power_mean)
    finite = grouped.band_power_fine(isfinite(grouped.band_power_fine(:, column)), column);
    if ~isempty(finite)
        power_mean(column) = mean(finite);
    end
end

figure('Name', 'Theta envelope spectrum', 'Color', 'w');
hold on
patch([2 3 3 2], [span(1) span(1) span(2) span(2)], [0.85 0.9 1], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
animal_line = plot(nan, nan, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Animal');
for animal_n = 1:size(grouped.band_power_fine, 1)
    plot(grouped.frequencies_fine, grouped.band_power_fine(animal_n, :), ...
        'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
end
mean_line = plot(grouped.frequencies_fine, power_mean, 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Mean');
hold off
xlim([grouped.frequencies_fine(1), grouped.frequencies_fine(end)])
xlabel('Frequency (Hz)')
ylabel('Envelope power')
title(sprintf('Envelope spectrum, confirmatory z p = %.4f', power_flip.p_greater))
legend([animal_line, mean_line], 'Location', 'northeast')
finite_peaks = grouped.peak_power_hz(isfinite(grouped.peak_power_hz));
if ~isempty(finite_peaks)
    text(0.02, 0.9, sprintf('mean peak %.2f Hz', mean(finite_peaks)), ...
        'Units', 'normalized')
end

%% Envelope-LFP coherence
y_values = grouped.band_coherence_fine(isfinite(grouped.band_coherence_fine));
if isempty(y_values)
    span = [-1 1];
else
    span = [min(y_values), max(y_values)];
    if span(1) == span(2)
        span = span + [-1 1];
    end
end
coherence_mean = nan(1, size(grouped.band_coherence_fine, 2));
for column = 1:numel(coherence_mean)
    finite = grouped.band_coherence_fine(isfinite(grouped.band_coherence_fine(:, column)), column);
    if ~isempty(finite)
        coherence_mean(column) = mean(finite);
    end
end

figure('Name', 'Theta envelope coherence', 'Color', 'w');
hold on
patch([2 3 3 2], [span(1) span(1) span(2) span(2)], [0.85 0.9 1], ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
animal_line = plot(nan, nan, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Animal');
for animal_n = 1:size(grouped.band_coherence_fine, 1)
    plot(grouped.frequencies_fine, grouped.band_coherence_fine(animal_n, :), ...
        'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
end
mean_line = plot(grouped.frequencies_fine, coherence_mean, 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Mean');
hold off
xlim([grouped.frequencies_fine(1), grouped.frequencies_fine(end)])
xlabel('Frequency (Hz)')
ylabel('Coherence')
title(sprintf('Envelope-LFP coherence, confirmatory z p = %.4f', coherence_flip.p_greater))
legend([animal_line, mean_line], 'Location', 'northeast')
finite_peaks = grouped.peak_coherence_hz(isfinite(grouped.peak_coherence_hz));
if ~isempty(finite_peaks)
    text(0.02, 0.9, sprintf('mean peak %.2f Hz', mean(finite_peaks)), ...
        'Units', 'normalized')
end

%% Save animal table
summary = table(grouped.animals, grouped.n_recordings, grouped.peak_power_hz, ...
    grouped.peak_coherence_hz, power_scalar, coherence_scalar, ...
    'VariableNames', {'animal', 'n_recordings', 'peak_power_hz', ...
    'peak_coherence_hz', 'confirmatory_power_z', 'confirmatory_coherence_z'});
[folder, ~, ~] = fileparts(result_file);
writetable(summary, fullfile(folder, 'theta_envelope_session_animals.csv'));
disp(summary)
