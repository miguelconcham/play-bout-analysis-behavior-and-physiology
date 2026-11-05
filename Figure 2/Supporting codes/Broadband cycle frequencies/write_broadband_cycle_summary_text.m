function write_broadband_cycle_summary_text(recordings, output_file)
%WRITE_BROADBAND_CYCLE_SUMMARY_TEXT  Plain-text cycle-frequency results.

file_id = fopen(output_file, 'w');
if file_id < 0
    error('write_broadband_cycle_summary_text:file', ...
        'Could not write %s.', output_file);
end
cleaner = onCleanup(@() fclose(file_id)); %#ok<NASGU>
parameters = recordings(1).parameters;

fprintf(file_id, 'Broadband cycle-frequency summary\n');
fprintf(file_id, '=================================\n\n');
fprintf(file_id, ['Question: on a %.4g Hz low-pass of mid-LPAG, do ' ...
    'trough-to-trough cycle frequencies form one broad cluster or ' ...
    'separated clusters?\n\n'], parameters.lowpass_hz);
fprintf(file_id, ['No delta or theta bandpass was used to place the cycles. ' ...
    'A non-sinusoidal slow cycle keeps the slow trough-to-trough period ' ...
    'even when its spectrum contains a harmonic.\n\n']);
fprintf(file_id, ['Accepted cycles lie between %.4g and %.4g Hz and pass ' ...
    'amplitude consistency >= %.3g, period consistency >= %.3g, ' ...
    'monotonicity >= %.3g, and amplitude at or above the %.4gth percentile ' ...
    'of that interval.\n\n'], ...
    parameters.frequency_limits_hz(1), parameters.frequency_limits_hz(2), ...
    parameters.amplitude_consistency_threshold, ...
    parameters.period_consistency_threshold, ...
    parameters.monotonicity_threshold, parameters.amplitude_percentile);

animals = unique({recordings.animal_id}, 'stable');
fprintf(file_id, 'Recordings: %d channels from %d animals.\n\n', ...
    numel(recordings), numel(animals));
fprintf(file_id, ['animal | recording | probe | channel | analyzed_s | ' ...
    'candidate_cycles | accepted_cycles | median_Hz | iqr_Hz | ' ...
    'fraction_0.5_to_2 | fraction_2_to_7 | fraction_7_to_12 | ' ...
    'fraction_12_to_20 | median_rise_decay_symmetry\n']);
for recording_n = 1:numel(recordings)
    recording = recordings(recording_n);
    frequency = accepted_frequency(recording);
    fprintf(file_id, ['%s | %s | %g | %g | %.4g | %g | %g | %.3g | %.3g | ' ...
        '%.3g | %.3g | %.3g | %.3g | %.3g\n'], ...
        recording.animal_id, recording.animal_name, ...
        recording.lpag_probe_number, recording.lpag_channel, ...
        recording.analyzed_seconds, recording.n_candidate_cycles, ...
        recording.n_accepted_cycles, finite_median(frequency), ...
        iqr_or_nan(frequency), band_fraction(frequency, [0.5 2]), ...
        band_fraction(frequency, [2 7]), band_fraction(frequency, [7 12]), ...
        band_fraction(frequency, [12 20]), ...
        finite_median(recording.rdsym(recording.accepted)));
end

edges = parameters.histogram_edges_hz;
[density, centers] = equal_animal_histogram(recordings, animals, edges);
fprintf(file_id, ['\nEqual-animal histogram. Each animal is normalized to ' ...
    'its own accepted-cycle count, then animals are averaged.\n']);
fprintf(file_id, 'bin_start_Hz | bin_end_Hz | mean_animal_fraction\n');
for bin_n = 1:numel(density)
    fprintf(file_id, '%.3g | %.3g | %.6g\n', edges(bin_n), edges(bin_n + 1), ...
        density(bin_n));
end
modes = histogram_modes(centers, density);
fprintf(file_id, '\nHistogram modes, local peaks at least half the maximum:\n');
if isempty(modes)
    fprintf(file_id, 'No local mode.\n');
else
    for mode_n = 1:size(modes, 1)
        fprintf(file_id, '%.3g Hz, mean animal fraction %.6g\n', ...
            modes(mode_n, 1), modes(mode_n, 2));
    end
end
if size(modes, 1) >= 2
    fprintf(file_id, ['%d separated local modes are present. Inspect the ' ...
        'valley between them before treating them as separate rhythms.\n'], ...
        size(modes, 1));
else
    fprintf(file_id, ['One dominant histogram mode is present. A broad ' ...
        'single mode favors one rhythm with a spectral shoulder rather ' ...
        'than two cycle populations.\n']);
end
fprintf(file_id, ['\nChannels from one animal are repeated measures. The ' ...
    'equal-animal histogram gives each animal one weight. This text does ' ...
    'not report a significance test.\n']);
end

function value = finite_median(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function value = finite_mean_rows(matrix)
value = nan(1, size(matrix, 2));
for column_n = 1:size(matrix, 2)
    column = matrix(:, column_n);
    column = column(isfinite(column));
    if ~isempty(column)
        value(column_n) = mean(column);
    end
end
end

function frequency = accepted_frequency(recording)
if isempty(recording.frequency_hz) || isempty(recording.accepted)
    frequency = zeros(0, 1);
    return
end
frequency = recording.frequency_hz(recording.accepted);
end

function value = iqr_or_nan(frequency)
if numel(frequency) < 2
    value = NaN;
    return
end
value = linear_percentile(frequency, 75) - linear_percentile(frequency, 25);
end

function fraction = band_fraction(frequency, band)
if isempty(frequency)
    fraction = NaN;
    return
end
fraction = mean(frequency >= band(1) & frequency < band(2));
end

function [density, centers] = equal_animal_histogram(recordings, animals, edges)
n_bins = numel(edges) - 1;
animal_density = nan(numel(animals), n_bins);
for animal_n = 1:numel(animals)
    frequency = zeros(0, 1);
    for recording_n = 1:numel(recordings)
        if strcmp(recordings(recording_n).animal_id, animals{animal_n})
            frequency = [frequency; accepted_frequency(recordings(recording_n))]; %#ok<AGROW>
        end
    end
    counts = histogram_counts(frequency, edges);
    if sum(counts) > 0
        animal_density(animal_n, :) = counts / sum(counts);
    end
end
density = finite_mean_rows(animal_density);
centers = (edges(1:end - 1) + edges(2:end)) / 2;
end

function counts = histogram_counts(values, edges)
counts = zeros(1, numel(edges) - 1);
for bin_n = 1:numel(counts)
    if bin_n < numel(counts)
        counts(bin_n) = nnz(values >= edges(bin_n) & values < edges(bin_n + 1));
    else
        counts(bin_n) = nnz(values >= edges(bin_n) & values <= edges(bin_n + 1));
    end
end
end

function modes = histogram_modes(centers, density)
modes = zeros(0, 2);
if numel(density) < 3 || all(~isfinite(density))
    return
end
peak_level = 0.5 * max(density);
for bin_n = 2:numel(density) - 1
    if density(bin_n) >= peak_level && ...
            density(bin_n) > density(bin_n - 1) && ...
            density(bin_n) >= density(bin_n + 1)
        modes(end + 1, :) = [centers(bin_n), density(bin_n)]; %#ok<AGROW>
    end
end
end
