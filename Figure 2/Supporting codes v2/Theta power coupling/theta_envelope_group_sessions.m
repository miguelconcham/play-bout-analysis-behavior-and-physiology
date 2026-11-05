function grouped = theta_envelope_group_sessions(recordings)
%THETA_ENVELOPE_GROUP_SESSIONS  Average session z-scores within each animal.
%
% Recordings that share animal_area_table_name are one animal. The peak is
% the 0.2 Hz center of the largest mean spectrum or coherence. The z-scores
% compare those values with the surrogate null.

grouped = empty_grouped();
if isempty(recordings)
    return
end
names = {recordings.animal_area_table_name};
animal_names = unique(names, 'stable');
n_animals = numel(animal_names);
n_test = numel(recordings(1).frequencies_test);
n_fine = numel(recordings(1).frequencies_fine);
grouped.animals = animal_names(:);
grouped.frequencies_test = recordings(1).frequencies_test;
grouped.frequencies_fine = recordings(1).frequencies_fine;
grouped.confirmatory_mask = recordings(1).confirmatory_mask;
grouped.power_z_test = nan(n_animals, n_test);
grouped.coherence_z_test = nan(n_animals, n_test);
grouped.power_z_fine = nan(n_animals, n_fine);
grouped.coherence_z_fine = nan(n_animals, n_fine);
grouped.band_power_fine = nan(n_animals, n_fine);
grouped.band_coherence_fine = nan(n_animals, n_fine);
grouped.peak_power_hz = nan(n_animals, 1);
grouped.peak_coherence_hz = nan(n_animals, 1);
grouped.n_recordings = zeros(n_animals, 1);
for animal_n = 1:n_animals
    rows = recordings(strcmp(names, animal_names{animal_n}));
    grouped.n_recordings(animal_n) = numel(rows);
    grouped.power_z_test(animal_n, :) = mean_rows(stack_field(rows, 'power_z_test'));
    grouped.coherence_z_test(animal_n, :) = mean_rows(stack_field(rows, 'coherence_z_test'));
    grouped.power_z_fine(animal_n, :) = mean_rows(stack_field(rows, 'power_z_fine'));
    grouped.coherence_z_fine(animal_n, :) = mean_rows(stack_field(rows, 'coherence_z_fine'));
    grouped.band_power_fine(animal_n, :) = mean_rows(stack_field(rows, 'band_power_fine'));
    grouped.band_coherence_fine(animal_n, :) = mean_rows(stack_field(rows, 'band_coherence_fine'));
    grouped.peak_power_hz(animal_n) = peak_hz( ...
        grouped.band_power_fine(animal_n, :), grouped.frequencies_fine);
    grouped.peak_coherence_hz(animal_n) = peak_hz( ...
        grouped.band_coherence_fine(animal_n, :), grouped.frequencies_fine);
end
end

function values = stack_field(rows, field_name)
values = [];
for row_n = 1:numel(rows)
    values = [values; rows(row_n).session.(field_name)(:).']; %#ok<AGROW>
end
end

function value = mean_rows(values)
value = nan(1, size(values, 2));
for column = 1:size(values, 2)
    finite = values(isfinite(values(:, column)), column);
    if ~isempty(finite)
        value(column) = mean(finite);
    end
end
end

function hz = peak_hz(z, centers)
if all(~isfinite(z))
    hz = NaN;
    return
end
z(~isfinite(z)) = -inf;
[~, idx] = max(z);
hz = centers(idx);
end

function grouped = empty_grouped()
grouped = struct( ...
    'animals', {{}}, ...
    'frequencies_test', [], ...
    'frequencies_fine', [], ...
    'confirmatory_mask', [], ...
    'power_z_test', [], ...
    'coherence_z_test', [], ...
    'power_z_fine', [], ...
    'coherence_z_fine', [], ...
    'band_power_fine', [], ...
    'band_coherence_fine', [], ...
    'peak_power_hz', [], ...
    'peak_coherence_hz', [], ...
    'n_recordings', []);
end
