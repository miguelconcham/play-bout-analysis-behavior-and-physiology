function [speed_z, reference] = within_animal_speed_z(animal, session_key, bin_id, speed)
%WITHIN_ANIMAL_SPEED_Z  Z-score speed inside each animal.
%
% The mean and SD use one speed per animal, session, and bin. Channels
% that share that bin are the same locomotion sample, so a second channel
% does not change the SD. Bins with a missing peak still define the
% animal's speed scale. An animal with fewer than two distinct finite
% speeds gets NaN, because that scale is undefined.
%
% session_key identifies the behavioral session without the channel.
% animal and session_key are cell arrays of character vectors.

animal = animal(:);
session_key = session_key(:);
bin_id = double(bin_id(:));
speed = double(speed(:));
n = numel(speed);
if ~iscell(animal) || ~iscell(session_key) || ...
        numel(animal) ~= n || numel(session_key) ~= n || numel(bin_id) ~= n
    error('within_animal_speed_z:size', ...
        ['animal and session_key must be cell arrays, and bin_id and ' ...
        'speed must have one entry per row.']);
end

speed_z = nan(n, 1);
animal_names = unique(animal, 'stable');
reference = repmat(struct( ...
    'animal', '', 'n_unique', 0, 'mean', NaN, 'sd', NaN), numel(animal_names), 1);
for animal_n = 1:numel(animal_names)
    reference(animal_n).animal = animal_names{animal_n};
    in_animal = strcmp(animal, animal_names{animal_n});
    finite_index = find(in_animal & isfinite(speed));
    if isempty(finite_index)
        continue
    end
    keys = cell(numel(finite_index), 1);
    for sample_n = 1:numel(finite_index)
        row = finite_index(sample_n);
        keys{sample_n} = sprintf('%s|%.10g', session_key{row}, bin_id(row));
    end
    unique_keys = unique(keys, 'stable');
    unique_speed = zeros(numel(unique_keys), 1);
    for key_n = 1:numel(unique_keys)
        unique_speed(key_n) = mean(speed(finite_index(strcmp(keys, unique_keys{key_n}))));
    end
    reference(animal_n).n_unique = numel(unique_speed);
    reference(animal_n).mean = mean(unique_speed);
    if numel(unique_speed) < 2
        continue
    end
    sd = std(unique_speed, 0);
    reference(animal_n).sd = sd;
    if sd > 0
        speed_z(in_animal) = (speed(in_animal) - reference(animal_n).mean) / sd;
        speed_z(in_animal & ~isfinite(speed)) = NaN;
    end
end
end
