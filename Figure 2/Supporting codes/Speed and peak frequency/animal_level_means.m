function summary = animal_level_means(differences, animals)
%ANIMAL_LEVEL_MEANS  One mean play-minus-non-play value per animal.
%
% Sessions from the same animal are averaged first, so an animal with more
% sessions does not contribute more than one number to the equal-animal mean.

differences = differences(:);
animals = animal_name_list(animals);
if numel(animals) ~= numel(differences)
    error('animal_level_means:Size', ...
        'Each difference needs one animal label.');
end
keep = isfinite(differences);
differences = differences(keep);
animals = animals(keep);

names = unique(animals, 'stable');
n_animals = numel(names);
means = nan(n_animals, 1);
counts = zeros(n_animals, 1);
for k = 1:n_animals
    selected = strcmp(animals, names{k});
    counts(k) = nnz(selected);
    means(k) = mean(differences(selected));
end

summary.animals = names;
summary.means = means;
summary.counts = counts;
summary.n_sessions = numel(differences);
summary.n_animals = n_animals;
if n_animals == 0
    summary.equal_animal_mean = NaN;
else
    summary.equal_animal_mean = mean(means);
end
end

function names = animal_name_list(animals)
animals = animals(:);
names = cell(numel(animals), 1);
for k = 1:numel(animals)
    item = animals{k};
    if exist('isstring', 'builtin') == 5 && isstring(item)
        item = char(item);
    end
    names{k} = char(item);
end
end
