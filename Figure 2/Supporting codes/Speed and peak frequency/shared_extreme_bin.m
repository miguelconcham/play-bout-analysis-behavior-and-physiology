function index = shared_extreme_bin(diff_mat, animals, direction)
%SHARED_EXTREME_BIN  Highest or lowest bin shared by two animals.
%
% A bin counts when at least two animals have a finite play-minus-non-play
% difference there. direction is 'high' or 'low'. The result is empty when
% no bin is shared.

if size(diff_mat, 1) ~= numel(animals)
    error('shared_extreme_bin:Size', ...
        'Difference rows (%d) do not match animals (%d).', ...
        size(diff_mat, 1), numel(animals));
end

n_bins = size(diff_mat, 2);
n_animals = zeros(1, n_bins);
for b = 1:n_bins
    keep = isfinite(diff_mat(:, b));
    if any(keep)
        n_animals(b) = numel(unique(animal_labels(animals(keep)), 'stable'));
    end
end
qualified = find(n_animals >= 2);
if isempty(qualified)
    index = [];
    return
end
if strcmp(direction, 'high')
    index = qualified(end);
elseif strcmp(direction, 'low')
    index = qualified(1);
else
    error('shared_extreme_bin:Direction', ...
        'direction must be ''high'' or ''low''.');
end
end

function names = animal_labels(animals)
names = cell(numel(animals), 1);
for k = 1:numel(animals)
    item = animals{k};
    if exist('isstring', 'builtin') == 5 && isstring(item)
        item = char(item);
    end
    names{k} = char(item);
end
end
