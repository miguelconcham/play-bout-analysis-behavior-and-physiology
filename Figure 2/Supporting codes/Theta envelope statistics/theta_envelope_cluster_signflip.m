function result = theta_envelope_cluster_signflip(animal_by_frequency, n_flip)
%THETA_ENVELOPE_CLUSTER_SIGNFLIP  Cluster sign-flip across the frequency grid.
%
% A cluster is a run of positive animal-mean values. Its statistic is the
% sum of that mean. The null is the largest cluster sum after flipping
% animal signs. The p of each observed cluster is the fraction of null maxima
% at least as large as that cluster, which is the correction for scanning
% the grid.

values = double(animal_by_frequency);
valid_animal = all(isfinite(values), 2);
values = values(valid_animal, :);
[n_animals, n_freq] = size(values);
result.n_animals = n_animals;
result.observed_mean = mean(values, 1);
result.clusters = positive_clusters(result.observed_mean);
result.p = NaN;
result.max_cluster_sum = 0;
result.null_max = zeros(0, 1);
if isempty(result.clusters)
    result.p = 1;
    return
end
result.max_cluster_sum = max([result.clusters.sum]);
if n_animals < 2 || n_flip < 1 || n_freq < 1
    return
end
signs = ones(n_flip, n_animals);
if n_flip > 1
    signs(2:end, :) = 2 * (rand(n_flip - 1, n_animals) > 0.5) - 1;
end
null_max = zeros(n_flip, 1);
for flip = 1:n_flip
    flipped = signs(flip, :).' .* values;
    clusters = positive_clusters(mean(flipped, 1));
    if ~isempty(clusters)
        null_max(flip) = max([clusters.sum]);
    end
end
result.null_max = null_max;
for c = 1:numel(result.clusters)
    result.clusters(c).p = mean(null_max >= result.clusters(c).sum);
end
result.p = mean(null_max >= result.max_cluster_sum);
end

function clusters = positive_clusters(mu)
mu = double(mu(:).');
above = mu > 0 & isfinite(mu);
edges = diff([false, above, false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
clusters = struct('start', {}, 'stop', {}, 'sum', {}, 'p', {});
for c = 1:numel(starts)
    clusters(c).start = starts(c); %#ok<AGROW>
    clusters(c).stop = stops(c);
    clusters(c).sum = sum(mu(starts(c):stops(c)));
    clusters(c).p = NaN;
end
end
