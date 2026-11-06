function result = paired_acg_clusters(low, high, centers, alpha, n_perm)
%PAIRED_ACG_CLUSTERS  Paired high-versus-low test along the autocorrelogram.
%
% Positive lags are tested. A cluster is a run of lags whose two-sided
% signed-rank p-value is below alpha. The cluster p-value comes from
% flipping the high/low label inside each neuron. The plotted bars use the
% uncorrected runs; this p-value says which of those runs survive.

centers = centers(:).';
pos = centers > 0;
lag = centers(pos);
difference = high(:, pos) - low(:, pos);
n_time = numel(lag);
n_neurons = size(difference, 1);
point_p = ones(1, n_time);
point_z = zeros(1, n_time);
valid = false(1, n_time);
prep = cell(1, n_time);
for t = 1:n_time
    [point_z(t), point_p(t), prep{t}] = column_wilcoxon(difference(:, t));
    valid(t) = prep{t}.valid;
end
sig = valid & point_p < alpha;
half = median(diff(centers)) / 2;
runs = true_runs(sig);
n_clusters = size(runs, 1);
clusters = repmat(empty_cluster(), n_clusters, 1);
if n_clusters > 0 && n_perm > 0
    stream = RandStream('mt19937ar', 'Seed', 1);
    neuron_flip = 2 * (rand(stream, n_neurons, n_perm) > 0.5) - 1;
    max_mass = zeros(1, n_perm);
    for perm = 1:n_perm
        [sig_perm, z_perm] = flipped_significance( ...
            difference, prep, neuron_flip(:, perm), alpha);
        max_mass(perm) = max_cluster_mass(sig_perm, z_perm);
    end
else
    max_mass = [];
end
for k = 1:n_clusters
    cols = runs(k, 1):runs(k, 2);
    clusters(k).x0 = lag(cols(1)) - half;
    clusters(k).x1 = lag(cols(end)) + half;
    clusters(k).mass = sum(abs(point_z(cols)));
    if isempty(max_mass)
        clusters(k).p = NaN;
    else
        clusters(k).p = (1 + sum(max_mass >= clusters(k).mass - 1e-12)) / (n_perm + 1);
    end
    clusters(k).survives = clusters(k).p <= alpha;
end

result.alpha = alpha;
result.n_perm = n_perm;
result.n = n_neurons;
result.lag = lag;
result.p = point_p;
result.clusters = clusters;
end


function [z, p, prep] = column_wilcoxon(x)
ok = isfinite(x) & x ~= 0;
xx = x(ok);
n = numel(xx);
prep.ok = ok;
prep.n = n;
prep.valid = false;
prep.R = [];
prep.sign = [];
prep.mu = NaN;
prep.sigma = NaN;
z = 0;
p = 1;
if n < 5
    return
end
R = tiedrank(abs(xx));
prep.R = R(:);
prep.sign = sign(xx(:));
prep.mu = n * (n + 1) / 4;
sigma2 = n * (n + 1) * (2 * n + 1) / 24 - tie_cubic(abs(xx)) / 48;
prep.sigma = sqrt(max(sigma2, 0));
if prep.sigma <= 0
    return
end
prep.valid = true;
W = sum(R(xx > 0));
z = wilcoxon_z(W, prep.mu, prep.sigma);
p = min(1, 2 * normcdf(-abs(z)));
end


function [sig, z] = flipped_significance(difference, prep, neuron_flip, alpha)
n_time = numel(prep);
sig = false(1, n_time);
z = zeros(1, n_time);
for t = 1:n_time
    if ~prep{t}.valid
        continue
    end
    flipped = neuron_flip(prep{t}.ok) .* prep{t}.sign;
    W = prep{t}.R.' * double(flipped > 0);
    z(t) = wilcoxon_z(W, prep{t}.mu, prep{t}.sigma);
    sig(t) = min(1, 2 * normcdf(-abs(z(t)))) < alpha;
end
end


function z = wilcoxon_z(W, mu, sigma)
delta = W - mu;
z = (delta - 0.5 * sign(delta)) / sigma;
z(delta == 0) = 0;
end


function s = tie_cubic(v)
sv = sort(v(:));
s = 0;
k = 1;
n = numel(sv);
while k <= n
    j = k;
    while j < n && sv(j + 1) == sv(k)
        j = j + 1;
    end
    t = j - k + 1;
    if t > 1
        s = s + t ^ 3 - t;
    end
    k = j + 1;
end
end


function runs = true_runs(sig)
d = diff([false, sig(:).', false]);
i1 = find(d == 1);
i2 = find(d == -1) - 1;
runs = [i1(:), i2(:)];
end


function mass = max_cluster_mass(sig, z)
runs = true_runs(sig);
mass = 0;
for k = 1:size(runs, 1)
    mass = max(mass, sum(abs(z(runs(k, 1):runs(k, 2)))));
end
end


function cluster = empty_cluster()
cluster = struct('x0', NaN, 'x1', NaN, 'mass', NaN, 'p', NaN, 'survives', false);
end
