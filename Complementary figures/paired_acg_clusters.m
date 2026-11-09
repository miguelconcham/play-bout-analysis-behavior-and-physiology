function result = paired_acg_clusters(low, high, centers, alpha, n_perm)
%PAIRED_ACG_CLUSTERS  Paired high-versus-low test along the autocorrelogram.
%
% Positive lags are tested. A cluster is a run of lags whose mean paired
% difference, high minus low, has a Student t interval that excludes zero.
% That is the same interval drawn on the difference. The cluster p-value
% comes from flipping the sign of each neuron's difference. The t statistic
% is only the score. Its null distribution is the sign-flip distribution,
% so the cluster decision does not require the differences to be normal.

centers = centers(:).';
pos = centers > 0;
lag = centers(pos);
difference = high(:, pos) - low(:, pos);
n_time = numel(lag);
n_neurons = size(difference, 1);
finite = isfinite(difference);
n = sum(finite, 1);
df = max(n - 1, 1);
tcrit = t_critical(df, alpha);
tcrit(n < 2) = Inf;
[point_t, point_p] = column_t(difference, finite, n, df);
sig = n >= 2 & abs(point_t) > tcrit;
half = median(diff(centers)) / 2;
if ~isfinite(half)
    half = 0;
end
runs = true_runs(sig);
n_clusters = size(runs, 1);
clusters = repmat(empty_cluster(), n_clusters, 1);
if n_clusters > 0 && n_perm > 0 && n_neurons > 0
    stream = RandStream('mt19937ar', 'Seed', 1);
    neuron_flip = 2 * (rand(stream, n_neurons, n_perm) > 0.5) - 1;
    max_mass = zeros(1, n_perm);
    for perm = 1:n_perm
        flipped = difference .* neuron_flip(:, perm);
        t_perm = column_t(flipped, finite, n, df);
        max_mass(perm) = max_cluster_mass(n >= 2 & abs(t_perm) > tcrit, t_perm);
    end
else
    max_mass = [];
end
for k = 1:n_clusters
    cols = runs(k, 1):runs(k, 2);
    clusters(k).x0 = lag(cols(1)) - half;
    clusters(k).x1 = lag(cols(end)) + half;
    clusters(k).mass = sum(abs(point_t(cols)));
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


function [tstat, p] = column_t(X, finite, n, df)
X0 = X;
X0(~finite) = 0;
mu = sum(X0, 1) ./ max(n, 1);
deviation = X0 - finite .* mu;
sem = sqrt(sum(deviation .^ 2, 1) ./ df) ./ sqrt(max(n, 1));
tstat = zeros(1, numel(n));
comparable = n >= 2 & sem > 0;
tstat(comparable) = mu(comparable) ./ sem(comparable);
perfect = n >= 2 & sem == 0 & mu ~= 0;
tstat(perfect) = sign(mu(perfect)) .* Inf;
if nargout < 2
    return
end
p = ones(1, numel(n));
p(n >= 2) = t_two_tailed(tstat(n >= 2), df(n >= 2));
end


function crit = t_critical(df, alpha)
if exist('tinv', 'file') == 2
    crit = tinv(1 - alpha / 2, df);
else
    crit = 1.96 * ones(size(df));
end
end


function p = t_two_tailed(tstat, dof)
p = nan(size(tstat));
ok = isfinite(tstat) & isfinite(dof) & dof > 0;
if ~any(ok)
    p(~isfinite(tstat) & dof > 0) = 0;
    return
end
ratio = dof(ok) ./ (dof(ok) + tstat(ok) .^ 2);
p(ok) = betainc(ratio, dof(ok) / 2, 0.5);
p(~isfinite(tstat) & dof > 0) = 0;
end


function runs = true_runs(sig)
d = diff([false, sig(:).', false]);
i1 = find(d == 1);
i2 = find(d == -1) - 1;
runs = [i1(:), i2(:)];
end


function mass = max_cluster_mass(sig, tstat)
runs = true_runs(sig);
mass = 0;
for k = 1:size(runs, 1)
    mass = max(mass, sum(abs(tstat(runs(k, 1):runs(k, 2)))));
end
end


function cluster = empty_cluster()
cluster = struct('x0', NaN, 'x1', NaN, 'mass', NaN, 'p', NaN, 'survives', false);
end
