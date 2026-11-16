function interval = acg_mean_interval(matrix)
%ACG_MEAN_INTERVAL  Mean and 95% t confidence interval across rows.
%
% Each column is one lag. The interval uses the neurons that are finite
% at that lag. Fewer than two finite values leaves the half-width NaN,
% and the low and high edges then sit on the mean.

matrix = double(matrix);
if isempty(matrix)
    interval = empty_mean_interval(0);
    return
end
if isvector(matrix)
    matrix = matrix(:).';
end
n_columns = size(matrix, 2);
n = zeros(1, n_columns);
mu = nan(1, n_columns);
sigma = nan(1, n_columns);
for j = 1:n_columns
    column = matrix(:, j);
    column = column(isfinite(column));
    n(j) = numel(column);
    if n(j) >= 1
        mu(j) = mean(column);
    end
    if n(j) >= 2
        sigma(j) = std(column, 0);
    end
end
half = nan(1, n_columns);
for j = 1:n_columns
    if n(j) >= 2 && isfinite(sigma(j)) && isfinite(mu(j))
        half(j) = mean_t_critical(0.025, n(j) - 1) * sigma(j) / sqrt(n(j));
    end
end
low = mu;
high = mu;
finite_half = isfinite(half);
low(finite_half) = mu(finite_half) - half(finite_half);
high(finite_half) = mu(finite_half) + half(finite_half);
interval.n = n;
interval.mean = mu;
interval.half = half;
interval.low = low;
interval.high = high;
end


function interval = empty_mean_interval(n_columns)
interval.n = zeros(1, n_columns);
interval.mean = nan(1, n_columns);
interval.half = nan(1, n_columns);
interval.low = nan(1, n_columns);
interval.high = nan(1, n_columns);
end


function tc = mean_t_critical(tail_p, df)
x = betaincinv(2 * tail_p, df / 2, 0.5);
tc = sqrt(df * (1 - x) / x);
end
