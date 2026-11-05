function correlation = theta_envelope_leave_one_out(epochs)
%THETA_ENVELOPE_LEAVE_ONE_OUT  Correlation of each behavior with the other mean.
%
% The template is the mean of the other finite rows from the same recording.
% A row that matches the up-and-down shared by the rest scores high.

n_rows = size(epochs, 1);
correlation = nan(n_rows, 1);
valid = all(isfinite(epochs), 2);
x = epochs(valid, :);
n = size(x, 1);
if n < 2
    return
end
x = x - mean(x, 2);
total = sum(x, 1);
scored = nan(n, 1);
for row = 1:n
    template = (total - x(row, :)) / (n - 1);
    scored(row) = one_correlation(x(row, :), template);
end
correlation(valid) = scored;
end

function value = one_correlation(a, b)
a = a - mean(a);
b = b - mean(b);
denominator = sqrt(sum(a .^ 2) * sum(b .^ 2));
if denominator == 0
    value = NaN;
else
    value = sum(a .* b) / denominator;
end
end
