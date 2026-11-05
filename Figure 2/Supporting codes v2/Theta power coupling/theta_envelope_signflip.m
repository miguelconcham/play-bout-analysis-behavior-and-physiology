function result = theta_envelope_signflip(animal_values, n_flip)
%THETA_ENVELOPE_SIGNFLIP  Right-tailed sign-flip test of an animal mean.
%
% animal_values is one number per animal. The first null draw keeps every
% sign, so the observed mean is inside the null. p is the fraction of draws
% whose mean is at least the observed mean.

values = double(animal_values(:));
values = values(isfinite(values));
n = numel(values);
result.n_animals = n;
result.mean = NaN;
result.p_greater = NaN;
result.null_mean = zeros(0, 1);
if n == 0
    return
end
result.mean = mean(values);
if n < 2 || n_flip < 1
    return
end
signs = ones(n_flip, n);
if n_flip > 1
    signs(2:end, :) = 2 * (rand(n_flip - 1, n) > 0.5) - 1;
end
result.null_mean = (signs * values) / n;
result.p_greater = mean(result.null_mean >= result.mean);
end
