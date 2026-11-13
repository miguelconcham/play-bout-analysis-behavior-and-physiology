function adjusted = holm_sidak_adjust(raw_p)
%HOLM_SIDAK_ADJUST  Step-down Holm-Sidak adjusted p-values.
%
% Non-finite raw p-values stay out of the family. For the rank-r test in
% a family of m, the Sidak value is 1 - (1 - p)^(m - r + 1). The adjusted
% values then cannot decrease as the raw p increases, and they are capped
% at 1. One finite p is returned unchanged.

adjusted = nan(size(raw_p));
usable = isfinite(raw_p);
p = raw_p(usable);
p = p(:);
n_tests = numel(p);
if n_tests == 0
    return
end
p = min(max(p, 0), 1);
[sorted_p, order] = sort(p);
sorted_adjusted = zeros(n_tests, 1);
running_max = 0;
for rank = 1:n_tests
    remaining = n_tests - rank + 1;
    this_adjusted = 1 - (1 - sorted_p(rank)) ^ remaining;
    if this_adjusted < running_max
        this_adjusted = running_max;
    end
    if this_adjusted > 1
        this_adjusted = 1;
    end
    running_max = this_adjusted;
    sorted_adjusted(rank) = this_adjusted;
end
unsorted = nan(n_tests, 1);
unsorted(order) = sorted_adjusted;
adjusted(usable) = unsorted;
end
