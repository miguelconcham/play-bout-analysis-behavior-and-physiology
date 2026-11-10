function pairs = match_power_periods(reference, candidate, caliper)
%MATCH_POWER_PERIODS  One-to-one match of two period lists on one power scale.
%
% reference and candidate are column vectors on the scale being matched
% (decibels in the play-power analysis). caliper is the largest absolute
% difference that may be paired. The pairing minimizes the sum of absolute
% differences. A reference period that cannot be paired costs twice the
% caliper, and a candidate may be left unused at no cost. A candidate is
% used at most once. This is the same decision as matchpairs with costs
% above the caliper removed.
%
% pairs is P x 2, [reference index, candidate index], sorted by the
% reference index. Indices point into the vectors that were passed in.

reference = double(reference(:));
candidate = double(candidate(:));
pairs = zeros(0, 2);
if ~isscalar(caliper) || ~isfinite(caliper) || caliper < 0
    error('match_power_periods:caliper', ...
        'Caliper must be a finite nonnegative scalar.');
end
valid_ref = find(isfinite(reference));
valid_cand = find(isfinite(candidate));
n = numel(valid_ref);
m = numel(valid_cand);
if n < 1 || m < 1
    return
end
[ref_v, ref_ord] = sort(reference(valid_ref), 'ascend');
[cand_v, cand_ord] = sort(candidate(valid_cand), 'ascend');
unmatched = 2 * caliper;
choice = zeros(n, m, 'uint8');
dp_prev = zeros(1, m + 1);
for i = 1:n
    dp_curr = inf(1, m + 1);
    dp_curr(1) = i * unmatched;
    for j = 1:m
        best = dp_curr(j);
        which = uint8(1);
        skip_ref = dp_prev(j + 1) + unmatched;
        if skip_ref < best
            best = skip_ref;
            which = uint8(2);
        end
        dist = abs(ref_v(i) - cand_v(j));
        if dist <= caliper + 1e-9
            match_cost = dp_prev(j) + dist;
            if match_cost <= best
                best = match_cost;
                which = uint8(3);
            end
        end
        dp_curr(j + 1) = best;
        choice(i, j) = which;
    end
    dp_prev = dp_curr;
end

i = n;
j = m;
ref_pos = zeros(0, 1);
cand_pos = zeros(0, 1);
while i > 0 && j > 0
    switch choice(i, j)
        case 1
            j = j - 1;
        case 2
            i = i - 1;
        case 3
            ref_pos(end + 1, 1) = i; %#ok<AGROW>
            cand_pos(end + 1, 1) = j; %#ok<AGROW>
            i = i - 1;
            j = j - 1;
        otherwise
            error('match_power_periods:path', 'The match path is incomplete.');
    end
end
if isempty(ref_pos)
    return
end
pairs = [valid_ref(ref_ord(ref_pos)), valid_cand(cand_ord(cand_pos))];
pairs = sortrows(pairs, 1);
end
