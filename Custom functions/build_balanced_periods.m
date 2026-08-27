function Q = build_balanced_periods(P)
%BUILD_BALANCED_PERIODS Generate a second set of non-overlapping periods Q
% that are before the given periods P, with similar total duration.
%
% Input:
%   P : n×2 array of [start, end] times (sorted, non-overlapping)
%
% Output:
%   Q : m×2 array of [start, end] times (non-overlapping with P)
%
% Logic:
%   - For each period in P, create a window before it (Q)
%     with length as close as possible to that of P(i).
%   - If total length(Q) < total length(P), extend Q backward
%     where possible, and if still not enough, add a final
%     window before the first period.

% --- Ensure P is sorted by start time
P = sortrows(P,1);

n = size(P,1);
Q = zeros(n,2);

% --- Step 1: Build initial Qs
for i = 1:n
    startP = P(i,1);
    endP   = P(i,2);
    lenP   = endP - startP;

    % previous P end time (0 if none)
    if i == 1
        prev_end = -inf;  % allow extending before first if needed later
    else
        prev_end = P(i-1,2);
    end

    % maximum gap available before current P
    available_gap = startP - prev_end;

    % choose Q interval length (limited by available gap)
    lenQ = min(lenP, available_gap);

    Q(i,2) = startP;           % Q ends where P starts
    Q(i,1) = startP - lenQ;    % start as far back as fits
end

% --- Step 2: Compute total durations
lenP_all = sum(P(:,2) - P(:,1));
lenQ_all = sum(Q(:,2) - Q(:,1));

% --- Step 3: Try to extend Q backward if total shorter
if lenQ_all < lenP_all
    extra_needed = lenP_all - lenQ_all;

    % Try extending backwards (without touching earlier P)
    for i = 1:n
        if i == 1
            limit = -inf; % can extend freely before first Q
        else
            limit = P(i-1,2);
        end

        max_extend = Q(i,1) - limit;
        actual_extend = min(extra_needed, max_extend);
        Q(i,1) = Q(i,1) - actual_extend;
        extra_needed = extra_needed - actual_extend;

        if extra_needed <= 0
            break;
        end
    end

    % --- Step 4: If still not enough, add one more window before first Q
    if extra_needed > 0
        new_window = [Q(1,1) - extra_needed, Q(1,1)];
        Q = [new_window; Q];
    end
end

end
