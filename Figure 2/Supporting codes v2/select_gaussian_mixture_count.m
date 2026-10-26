function best_count = select_gaussian_mixture_count( ...
        n_components, criterion_value, eligible)
%SELECT_GAUSSIAN_MIXTURE_COUNT  Smallest information criterion.
%
% Eligible components with a finite criterion are compared. Lower is
% better. A tie keeps the smallest component count, because the counts
% are stored in ascending order and min returns the first minimum.

n_components = n_components(:);
criterion_value = criterion_value(:);
eligible = logical(eligible(:));
if numel(n_components) ~= numel(criterion_value) || ...
        numel(n_components) ~= numel(eligible)
    error('regionGaussianCenters:selection', ...
        'Component count, criterion, and eligibility must match.');
end

keep = eligible & isfinite(criterion_value);
if ~any(keep)
    best_count = NaN;
    return
end
[~, best_index] = min(criterion_value(keep));
kept_counts = n_components(keep);
best_count = kept_counts(best_index);
end
