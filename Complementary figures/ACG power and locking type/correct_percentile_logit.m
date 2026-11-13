function logit = correct_percentile_logit(percentile, n_shuffles)
%CORRECT_PERCENTILE_LOGIT  log10(p/(1-p)) with p kept off 0 and 1.
%
% percentile is on a 0-100 scale. A value of 0 or 100 would make the
% logit infinite, so each value is pulled in by half a count out of
% n_shuffles + 1 before it is turned into a fraction. NaNs stay NaN.
% The logit of the corrected fraction is log10(p / (1-p)).

percentile = double(percentile);
if nargin < 2 || isempty(n_shuffles) || ~isfinite(n_shuffles) || n_shuffles < 1
    error('correct_percentile_logit:shuffles', ...
        'n_shuffles must be a positive finite number.');
end
edge = 50 / (double(n_shuffles) + 1);
logit = nan(size(percentile));
finite = isfinite(percentile);
corrected = min(max(percentile(finite), edge), 100 - edge);
p = corrected / 100;
logit(finite) = log10(p ./ (1 - p));
end
