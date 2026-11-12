function exponent = acg_rate_exponent(method)
%ACG_RATE_EXPONENT  Power of firing rate in the Poisson ratio.
%
% Eligible-spike curves are partners per spike, so the Poisson baseline
% is rate times bin width (exponent 1). Overlap-weighted curves are pair
% counts per second, so the baseline is rate squared times bin width
% (exponent 2). A missing method is the overlap exponent, which keeps
% summaries saved before the selector on the previous scale.

if nargin < 1 || isempty(method)
    exponent = 2;
    return
end
if strcmp(canonical_acg_method(method), 'eligible')
    exponent = 1;
else
    exponent = 2;
end
end
