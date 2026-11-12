function [row, mass] = sum_to_one_acg(rate, centers, refractory_s)
%SUM_TO_ONE_ACG  Blank the refractory bins, then scale the row to sum to 1.
%
% Bins whose lag interval overlaps +/- refractory_s are set to NaN before
% the sum, so the zero-lag bin does not absorb the mass. A row with no
% finite positive mass becomes NaN. mass is the sum of the blanked rate
% and is NaN when the row cannot be scaled.

centers = double(centers(:))';
rate = double(rate(:))';
n_bins = numel(centers);
if numel(rate) ~= n_bins
    error('sum_to_one_acg:width', ...
        'The autocorrelogram has %d bins and the centers have %d.', numel(rate), n_bins);
end
blanked = blank_refractory_bins(rate, centers, refractory_s);
mass = sum(blanked, 'omitnan');
if ~isfinite(mass) || mass <= 0
    mass = NaN;
end
row = normalize_acg_rows(blanked);
if isempty(row)
    row = nan(1, n_bins);
end
end


function blanked = blank_refractory_bins(rate, centers, refractory_s)
blanked = rate;
if isempty(centers) || ~isfinite(refractory_s) || refractory_s < 0
    return
end
gaps = diff(centers);
if isempty(gaps) || ~isfinite(gaps(1)) || gaps(1) <= 0
    half = refractory_s;
else
    half = abs(median(gaps)) / 2;
end
bin_start = centers - half;
bin_end = centers + half;
blanked(bin_start < refractory_s & bin_end > -refractory_s) = NaN;
end
