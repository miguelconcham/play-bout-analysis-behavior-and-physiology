function counts = crosscorr_acg_counts(spikes, t0, t1, bin_edges, closed_end)
%CROSSCORR_ACG_COUNTS  Autocorrelogram from crosscorrelogram.
%
% Reference spikes are inset from each end by the longest lag in
% bin_edges. Partner spikes are every spike in the interval. For a window
% [0, 10] and lags of +/- 1 s, the references are the spikes in [1, 9] and
% the partners are the spikes in [0, 10]. Each reference then has the full
% lag range inside the interval, so the histogram is not divided by
% (T-|lag|).
%
% crosscorrelogram starts at the second center time and stops its partner
% scan one sample short of the last timestamp. One time is placed before
% the first reference and one after the last partner so those samples stay
% in the correlogram. A zero offset is the spike paired with itself and is
% removed.

if nargin < 5 || isempty(closed_end)
    closed_end = true;
end
bin_edges = double(bin_edges(:))';
n_bins = max(numel(bin_edges) - 1, 0);
counts = nan(1, n_bins);
if n_bins < 1 || ~isfinite(t0) || ~isfinite(t1) || t1 <= t0
    return
end
margin = max(abs(bin_edges(1)), abs(bin_edges(end)));
if ~(margin > 0) || (t1 - t0) <= 2 * margin
    return
end
spikes = double(spikes(:));
spikes = sort(spikes(isfinite(spikes)));
if closed_end
    partners = spikes(spikes >= t0 & spikes <= t1);
else
    partners = spikes(spikes >= t0 & spikes < t1);
end
references = partners(partners >= (t0 + margin) & partners <= (t1 - margin));
if isempty(references) || isempty(partners)
    return
end
references = [references(1) - (margin + 1); references];
partners = [partners; partners(end) + (margin + 1)];
offsets = crosscorrelogram(references, partners, [bin_edges(1), bin_edges(end)]);
offsets = double(offsets(:));
offsets = offsets(isfinite(offsets) & offsets ~= 0);
if isempty(offsets)
    counts = zeros(1, n_bins);
    return
end
counts = histcounts(offsets, bin_edges);
end
