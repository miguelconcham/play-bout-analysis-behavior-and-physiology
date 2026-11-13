function acg = run_window_acgs(spike_times, run_starts, run_ends, bin_edges, method, closed_end)
%RUN_WINDOW_ACGS  Autocorrelogram of each firing run or epoch.
%
%   acg = run_window_acgs(spike_times, run_starts, run_ends, bin_edges)
%   acg = run_window_acgs(..., method)
%   acg = run_window_acgs(..., method, closed_end)
%
% method chooses how the pairs are counted:
%   'pairs'      autocorrelogram of every spike in the window. This is
%                the default. Counts are divided by the overlap (T-|lag|).
%   'fft'        the same overlap correction, with counts from the FFT of
%                the binned window.
%   'crosscorr'  crosscorrelogram. Reference spikes are inset from each
%                edge by the lag limit. Partner spikes are every spike in
%                the window. A run from 0 to 10 s with lags of +/- 1 s
%                uses references in [1 9] and partners in [0 10]. Every
%                reference spike has the full lag window, so the count is
%                divided by that inner duration.
%
% closed_end true keeps a spike at the window end. That is the firing-run
% rule and the default. closed_end false uses [start, end), which is the
% play-epoch rule. Bins that the chosen count cannot form are NaN.

if nargin < 5 || isempty(method)
    method = 'pairs';
end
if nargin < 6 || isempty(closed_end)
    closed_end = true;
end
method = canonical_count_method(method);
closed_end = logical(closed_end);
spike_times = double(spike_times(:));
run_starts = double(run_starts(:));
run_ends = double(run_ends(:));
bin_edges = bin_edges(:)';
n_runs = numel(run_starts);
n_bins = max(numel(bin_edges) - 1, 0);
acg = nan(n_runs, n_bins);
if n_runs < 1 || n_bins < 1 || numel(run_ends) ~= n_runs
    return
end
range_limits = [bin_edges(1), bin_edges(end)];
for j = 1:n_runs
    if ~isfinite(run_starts(j)) || ~isfinite(run_ends(j)) || run_ends(j) <= run_starts(j)
        continue
    end
    in_run = spikes_in_window(spike_times, run_starts(j), run_ends(j), closed_end);
    switch method
        case 'pairs'
            lags = autocorrelogram(spike_times(in_run), range_limits);
            counts = lag_bin_counts(lags, bin_edges);
            acg(j, :) = edge_correct_acg(counts, bin_edges, run_ends(j) - run_starts(j));
        case 'fft'
            counts = fft_acg_counts(spike_times(in_run), run_starts(j), run_ends(j), bin_edges);
            acg(j, :) = edge_correct_acg(counts, bin_edges, run_ends(j) - run_starts(j));
        case 'crosscorr'
            acg(j, :) = crosscorr_acg(spike_times(in_run), run_starts(j), run_ends(j), ...
                bin_edges, closed_end);
    end
end
end


function inside = spikes_in_window(spike_times, t0, t1, closed_end)
if closed_end
    inside = spike_times >= t0 & spike_times <= t1;
else
    inside = spike_times >= t0 & spike_times < t1;
end
end


function rate = crosscorr_acg(spikes, t0, t1, bin_edges, closed_end)
% Partners from the whole window, references from the window inset by the lag.
% On [0 10] with lags +/- 1 this is partners in [0 10] and references in [1 9].
n_bins = numel(bin_edges) - 1;
rate = nan(1, n_bins);
margin = max(abs(bin_edges(1)), abs(bin_edges(end)));
inner0 = t0 + margin;
inner1 = t1 - margin;
if ~(inner1 > inner0) || ~isfinite(margin)
    return
end
spikes = sort(spikes(isfinite(spikes)));
ts2 = spikes(spikes_in_window(spikes, t0, t1, closed_end));
ts1 = spikes(spikes_in_window(spikes, inner0, inner1, closed_end));
if isempty(ts1) || isempty(ts2)
    rate(:) = 0;
    return
end
% crosscorrelogram starts at the second reference and stops short of the
% last partner index. A leading reference outside the lag window and a
% trailing partner past the window put every real spike where the function
% will see it.
ts1_call = [inner0 - margin - 1; ts1(:)];
ts2_call = [ts2(:); t1 + margin + 1];
offsets = crosscorrelogram(ts1_call, ts2_call, [bin_edges(1), bin_edges(end)]);
offsets = offsets(isfinite(offsets) & offsets ~= 0);
counts = lag_bin_counts(offsets, bin_edges);
rate = counts / (inner1 - inner0);
end


function counts = fft_acg_counts(spikes, t0, t1, bin_edges)
dt = median(diff(bin_edges));
n_bins = numel(bin_edges) - 1;
counts = zeros(1, n_bins);
if ~isfinite(dt) || dt <= 0 || ~(t1 > t0)
    return
end
n = floor((t1 - t0) / dt);
if n < 2
    return
end
x = zeros(1, n);
spikes = spikes(isfinite(spikes));
for s = 1:numel(spikes)
    if spikes(s) < t0 || spikes(s) >= t0 + n * dt
        continue
    end
    bin = floor((spikes(s) - t0) / dt) + 1;
    if bin >= 1 && bin <= n
        x(bin) = x(bin) + 1;
    end
end
nfft = 2 ^ nextpow2(2 * n - 1);
spectrum = fft(x, nfft);
corr = real(ifft(spectrum .* conj(spectrum)));
corr(1) = corr(1) - sum(x);
corr(abs(corr) < 1e-8) = 0;
for k = 1:n - 1
    lag = k * dt;
    counts = add_lag_count(counts, bin_edges, lag, corr(k + 1));
    counts = add_lag_count(counts, bin_edges, -lag, corr(k + 1));
end
end


function counts = lag_bin_counts(lags, edges)
% Same edges as histcounts: every bin is [left, right), and the last bin
% also keeps a lag that lands on the right edge.
lags = double(lags(:));
edges = double(edges(:)');
n_bins = numel(edges) - 1;
counts = zeros(1, n_bins);
if n_bins < 1
    return
end
lags = lags(isfinite(lags));
for b = 1:n_bins
    if b < n_bins
        in_bin = lags >= edges(b) & lags < edges(b + 1);
    else
        in_bin = lags >= edges(b) & lags <= edges(b + 1);
    end
    counts(b) = sum(in_bin);
end
end


function counts = add_lag_count(counts, edges, lag, value)
if ~isfinite(value) || value == 0
    return
end
n_bins = numel(edges) - 1;
if lag < edges(1) || lag > edges(end)
    return
end
bin = find(lag >= edges(1:end - 1) & lag < edges(2:end), 1, 'first');
if isempty(bin) && lag == edges(end)
    bin = n_bins;
end
if ~isempty(bin)
    counts(bin) = counts(bin) + value;
end
end
