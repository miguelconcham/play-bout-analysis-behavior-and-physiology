function summary = component_phase_summary(phases, n_bins)
%COMPONENT_PHASE_SUMMARY  Phase lock of one spike sample.
%
% phases are radians in the same circularly uniformized Hilbert coordinate
% used for the delta-lock labels: 0 is the peak side, +/- pi is the trough
% side. proj_peak is the mean of cos(phase). proj_trough is the mean of
% cos(phase - pi), which is -proj_peak. frac_peak is the fraction of spikes
% inside [-pi/2, pi/2], the same half-circle as delta_lock_group.
%
% p_peak and p_trough are one-sided V-test p-values from the normal
% approximation. They describe one neuron. The group test is in Analyze.

if nargin < 2 || isempty(n_bins)
    n_bins = 32;
end
n_bins = max(1, round(n_bins));
phases = double(phases(:));
phases = phases(isfinite(phases));
phases = mod(phases + pi, 2 * pi) - pi;

summary.n = numel(phases);
summary.angle = NaN;
summary.mvl = NaN;
summary.ppc = NaN;
summary.proj_peak = NaN;
summary.proj_trough = NaN;
summary.p_peak = NaN;
summary.p_trough = NaN;
summary.frac_peak = NaN;
summary.edges = linspace(-pi, pi, n_bins + 1);
summary.hist = zeros(1, n_bins);
if summary.n < 1
    return
end

summary.hist = phase_hist_counts(phases, summary.edges);
z = exp(1i * phases);
resultant = sum(z);
summary.mvl = abs(resultant) / summary.n;
summary.angle = atan2(imag(resultant), real(resultant));
summary.proj_peak = mean(cos(phases));
summary.proj_trough = mean(cos(phases - pi));
summary.frac_peak = mean(phases >= -pi / 2 & phases <= pi / 2);
if summary.n >= 2
    summary.ppc = (abs(resultant) ^ 2 - summary.n) / (summary.n * (summary.n - 1));
    summary.p_peak = v_test_p(summary.proj_peak, summary.n);
    summary.p_trough = v_test_p(summary.proj_trough, summary.n);
end
end


function counts = phase_hist_counts(phases, edges)
% Same edges as histcounts: every bin is [left, right), and the last bin
% also includes the right edge.
n_bins = numel(edges) - 1;
counts = zeros(1, n_bins);
if isempty(phases) || n_bins < 1
    return
end
for k = 1:n_bins - 1
    counts(k) = sum(phases >= edges(k) & phases < edges(k + 1));
end
counts(n_bins) = sum(phases >= edges(n_bins) & phases <= edges(end));
end


function p = v_test_p(projection, n)
% projection is mean cos(phase - target). Under uniformity the sum of
% those cosines is normal with variance n/2.
u = projection * sqrt(2 * n);
p = 0.5 * erfc(u / sqrt(2));
end
