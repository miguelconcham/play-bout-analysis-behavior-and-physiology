function snippets = center_waveform_snippets(snippets, time_samples)
%CENTER_WAVEFORM_SNIPPETS  Subtract the pre-spike baseline on each channel.
%
% snippets is n_spikes x n_channels x n_times. For every spike and channel
% the baseline is the mean of finite samples with time_samples <= -10.
% Those samples sit at or before 10 samples prior to the spike (1/3 ms at
% 30 kHz when the spike is sample 0).

time_samples = double(time_samples(:)');
pre = time_samples <= -10;
if ~any(pre)
    pre = false(1, numel(time_samples));
    pre(1) = true;
end

n_spikes = size(snippets, 1);
n_channels = size(snippets, 2);
n_times = size(snippets, 3);
for channel_index = 1:n_channels
    traces = reshape(snippets(:, channel_index, :), n_spikes, n_times);
    pre_values = traces(:, pre);
    finite_pre = isfinite(pre_values);
    pre_values(~finite_pre) = 0;
    counts = sum(finite_pre, 2);
    baseline = zeros(n_spikes, 1);
    has_pre = counts > 0;
    baseline(has_pre) = sum(pre_values(has_pre, :), 2) ./ counts(has_pre);
    traces = bsxfun(@minus, traces, baseline);
    snippets(:, channel_index, :) = reshape(traces, n_spikes, 1, n_times);
end
end
