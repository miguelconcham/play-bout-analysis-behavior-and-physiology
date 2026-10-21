function samples = template_spike_samples(spike_times, spike_samples, ...
    bout_starts, hist_range, max_spikes, min_in_window)
%TEMPLATE_SPIKE_SAMPLES  Spikes used for a mean waveform outside the pairs.
%
% Spikes inside the play-bout windows are kept when there are at least
% min_in_window of them. Otherwise the whole train is used. When more than
% max_spikes remain, they are taken at even spacing through that train.

spike_times = double(spike_times(:));
spike_samples = double(spike_samples(:));
keep = isfinite(spike_times) & isfinite(spike_samples);
spike_times = spike_times(keep);
spike_samples = spike_samples(keep);
if isempty(spike_samples)
    samples = zeros(0, 1);
    return
end
if nargin < 5 || isempty(max_spikes) || max_spikes < 1
    max_spikes = numel(spike_samples);
end
if nargin < 6 || isempty(min_in_window)
    min_in_window = 20;
end
in_window = spikes_in_bouts(spike_times, bout_starts, hist_range);
if sum(in_window) >= min_in_window
    spike_samples = spike_samples(in_window);
end
samples = evenly_spaced(spike_samples, round(max_spikes));
end

function in_window = spikes_in_bouts(spike_times, bout_starts, hist_range)
in_window = false(size(spike_times));
if isempty(bout_starts) || numel(hist_range) < 2
    return
end
starts = bout_starts(:) + hist_range(1);
stops = bout_starts(:) + hist_range(2);
for index = 1:numel(starts)
    in_window = in_window | (spike_times >= starts(index) & spike_times <= stops(index));
end
end

function samples = evenly_spaced(samples, max_spikes)
if numel(samples) <= max_spikes
    return
end
index = unique(round(linspace(1, numel(samples), max_spikes)));
samples = samples(index);
end
