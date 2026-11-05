function [epochs, time] = theta_envelope_epochs(envelope, rate_out, onsets, window, rate_in)
%THETA_ENVELOPE_EPOCHS  Envelope rows aligned to onsets.
%
% A row is NaN when its window leaves the trace or crosses a non-finite
% sample. `time` is seconds relative to onset.

if nargin < 5 || isempty(rate_in)
    rate_in = rate_out;
end
step = 1 / rate_out;
time = (window(1):step:window(2)).';
if ~isempty(time) && time(end) > window(2) + 0.1 * step
    time(end) = [];
end
offsets = round(time * rate_out);
n_events = numel(onsets);
epochs = nan(n_events, numel(time));
n_samples = numel(envelope);
for event_n = 1:n_events
    if ~isfinite(onsets(event_n))
        continue
    end
    center = theta_envelope_time_index(onsets(event_n), rate_out, rate_in);
    idx = center + offsets(:).';
    if any(idx < 1 | idx > n_samples)
        continue
    end
    snippet = envelope(idx);
    if any(~isfinite(snippet))
        continue
    end
    epochs(event_n, :) = snippet(:).';
end
end
