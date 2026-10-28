function [traces, time_s, complete] = cut_event_lfp(lfp, event_times_s, window_s, sampling_rate)
%CUT_EVENT_LFP  Raw LFP epoch around each event time.
%
% Sample i of a recording is at time i/sampling_rate seconds. That is the
% clock used by GENERATE_PSTH_PLAY_BOUT. The sample nearest an event is
% kept, and an exact tie goes to the earlier sample.
%
% window_s is [start end] in seconds relative to the event and must contain
% the onset in its interior. Every returned row has the same samples. A row
% whose window is not entirely inside the recording is NaN, with complete
% false. Traces are the stored sample values, not power and not filtered.

if ~isvector(lfp)
    error('cut_event_lfp:trace', 'Pass one LFP channel as a vector.');
end
if nargin < 4 || ~isscalar(sampling_rate) || ~isfinite(sampling_rate) || sampling_rate <= 0
    error('cut_event_lfp:rate', 'sampling_rate must be a positive finite scalar.');
end
if numel(window_s) ~= 2 || any(~isfinite(window_s)) || window_s(1) >= window_s(2)
    error('cut_event_lfp:window', 'window_s must be [start end] in seconds.');
end
if ~(window_s(1) < 0 && window_s(2) > 0)
    error('cut_event_lfp:window', ...
        'The window must contain the onset, with a pre-onset and a post-onset side.');
end

lfp = double(lfp(:)).';
event_times_s = double(event_times_s(:));
sample_index = round(window_s(1) * sampling_rate):round(window_s(2) * sampling_rate);
if ~any(sample_index == 0)
    error('cut_event_lfp:window', ...
        'The window does not land on the onset sample at this sampling rate.');
end

time_s = sample_index / sampling_rate;
n_samples = numel(lfp);
n_events = numel(event_times_s);
n_epoch = numel(sample_index);
traces = nan(n_events, n_epoch);
complete = false(n_events, 1);

for event_n = 1:n_events
    center = nearest_lfp_sample(event_times_s(event_n), sampling_rate);
    if ~isscalar(center) || ~isfinite(center)
        continue
    end
    idx = center + sample_index;
    if idx(1) < 1 || idx(end) > n_samples
        continue
    end
    traces(event_n, :) = lfp(idx);
    complete(event_n) = true;
end
end

function index = nearest_lfp_sample(time_s, sampling_rate)
scaled = double(time_s) * sampling_rate;
index = floor(scaled);
fraction = scaled - index;
if fraction > 0.5
    index = index + 1;
end
if ~isfinite(index)
    index = NaN;
end
end
