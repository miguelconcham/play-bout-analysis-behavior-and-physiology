function erp = mean_event_erp(traces)
%MEAN_EVENT_ERP  Mean waveform across trials that are finite at every sample.
%
% Incomplete rows (any NaN) are left out. An empty set of finite trials
% returns a NaN waveform with one column per input sample.

if isempty(traces)
    erp = nan(1, size(traces, 2));
    return
end
if ~ismatrix(traces)
    error('mean_event_erp:traces', 'traces must be trials by time.');
end

finite_trial = all(isfinite(traces), 2);
if ~any(finite_trial)
    erp = nan(1, size(traces, 2));
    return
end
erp = mean(traces(finite_trial, :), 1);
end
