function corrected = baseline_correct_trials(traces, time_s, baseline_window)
%BASELINE_CORRECT_TRIALS  Subtract each trial's mean on the baseline window.
%
% baseline_window is [start end) in the same seconds as time_s. An empty
% window returns the traces unchanged. The end of the window is excluded,
% so the default [-2, 0) does not include the onset sample.

if isempty(baseline_window)
    corrected = traces;
    return
end
if numel(baseline_window) ~= 2 || any(~isfinite(baseline_window)) || ...
        baseline_window(1) >= baseline_window(2)
    error('baseline_correct_trials:window', ...
        'baseline_window must be [start end) in seconds.');
end
time_s = time_s(:).';
if numel(time_s) ~= size(traces, 2)
    error('baseline_correct_trials:time', ...
        'time_s must have one entry per trace column.');
end

baseline_index = time_s >= baseline_window(1) & time_s < baseline_window(2);
if ~any(baseline_index)
    error('baseline_correct_trials:window', ...
        'No samples fall inside [%.4g, %.4g) s.', ...
        baseline_window(1), baseline_window(2));
end
if size(traces, 1) == 0
    corrected = traces;
    return
end

baseline = mean(traces(:, baseline_index), 2, 'omitnan');
corrected = bsxfun(@minus, traces, baseline);
end
