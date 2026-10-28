function [erp, info] = clean_session_erp(traces, time_s, saturation_limit, baseline_window)
%CLEAN_SESSION_ERP  Drop saturated trials, baseline-correct, then average.
%
% A trial enters the session ERP only when every sample is finite and none
% of them reach saturation_limit. Each kept trial then loses its own mean
% on baseline_window. The session ERP is the mean of those trials.

if nargin < 3 || isempty(saturation_limit)
    saturation_limit = playbout_erp_defaults().saturation_limit;
end
if nargin < 4
    baseline_window = playbout_erp_defaults().baseline_window;
end

saturated = saturated_lfp_trials(traces, saturation_limit);
complete = all(isfinite(traces), 2);
if isempty(traces)
    complete = false(0, 1);
end
keep = complete & ~saturated;

kept_traces = traces;
if ~isempty(kept_traces)
    kept_traces(~keep, :) = NaN;
end
corrected = baseline_correct_trials(kept_traces, time_s, baseline_window);
erp = mean_event_erp(corrected);

info = struct();
info.n_complete = sum(complete);
info.n_saturated = sum(complete & saturated);
info.n_kept = sum(keep);
end
