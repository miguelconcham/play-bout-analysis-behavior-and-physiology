function mask = saturated_lfp_trials(traces, saturation_limit)
%SATURATED_LFP_TRIALS  Trials that hit the LFP acquisition rails.
%
% LFP_PAG is stored as int16. A trial is saturated when any sample reaches
% either rail: abs(sample) >= saturation_limit. The default limit is 32767,
% so both 32767 and -32768 are rejected. A trial with no finite samples is
% incomplete, not saturated.

if nargin < 2 || isempty(saturation_limit)
    saturation_limit = playbout_erp_defaults().saturation_limit;
end
if ~isscalar(saturation_limit) || ~isfinite(saturation_limit) || saturation_limit <= 0
    error('saturated_lfp_trials:limit', ...
        'saturation_limit must be a positive finite scalar.');
end
if isempty(traces)
    mask = false(size(traces, 1), 1);
    return
end
if ~ismatrix(traces)
    error('saturated_lfp_trials:traces', 'traces must be trials by time.');
end

mask = any(abs(double(traces)) >= saturation_limit, 2);
end
