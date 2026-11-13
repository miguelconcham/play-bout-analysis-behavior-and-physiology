function names = session_acg_names(kind)
%SESSION_ACG_NAMES  Curve names or lock-group names used by the estimate.
%
% kind is 'curves' or 'groups'. Delta-locked is a group name only. It is
% the trough-locked and peak-locked neurons together.

if strcmp(kind, 'groups')
    names = {'trough', 'peak', 'delta', 'nonentrained'};
else
    names = {'session', 'fooof_low', 'fooof_high', 'hilbert_low', 'hilbert_high'};
end
end
